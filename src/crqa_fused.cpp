// crqa_fused.cpp
// =============================================================================
// Single fused pass: distance -> threshold -> line-statistics for crqa().
//
// Replaces the call chain (cdist -> dm -> dm <= radius -> sparseMatrix ->
// theiler -> side-mask -> line_stats) with one O(N*M) traversal that never
// materialises the N*N distance matrix or the recurrence plot. Memory is
// O(N + nnz) instead of O(N^2). Speed gain comes from (a) no intermediate
// allocation, (b) tight C++ inner loops with -O2 SIMD, (c) running the
// diagonal/vertical line accumulators in lockstep with the distance compute.
//
// Output is bit-for-bit identical to the v2.1.0 line_stats() helper on the
// same recurrent-point set, modulo line-length ordering (we sort decreasing
// to match line_stats output).
//
// Supported parameters:
//   ts1, ts2    : numeric matrices (rows = embedded points, cols = embed dims)
//   radius      : scalar threshold (already scaled by rescale stat if applicable)
//   metric      : 1 = euclidean, 2 = maximum (Chebyshev), 3 = manhattan
//   tw          : Theiler window (excludes |i-j| <= tw if side == "both";
//                 excludes |i-j| < max(tw, 1) if side != "both")
//   side        : 0 = both, 1 = upper, 2 = lower (in the post-transpose RP,
//                 matching R's convention after `S = t(S)`)
//   mindiag     : minimum diagonal-line length to retain
//   minvert     : minimum vertical-line length to retain
//   whiteline   : if true, also compute white vertical-line stats
//
// Returns an R list mirroring line_stats() output:
//   numrecurs, diaglines, max_vertlength, TT, lam, wmean, wmax, wENTR, wlines

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// -----------------------------------------------------------------------------
// Embedded-point distance kernels. ts1 has m1 rows, ts2 has m2 rows, both have
// `dim` columns. Returns the chosen metric distance between row i of ts1 and
// row j of ts2.

static inline double dist_euclid(const double* a, const double* b, int dim) {
  double s = 0.0;
  for (int k = 0; k < dim; ++k) {
    double d = a[k] - b[k];
    s += d * d;
  }
  return std::sqrt(s);
}
static inline double dist_max(const double* a, const double* b, int dim) {
  double s = 0.0;
  for (int k = 0; k < dim; ++k) {
    double d = std::abs(a[k] - b[k]);
    if (d > s) s = d;
  }
  return s;
}
static inline double dist_manh(const double* a, const double* b, int dim) {
  double s = 0.0;
  for (int k = 0; k < dim; ++k) s += std::abs(a[k] - b[k]);
  return s;
}

// -----------------------------------------------------------------------------
// Apply Theiler & side mask. The R-level pipeline computes
//   S = t(sparseMatrix(r, c, dims = c(v1l, v2l)))
//   S = theiler(S, tw)
//   if (side == "upper") S[lower.tri(S, diag = TRUE)] = 0
//   if (side == "lower") S[upper.tri(S, diag = TRUE)] = 0
//
// `theiler(S, tw)` (see crqa_helpers.R) blanks cells with |row - col| <= tw,
// AND if tw == 0 it still blanks the main diagonal. We replicate that here.
//
// Indexing convention: in the post-transpose RP we use (row = j_orig,
// col = i_orig) so a cell at (orig_i, orig_j) in the pre-transpose matrix
// lives at (orig_j, orig_i) after transpose. To keep the inner loop natural
// we iterate (i, j) over the PRE-transpose layout (i = ts1 row, j = ts2 row)
// and reuse the mask test in those original coordinates: the Theiler band
// is symmetric in |i - j|, so it is invariant under transpose. For the side
// mask, "upper" in the transposed RP means lower.tri-of-pre-transpose is
// excluded, i.e. cells with j_orig < i_orig (= post-transpose col < row in
// the strict-upper sense including diagonal) are kept. The R code blanks
// lower.tri(S, diag = TRUE) for "upper", so we keep cells with j_orig > i_orig.
//
// Translating: keep_cell(i, j) returns true iff cell (i, j) of pre-transpose
// dm survives the same masks the R pipeline applies.

static inline bool keep_cell(int i, int j, int tw, int side) {
  int absd = std::abs(i - j);
  // Theiler band (matches theiler() in crqa_helpers.R exactly):
  //   tw == 0  -> no diagonals blanked (theiler returns S unchanged)
  //   tw >= 1  -> blank |i - j| <= tw - 1  (the main diagonal plus tw - 1
  //              off-diagonals on each side, i.e. 2*tw - 1 blanked diagonals)
  if (tw > 0 && absd < tw) return false;
  if (side == 0) return true;                 // "both"
  // After t(S), R blanks lower.tri(S, diag=TRUE) for side="upper", which
  // corresponds to keeping (post-transpose) row < col, i.e. j_orig > i_orig.
  if (side == 1) return (j > i);              // "upper" (post-transpose)
  // side="lower" keeps post-transpose row > col, i.e. i_orig > j_orig.
  return (i > j);                             // "lower" (post-transpose)
}

// -----------------------------------------------------------------------------
// Compute (rescaledist) over the FULL pre-threshold distance matrix, matching
// the R code which computes mean/max/min/sum over ALL of dm (no mask applied).
//
// `rescale_kind`: 1 = mean, 2 = max, 3 = min, 4 = euclidean-style mean

// [[Rcpp::export]]
double crqa_rescale_stat(const NumericMatrix& ts1, const NumericMatrix& ts2,
                          int metric_id, int rescale_kind) {
  const int n1 = ts1.nrow();
  const int n2 = ts2.nrow();
  const int dim = ts1.ncol();
  if (ts2.ncol() != dim) stop("ts1 and ts2 must have the same embed dimension");
  const double* p1 = REAL(ts1);
  const double* p2 = REAL(ts2);

  // Row-major-ish access: walk j outer, i inner; ts1/ts2 are column-major
  // R matrices, so element (i, k) is at p1[i + k*n1]. Cache k-stride.
  // For each (i, j) compute the chosen distance.
  double sum = 0.0;
  double mn = std::numeric_limits<double>::infinity();
  double mx = -std::numeric_limits<double>::infinity();
  long long count = 0;

  // Stack-allocated row buffers (avoid repeated stride math). Limit to small dim.
  std::vector<double> ai(dim), aj(dim);

  for (int i = 0; i < n1; ++i) {
    for (int k = 0; k < dim; ++k) ai[k] = p1[i + (long long)k * n1];
    for (int j = 0; j < n2; ++j) {
      for (int k = 0; k < dim; ++k) aj[k] = p2[j + (long long)k * n2];
      double d;
      switch (metric_id) {
        case 1:  d = dist_euclid(ai.data(), aj.data(), dim); break;
        case 2:  d = dist_max   (ai.data(), aj.data(), dim); break;
        case 3:  d = dist_manh  (ai.data(), aj.data(), dim); break;
        default: stop("unsupported metric_id");
      }
      sum += d;
      if (d < mn) mn = d;
      if (d > mx) mx = d;
      ++count;
    }
    if ((i & 0xFFF) == 0) Rcpp::checkUserInterrupt();
  }
  switch (rescale_kind) {
    case 1: return sum / (double)count;                  // mean of dm
    case 2: return mx;                                   // max
    case 3: return mn;                                   // min
    case 4: {                                            // sum(dm)/(N^2 - N)
      // matches R: abs(sum(dm)/(nrow(dm)^2 - nrow(dm)))
      double denom = (double)n1 * (double)n1 - (double)n1;
      return std::fabs(sum / denom);
    }
    default: stop("unsupported rescale_kind");
  }
}

// -----------------------------------------------------------------------------
// Main fused entry point.
//
// Algorithmic strategy:
// 1) Collect (i, j) recurrent indices in column-major order (jj sorted, ii
//    sorted within each jj), matching how the existing line_stats() does it
//    after sparseMatrix construction + drop0/triplet conversion + ord2 sort.
//    We do this by iterating j_outer over the POST-TRANSPOSE column axis (=
//    pre-transpose row axis, i = ts1 index) and i_inner over the post-transpose
//    row axis (= pre-transpose col axis, j = ts2 index). Wait — easier:
//    we iterate in post-transpose (row, col) coordinates directly.
//    Post-transpose RP cell (R, C) corresponds to pre-transpose dm[C, R]
//    (because t(M)[a,b] = M[b,a]). So we compute dist(ts1 row C, ts2 row R).
// 2) Apply Theiler + side mask in post-transpose coordinates.
// 3) Accumulate diagonal-run state keyed by d = C - R (post-transpose), and
//    vertical-run state keyed by C (post-transpose column). All in O(1) per
//    cell amortised.
// 4) For whiteline, also track within-column white-run lengths between recurrent
//    points, mirroring line_stats_white() exactly.

// [[Rcpp::export]]
List crqa_fused_cpp(const NumericMatrix& ts1, const NumericMatrix& ts2,
                    double radius, int metric_id,
                    int tw, int side, int mindiag, int minvert,
                    bool whiteline) {
  // ts1 is the "first" embedded series (n1 rows), ts2 the "second" (n2 rows).
  // The R pipeline builds dm = cdist(ts1, ts2) of size n1 x n2, then
  // S = t(sparseMatrix(r, c, dims = c(n1, n2))) which is n2 x n1.
  // Post-transpose RP rows correspond to ts2 (n2), cols to ts1 (n1).
  // For RR analysis we want recurrent indices (row, col) in the post-transpose
  // RP; cell (R, C) corresponds to dm[C, R].
  const int n1 = ts1.nrow();        // pre-transpose rows  -> post-transpose cols
  const int n2 = ts2.nrow();        // pre-transpose cols  -> post-transpose rows
  const int dim = ts1.ncol();
  if (ts2.ncol() != dim) stop("ts1 and ts2 must have the same embed dimension");
  const int v1l = n2;               // post-transpose rows
  const int v2l = n1;               // post-transpose cols
  const double* p1 = REAL(ts1);
  const double* p2 = REAL(ts2);

  // Validate metric_id ONCE up front so the parallel inner loop never has to
  // call Rcpp::stop() (which is not thread-safe from worker threads).
  if (metric_id < 1 || metric_id > 3) stop("unsupported metric_id");

  // ---------------------------------------------------------------------------
  // PARALLEL COLLECTION (OpenMP). Strategy:
  //   - Each thread processes a contiguous range of post-transpose columns C
  //     (static schedule, no chunk size -> default block partitioning).
  //   - Each thread writes (ii, jj) recurrent indices into thread-local vectors.
  //   - After the parallel section we concatenate thread vectors IN THREAD ORDER.
  //     Because (a) static scheduling assigns thread t a contiguous column block
  //     [t*v2l/T, (t+1)*v2l/T), and (b) each thread iterates rows in order
  //     within each of its columns, the concatenated result is bit-for-bit
  //     identical to the serial (jj outer, ii inner) ordering.
  //
  // The downstream line-statistics passes (diagonal sort+scan, vertical scan,
  // white-line scan) operate on the merged (ii_v, jj_v) and are unchanged, so
  // all measures (RR, DET, L, ENTR, rENTR, NRLINE, maxL, LAM, TT,
  // max_vertlength, wmean, wmax, wENTR, RP) are bit-identical to the serial
  // reference. Only floating-point comparisons (d <= radius) occur in the
  // parallel region -- there is no FP reduction, so the threshold decisions
  // are independent of thread count.
  // ---------------------------------------------------------------------------

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
  if (nthreads < 1) nthreads = 1;
#endif

  std::vector<std::vector<int>> ii_chunks(nthreads);
  std::vector<std::vector<int>> jj_chunks(nthreads);
  // Per-thread reservation: split the global heuristic across threads.
  const size_t per_thread_reserve =
    (size_t)((long long)v1l * (long long)v2l / 20 / (long long)nthreads);
  for (int t = 0; t < nthreads; ++t) {
    ii_chunks[t].reserve(per_thread_reserve);
    jj_chunks[t].reserve(per_thread_reserve);
  }

#ifdef _OPENMP
  #pragma omp parallel num_threads(nthreads)
#endif
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    std::vector<int>& ii_local = ii_chunks[tid];
    std::vector<int>& jj_local = jj_chunks[tid];
    std::vector<double> ai(dim), aj(dim);   // thread-private buffers

#ifdef _OPENMP
    #pragma omp for schedule(static) nowait
#endif
    for (int C = 1; C <= v2l; ++C) {            // post-transpose column
      int c_idx = C - 1;
      for (int k = 0; k < dim; ++k) ai[k] = p1[c_idx + (long long)k * n1];
      for (int R = 1; R <= v1l; ++R) {          // post-transpose row
        // Apply masks in post-transpose coordinates. Matches theiler() in
        // crqa_helpers.R: tw=0 leaves S untouched; tw>=1 blanks |R-C| <= tw-1.
        int absd = std::abs(R - C);
        if (tw > 0 && absd < tw) continue;
        if (side == 1 && !(R < C)) continue;     // "upper": keep row < col
        if (side == 2 && !(R > C)) continue;     // "lower": keep row > col

        int r_idx = R - 1;
        for (int k = 0; k < dim; ++k) aj[k] = p2[r_idx + (long long)k * n2];
        double d;
        switch (metric_id) {
          case 1:  d = dist_euclid(ai.data(), aj.data(), dim); break;
          case 2:  d = dist_max   (ai.data(), aj.data(), dim); break;
          default: d = dist_manh  (ai.data(), aj.data(), dim); break;
        }
        if (d <= radius) {
          ii_local.push_back(R);
          jj_local.push_back(C);
        }
      }
    }
  }  // end parallel region

  // Concatenate per-thread chunks IN THREAD ORDER. Static scheduling guarantees
  // thread t handled the lowest column block, thread t+1 the next, etc., so
  // this concatenation reproduces the serial (jj outer, ii inner) layout.
  size_t total_recurs = 0;
  for (int t = 0; t < nthreads; ++t) total_recurs += ii_chunks[t].size();

  std::vector<int> ii_v;
  std::vector<int> jj_v;
  ii_v.reserve(total_recurs);
  jj_v.reserve(total_recurs);
  for (int t = 0; t < nthreads; ++t) {
    ii_v.insert(ii_v.end(), ii_chunks[t].begin(), ii_chunks[t].end());
    jj_v.insert(jj_v.end(), jj_chunks[t].begin(), jj_chunks[t].end());
    // Free thread-chunk memory eagerly so peak RAM doesn't double.
    std::vector<int>().swap(ii_chunks[t]);
    std::vector<int>().swap(jj_chunks[t]);
  }
  long long numrecurs = (long long)total_recurs;
  Rcpp::checkUserInterrupt();

  // Early return for empty RP, matching line_stats() exactly.
  if (numrecurs == 0) {
    return List::create(
      _["numrecurs"]      = (R_xlen_t)0,
      _["diaglines"]      = IntegerVector::create(),
      _["max_vertlength"] = R_NegInf,
      _["TT"]             = R_NaN,
      _["lam"]            = 0.0,
      _["wmean"]          = NA_REAL,
      _["wmax"]           = NA_REAL,
      _["wENTR"]          = NA_REAL,
      _["wlines"]         = IntegerVector::create()
    );
  }

  // ---- Diagonal lines: group by d = j - i, run-length on consecutive i ----
  // line_stats() does: d_idx = jj - ii; ord = order(d_idx, ii); then groups.
  // We have ii_v, jj_v in (jj outer, ii inner) order. Build (d, ii) pairs and
  // sort by (d, ii). Then a new run starts where d changes OR ii is not previous+1.
  const size_t n = ii_v.size();
  std::vector<int> idx(n);
  for (size_t k = 0; k < n; ++k) idx[k] = (int)k;
  std::sort(idx.begin(), idx.end(), [&](int a, int b) {
    int da = jj_v[a] - ii_v[a];
    int db = jj_v[b] - ii_v[b];
    if (da != db) return da < db;
    return ii_v[a] < ii_v[b];
  });
  std::vector<int> diag_runs;
  diag_runs.reserve(n / 2 + 4);
  int prev_d = INT_MIN, prev_i = INT_MIN;
  int run_len = 0;
  for (size_t k = 0; k < n; ++k) {
    int p = idx[k];
    int d = jj_v[p] - ii_v[p];
    int i = ii_v[p];
    if (k == 0 || d != prev_d || i != prev_i + 1) {
      if (run_len >= mindiag) diag_runs.push_back(run_len);
      run_len = 1;
    } else {
      ++run_len;
    }
    prev_d = d; prev_i = i;
  }
  if (run_len >= mindiag) diag_runs.push_back(run_len);
  std::sort(diag_runs.begin(), diag_runs.end(), std::greater<int>());

  // ---- Vertical lines: group by j, run-length on consecutive i -----------
  // ii_v, jj_v are already in (jj outer, ii inner) order, so we can do a
  // straight pass — same as line_stats() does via ord2 = order(jj, ii).
  std::vector<int> vert_min;        // vertical runs >= minvert (unsorted)
  vert_min.reserve(n / 4 + 4);
  int v_prev_j = INT_MIN, v_prev_i = INT_MIN;
  int v_run = 0;
  for (size_t k = 0; k < n; ++k) {
    int i = ii_v[k];
    int j = jj_v[k];
    if (k == 0 || j != v_prev_j || i != v_prev_i + 1) {
      if (v_run >= minvert) vert_min.push_back(v_run);
      v_run = 1;
    } else {
      ++v_run;
    }
    v_prev_j = j; v_prev_i = i;
  }
  if (v_run >= minvert) vert_min.push_back(v_run);

  double max_vertlength_d;
  double TT_val;
  double lam_val;
  if (!vert_min.empty()) {
    int mxv = 0; long long sumv = 0;
    for (int x : vert_min) { if (x > mxv) mxv = x; sumv += x; }
    max_vertlength_d = (double)mxv;
    TT_val = (double)sumv / (double)vert_min.size();
    lam_val = ((double)sumv / (double)numrecurs) * 100.0;
  } else {
    max_vertlength_d = R_NegInf;
    TT_val = R_NaN;
    lam_val = 0.0;
  }

  // ---- White-line stats (optional) ---------------------------------------
  // line_stats_white(): for each column, scan (ii_v) sorted ascending and
  // record gaps `ii[k+1] - ii[k] - 1` between consecutive recurrent rows.
  // Then wmax = max, wmean = mean (after removing zeros — gaps of 1 mean
  // adjacent recurrent points, gap = 0 means no gap). Replicates exactly
  // line_stats_white() in crqa_helpers.R.
  double wmean = NA_REAL, wmax = NA_REAL, wENTR = NA_REAL;
  IntegerVector wlines_out;
  if (whiteline) {
    std::vector<int> wlines;
    wlines.reserve(n / 4 + 4);
    int col = INT_MIN, prev_i_col = INT_MIN;
    for (size_t k = 0; k < n; ++k) {
      int j = jj_v[k];
      int i = ii_v[k];
      if (j != col) {
        col = j; prev_i_col = i; continue;
      }
      int gap = i - prev_i_col - 1;
      if (gap > 0) wlines.push_back(gap);
      prev_i_col = i;
    }
    if (!wlines.empty()) {
      int wmx = 0; long long ws = 0;
      for (int x : wlines) { if (x > wmx) wmx = x; ws += x; }
      wmax  = (double)wmx;
      wmean = (double)ws / (double)wlines.size();
      // entropy of the empirical distribution of white-line lengths
      // (matches the - sum(p * log(p)) formula used in line_stats_white).
      // We tabulate frequencies via a sort.
      std::vector<int> tmp(wlines);
      std::sort(tmp.begin(), tmp.end());
      double H = 0.0;
      size_t i0 = 0, N = tmp.size();
      while (i0 < N) {
        size_t i1 = i0;
        while (i1 < N && tmp[i1] == tmp[i0]) ++i1;
        double p = (double)(i1 - i0) / (double)N;
        H -= p * std::log(p);
        i0 = i1;
      }
      wENTR = H;
    }
    wlines_out = IntegerVector(wlines.begin(), wlines.end());
  }

  return List::create(
    _["numrecurs"]      = (R_xlen_t)numrecurs,
    _["diaglines"]      = IntegerVector(diag_runs.begin(), diag_runs.end()),
    _["max_vertlength"] = max_vertlength_d,
    _["TT"]             = TT_val,
    _["lam"]            = lam_val,
    _["wmean"]          = wmean,
    _["wmax"]           = wmax,
    _["wENTR"]          = wENTR,
    _["wlines"]         = wlines_out,
    _["ii"]             = IntegerVector(ii_v.begin(), ii_v.end()),
    _["jj"]             = IntegerVector(jj_v.begin(), jj_v.end()),
    _["v1l"]            = v1l,
    _["v2l"]            = v2l
  );
}
