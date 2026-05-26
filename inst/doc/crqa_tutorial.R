## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  eval      = TRUE
)
library(crqa)

## ----benchmark-fig, echo = FALSE, out.width = "100%", fig.cap = "Figure 1. Wall time per call for crqa v2.1.0 vs. other RQA implementations on the Rössler-system benchmark (Marwan & Kraemer 2023). Filled markers: measured on this hardware (4-core Ryzen-class CPU). Open markers: reference levels from Marwan & Kraemer 2023 Fig 3A. crqa v2.1.0 uses the fused C++ kernel with OpenMP enabled (default); thread count = all available cores. AccRQA 1.0.1 (Adámek et al. 2026) measured single-threaded — their library also offers OpenMP and CUDA paths not benchmarked here. The aRQA path is a fundamentally different algorithm (phase-space histogram) with O(N) memory and scales to N = 200 000 in seconds."----
knitr::include_graphics("figure3_speed_memory.png")

## ----install, eval = FALSE----------------------------------------------------
# install.packages("crqa")          # CRAN release
# # devtools::install_github("morenococo/crqa")  # development version

## ----rqa-cat------------------------------------------------------------------
data(crqa)

ans_rqa <- crqa(text, text,
                delay = 1, embed = 1, rescale = 0, radius = 0.0001,
                normalize = 0, mindiagline = 2, minvertline = 2,
                tw = 1, whiteline = FALSE, recpt = FALSE,
                side = "both", method = "rqa",
                metric = "euclidean", datatype = "categorical")

# Print all scalar measures (exclude the RP itself)
print(ans_rqa[!names(ans_rqa) %in% "RP"])

## ----rqa-plot, fig.width = 5, fig.height = 5, out.width = "60%"---------------
plot_rp(ans_rqa$RP,
        title = "Auto-RP: Wheels on the Bus",
        xlabel = "Word index", ylabel = "Word index")

## ----crqa-cat-----------------------------------------------------------------
listener <- eyemovement$listener
narrator <- eyemovement$narrator

ans_crqa <- crqa(narrator, listener,
                 delay = 1, embed = 1, rescale = 0, radius = 0.01,
                 normalize = 0, mindiagline = 2, minvertline = 2,
                 tw = 0, whiteline = FALSE, recpt = FALSE,
                 side = "both", method = "crqa",
                 metric = "euclidean", datatype = "categorical")

cat(sprintf("RR = %.2f%%   DET = %.2f%%   LAM = %.2f%%\n",
            ans_crqa$RR, ans_crqa$DET, ans_crqa$LAM))

## ----whiteline----------------------------------------------------------------
ans_wl <- crqa(narrator, narrator,
               delay = 1, embed = 1, rescale = 0, radius = 0.001,
               normalize = 0, mindiagline = 2, minvertline = 2,
               tw = 0, whiteline = TRUE, recpt = FALSE,
               side = "both", method = "rqa",
               metric = "euclidean", datatype = "categorical")

cat(sprintf("wmean = %.2f   wmax = %d   wENTR = %.3f\n",
            ans_wl$wmean, ans_wl$wmax, ans_wl$wENTR))

## ----mdcrqa-------------------------------------------------------------------

P1 <- cbind(handmovement$P1_TT_d, handmovement$P1_TT_n)
P2 <- cbind(handmovement$P2_TT_d, handmovement$P2_TT_n)

ans_md <- crqa(P1, P2,
               delay = 5, embed = 2, rescale = 0, radius = 0.3,
               normalize = 0, mindiagline = 2, minvertline = 2,
               tw = 0, whiteline = FALSE, recpt = FALSE,
               side = "both", method = "mdcrqa",
               metric = "euclidean", datatype = "continuous")

cat(sprintf("RR = %.2f%%   DET = %.2f%%   LAM = %.2f%%\n",
            ans_md$RR, ans_md$DET, ans_md$LAM))

## ----drp, fig.width = 6, fig.height = 4, out.width = "75%"--------------------
res <- drpfromts(narrator, listener,
                 windowsize = 100,
                 radius = 0.001, delay = 1, embed = 1,
                 rescale = 0, normalize = 0,
                 mindiagline = 2, minvertline = 2,
                 tw = 0, whiteline = FALSE, recpt = FALSE,
                 side = "both", method = "crqa",
                 metric = "euclidean", datatype = "categorical")

timecourse <- seq_along(res$profile) - ceiling(length(res$profile) / 2)
plot(timecourse, res$profile * 100,
     type = "l", lwd = 2,
     xlab = "Lag (samples)", ylab = "Recurrence rate (%)",
     main = "Diagonal cross-recurrence profile")
abline(v = res$maxlag - ceiling(length(res$profile) / 2),
       col = "red", lty = 2)

## ----wincrqa, fig.width = 6, fig.height = 4, out.width = "75%"----------------
win_res <- wincrqa(narrator, listener,
                   windowstep = 20, windowsize = 80,
                   delay = 1, embed = 1,
                   radius = 0.001, rescale = 0, normalize = 0,
                   mindiagline = 2, minvertline = 2,
                   tw = 0, whiteline = FALSE, recpt = FALSE,
                   side = "both", method = "crqa",
                   metric = "euclidean", datatype = "categorical",
                   trend = FALSE, workers = 1L)

plot(win_res$win, win_res$RR,
     type = "l", lwd = 2,
     xlab = "Window", ylab = "RR (%)",
     main = "Windowed RR (narrator-listener gaze coupling)")

## ----wincrqa-parallel, eval = FALSE-------------------------------------------
# win_par <- wincrqa(narrator, listener,
#                    windowstep = 20, windowsize = 80,
#                    delay = 1, embed = 1, radius = 0.001,
#                    workers = 4L)   # spread windows across 4 CPU cores

## ----arqa---------------------------------------------------------------------
# Simulate 5000 points of AR(1) data as a stand-in for a long series
set.seed(42)
x <- as.numeric(arima.sim(model = list(ar = 0.7), n = 5000))

res_arqa <- crqa(x, x,
                 delay = 1, embed = 3, radius = 0.5,
                 normalize = 2, mindiagline = 2,
                 method = "aRQA")

cat(sprintf("aRQA:  RR = %.3f%%   DET = %.2f%%   L = %.2f\n",
            res_arqa$RR, res_arqa$DET, res_arqa$L))

## ----arqa-compare, eval = FALSE-----------------------------------------------
# # Exact crqa for the same series at N = 1000 (comparable parameter set)
# x_short <- x[1:1000]
# 
# res_exact <- crqa(x_short, x_short,
#                   delay = 1, embed = 3, radius = 0.5,
#                   normalize = 2, mindiagline = 2,
#                   method = "rqa")
# 
# res_approx <- crqa(x_short, x_short,
#                    delay = 1, embed = 3, radius = 0.5,
#                    normalize = 2, mindiagline = 2,
#                    method = "aRQA")
# 
# cat(sprintf("Exact:  RR = %.3f%%  DET = %.2f%%\n", res_exact$RR,  res_exact$DET))
# cat(sprintf("aRQA:   RR = %.3f%%  DET = %.2f%%\n", res_approx$RR, res_approx$DET))

## ----optim, eval = FALSE------------------------------------------------------
# data(crqa)
# 
# P1 <- cbind(handmovement$P1_TT_d, handmovement$P1_TT_n)
# P2 <- cbind(handmovement$P2_TT_d, handmovement$P2_TT_n)
# 
# par <- list(method = "mdcrqa", metric = "euclidean",
#             maxlag = 20, radiusspan = 100, normalize = 0,
#             rescale = 4, mindiagline = 10, minvertline = 10,
#             tw = 0, whiteline = FALSE, recpt = FALSE,
#             side = "both", datatype = "continuous",
#             fnnpercent = NA, typeami = NA,
#             nbins = 50, criterion = "firstBelow",
#             threshold = 1, maxEmb = 20, numSamples = 500,
#             Rtol = 10, Atol = 2)
# 
# results <- optimizeParam(P1, P2, par, min.rec = 2, max.rec = 5)
# print(unlist(results))

## ----rr-denom-----------------------------------------------------------------
r_valid <- crqa(narrator, listener,
                delay = 1, embed = 1, rescale = 0, radius = 0.001,
                normalize = 0, mindiagline = 2, minvertline = 2,
                tw = 2, side = "both", method = "crqa",
                metric = "euclidean", datatype = "categorical",
                rr_denom = "valid")

r_full <- crqa(narrator, listener,
               delay = 1, embed = 1, rescale = 0, radius = 0.001,
               normalize = 0, mindiagline = 2, minvertline = 2,
               tw = 2, side = "both", method = "crqa",
               metric = "euclidean", datatype = "categorical",
               rr_denom = "full")

cat(sprintf("rr_denom='valid': RR = %.4f%%\n", r_valid$RR))
cat(sprintf("rr_denom='full':  RR = %.4f%%  (v2.0.7 convention)\n", r_full$RR))

