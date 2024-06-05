#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_SUB(jspd)(int *m, int *n, int *k, double *Adata, int *jb, double *Bdata, int *d, double *tv, int *na, int *nb, int *nd);

static const R_FortranMethodDef FortranEntries[] = {
    {"jspd", (DL_FUNC) &F77_SUB(jspd), 11},
    {NULL, NULL, 0}
};

void R_init_crqa(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
