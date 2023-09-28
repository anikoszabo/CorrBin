#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP makeSmatrix(SEXP, SEXP);
extern SEXP MixReprodQ(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ReprodEstimates(SEXP, SEXP, SEXP);
extern SEXP ReprodISDM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"makeSmatrix",     (DL_FUNC) &makeSmatrix,     2},
    {"MixReprodQ",      (DL_FUNC) &MixReprodQ,      6},
    {"ReprodEstimates", (DL_FUNC) &ReprodEstimates, 3},
    {"ReprodISDM",      (DL_FUNC) &ReprodISDM,      7},
    {NULL, NULL, 0}
};

void R_init_CorrBin(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
