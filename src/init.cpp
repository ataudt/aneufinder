#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "R_interface.h"


static const R_CMethodDef CEntries[]  = {
    {"C_univariate_hmm", (DL_FUNC) &univariate_hmm, 24, (R_NativePrimitiveArgType[24]) {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, LGLSXP, INTSXP, INTSXP, INTSXP, INTSXP}},
    {"C_multivariate_hmm", (DL_FUNC) &multivariate_hmm, 18, (R_NativePrimitiveArgType[18]) {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, LGLSXP, INTSXP, INTSXP, INTSXP}},
    {"C_univariate_cleanup", (DL_FUNC) &univariate_cleanup, 0, NULL},
    {"C_multivariate_cleanup", (DL_FUNC) &multivariate_cleanup, 1, (R_NativePrimitiveArgType[1]) {INTSXP}},
    {NULL, NULL, 0, NULL}
};


extern "C" {
void R_init_AneuFinder(DllInfo *dll)
{
	R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
// 	R_forceSymbols(dll, TRUE);
}
}
