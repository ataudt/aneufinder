#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "R_interface.h"


R_NativePrimitiveArgType arg1[] = {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, LGLSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP};
R_NativePrimitiveArgType arg2[] = {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, LGLSXP, INTSXP, INTSXP, INTSXP, INTSXP};
R_NativePrimitiveArgType arg4[] = {INTSXP};
R_NativePrimitiveArgType arg5[] = {REALSXP, INTSXP, INTSXP, REALSXP};

static const R_CMethodDef CEntries[]  = {
    {"C_univariate_hmm", (DL_FUNC) &univariate_hmm, 26, arg1},
    {"C_multivariate_hmm", (DL_FUNC) &multivariate_hmm, 19, arg2},
    {"C_univariate_cleanup", (DL_FUNC) &univariate_cleanup, 0, NULL},
    {"C_multivariate_cleanup", (DL_FUNC) &multivariate_cleanup, 1, arg4},
    {"C_array2D_which_max", (DL_FUNC) &array2D_which_max, 4, arg5},
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
