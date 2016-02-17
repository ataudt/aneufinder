// #include <Rinternals.h>
// #include <R_ext/Rdynload.h>
// #include <R_ext/Visibility.h>
// #include "R_interface.h"
// 
// 
// static const R_CMethodDef CEntries[]  = {
//     {"R_univariate_hmm", (DL_FUNC) &R_univariate_hmm, 24},
//     {"R_multivariate_hmm", (DL_FUNC) &R_multivariate_hmm, 18},
//     {"R_univariate_cleanup", (DL_FUNC) &R_univariate_cleanup, 0},
//     {"R_multivariate_cleanup", (DL_FUNC) &R_multivariate_cleanup, 1},
//     {NULL, NULL, 0}
// };
// 
// 
// void attribute_visible R_init_aneufinder(DllInfo *dll)
// {
//     R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
//     R_useDynamicSymbols(dll, FALSE);
//     R_forceSymbols(dll, TRUE);
// }
