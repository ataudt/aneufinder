


#include "utility.h"
#include "scalehmm.h"
#include "loghmm.h"
#include <string> // strcmp

#if defined TARGET_OS_MAC || defined __APPLE__
#include <libiomp/omp.h> // parallelization options on mac
#elif defined __linux__ || defined _WIN32 || defined _WIN64
#include <omp.h> // parallelization options
#endif

// void R_univariate_hmm(int* O, int* T, int* N, int* state_labels, double* size, double* prob, int* maxiter, int* maxtime, double* eps, int* states, double* A, double* proba, double* loglik, double* weights, int* distr_type, double* initial_size, double* initial_prob, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error, int* read_cutoff, int* algorithm);
// 
// void R_multivariate_hmm(double* D, int* T, int* N, int *Nmod, int* comb_states, int* maxiter, int* maxtime, double* eps, int* states, double* A, double* proba, double* loglik, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error, int* algorithm);
// 
// void R_univariate_cleanup();
// 
// void R_multivariate_cleanup(int* N);
// 
