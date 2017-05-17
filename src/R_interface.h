


#include "utility.h"
#include "scalehmm.h"
#include "loghmm.h"
#include <string> // strcmp

// #if defined TARGET_OS_MAC || defined __APPLE__
// #include <libiomp/omp.h> // parallelization options on mac
// #elif defined __linux__ || defined _WIN32 || defined _WIN64
// #include <omp.h> // parallelization options
// #endif

extern "C"
void univariate_hmm(int* O, int* T, int* N, int* state_labels, double* size, double* prob, int* maxiter, int* maxtime, double* eps, double* maxPosterior, int* states, double* A, double* proba, double* loglik, double* weights, int* distr_type, double* initial_size, double* initial_prob, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error, int* read_cutoff, int* algorithm, int* verbosity);

extern "C"
void multivariate_hmm(double* D, int* T, int* N, int *Nmod, int* comb_states, int* maxiter, int* maxtime, double* eps, int* states, double* A, double* proba, double* loglik, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error, int* algorithm, int* verbosity);

extern "C"
void univariate_cleanup();

extern "C"
void multivariate_cleanup(int* N);

extern "C"
void array2D_which_max(double* array2D, int* dim, int* ind_max, double* value_max);
