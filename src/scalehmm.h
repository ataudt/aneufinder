


#ifndef SCALEHMM_H
#define SCALEHMM_H

#include "utility.h"
#include "densities.h"
#include <cmath>
#include <R.h> // R_CheckUserInterrupt()
#include <vector> // storing density functions
#include <time.h> // time(), difftime()
#include <string> // strcmp

#if defined TARGET_OS_MAC || defined __APPLE__
#include <libiomp/omp.h> // parallelization options on mac
#elif defined __linux__ || defined _WIN32 || defined _WIN64
#include <omp.h> // parallelization options
#endif

class ScaleHMM  {

	public:
		// Constructor and Destructor
		ScaleHMM(int T, int N);
		ScaleHMM(int T, int N, int Nmod, double** densities);
		~ScaleHMM();

		// Member variables
		std::vector<Density*> densityFunctions; ///< density functions for each state

		// Methods
		void initialize_transition_probs(double* initial_A, bool use_initial_params);
		void initialize_proba(double* initial_proba, bool use_initial_params);
		void baumWelch();
		void EM(int* maxiter, int* maxtime, double* eps);
		void check_for_state_swap();
		std::vector<double> calc_weights();
		void calc_weights(double* weights);

		// Getters and Setters
		void get_posteriors(double** post);
		double get_posterior(int iN, int t);
		double get_proba(int i);
		double get_A(int i, int j);
		double get_logP();
		void set_cutoff(int cutoff);

	private:
		// Member variables
		int T; ///< length of observed sequence
		int* obs; ///< vector [T] of observations
		int N; ///< number of states
		int Nmod; ///< number of modifications / marks
		int cutoff; ///< a cutoff for observations
		double* sumgamma; ///< vector[N] of sum of posteriors (gamma values)
		double** sumxi; ///< matrix[N x N] of xi values
		double** gamma; ///< matrix[N x T] of posteriors
		double logP; ///< loglikelihood
		double dlogP; ///< difference in loglikelihood from one iteration to the next
		double** A; ///< matrix [N x N] of transition probabilities
		double* proba; ///< initial probabilities (length N)
		double* scalefactoralpha; ///< vector[T] of scaling factors
		double** scalealpha; ///< matrix [T x N] of forward probabilities
		double** scalebeta; ///<  matrix [T x N] of backward probabilities
		double** densities; ///< matrix [N x T] of density values
// 		double** tdensities; ///< matrix [T x N] of density values, for use in multivariate !increases speed, but on cost of RAM usage and that seems to be limiting
		time_t EMStartTime_sec; ///< start time of the EM in sec
		int EMTime_real; ///< elapsed time from start of the 0th iteration
		int sumdiff_state_last; ///< sum of the difference in the state 1 assignments from one iteration to the next
		double sumdiff_posterior; ///< sum of the difference in posterior (gamma) values from one iteration to the next
// 		bool use_tdens; ///< switch for using the tdensities in the calculations
		whichvariate xvariate; ///< enum which stores if UNIVARIATE or MULTIVARIATE

		// Methods
		void forward(); ///< calculate forward variables (alpha)
		void backward(); ///< calculate backward variables (beta)
		void calc_sumgamma();
		void calc_sumxi();
		void calc_loglikelihood();
		void calc_densities();
		void print_uni_iteration(int iteration);
		void print_multi_iteration(int iteration);
		void print_uni_params();
};

#endif
