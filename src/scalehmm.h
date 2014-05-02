#ifndef SCALEHMM_H
#define SCALEHMM_H

#include "utility.h"
#include "densities.h"
#include "logging.h"
#include <omp.h>

class ScaleHMM  {

	public:
		int T; ///< length of observed sequence
		int N; ///< number of states
		int Nmod; ///< number of modifications / marks
		double** A; ///< matrix [N x N] of transition probabilities
		double* proba; ///< initial probabilities (length N)
		double logP; ///< loglikelihood
		int* O; ///< vector [T] of observations
		int** multiO; ///< matrix [Nmod x T] of observations
		vector<Density*> densityFunctions; ///< density functions for each state
		
		ScaleHMM();
		ScaleHMM(int T, int N);
		ScaleHMM(int T, int N, int Nmod);
		~ScaleHMM();
		void initialize_transition_probs(double* initial_A, bool use_initial_params);
		void initialize_proba(double* initial_proba, bool use_initial_params);
		void get_posteriors(double** post);
		void baumWelch(int* maxiter, int* maxtime, double* eps);
		void check_for_state_swap();
		void calc_weights(double* weights);
// 		void viterbi(int* path, int recompute);
		double get_proba(int i);

	private:
		void forward();
		void backward();
		void calc_sumgamma();
		void calc_sumxi();
		void calc_loglikelihood();
		void computeDensities();
		void print_uni_params();
		void print_multi_params();
		void print_uni_iteration(int iteration);
		void print_multi_iteration(int iteration);
		double* scalefactoralpha; ///< vector[T] of scaling factors
		double** scalealpha; ///< matrix [T x N] of forward probabilities
		double** scalebeta; ///<  matrix [T x N] of backward probabilities
		double** densities; ///< matrix [N x T] of density values
		double** tdensities; ///< matrix [T x N] of density values, for use in multivariate
		double* sumgamma; ///< vector[N] of sum of posteriors (gamma values)
		double** sumxi; ///< matrix[N x N] of xi values
		double** gamma; ///< matrix[N x T] of posteriors
		double dlogP; ///< difference in loglikelihood from one iteration to the next
		time_t baumWelchStartTime_sec; ///< start time of the Baum-Welch in sec
		int baumWelchTime_real; ///< elapsed time from start of the 0th iteration
		int sumdiff_state1; ///< sum of the difference in the state 1 assignments from one iteration to the next
		double sumdiff_posterior; ///< sum of the difference in posterior (gamma) values from one iteration to the next
		whichvariate xvariate; ///< enum which stores if UNIVARIATE or MULTIVARIATE
};

#endif
