//aneufinder - An R-package for CNV detection in whole-genome single cell sequencing data
//Copyright (C) 2015  Aaron Taudt
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.



#ifndef LogHMM_H
#define LogHMM_H

// #include <omp.h> // #pragma omp parallel
#include <R.h> // R_CheckUserInterrupt()
#include <vector> // storing density functions
#include <time.h> // time(), difftime()
#include "utility.h"
#include "densities.h"

class LogHMM  {

	public:
		// Constructor and Destructor
		LogHMM(int T, int N);
		~LogHMM();

		// Member variables
		std::vector<Density*> densityFunctions; ///< density functions for each state

		// Methods
		void initialize_transition_probs(double* initial_A, bool use_initial_params);
		void initialize_proba(double* initial_proba, bool use_initial_params);
		void baumWelch(int* maxiter, int* maxtime, double* eps);
		void check_for_state_swap();
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
		int cutoff; ///< a cutoff for observations
		double* sumgamma; ///< vector[N] of sum of posteriors (gamma values)
		double** sumxi; ///< matrix[N x N] of xi values
		double** gamma; ///< matrix[N x T] of posteriors
		double logP; ///< loglikelihood
		double dlogP; ///< difference in loglikelihood from one iteration to the next
		// Scaling approach
			double** A; ///< matrix [N x N] of transition probabilities
			double* proba; ///< initial probabilities (length N)
			double* scalefactoralpha; ///< vector[T] of scaling factors
			double** scalealpha; ///< matrix [T x N] of forward probabilities
			double** scalebeta; ///<  matrix [T x N] of backward probabilities
			double** densities; ///< matrix [N x T] of density values
		// Log approach
			double** logA; ///< matrix [N x N] of transition probabilities
			double* logproba; ///< initial probabilities (length N)
			double** logalpha; ///< matrix [T x N] of forward probabilities
			double** logbeta; ///<  matrix [T x N] of backward probabilities
			double** logdensities; ///< matrix [N x T] of density values
		// Miscellany
			time_t baumWelchStartTime_sec; ///< start time of the Baum-Welch in sec
			int baumWelchTime_real; ///< elapsed time from start of the 0th iteration
			int sumdiff_state_last; ///< sum of the difference in the state 1 assignments from one iteration to the next
			double sumdiff_posterior; ///< sum of the difference in posterior (gamma) values from one iteration to the next
	// 		bool use_tdens; ///< switch for using the tdensities in the calculations

		// Methods
		void forward(); ///< calculate forward variables (alpha)
		void backward(); ///< calculate backward variables (beta)
		void calc_sumgamma();
		void calc_sumxi();
		void calc_loglikelihood();
		void calc_densities();
		void print_uni_iteration(int iteration);
		void print_uni_params();
};

#endif
