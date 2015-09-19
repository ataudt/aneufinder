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



#include <Rmath.h> // dnorm(), dnbinom() and digamma() etc.
#include <vector> // storing density functions in MVCopula
#include "utility.h" // //FILE_LOG(), intMax()

#ifndef DENSITIES_H
#define DENSITIES_H

enum whichvariate {UNIVARIATE, MULTIVARIATE};
enum DensityName {ZERO_INFLATION, NORMAL, NEGATIVE_BINOMIAL, GEOMETRIC, POISSON, BINOMIAL, OTHER};

class Density
{
	public:
		// Constructor and Destructor
		virtual ~Density() {};
		// Methods
		virtual void calc_logdensities(double*) {};
		virtual void calc_densities(double*) {};
		virtual void update(double*) {}; 
		virtual void update_constrained(double**, int, int) {};
		// Getter and Setter
		virtual DensityName get_name() { return(OTHER); };
		virtual void set_name(DensityName) {};
		virtual double get_mean() { return(0); };
		virtual void set_mean(double) {};
		virtual double get_variance() { return(0); };
		virtual void set_variance (double) {};

};  


class Normal : public Density
{
	public:
		// Constructor and Destructor
		Normal(int* observations, int T, double mean, double variance);
		~Normal();

		// Methods
		void calc_densities(double* density);
		void calc_logdensities(double* logdensity);
		void update(double* weights);

		// Getter and Setter
		DensityName get_name();
		void set_name(DensityName name);
		double get_mean();
		void set_mean(double mean);
		double get_variance();
		void set_variance(double variance);
		double get_stdev();
		void set_stdev(double stdev);

	private:
		DensityName name; ///< name of the distribution
		int T; ///< length of observation vector
		int* obs; ///< vector [T] of observations
		double mean; ///< mean of the normal
		double variance; ///< variance of the normal
		double sd; ///< standard deviation of the normal
};


class Poisson : public Density
{
	public:
		// Constructor and Destructor
		Poisson(int* observations, int T, double lambda);
		~Poisson();

		// Methods
		DensityName get_name();
		void set_name(DensityName name);
		void calc_densities(double* density);
		void calc_logdensities(double* logdensity);
		void update(double* weights);
		void update_constrained(double** weights, int fromState, int toState);

		// Getter and Setter
		double get_mean();
		void set_mean(double mean);
		double get_variance();
		void set_variance(double variance);
		double get_lambda();

	private:
		// Member variables
		DensityName name; ///< name of the distribution
		int T; ///< length of observation vector
		int* obs; ///< vector [T] of observations
		double lambda; ///< lambda parameter of the poisson
		double mean; ///< mean of the poisson
		double variance; ///< variance of the poisson
		int max_obs; ///< maximum observation
		double* lxfactorials; ///< vector of precomputed factorials log(x!)

};


class NegativeBinomial : public Density
{
	public:
		// Constructor and Destructor
		NegativeBinomial(int* observations, int T, double size, double prob);
		~NegativeBinomial();

		// Methods
		DensityName get_name();
		void set_name(DensityName name);
		void calc_densities(double* density);
		void calc_logdensities(double* logdensity);
		void update(double* weights);
		void update_constrained(double** weights, int fromState, int toState);
		double fsize(double mean, double variance);
		double fprob(double mean, double variance);
		double fmean(double size, double prob);
		double fvariance(double size, double variance);

		// Getter and Setter
		double get_mean();
		void set_mean(double mean);
		double get_variance();
		void set_variance(double variance);
		double get_size();
		double get_prob();

	private:
		// Member variables
		DensityName name; ///< name of the distribution
		int T; ///< length of observation vector
		int* obs; ///< vector [T] of observations
		double size; ///< size parameter of the negative binomial
		double prob; ///< probability parameter of the negative binomial
		double mean; ///< mean of the negative binomial
		double variance; ///< variance of the negative binomial
		int max_obs; ///< maximum observation
		double* lxfactorials; ///< vector of precomputed factorials log(x!)

};


class Binomial : public Density
{
	public:
		// Constructor and Destructor
		Binomial(int* observations, int T, double size, double prob);
		~Binomial();

		// Methods
		DensityName get_name();
		void set_name(DensityName name);
		void calc_densities(double* density);
		void calc_logdensities(double* logdensity);
		void update(double* weights);
		void update_constrained(double** weights, int fromState, int toState);
		double fsize(double mean, double variance);
		double fprob(double mean, double variance);
		double fmean(double size, double prob);
		double fvariance(double size, double variance);

		// Getter and Setter
		double get_mean();
		void set_mean(double mean);
		double get_variance();
		void set_variance(double variance);
		double get_size();
		double get_prob();

	private:
		// Member variables
		DensityName name; ///< name of the distribution
		int T; ///< length of observation vector
		int* obs; ///< vector [T] of observations
		double size; ///< size parameter of the  binomial
		double prob; ///< probability parameter of the  binomial
		double mean; ///< mean of the  binomial
		double variance; ///< variance of the  binomial
		int max_obs; ///< maximum observation
		double* lxfactorials; ///< vector of precomputed factorials log(x!)

};


class ZeroInflation : public Density
{
	public:
		// Constructor and Destructor
		ZeroInflation(int* observations, int T);
		~ZeroInflation();

		// Methods
		DensityName get_name();
		void set_name(DensityName name);
		void calc_densities(double* density);
		void calc_logdensities(double* logdensity);
		void update(double* weights);

		// Getter and Setter
		double get_mean();
		void set_mean(double mean);
		double get_variance();
		void set_variance(double variance);

	private:
		DensityName name; ///< name of the distribution
		int T; ///< length of observation vector
		int* obs; ///< vector [T] of observations
};


class Geometric : public Density
{
	public:
		// Constructor and Destructor
		Geometric(int* observations, int T, double prob);
		~Geometric();

		// Methods
		DensityName get_name();
		void set_name(DensityName name);
		void calc_densities(double* density);
		void calc_logdensities(double* logdensity);
		void update(double* weights);
		double fprob(double mean, double variance);
		double fmean(double prob);
		double fvariance(double prob);

		// Getter and Setter
		double get_mean();
		void set_mean(double mean);
		double get_variance();
		void set_variance(double variance);
		double get_prob();

	private:
		// Member variables
		DensityName name; ///< name of the distribution
		int T; ///< length of observation vector
		int* obs; ///< vector [T] of observations
		int max_obs; ///< maximum observation
		double prob; ///< probability parameter of the geometric distribution
		double mean; ///< mean of the geometric distribution
		double variance; ///< variance of the geometric distribution

};


class MVCopulaApproximation : public Density {
	public:
		// Constructor and Destructor
		MVCopulaApproximation(int** multiobservations, int T, std::vector<Density*> marginals, double* cor_matrix_inv, double cor_matrix_determinant);
		~MVCopulaApproximation();
	
		// Methods

		// Getters and Setters
		DensityName get_name();
		void set_name(DensityName name);

	private:
		// Member variables
		DensityName name; ///< name of the distribution
		int Nmod; ///< number of modifications
		int** multi_obs; ///< matrix [Nmod x T] of observations
		int T; ///< length of observation vector
		std::vector<Density*> marginals; ///< vector [Nmod] of marginal distributions
		double* cor_matrix_inv; ///< vector with elements of the inverse of the correlation matrix
		double cor_matrix_determinant; ///< determinant of the correlation matrix
};


#endif
