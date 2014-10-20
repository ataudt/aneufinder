#include <Rmath.h> // dnorm(), dnbinom() and digamma() etc.
#include "utility.h" // FILE_LOG(), intMax()

#ifndef DENSITIES_H
#define DENSITIES_H

enum DensityName {ZERO_INFLATION, NORMAL, NEGATIVE_BINOMIAL, GEOMETRIC, POISSON, BINOMIAL};

class Density
{
	public:
		// Constructor and Destructor
		virtual ~Density() {};
		// Methods
		virtual void calc_logdensities(double* logdensity) {};
		virtual void calc_densities(double* density) {};
		virtual void update(double* weight) {}; 
		virtual void update_constrained(double** weights, int fromState, int toState) {};
		// Getter and Setter
		virtual DensityName get_name() {};
		virtual double get_mean() {};
		virtual void set_mean(double mean) {};
		virtual double get_variance() {};
		virtual void set_variance (double variance) {};

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
		double get_mean();
		void set_mean(double mean);
		double get_variance();
		void set_variance(double variance);
		double get_stdev();
		void set_stdev(double stdev);

	private:
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
		void calc_densities(double* density);
		void calc_logdensities(double* logdensity);
		void update(double* weights);

		// Getter and Setter
		double get_mean();
		void set_mean(double mean);
		double get_variance();
		void set_variance(double variance);

	private:
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
		int T; ///< length of observation vector
		int* obs; ///< vector [T] of observations
		int max_obs; ///< maximum observation
		double prob; ///< probability parameter of the geometric distribution
		double mean; ///< mean of the geometric distribution
		double variance; ///< variance of the geometric distribution

};


#endif
