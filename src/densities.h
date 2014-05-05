#include <Rmath.h> // digamma() and qnorm()
#include "logging.h" // FILE_LOG() capability
#include "utility.h"

#ifndef DENSITIES_H
#define DENSITIES_H

class Density
{
	public:
		// Constructor and Destructor
		virtual ~Density() {};
		// Methods
		virtual void calc_logdensities(double* logdensity) {};
		virtual void calc_densities(double* density) {};
		virtual void update(double* weight) {}; 
		// Getter and Setter
		virtual double get_mean() {};
		virtual void set_mean() {};
		virtual double get_variance() {};
		virtual void set_variance () {};

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
		void update(double* weight);

		// Getter and Setter
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
		double var; ///< variance of the normal
		double sd; ///< standard deviation of the normal
};

#endif
