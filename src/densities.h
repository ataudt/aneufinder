#include <vector> // needed for storing the density functions
using std::vector;
#include <R.h> // M_PI and srand()
#include <Rmath.h> // digamma() and qnorm()
#include <gsl/gsl_sf_hyperg.h> //needed for the hypergeometric function
#include "logging.h" // FILE_LOG() capability

const double pi = M_PI;

#ifndef DENSITIES_H
#define DENSITIES_H

enum whichvariate {UNIVARIATE, MULTIVARIATE};
enum DensityName {ZI, NB, ZINB, Other};

class Density {
	public:
		virtual ~Density() {};
		virtual void calc_logdensities(double* logdensity) {};
		virtual void calc_densities(double* density) {};
		virtual void calc_logCDFs(double* logCDF) {};
		virtual void calc_CDFs(double* CDF) {};
		virtual void update(double* weight) {}; 
		virtual void copy(Density* other) {};
		virtual DensityName getType() {};
		virtual double getMean() {};
		virtual double getVariance() {};
		virtual double getLogDensityAt10Variance() {};

};  


// Maria's univariate densities
class ZiNB : public Density {
	public:
		ZiNB();
		ZiNB(int* observations, int T, double r, double p, double w);
		~ZiNB();
		void calc_logdensities(double* logdensity);
		void calc_densities(double* logdensity);
		void calc_logCDFs(double* logCDF);
		void calc_CDFs(double* CDF);
		double getMean();
		double getVariance();
		double logLikelihood();
		DensityName getType();
		void setR (double newR);
		double getR();
		void setP (double newP);
		double getP();
		void setW (double newW);
		double getW();
		void copy(Density* other);
		double p;
		double r;
		double w;
		double getLogDensityAt10Variance();
	private:
		int* O; // observations
		int T;//size of observations
		double* weight; // temp storage for weights in update
		int max_O;
		double* lxfactorials;	//precomputed factorials
};

class NegativeBinomial : public Density {
	public:
		NegativeBinomial();
		NegativeBinomial(int* observations, int T, double r, double p);
		~NegativeBinomial();
		void calc_logdensities(double* logdensity);
		void calc_densities(double* density);
		void calc_logCDFs(double* logCDF);
		void calc_CDFs(double* CDF);
		void update(double* weight);
		double getMean();
		double getVariance();
		DensityName getType();
		void setR (double newR);
		double getR();
		void setP (double newP);
		double getP();
		void copy(Density* other);
		double getLogDensityAt10Variance();
	private:
		int* O; // observations
		int T;
		double p;
		double r;
		int max_O;	// maximum number of reads
		double* lxfactorials;	//precomputed factorials
		vector<double>* cdf;
};



class OnlyZeros : public Density {
	public:
		OnlyZeros();
		OnlyZeros(int* observations, int T);
		~OnlyZeros();
		void calc_logdensities(double* logdensity);
		void calc_densities(double* density);
		void update(double* weight);
		double getMean();
		double getVariance();
		DensityName getType();
		void copy(Density* other);
		double getLogDensityAt10Variance();
	private:
		int* O; // observations
		int T;
};


class MVCopulaApproximation : public Density {
	public:
		MVCopulaApproximation(int** multiobservations, int T, vector<Density*> marginals, double* cor_matrix_inv, double cor_matrix_determinant);
		~MVCopulaApproximation();
		void calc_logdensities(double* logdensity);
		void calc_densities(double* density);
		DensityName getType();
	private:
		int** multiO;
		int T;
		vector<Density*> marginals;
		double* cor_matrix_inv;
		double cor_matrix_determinant;
		int Nmod;
};


class BernoulliProduct : public Density {
	public:
		BernoulliProduct(double** multiobservations, bool* binary_states, int T, int Nmod);
		~BernoulliProduct();
		void calc_logdensities(double* logdens);
		DensityName getType();
	private:
		double** multiO;
		bool* binary_states;
		int T;
		int Nmod;
};


#endif
