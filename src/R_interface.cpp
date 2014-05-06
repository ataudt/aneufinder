#include "scalehmm.h"

// ---------------------------------------------------------------
// void R_univariate_hmm()
// This function takes parameters from R, creates a univariate HMM object, creates the distributions, runs the Baum-Welch and returns the result to R.
// ---------------------------------------------------------------
extern "C" {
void R_univariate_hmm(int* O, int* T, int* N, double* means, double* variances, int* maxiter, int* maxtime, double* eps, double* posteriors, double* A, double* proba, double* loglik, double* weights, double* initial_means, double* initial_variances, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads) {

	// Define logging level
// 	FILE* pFile = fopen("aneufinder.log", "w");
// 	Output2FILE::Stream() = pFile;
//  	FILELog::ReportingLevel() = FILELog::FromString("INFO");
 	FILELog::ReportingLevel() = FILELog::FromString("ITERATION");
//  	FILELog::ReportingLevel() = FILELog::FromString("DEBUG1");

	// Parallelization settings
	omp_set_num_threads(*num_threads);

	// Print some information
	FILE_LOG(logINFO) << "number of states = " << *N;
	FILE_LOG(logINFO) << "number of bins = " << *T;
	if (*maxiter < 0)
	{
		FILE_LOG(logINFO) << "maximum number of iterations = none";
	} else {
		FILE_LOG(logINFO) << "maximum number of iterations = " << *maxiter;
	}
	if (*maxtime < 0)
	{
		FILE_LOG(logINFO) << "maximum running time = none";
	} else {
		FILE_LOG(logINFO) << "maximum running time = " << *maxtime << " sec";
	}
	FILE_LOG(logINFO) << "epsilon = " << *eps;
	FILE_LOG(logDEBUG3) << "observation vector";
	for (int t=0; t<50; t++) {
		FILE_LOG(logDEBUG3) << "O["<<t<<"] = " << O[t];
	}

	// Create the HMM
	FILE_LOG(logDEBUG1) << "Creating a univariate HMM";
	ScaleHMM* hmm = new ScaleHMM(*T, *N);
	// Initialize the transition probabilities and proba
	hmm->initialize_transition_probs(initial_A, *use_initial_params);
	hmm->initialize_proba(initial_proba, *use_initial_params);
    
	// Calculate mean and variance of data
	double mean = 0, variance = 0;
	for(int t=0; t<*T; t++)
	{
		mean+= O[t];
	}
	mean = mean / *T;
	for(int t=0; t<*T; t++)
	{
		variance+= pow(O[t] - mean, 2);
	}
	variance = variance / *T;
	FILE_LOG(logINFO) << "data mean = " << mean << ", data variance = " << variance;		
	
	// Go through all states of the hmm and assign the density functions
	srand (clock());
	int rand1, rand2;
	double imean, ivariance;
	for (int istate=0; istate<*N; istate++)
	{

		if (*use_initial_params) {
			FILE_LOG(logINFO) << "Using given parameters for mean and variance";
			imean = initial_means[istate];
			ivariance = initial_variances[istate];
		} else {
			// Simple initialization
			imean = mean * pow(2, istate-2);
			ivariance = variance * pow(2, istate-2);
			initial_means[istate] = imean;
			initial_variances[istate] = ivariance;
		}
		FILE_LOG(logDEBUG3) << "imean = " << imean;
		FILE_LOG(logDEBUG3) << "ivariance = " << ivariance;

		if (istate == 0)
		{
			ZeroInflation *d = new ZeroInflation(O, *T); // delete is done inside ~ScaleHMM()
			FILE_LOG(logDEBUG1) << "Using "<< d->get_name() <<" for state " << istate;
			hmm->densityFunctions.push_back(d);
		}
		else
		{
			NegativeBinomial *d = new NegativeBinomial(O, *T, imean, ivariance); // delete is done inside ~ScaleHMM()
			FILE_LOG(logDEBUG1) << "Using "<< d->get_name() <<" for state " << istate;
			hmm->densityFunctions.push_back(d);
		}

	}

	// Do the Baum-Welch to estimate the parameters
	FILE_LOG(logDEBUG1) << "Starting Baum-Welch estimation";
	hmm->baumWelch(maxiter, maxtime, eps);
	FILE_LOG(logDEBUG1) << "Finished with Baum-Welch estimation";
	// Compute the posteriors and save results directly to the R pointer
	double** post = allocDoubleMatrix(*N, *T);
	hmm->get_posteriors(post);
	FILE_LOG(logDEBUG1) << "Recode posteriors into column representation";
	#pragma omp parallel for
	for (int iN=0; iN<*N; iN++)
	{
		for (int t=0; t<*T; t++)
		{
			posteriors[t + iN * (*T)] = post[iN][t];
		}
	}
	freeDoubleMatrix(post, *N);

	FILE_LOG(logDEBUG1) << "Return parameters";
	// also return the estimated transition matrix and the initial probs
	for (int i=0; i<*N; i++)
	{
		proba[i] = hmm->get_proba(i);
		for (int j=0; j<*N; j++)
		{
			A[i * (*N) + j] = hmm->get_A(i,j);
		}
	}

	// copy the estimated distribution params
	for (int i=0; i<*N; i++)
	{
		means[i] = hmm->densityFunctions[i]->get_mean();
		variances[i] = hmm->densityFunctions[i]->get_variance();
	}
	*loglik = hmm->get_logP();
	hmm->calc_weights(weights);
	
	FILE_LOG(logDEBUG1) << "Deleting the hmm";
	delete hmm;
}
} // extern C

