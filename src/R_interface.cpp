#include "utility.h"
#include "scalehmm.h"
#include "loghmm.h"
#include <omp.h> // parallelization options
#include <string> // strcmp

static ScaleHMM* hmm; // declare as static outside the function because we only need one and this enables memory-cleanup on R_CheckUserInterrupt()
static double** multiD;

// ===================================================================================================================================================
// This function takes parameters from R, creates a univariate HMM object, creates the distributions, runs the Baum-Welch and returns the result to R.
// ===================================================================================================================================================
extern "C" {
void R_univariate_hmm(int* O, int* T, int* N, int* state_labels, double* size, double* prob, int* maxiter, int* maxtime, double* eps, int* states, double* A, double* proba, double* loglik, double* weights, int* distr_type, double* initial_size, double* initial_prob, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error, int* read_cutoff)
{

	// Define logging level
// 	FILE* pFile = fopen("chromStar.log", "w");
// 	Output2FILE::Stream() = pFile;
//  	FILELog::ReportingLevel() = FILELog::FromString("NONE");
//  	FILELog::ReportingLevel() = FILELog::FromString("DEBUG2");

// 	// Parallelization settings
	omp_set_num_threads(*num_threads);

	// Print some information
	//FILE_LOG(logINFO) << "number of states = " << *N;
	Rprintf("number of states = %d\n", *N);
	//FILE_LOG(logINFO) << "number of bins = " << *T;
	Rprintf("number of bins = %d\n", *T);
	if (*maxiter < 0)
	{
		//FILE_LOG(logINFO) << "maximum number of iterations = none";
		Rprintf("maximum number of iterations = none\n");
	} else {
		//FILE_LOG(logINFO) << "maximum number of iterations = " << *maxiter;
		Rprintf("maximum number of iterations = %d\n", *maxiter);
	}
	if (*maxtime < 0)
	{
		//FILE_LOG(logINFO) << "maximum running time = none";
		Rprintf("maximum running time = none\n");
	} else {
		//FILE_LOG(logINFO) << "maximum running time = " << *maxtime << " sec";
		Rprintf("maximum running time = %d sec\n", *maxtime);
	}
	//FILE_LOG(logINFO) << "epsilon = " << *eps;
	Rprintf("epsilon = %g\n", *eps);

	//FILE_LOG(logDEBUG3) << "observation vector";
	for (int t=0; t<50; t++) {
		//FILE_LOG(logDEBUG3) << "O["<<t<<"] = " << O[t];
	}

	// Flush Rprintf statements to console
	R_FlushConsole();

	// Create the HMM
	//FILE_LOG(logDEBUG1) << "Creating a univariate HMM";
	hmm = new ScaleHMM(*T, *N);
// 	LogHMM* hmm = new LogHMM(*T, *N);
	hmm->set_cutoff(*read_cutoff);
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
	//FILE_LOG(logINFO) << "data mean = " << mean << ", data variance = " << variance;		
	Rprintf("data mean = %g, data variance = %g\n", mean, variance);		
	
	// Create the emission densities and initialize
	for (int i_state=0; i_state<*N; i_state++)
	{
		if (distr_type[i_state] == 1)
		{
			//FILE_LOG(logDEBUG1) << "Using delta distribution for state " << i_state;
			ZeroInflation *d = new ZeroInflation(O, *T); // delete is done inside ~ScaleHMM()
			hmm->densityFunctions.push_back(d);
		}
		else if (distr_type[i_state] == 2)
		{
			//FILE_LOG(logDEBUG1) << "Using geometric distribution for state " << i_state;
			Geometric *d = new Geometric(O, *T, initial_prob[i_state]); // delete is done inside ~ScaleHMM()
			hmm->densityFunctions.push_back(d);
		}
		else if (distr_type[i_state] == 3)
		{
			//FILE_LOG(logDEBUG1) << "Using negative binomial for state " << i_state;
			NegativeBinomial *d = new NegativeBinomial(O, *T, initial_size[i_state], initial_prob[i_state]); // delete is done inside ~ScaleHMM()
			hmm->densityFunctions.push_back(d);
		}
		else if (distr_type[i_state] == 4)
		{
			//FILE_LOG(logDEBUG1) << "Using binomial for state " << i_state;
			NegativeBinomial *d = new NegativeBinomial(O, *T, initial_size[i_state], initial_prob[i_state]); // delete is done inside ~ScaleHMM()
			hmm->densityFunctions.push_back(d);
		}
		else
		{
			//FILE_LOG(logWARNING) << "Density not specified, using default negative binomial for state " << i_state;
			NegativeBinomial *d = new NegativeBinomial(O, *T, initial_size[i_state], initial_prob[i_state]);
			hmm->densityFunctions.push_back(d);
		}
	}

	// Flush Rprintf statements to console
	R_FlushConsole();

	// Do the Baum-Welch to estimate the parameters
	//FILE_LOG(logDEBUG1) << "Starting Baum-Welch estimation";
	try
	{
		hmm->baumWelch(maxiter, maxtime, eps);
	}
	catch (std::exception& e)
	{
		//FILE_LOG(logERROR) << "Error in Baum-Welch: " << e.what();
		Rprintf("Error in Baum-Welch: %s\n", e.what());
		if (strcmp(e.what(),"nan detected")==0) { *error = 1; }
		else { *error = 2; }
	}
		
	//FILE_LOG(logDEBUG1) << "Finished with Baum-Welch estimation";

// 	// Compute the posteriors and save results directly to the R pointer
// 	//FILE_LOG(logDEBUG1) << "Recode posteriors into column representation";
// 	#pragma omp parallel for
// 	for (int iN=0; iN<*N; iN++)
// 	{
// 		for (int t=0; t<*T; t++)
// 		{
// 			posteriors[t + iN * (*T)] = hmm->get_posterior(iN, t);
// 		}
// 	}

	// Compute the states from posteriors
	//FILE_LOG(logDEBUG1) << "Computing states from posteriors";
	int ind_max;
	std::vector<double> posterior_per_t(*N);
	for (int t=0; t<*T; t++)
	{
		for (int iN=0; iN<*N; iN++)
		{
			posterior_per_t[iN] = hmm->get_posterior(iN, t);
		}
		ind_max = std::distance(posterior_per_t.begin(), std::max_element(posterior_per_t.begin(), posterior_per_t.end()));
		states[t] = state_labels[ind_max];
	}

	//FILE_LOG(logDEBUG1) << "Return parameters";
	// also return the estimated transition matrix and the initial probs
	for (int i=0; i<*N; i++)
	{
		proba[i] = hmm->get_proba(i);
		for (int j=0; j<*N; j++)
		{
			A[i * (*N) + j] = hmm->get_A(j,i);
		}
	}

	// copy the estimated distribution params
	for (int i=0; i<*N; i++)
	{
		if (hmm->densityFunctions[i]->get_name() == NEGATIVE_BINOMIAL) 
		{
			NegativeBinomial* d = (NegativeBinomial*)(hmm->densityFunctions[i]);
			size[i] = d->get_size();
			prob[i] = d->get_prob();
		}
		else if (hmm->densityFunctions[i]->get_name() == GEOMETRIC)
		{
			Geometric* d = (Geometric*)(hmm->densityFunctions[i]);
			size[i] = 0;
			prob[i] = d->get_prob();
		}
		else if (hmm->densityFunctions[i]->get_name() == ZERO_INFLATION)
		{
			// These values for a Negative Binomial define a zero-inflation (delta distribution)
			size[i] = 0;
			prob[i] = 1;
		}
		else if (hmm->densityFunctions[i]->get_name() == BINOMIAL) 
		{
			Binomial* d = (Binomial*)(hmm->densityFunctions[i]);
			size[i] = d->get_size();
			prob[i] = d->get_prob();
		}
	}
	*loglik = hmm->get_logP();
	hmm->calc_weights(weights);
	
	//FILE_LOG(logDEBUG1) << "Deleting the hmm";
	delete hmm;
	hmm = NULL; // assign NULL to defuse the additional delete in on.exit() call
}
} // extern C


// =====================================================================================================================================================
// This function takes parameters from R, creates a multivariate HMM object, runs the Baum-Welch and returns the result to R.
// =====================================================================================================================================================
extern "C" {
void R_multivariate_hmm(double* D, int* T, int* N, int *Nmod, int* comb_states, int* maxiter, int* maxtime, double* eps, int* states, double* A, double* proba, double* loglik, double* initial_A, double* initial_proba, bool* use_initial_params, int* num_threads, int* error)
{

	// Define logging level {"ERROR", "WARNING", "INFO", "ITERATION", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4"}
// 	FILE* pFile = fopen("chromStar.log", "w");
// 	Output2FILE::Stream() = pFile;
//  	FILELog::ReportingLevel() = FILELog::FromString("ITERATION");
//  	FILELog::ReportingLevel() = FILELog::FromString("DEBUG2");
//  	FILELog::ReportingLevel() = FILELog::FromString("ERROR");

	// Parallelization settings
	omp_set_num_threads(*num_threads);

	// Print some information
	//FILE_LOG(logINFO) << "number of states = " << *N;
	Rprintf("number of states = %d\n", *N);
	//FILE_LOG(logINFO) << "number of bins = " << *T;
	Rprintf("number of bins = %d\n", *T);
	if (*maxiter < 0)
	{
		//FILE_LOG(logINFO) << "maximum number of iterations = none";
		Rprintf("maximum number of iterations = none\n");
	} else {
		//FILE_LOG(logINFO) << "maximum number of iterations = " << *maxiter;
		Rprintf("maximum number of iterations = %d\n", *maxiter);
	}
	if (*maxtime < 0)
	{
		//FILE_LOG(logINFO) << "maximum running time = none";
		Rprintf("maximum running time = none\n");
	} else {
		//FILE_LOG(logINFO) << "maximum running time = " << *maxtime << " sec";
		Rprintf("maximum running time = %d sec\n", *maxtime);
	}
	//FILE_LOG(logINFO) << "epsilon = " << *eps;
	Rprintf("epsilon = %g\n", *eps);
	//FILE_LOG(logINFO) << "number of modifications = " << *Nmod;
	Rprintf("number of modifications = %d\n", *Nmod);

	// Flush Rprintf statements to console
	R_FlushConsole();

	// Recode the densities vector to matrix representation
// 	clock_t clocktime = clock(), dtime;
	multiD = CallocDoubleMatrix(*N, *T);
	for (int iN=0; iN<*N; iN++)
	{
		for (int t=0; t<*T; t++)
		{
			multiD[iN][t] = D[iN*(*T)+t];
		}
	}
// 	dtime = clock() - clocktime;
// 	//FILE_LOG(logDEBUG1) << "recoding densities vector to matrix representation: " << dtime << " clicks";

	// Create the HMM
	//FILE_LOG(logDEBUG1) << "Creating the multivariate HMM";
	hmm = new ScaleHMM(*T, *N, *Nmod, multiD);
	// Initialize the transition probabilities and proba
	hmm->initialize_transition_probs(initial_A, *use_initial_params);
	hmm->initialize_proba(initial_proba, *use_initial_params);
	
	// Print logproba and A
// 	for (int iN=0; iN<*N; iN++)
// 	{
// 		//FILE_LOG(logDEBUG) << "proba["<<iN<<"] = " <<exp(hmm->logproba[iN]);
// 		for (int jN=0; jN<*N; jN++)
// 		{
// 			//FILE_LOG(logDEBUG) << "A["<<iN<<"]["<<jN<<"] = " << hmm->A[iN][jN];
// 		}
// 	}

	// Estimate the parameters
	//FILE_LOG(logDEBUG1) << "Starting Baum-Welch estimation";
	try
	{
		hmm->baumWelch(maxiter, maxtime, eps);
	}
	catch (std::exception& e)
	{
		//FILE_LOG(logERROR) << "Error in Baum-Welch: " << e.what();
		Rprintf("Error in Baum-Welch: %s\n", e.what());
		if (strcmp(e.what(),"nan detected")==0) { *error = 1; }
		else { *error = 2; }
	}
	//FILE_LOG(logDEBUG1) << "Finished with Baum-Welch estimation";
	
// 	// Compute the posteriors and save results directly to the R pointer
// 	//FILE_LOG(logDEBUG1) << "Recode posteriors into column representation";
// 	for (int iN=0; iN<*N; iN++)
// 	{
// 		for (int t=0; t<*T; t++)
// 		{
// 			posteriors[t + iN * (*T)] = hmm->get_posterior(iN, t);
// 		}
// 	}

	// Compute the states from posteriors
	//FILE_LOG(logDEBUG1) << "Computing states from posteriors";
	int ind_max;
	std::vector<double> posterior_per_t(*N);
	for (int t=0; t<*T; t++)
	{
		for (int iN=0; iN<*N; iN++)
		{
			posterior_per_t[iN] = hmm->get_posterior(iN, t);
		}
		ind_max = std::distance(posterior_per_t.begin(), std::max_element(posterior_per_t.begin(), posterior_per_t.end()));
		states[t] = comb_states[ind_max];
	}
	
	//FILE_LOG(logDEBUG1) << "Return parameters";
	// also return the estimated transition matrix and the initial probs
	for (int i=0; i<*N; i++)
	{
		proba[i] = hmm->get_proba(i);
		for (int j=0; j<*N; j++)
		{
				A[i * (*N) + j] = hmm->get_A(j,i);
		}
	}
	*loglik = hmm->get_logP();

	//FILE_LOG(logDEBUG1) << "Deleting the hmm";
	delete hmm;
	hmm = NULL; // assign NULL to defuse the additional delete in on.exit() call
// 	FreeDoubleMatrix(multiD, *N);
}
} // extern C


// =======================================================
// This function make a cleanup if anything was left over
// =======================================================
extern "C" {
void R_univariate_cleanup()
{
// 	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__; // This message will be shown if interrupt happens before start of C-code
	delete hmm;
}
}

extern "C" {
void R_multivariate_cleanup(int* N)
{
	delete hmm;
	FreeDoubleMatrix(multiD, *N);
}
}

