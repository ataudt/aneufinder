#include "scalehmm.h"

// ============================================================
// Hidden Markov Model implemented with scaling strategy
// ============================================================

// Public =====================================================

// Constructor and Destructor ---------------------------------
ScaleHMM::ScaleHMM(int T, int N)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	//FILE_LOG(logDEBUG2) << "Initializing univariate ScaleHMM";
	this->xvariate = UNIVARIATE;
	this->T = T;
	this->N = N;
	this->A = CallocDoubleMatrix(N, N);
	this->scalefactoralpha = (double*) Calloc(T, double);
	this->scalealpha = CallocDoubleMatrix(T, N);
	this->scalebeta = CallocDoubleMatrix(T, N);
	this->densities = CallocDoubleMatrix(N, T);
// 	this->tdensities = CallocDoubleMatrix(T, N);
	this->proba = (double*) Calloc(N, double);
	this->gamma = CallocDoubleMatrix(N, T);
	this->sumgamma = (double*) Calloc(N, double);
	this->sumxi = CallocDoubleMatrix(N, N);
	this->logP = -INFINITY;
	this->dlogP = INFINITY;
	this->sumdiff_state_last = 0;
	this->sumdiff_posterior = 0.0;
// 	this->use_tdens = false;

}


ScaleHMM::ScaleHMM(int T, int N, int Nmod, double** densities)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	//FILE_LOG(logDEBUG2) << "Initializing multivariate ScaleHMM";
	this->xvariate = MULTIVARIATE;
	this->T = T;
	this->N = N;
	this->A = CallocDoubleMatrix(N, N);
	this->scalefactoralpha = (double*) Calloc(T, double);
	this->scalealpha = CallocDoubleMatrix(T, N);
	this->scalebeta = CallocDoubleMatrix(T, N);
	this->densities = densities;
	this->proba = (double*) Calloc(N, double);
	this->gamma = CallocDoubleMatrix(N, T);
	this->sumgamma = (double*) Calloc(N, double);
	this->sumxi = CallocDoubleMatrix(N, N);
	this->logP = -INFINITY;
	this->dlogP = INFINITY;
	this->Nmod = Nmod;

}

ScaleHMM::~ScaleHMM()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	FreeDoubleMatrix(this->A, this->N);
	Free(this->scalefactoralpha);
	FreeDoubleMatrix(this->scalealpha, this->T);
	FreeDoubleMatrix(this->scalebeta, this->T);
// 	FreeDoubleMatrix(this->tdensities, this->T);
	FreeDoubleMatrix(this->gamma, this->N);
	FreeDoubleMatrix(this->sumxi, this->N);
	Free(this->proba);
	Free(this->sumgamma);
	if (this->xvariate == UNIVARIATE)
	{
		FreeDoubleMatrix(this->densities, this->N);
		for (int iN=0; iN<this->N; iN++)
		{
			//FILE_LOG(logDEBUG1) << "Deleting density functions"; 
			delete this->densityFunctions[iN];
		}
	}
}

// Methods ----------------------------------------------------
void ScaleHMM::initialize_transition_probs(double* initial_A, bool use_initial_params)
{

	if (use_initial_params)
	{
		for (int iN=0; iN<this->N; iN++)
		{
			for (int jN=0; jN<this->N; jN++)
			{
				// convert from vector to matrix representation
				this->A[jN][iN] = initial_A[iN*this->N + jN];
			}
		}
	}
	else
	{
		double self = 0.9;
// 		self = 1.0 / this->N; // set to uniform
		double other = (1.0 - self) / (this->N - 1.0);
		for (int iN=0; iN<this->N; iN++)
		{
			for (int jN=0; jN<this->N; jN++)
			{
				if (iN == jN)
					this->A[iN][jN] = self;
				else
					this->A[iN][jN] = other;
				// Save value to initial A
				initial_A[jN*this->N + iN] = this->A[iN][jN];
			}
		}
	}
	
}

void ScaleHMM::initialize_proba(double* initial_proba, bool use_initial_params)
{

	if (use_initial_params)
	{
		for (int iN=0; iN<this->N; iN++)
		{
			this->proba[iN] = initial_proba[iN];
		}
	}
	else
	{
		for (int iN=0; iN<this->N; iN++)
		{
			this->proba[iN] = (double)1/this->N;
			// Save value to initial proba
			initial_proba[iN] = this->proba[iN];
		}
	}

}

void ScaleHMM::baumWelch(int* maxiter, int* maxtime, double* eps)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;

	double logPold = -INFINITY;
	double logPnew;
	double** gammaold = CallocDoubleMatrix(this->N, this->T);

	// Parallelization settings
// 	omp_set_nested(1);
	
	// measuring the time
	this->baumWelchStartTime_sec = time(NULL);

	// Print some initial information
	if (this->xvariate == UNIVARIATE)
	{
		//FILE_LOG(logINFO) << "";
		//FILE_LOG(logINFO) << "INITIAL PARAMETERS";
// 		this->print_uni_params();
		this->print_uni_iteration(0);
	}
	else if (this->xvariate == MULTIVARIATE)
	{
		this->print_multi_iteration(0);
	}

	R_CheckUserInterrupt();

	// Do the Baum-Welch
	int iteration = 0;
	while (((this->baumWelchTime_real < *maxtime) or (*maxtime < 0)) and ((iteration < *maxiter) or (*maxiter < 0)))
	{

		iteration++;
		
		if (this->xvariate == UNIVARIATE)
		{
			//FILE_LOG(logDEBUG1) << "Calling calc_densities() from baumWelch()";
			try { this->calc_densities(); } catch(...) { throw; }
			R_CheckUserInterrupt();
		}

		//FILE_LOG(logDEBUG1) << "Calling forward() from baumWelch()";
		try { this->forward(); } catch(...) { throw; }
		R_CheckUserInterrupt();

		//FILE_LOG(logDEBUG1) << "Calling backward() from baumWelch()";
		try { this->backward(); } catch(...) { throw; }
		R_CheckUserInterrupt();

		//FILE_LOG(logDEBUG1) << "Calling calc_loglikelihood() from baumWelch()";
		this->calc_loglikelihood();
		logPnew = this->logP;
		if(isnan(logPnew))
		{
			//FILE_LOG(logERROR) << "logPnew = " << logPnew;
			throw nan_detected;
		}
		this->dlogP = logPnew - logPold;

		//FILE_LOG(logDEBUG1) << "Calling calc_sumxi() from baumWelch()";
		this->calc_sumxi();
		R_CheckUserInterrupt();

		//FILE_LOG(logDEBUG1) << "Calling calc_sumgamma() from baumWelch()";
		this->calc_sumgamma();
		R_CheckUserInterrupt();

		if (this->xvariate == UNIVARIATE)
		{
// 			clock_t clocktime = clock(), dtime;
			// difference in posterior
			//FILE_LOG(logDEBUG1) << "Calculating differences in posterior in baumWelch()";
			double postsum = 0.0;
			for (int t=0; t<this->T; t++)
			{
				for (int iN=0; iN<this->N; iN++)
				{
					postsum += fabs(this->gamma[iN][t] - gammaold[iN][t]);
					gammaold[iN][t] = this->gamma[iN][t];
				}
			}
			this->sumdiff_posterior = postsum;
// 			dtime = clock() - clocktime;
// 			//FILE_LOG(logDEBUG) << "differences in posterior: " << dtime << " clicks";
		}

		R_CheckUserInterrupt();

		// Print information about current iteration
		if (this->xvariate == UNIVARIATE)
		{
			this->print_uni_iteration(iteration);
		}
		else if (this->xvariate == MULTIVARIATE)
		{
			this->print_multi_iteration(iteration);
		}

		// Check convergence
		if(this->dlogP < *eps) //it has converged
		{
			//FILE_LOG(logINFO) << "Convergence reached!\n";
			Rprintf("Convergence reached!\n");
			this->check_for_state_swap();
			break;
		}
		else
		{ // not converged
			this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
			if (iteration == *maxiter)
			{
				//FILE_LOG(logINFO) << "Maximum number of iterations reached!";
				Rprintf("Maximum number of iterations reached!\n");
				this->check_for_state_swap();
			}
			else if ((this->baumWelchTime_real >= *maxtime) and (*maxtime >= 0))
			{
				//FILE_LOG(logINFO) << "Exceeded maximum time!";
				Rprintf("Exceeded maximum time!\n");
				this->check_for_state_swap();
			}
			logPold = logPnew;
		}
		
// 		// Check weights
// 		double weights [this->N];
// 		double sumweights=0;
// 		this->calc_weights(weights);
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			//FILE_LOG(logERROR) << "weights["<<iN<<"] = " << weights[iN];
// 			sumweights += weights[iN];
// 		}
// 		//FILE_LOG(logERROR) << "sumweights = " << sumweights;
// 		

		// Updating initial probabilities proba and transition matrix A
		for (int iN=0; iN<this->N; iN++)
		{
			this->proba[iN] = this->gamma[iN][0];
			//FILE_LOG(logDEBUG4) << "sumgamma["<<iN<<"] = " << sumgamma[iN];
			if (this->sumgamma[iN] == 0)
			{
				//FILE_LOG(logINFO) << "Not reestimating A["<<iN<<"][x] because sumgamma["<<iN<<"] = 0";
				Rprintf("Not reestimating A[%d][x] because sumgamma[%d] = 0\n", iN, iN);
			}
			else
			{
				for (int jN=0; jN<this->N; jN++)
				{
					//FILE_LOG(logDEBUG4) << "sumxi["<<iN<<"]["<<jN<<"] = " << sumxi[iN][jN];
					this->A[iN][jN] = this->sumxi[iN][jN] / this->sumgamma[iN];
					if (isnan(this->A[iN][jN]))
					{
						//FILE_LOG(logERROR) << "updating transition probabilities";
						//FILE_LOG(logERROR) << "A["<<iN<<"]["<<jN<<"] = " << A[iN][jN];
						//FILE_LOG(logERROR) << "sumxi["<<iN<<"]["<<jN<<"] = " << sumxi[iN][jN];
						//FILE_LOG(logERROR) << "sumgamma["<<iN<<"] = " << sumgamma[iN];
						throw nan_detected;
					}
				}
			}
		}

		if (this->xvariate == UNIVARIATE)
		{
// 			clock_t clocktime = clock(), dtime;
// 
// 			// Update all distributions independantly
// 			for (int iN=0; iN<this->N; iN++)
// 			{
// 				this->densityFunctions[iN]->update(this->gamma[iN]);
// 			}

			// Update distribution of state 'null-mixed' and 'monosomic', set others as multiples of 'monosomic'
			// This loop assumes that the negative binomial states come last and are consecutive
			int xsomy = 1;
			for (int iN=0; iN<this->N; iN++)
			{
				if (this->densityFunctions[iN]->get_name() == ZERO_INFLATION) {}
				if (this->densityFunctions[iN]->get_name() == GEOMETRIC)
				{
					this->densityFunctions[iN]->update(this->gamma[iN]);
				}
				if (this->densityFunctions[iN]->get_name() == NEGATIVE_BINOMIAL)
				{
					if (xsomy==1)
					{
						//FILE_LOG(logDEBUG1) << "mean(state="<<iN<<") = " << this->densityFunctions[iN]->get_mean() << ", var(state="<<iN<<") = " << this->densityFunctions[iN]->get_variance();
						this->densityFunctions[iN]->update_constrained(this->gamma, iN, this->N);
						double mean1 = this->densityFunctions[iN]->get_mean();
						double variance1 = this->densityFunctions[iN]->get_variance();
						//FILE_LOG(logDEBUG1) << "mean(state="<<iN<<") = " << this->densityFunctions[iN]->get_mean() << ", var(state="<<iN<<") = " << this->densityFunctions[iN]->get_variance();
						// Set others as multiples
						for (int jN=iN+1; jN<this->N; jN++)
						{
							this->densityFunctions[jN]->set_mean(mean1 * (jN-iN+1));
							this->densityFunctions[jN]->set_variance(variance1 * (jN-iN+1));
							//FILE_LOG(logDEBUG1) << "mean(state="<<jN<<") = " << this->densityFunctions[jN]->get_mean() << ", var(state="<<jN<<") = " << this->densityFunctions[jN]->get_variance();
						}
						break;
					}
					xsomy++;
				}
			}
// 			dtime = clock() - clocktime;
// 			//FILE_LOG(logDEBUG) << "updating distributions: " << dtime << " clicks";
			R_CheckUserInterrupt();
		}

	} /* main loop end */
    
    
	//Print the last results
	//FILE_LOG(logINFO) << "";
	//FILE_LOG(logINFO) << "FINAL ESTIMATION RESULTS";
// 	this->print_uni_params();

	// free memory
	FreeDoubleMatrix(gammaold, this->N);

	// Return values
	*maxiter = iteration;
	*eps = this->dlogP;
	this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
	*maxtime = this->baumWelchTime_real;
}

void ScaleHMM::check_for_state_swap()
{
// TODO
}

std::vector<double> ScaleHMM::calc_weights()
{
	std::vector<double> weights(this->N);
	#pragma omp parallel for
	for (int iN=0; iN<this->N; iN++)
	{
		// Do not use weights[iN] = ( this->sumgamma[iN] + this->gamma[iN][T-1] ) / this->T; here, since states are swapped and gammas not
		double sum_over_gammas_per_state = 0;
		for (int t=0; t<this->T; t++)
		{
			sum_over_gammas_per_state += this->gamma[iN][t];
		}
		weights[iN] = sum_over_gammas_per_state / this->T;
	}
	return(weights);
}

void ScaleHMM::calc_weights(double* weights)
{
	#pragma omp parallel for
	for (int iN=0; iN<this->N; iN++)
	{
		// Do not use weights[iN] = ( this->sumgamma[iN] + this->gamma[iN][T-1] ) / this->T; here, since states are swapped and gammas not
		double sum_over_gammas_per_state = 0;
		for (int t=0; t<this->T; t++)
		{
			sum_over_gammas_per_state += this->gamma[iN][t];
		}
		weights[iN] = sum_over_gammas_per_state / this->T;
	}
}

// Getters and Setters ----------------------------------------
void ScaleHMM::get_posteriors(double** post)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int iN=0; iN<this->N; iN++)
	{
		for (int t=0; t<this->T; t++)
		{
			post[iN][t] = this->gamma[iN][t];
		}
	}
}

double ScaleHMM::get_posterior(int iN, int t)
{
	//FILE_LOG(logDEBUG3) << __PRETTY_FUNCTION__;
	return(this->gamma[iN][t]);
}

double ScaleHMM::get_proba(int i)
{
	return( this->proba[i] );
}

double ScaleHMM::get_A(int i, int j)
{
	return( this->A[i][j] );
}

double ScaleHMM::get_logP()
{
	return( this->logP );
}

void ScaleHMM::set_cutoff(int cutoff)
{
	this->cutoff = cutoff;
}

// Private ====================================================
// Methods ----------------------------------------------------
void ScaleHMM::forward()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
//	clock_t time = clock(), dtime;

// 	if (not this->use_tdens)
// 	{

		std::vector<double> alpha(this->N);
		// Initialization
		this->scalefactoralpha[0] = 0.0;
		for (int iN=0; iN<this->N; iN++)
		{
			alpha[iN] = this->proba[iN] * this->densities[iN][0];
			//FILE_LOG(logDEBUG4) << "alpha["<<iN<<"] = " << alpha[iN];
			this->scalefactoralpha[0] += alpha[iN];
		}
		//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<0<<"] = " << scalefactoralpha[0];
		for (int iN=0; iN<this->N; iN++)
		{
			this->scalealpha[0][iN] = alpha[iN] / this->scalefactoralpha[0];
			//FILE_LOG(logDEBUG4) << "scalealpha["<<0<<"]["<<iN<<"] = " << scalealpha[0][iN];
		}
		// Induction
		for (int t=1; t<this->T; t++)
		{
			this->scalefactoralpha[t] = 0.0;
			for (int iN=0; iN<this->N; iN++)
			{
				double helpsum = 0.0;
				for (int jN=0; jN<this->N; jN++)
				{
					helpsum += this->scalealpha[t-1][jN] * this->A[jN][iN];
				}
				alpha[iN] = helpsum * this->densities[iN][t];
				//FILE_LOG(logDEBUG4) << "alpha["<<iN<<"] = " << alpha[iN];
				this->scalefactoralpha[t] += alpha[iN];
			}
			//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<t<<"] = " << scalefactoralpha[t];
			for (int iN=0; iN<this->N; iN++)
			{
				this->scalealpha[t][iN] = alpha[iN] / this->scalefactoralpha[t];
				//FILE_LOG(logDEBUG4) << "scalealpha["<<t<<"]["<<iN<<"] = " << scalealpha[t][iN];
				if(isnan(this->scalealpha[t][iN]))
				{
					//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
					for (int jN=0; jN<this->N; jN++)
					{
						//FILE_LOG(logERROR) << "scalealpha["<<t-1<<"]["<<jN<<"] = " << scalealpha[t-1][jN];
						//FILE_LOG(logERROR) << "A["<<iN<<"]["<<jN<<"] = " << A[iN][jN];
					}
					//FILE_LOG(logERROR) << "scalefactoralpha["<<t<<"] = "<<scalefactoralpha[t] << ", densities = "<<densities[iN][t];
					//FILE_LOG(logERROR) << "scalealpha["<<t<<"]["<<iN<<"] = " << scalealpha[t][iN];
					throw nan_detected;
				}
			}
		}

// 	}
// 	else if (this->use_tdens)
// 	{
// 
// 		double alpha [this->N];
// 		// Initialization
// 		this->scalefactoralpha[0] = 0.0;
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			alpha[iN] = this->proba[iN] * this->tdensities[0][iN];
// 			//FILE_LOG(logDEBUG4) << "alpha["<<iN<<"] = " << alpha[iN];
// 			this->scalefactoralpha[0] += alpha[iN];
// 		}
// 		//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<0<<"] = " << scalefactoralpha[0];
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			this->scalealpha[0][iN] = alpha[iN] / this->scalefactoralpha[0];
// 			//FILE_LOG(logDEBUG4) << "scalealpha["<<0<<"]["<<iN<<"] = " << scalealpha[0][iN];
// 		}
// 		// Induction
// 		for (int t=1; t<this->T; t++)
// 		{
// 			this->scalefactoralpha[t] = 0.0;
// 			for (int iN=0; iN<this->N; iN++)
// 			{
// 				double helpsum = 0.0;
// 				for (int jN=0; jN<this->N; jN++)
// 				{
// 					helpsum += this->scalealpha[t-1][jN] * this->A[jN][iN];
// 				}
// 				alpha[iN] = helpsum * this->tdensities[t][iN];
// 				//FILE_LOG(logDEBUG4) << "alpha["<<iN<<"] = " << alpha[iN];
// 				this->scalefactoralpha[t] += alpha[iN];
// 			}
// 			//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<t<<"] = " << scalefactoralpha[t];
// 			for (int iN=0; iN<this->N; iN++)
// 			{
// 				this->scalealpha[t][iN] = alpha[iN] / this->scalefactoralpha[t];
// 				//FILE_LOG(logDEBUG4) << "scalealpha["<<t<<"]["<<iN<<"] = " << scalealpha[t][iN];
// 				if(isnan(this->scalealpha[t][iN]))
// 				{
// 					//FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
// 					for (int jN=0; jN<this->N; jN++)
// 					{
// 						//FILE_LOG(logWARNING) << "scalealpha["<<t-1<<"]["<<jN<<"] = " << scalealpha[t-1][jN];
// 						//FILE_LOG(logWARNING) << "A["<<iN<<"]["<<jN<<"] = " << A[iN][jN];
// 					}
// 					//FILE_LOG(logWARNING) << "scalefactoralpha["<<t<<"] = "<<scalefactoralpha[t] << ", tdensities = "<<tdensities[t][iN];
// 					//FILE_LOG(logWARNING) << "scalealpha["<<t<<"]["<<iN<<"] = " << scalealpha[t][iN];
// 					exit(1);
// 				}
// 			}
// 		}
// 
// 	}

//	dtime = clock() - time;
//	//FILE_LOG(logDEBUG) << "forward(): " << dtime << " clicks";
}

void ScaleHMM::backward()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
//	clock_t time = clock(), dtime;

// 	if (not this->use_tdens)
// 	{

		std::vector<double> beta(this->N);
		// Initialization
		for (int iN=0; iN<this->N; iN++)
		{
			beta[iN] = 1.0;
			//FILE_LOG(logDEBUG4) << "beta["<<iN<<"] = " << beta[iN];
		}
		//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<T-1<<"] = " << scalefactoralpha[T-1];
		for (int iN=0; iN<this->N; iN++)
		{
			this->scalebeta[T-1][iN] = beta[iN] / this->scalefactoralpha[T-1];
			//FILE_LOG(logDEBUG4) << "scalebeta["<<T-1<<"]["<<iN<<"] = " << scalebeta[T-1][iN];
		}
		// Induction
		for (int t=this->T-2; t>=0; t--)
		{
			for (int iN=0; iN<this->N; iN++)
			{
				beta[iN] = 0.0;
				for(int jN=0; jN<this->N; jN++)
				{
					beta[iN] += this->A[iN][jN] * this->densities[jN][t+1] * this->scalebeta[t+1][jN];
				}
				//FILE_LOG(logDEBUG4) << "beta["<<iN<<"] = " << beta[iN];
			}
			//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<t<<"] = " << scalefactoralpha[t];
			for (int iN=0; iN<this->N; iN++)
			{
				this->scalebeta[t][iN] = beta[iN] / this->scalefactoralpha[t];
				//FILE_LOG(logDEBUG4) << "scalebeta["<<t<<"]["<<iN<<"] = " << scalebeta[t][iN];
				if (isnan(this->scalebeta[t][iN]))
				{
					//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
					for (int jN=0; jN<this->N; jN++)
					{
						//FILE_LOG(logERROR) << "scalebeta["<<jN<<"]["<<t+1<<"] = " << scalebeta[t+1][jN];
						//FILE_LOG(logERROR) << "A["<<iN<<"]["<<jN<<"] = " << A[iN][jN];
						//FILE_LOG(logERROR) << "densities["<<jN<<"]["<<t+1<<"] = " << densities[jN][t+1];
					}
					//FILE_LOG(logERROR) << "this->scalefactoralpha[t]["<<t<<"] = "<<this->scalefactoralpha[t] << ", densities = "<<densities[iN][t];
					//FILE_LOG(logERROR) << "scalebeta["<<iN<<"]["<<t<<"] = " << scalebeta[t][iN];
					throw nan_detected;
				}
			}
		}

// 	}
// 	else if (this->use_tdens)
// 	{
// 
// 		double beta [this->N];
// 		// Initialization
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			beta[iN] = 1.0;
// 			//FILE_LOG(logDEBUG4) << "beta["<<iN<<"] = " << beta[iN];
// 		}
// 		//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<T-1<<"] = " << scalefactoralpha[T-1];
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			this->scalebeta[T-1][iN] = beta[iN] / this->scalefactoralpha[T-1];
// 			//FILE_LOG(logDEBUG4) << "scalebeta["<<T-1<<"]["<<iN<<"] = " << scalebeta[T-1][iN];
// 		}
// 		// Induction
// 		for (int t=this->T-2; t>=0; t--)
// 		{
// 			for (int iN=0; iN<this->N; iN++)
// 			{
// 				//FILE_LOG(logDEBUG4) << "Calculating backward variable for state " << iN;
// 				beta[iN] = 0.0;
// 				for(int jN=0; jN<this->N; jN++)
// 				{
// 					beta[iN] += this->A[iN][jN] * this->tdensities[t+1][jN] * this->scalebeta[t+1][jN];
// 				}
// 			}
// 			//FILE_LOG(logDEBUG4) << "scalefactoralpha["<<t<<"] = " << scalefactoralpha[t];
// 			for (int iN=0; iN<this->N; iN++)
// 			{
// 				this->scalebeta[t][iN] = beta[iN] / this->scalefactoralpha[t];
// 				//FILE_LOG(logDEBUG4) << "scalebeta["<<t<<"]["<<iN<<"] = " << scalebeta[t][iN];
// 				if (isnan(this->scalebeta[t][iN]))
// 				{
// 					//FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
// 					for (int jN=0; jN<this->N; jN++)
// 					{
// 						//FILE_LOG(logWARNING) << "scalebeta["<<jN<<"]["<<t+1<<"] = " << scalebeta[t+1][jN];
// 						//FILE_LOG(logWARNING) << "A["<<iN<<"]["<<jN<<"] = " << A[iN][jN];
// 						//FILE_LOG(logWARNING) << "tdensities["<<t+1<<"]["<<jN<<"] = " << tdensities[t+1][jN];
// 					}
// 					//FILE_LOG(logWARNING) << "this->scalefactoralpha[t]["<<t<<"] = "<<this->scalefactoralpha[t] << ", tdensities = "<<densities[t][iN];
// 					//FILE_LOG(logWARNING) << "scalebeta["<<iN<<"]["<<t<<"] = " << scalebeta[t][iN];
// 					exit(1);
// 				}
// 			}
// 		}
// 
// 	}

//	dtime = clock() - time;
//	//FILE_LOG(logDEBUG) << "backward(): " << dtime << " clicks";
}

void ScaleHMM::calc_sumgamma()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
//	clock_t time = clock(), dtime;

	// Initialize the sumgamma
	for (int iN=0; iN<this->N; iN++)
	{
		this->sumgamma[iN] = 0.0;
	}

	// Compute the gammas (posteriors) and sumgamma
	#pragma omp parallel for
	for (int iN=0; iN<this->N; iN++)
	{
		for (int t=0; t<this->T; t++)
		{
			this->gamma[iN][t] = this->scalealpha[t][iN] * this->scalebeta[t][iN] * this->scalefactoralpha[t];
			this->sumgamma[iN] += this->gamma[iN][t];
		}
	}
	// Subtract the last value because sumgamma goes only until T-1 and we computed until T to get also loggamma at T
	for (int iN=0; iN<this->N; iN++)
	{
		this->sumgamma[iN] -= this->gamma[iN][T-1];
	}

//	dtime = clock() - time;
//	//FILE_LOG(logDEBUG) << "calc_sumgamma(): " << dtime << " clicks";
}

void ScaleHMM::calc_sumxi()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
//	clock_t time = clock(), dtime;

	double xi;
	// Initialize the sumxi
	for (int iN=0; iN<this->N; iN++)
	{
		for (int jN=0; jN<this->N; jN++)
		{
			this->sumxi[iN][jN] = 0.0;
		}
	}	

	// Compute the sumxi
// 	if (not this->use_tdens)
// 	{

		#pragma omp parallel for
		for (int iN=0; iN<this->N; iN++)
		{
			//FILE_LOG(logDEBUG3) << "Calculating sumxi["<<iN<<"][jN]";
			for (int t=0; t<this->T-1; t++)
			{
				for (int jN=0; jN<this->N; jN++)
				{
					xi = this->scalealpha[t][iN] * this->A[iN][jN] * this->densities[jN][t+1] * this->scalebeta[t+1][jN];
					this->sumxi[iN][jN] += xi;
				}
			}
		}

// 	}
// 	else if (this->use_tdens)
// 	{
// 
// 		#pragma omp parallel for
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			//FILE_LOG(logDEBUG3) << "Calculating sumxi["<<iN<<"][jN]";
// 			for (int t=0; t<this->T-1; t++)
// 			{
// 				for (int jN=0; jN<this->N; jN++)
// 				{
// 					xi = this->scalealpha[t][iN] * this->A[iN][jN] * this->tdensities[t+1][jN] * this->scalebeta[t+1][jN];
// 					this->sumxi[iN][jN] += xi;
// 				}
// 			}
// 		}
// 
// 	}

//	dtime = clock() - time;
//	//FILE_LOG(logDEBUG) << "calc_sumxi(): " << dtime << " clicks";
}

void ScaleHMM::calc_loglikelihood()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
//	clock_t time = clock(), dtime;

	this->logP = 0.0;
	for (int t=0; t<this->T; t++)
	{
		this->logP += log(this->scalefactoralpha[t]);
	}

//	dtime = clock() - time;
//	//FILE_LOG(logDEBUG) << "calc_loglikelihood(): " << dtime << " clicks";
}

void ScaleHMM::calc_densities()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
//	clock_t time = clock(), dtime;
	// Errors thrown inside a #pragma must be handled inside the thread
	std::vector<bool> nan_encountered(this->N);
	#pragma omp parallel for
	for (int iN=0; iN<this->N; iN++)
	{
		//FILE_LOG(logDEBUG3) << "Calculating densities for state " << iN;
		try
		{
			this->densityFunctions[iN]->calc_densities(this->densities[iN]);
		}
		catch(std::exception& e)
		{
			if (strcmp(e.what(),"nan detected")==0) { nan_encountered[iN]=true; }
			else { throw; }
		}
	}
	for (int iN=0; iN<this->N; iN++)
	{
		if (nan_encountered[iN]==true)
		{
			throw nan_detected;
		}
	}

	// Check if the density for all states is numerically zero and correct to prevent NaNs
	std::vector<double> temp(this->N);
	// t=0
	for (int iN=0; iN<this->N; iN++)
	{
		temp[iN] = this->densities[iN][0];
	}
	if (*std::max_element(temp.begin(), temp.end()) == 0.0)
	{
		for (int iN=0; iN<this->N; iN++)
		{
			this->densities[iN][0] = 0.00000000001;
		}
	}
	// t>0
	for (int t=1; t<this->T; t++)
	{
		for (int iN=0; iN<this->N; iN++)
		{
			temp[iN] = this->densities[iN][t];
		}
		if (*std::max_element(temp.begin(), temp.end()) == 0.0)
		{
			for (int iN=0; iN<this->N; iN++)
			{
				this->densities[iN][t] = this->densities[iN][t-1];
			}
		}
	}

//	dtime = clock() - time;
//	//FILE_LOG(logDEBUG) << "calc_densities(): " << dtime << " clicks";
}

void ScaleHMM::print_uni_iteration(int iteration)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
	int bs = 106;
	char buffer [106];
	if (iteration % 20 == 0)
	{
		snprintf(buffer, bs, "%10s%20s%20s%20s%15s", "Iteration", "log(P)", "dlog(P)", "Diff in posterior", "Time in sec");
		//FILE_LOG(logITERATION) << buffer;
		Rprintf("%s\n", buffer);
	}
	if (iteration == 0)
	{
		snprintf(buffer, bs, "%10s%20s%20s%20s%*d", "0", "-inf", "-", "-", 15, this->baumWelchTime_real);
	}
	else if (iteration == 1)
	{
		snprintf(buffer, bs, "%*d%*f%20s%*f%*d", 10, iteration, 20, this->logP, "inf", 20, this->sumdiff_posterior, 15, this->baumWelchTime_real);
	}
	else
	{
		snprintf(buffer, bs, "%*d%*f%*f%*f%*d", 10, iteration, 20, this->logP, 20, this->dlogP, 20, this->sumdiff_posterior, 15, this->baumWelchTime_real);
	}
	//FILE_LOG(logITERATION) << buffer;
	Rprintf("%s\n", buffer);

	// Flush Rprintf statements to R console
	R_FlushConsole();
}

void ScaleHMM::print_multi_iteration(int iteration)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
	int bs = 86;
	char buffer [86];
	if (iteration % 20 == 0)
	{
		snprintf(buffer, bs, "%10s%20s%20s%15s", "Iteration", "log(P)", "dlog(P)", "Time in sec");
		//FILE_LOG(logITERATION) << buffer;
		Rprintf("%s\n", buffer);
	}
	if (iteration == 0)
	{
		snprintf(buffer, bs, "%10s%20s%20s%*d", "0", "-inf", "-", 15, this->baumWelchTime_real);
	}
	else if (iteration == 1)
	{
		snprintf(buffer, bs, "%*d%*f%20s%*d", 10, iteration, 20, this->logP, "inf", 15, this->baumWelchTime_real);
	}
	else
	{
		snprintf(buffer, bs, "%*d%*f%*f%*d", 10, iteration, 20, this->logP, 20, this->dlogP, 15, this->baumWelchTime_real);
	}
	//FILE_LOG(logITERATION) << buffer;
	Rprintf("%s\n", buffer);

	// Flush Rprintf statements to R console
	R_FlushConsole();
}

void ScaleHMM::print_uni_params()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	int bs = 82;
	char buffer [82];
	int cx;
	snprintf(buffer, bs, " -------------------------------------------------------------------------------");
	//FILE_LOG(logINFO) << buffer;
	Rprintf("%s\n", buffer);
	snprintf(buffer, bs, "|%80s", "|");
	//FILE_LOG(logINFO) << buffer;
	Rprintf("%s\n", buffer);
	// print loglik
	snprintf(buffer, bs, "| log(P) = %*.6f%54s", 16, this->logP, "|");
	//FILE_LOG(logINFO) << buffer;
	Rprintf("%s\n", buffer);
	snprintf(buffer, bs, "|%80s", "|");
	//FILE_LOG(logINFO) << buffer;
	Rprintf("%s\n", buffer);
	// print initial probabilities
	cx = snprintf(buffer, bs, "|%7s", "");
	for (int iN=0; iN<this->N; iN++)
	{
		cx += snprintf(buffer+cx, bs-cx, "proba[%d] = %.6f    ", iN, this->proba[iN]);
	}
	cx += snprintf(buffer+cx, bs-cx, "   |");
	//FILE_LOG(logINFO) << buffer;
	Rprintf("%s\n", buffer);
	snprintf(buffer, bs, "|%80s", "|");
	//FILE_LOG(logINFO) << buffer;
	Rprintf("%s\n", buffer);
	// print transition probabilities
	for (int iN=0; iN<this->N; iN++)
	{
		cx = snprintf(buffer, bs, "|%7s", "");
		for (int jN=0; jN<this->N; jN++)
		{
			cx += snprintf(buffer+cx, bs-cx, "A[%d][%d] = %.6f    ", iN, jN, this->A[iN][jN]);
		}
		cx += snprintf(buffer+cx, bs-cx, "      |");
		//FILE_LOG(logINFO) << buffer;
		Rprintf("%s\n", buffer);
	}
	// print emission parameters
	snprintf(buffer, bs, "|%80s", "|");
	//FILE_LOG(logINFO) << buffer;
	Rprintf("%s\n", buffer);
	for (int iN=0; iN<this->N; iN++)
	{
		if (iN == 1)
		{
			snprintf(buffer, bs, "| unmodified component%59s", "|");
			//FILE_LOG(logINFO) << buffer;
			Rprintf("%s\n", buffer);
		}
		if (iN == 2)
		{
			snprintf(buffer, bs, "| modified component%61s", "|");
			//FILE_LOG(logINFO) << buffer;
			Rprintf("%s\n", buffer);
		}
		double curMean = this->densityFunctions[iN]->get_mean();
		double curVar = this->densityFunctions[iN]->get_variance();
		snprintf(buffer, bs, "| mean = %*.2f, var = %*.2f%20s", 6, curMean, 8, curVar, "|");
		//FILE_LOG(logINFO) << buffer;
		Rprintf("%s\n", buffer);
	}
	
	snprintf(buffer, bs, "|%80s", "|");
	//FILE_LOG(logINFO) << buffer;
	Rprintf("%s\n", buffer);
	snprintf(buffer, bs, " -------------------------------------------------------------------------------");
	//FILE_LOG(logINFO) << buffer;
	Rprintf("%s\n", buffer);
	//FILE_LOG(logINFO) << "";
	Rprintf("\n");

	// Flush Rprintf statements to R console
	R_FlushConsole();
}


