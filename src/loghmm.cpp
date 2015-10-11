


#include "loghmm.h"

// ============================================================
// Hidden Markov Model implemented with scaling strategy
// ============================================================

// Public =====================================================

// Constructor and Destructor ---------------------------------
LogHMM::LogHMM(int T, int N)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->T = T;
	this->N = N;
	this->A = CallocDoubleMatrix(N, N);
	this->logA = CallocDoubleMatrix(N, N);
	this->logalpha = CallocDoubleMatrix(T, N);
	this->logbeta = CallocDoubleMatrix(T, N);
	this->logdensities = CallocDoubleMatrix(N, T);
// 	this->tdensities = CallocDoubleMatrix(T, N);
	this->proba = (double*) Calloc(N, double);
	this->logproba = (double*) Calloc(N, double);
	this->gamma = CallocDoubleMatrix(N, T);
	this->sumgamma = (double*) Calloc(N, double);
	this->sumxi = CallocDoubleMatrix(N, N);
	this->logP = -INFINITY;
	this->dlogP = INFINITY;
	this->sumdiff_state_last = 0;
	this->sumdiff_posterior = 0.0;
// 	this->use_tdens = false;

}

LogHMM::~LogHMM()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	FreeDoubleMatrix(this->A, this->N);
	FreeDoubleMatrix(this->logA, this->N);
	FreeDoubleMatrix(this->logalpha, this->T);
	FreeDoubleMatrix(this->logbeta, this->T);
	FreeDoubleMatrix(this->logdensities, this->N);
// 	FreeDoubleMatrix(this->tdensities, this->T);
	FreeDoubleMatrix(this->gamma, this->N);
	FreeDoubleMatrix(this->sumxi, this->N);
	Free(this->proba);
	Free(this->logproba);
	Free(this->sumgamma);
	for (int iN=0; iN<this->N; iN++)
	{
		//FILE_LOG(logDEBUG1) << "Deleting density functions"; 
		delete this->densityFunctions[iN];
	}
}

// Methods ----------------------------------------------------
void LogHMM::initialize_transition_probs(double* initial_A, bool use_initial_params)
{

	if (use_initial_params)
	{
		for (int iN=0; iN<this->N; iN++)
		{
			for (int jN=0; jN<this->N; jN++)
			{
				// convert from vector to matrix representation
				this->A[jN][iN] = initial_A[iN*this->N + jN];
				this->logA[jN][iN] = log(this->A[jN][iN]);
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
				{
					this->A[iN][jN] = self;
					this->logA[iN][jN] = log(this->A[iN][jN]);
				}
				else
				{
					this->A[iN][jN] = other;
					this->logA[iN][jN] = log(this->A[iN][jN]);
				}
				// Save value to initial A
				initial_A[iN*this->N + jN] = this->A[iN][jN];
			}
		}
	}
	
}

void LogHMM::initialize_proba(double* initial_proba, bool use_initial_params)
{

	if (use_initial_params)
	{
		for (int iN=0; iN<this->N; iN++)
		{
			this->proba[iN] = initial_proba[iN];
			this->proba[iN] = log(this->proba[iN]);
		}
	}
	else
	{
		for (int iN=0; iN<this->N; iN++)
		{
			this->proba[iN] = (double)1/this->N;
			this->proba[iN] = log(this->proba[iN]);
			// Save value to initial proba
			initial_proba[iN] = this->proba[iN];
		}
	}

}

void LogHMM::EM(int* maxiter, int* maxtime, double* eps)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;

	double logPold = -INFINITY;
	double logPnew;
	double** gammaold = CallocDoubleMatrix(this->N, this->T);

	// Parallelization settings
// 	omp_set_nested(1);
	
	// measuring the time
	this->EMStartTime_sec = time(NULL);

	// Print some initial information
	//FILE_LOG(logINFO) << "";
	//FILE_LOG(logINFO) << "INITIAL PARAMETERS";
// 	this->print_uni_params();
	this->print_uni_iteration(0);

	R_CheckUserInterrupt();

	// Do the Baum-Welch and updates
	int iteration = 0;
	while (((this->EMTime_real < *maxtime) or (*maxtime < 0)) and ((iteration < *maxiter) or (*maxiter < 0)))
	{

		iteration++;
		
		//FILE_LOG(logDEBUG1) << "Calling calc_densities() from EM()";
		try { this->calc_densities(); } catch(...) { throw; }
		R_CheckUserInterrupt();

		//FILE_LOG(logDEBUG1) << "Calling forward() from EM()";
		try { this->forward(); } catch(...) { throw; }
		R_CheckUserInterrupt();

		//FILE_LOG(logDEBUG1) << "Calling backward() from EM()";
		try { this->backward(); } catch(...) { throw; }
		R_CheckUserInterrupt();

		//FILE_LOG(logDEBUG1) << "Calling calc_loglikelihood() from EM()";
		this->calc_loglikelihood();
		logPnew = this->logP;
		if(isnan(logPnew))
		{
			//FILE_LOG(logERROR) << "logPnew = " << logPnew;
			throw nan_detected;
		}
		this->dlogP = logPnew - logPold;

		//FILE_LOG(logDEBUG1) << "Calling calc_sumxi() from EM()";
		this->calc_sumxi();
		R_CheckUserInterrupt();

		//FILE_LOG(logDEBUG1) << "Calling calc_sumgamma() from EM()";
		this->calc_sumgamma();
		R_CheckUserInterrupt();

// 		clock_t clocktime = clock(), dtime;
		// difference in posterior
		//FILE_LOG(logDEBUG1) << "Calculating differences in posterior in EM()";
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
// 		dtime = clock() - clocktime;
// 		//FILE_LOG(logDEBUG) << "differences in posterior: " << dtime << " clicks";

		R_CheckUserInterrupt();

		// Print information about current iteration
		this->print_uni_iteration(iteration);

		// Check convergence
		if(this->dlogP < *eps) //it has converged
		{
			//FILE_LOG(logINFO) << "\nConvergence reached!\n";
			Rprintf("\nConvergence reached!\n\n");
			this->check_for_state_swap();
			break;
		}
		else
		{ // not converged
			this->EMTime_real = difftime(time(NULL),this->EMStartTime_sec);
			if (iteration == *maxiter)
			{
				//FILE_LOG(logINFO) << "Maximum number of iterations reached!";
				Rprintf("Maximum number of iterations reached!\n");
				this->check_for_state_swap();
			}
			else if ((this->EMTime_real >= *maxtime) and (*maxtime >= 0))
			{
				//FILE_LOG(logINFO) << "Exceeded maximum time!";
				Rprintf("Exceeded maximum time!\n");
				this->check_for_state_swap();
			}
			logPold = logPnew;
		}
		
		// Updating initial probabilities proba and transition matrix A
		for (int iN=0; iN<this->N; iN++)
		{
			this->logproba[iN] = log( this->gamma[iN][0] );
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
					this->logA[iN][jN] = log( this->A[iN][jN] );
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

		// Update distribution of state 1 (null-mixed)
		this->densityFunctions[1]->update(this->gamma[1]);
		// Update distribution of state 2 (monosomic)
		int iN=2;
		this->densityFunctions[iN]->update_constrained(this->gamma, iN, this->N);
		double mean1 = this->densityFunctions[iN]->get_mean();
		double variance1 = this->densityFunctions[iN]->get_variance();
		// Set others as multiples
		for (int iN=0; iN<this->N; iN++)
		{
			if (iN!=1 and iN!=2)
			{
				this->densityFunctions[iN]->set_mean(mean1 * (iN-1));
				this->densityFunctions[iN]->set_variance(variance1 * (iN-1));
			}
			//FILE_LOG(logDEBUG1) << "mean(state="<<iN<<") = " << this->densityFunctions[iN]->get_mean();
		}
		
// 		// Update distribution of state 3 (disomic)
// // 		clock_t clocktime = clock(), dtime;
// 		this->densityFunctions[3]->update(this->gamma[3]);
// 		double mean2 = this->densityFunctions[3]->get_mean();
// 		double variance2 = this->densityFunctions[3]->get_variance();
// 		// Update distribution of state 1 (null-mixed)
// 		this->densityFunctions[1]->update(this->gamma[1]);
// // 		// Set mean of state 1 (null-mixed) to half of monosomic
// // 		this->densityFunctions[1]->set_mean(mean2/2 / 2);
// // 		this->densityFunctions[1]->set_variance(variance2/2 / 2);
// 		// Set the others as multiples of disomic
// 		for (int iN=0; iN<this->N; iN++)
// 		{
// 			if (iN!=1 and iN!=3)
// 			{
// 				this->densityFunctions[iN]->set_mean(mean2/2 * (iN-1));
// 				this->densityFunctions[iN]->set_variance(variance2/2 * (iN-1));
// 			}
// 		}
// // 		dtime = clock() - clocktime;
// // 	 	//FILE_LOG(logDEBUG) << "updating distributions: " << dtime << " clicks";
		R_CheckUserInterrupt();

	} /* main loop end */
    
    
	//Print the last results
	//FILE_LOG(logINFO) << "";
	//FILE_LOG(logINFO) << "FINAL ESTIMATION RESULTS";
// 	this->print_uni_params();

	/* free memory */
	FreeDoubleMatrix(gammaold, this->N);

	// Return values
	*maxiter = iteration;
	*eps = this->dlogP;
	this->EMTime_real = difftime(time(NULL),this->EMStartTime_sec);
	*maxtime = this->EMTime_real;
}

void LogHMM::check_for_state_swap()
{
// TODO
}

void LogHMM::calc_weights(double* weights)
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
void LogHMM::get_posteriors(double** post)
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

double LogHMM::get_posterior(int iN, int t)
{
	//FILE_LOG(logDEBUG4) << __PRETTY_FUNCTION__;
	return(this->gamma[iN][t]);
}

double LogHMM::get_proba(int i)
{
	return( this->proba[i] );
}

double LogHMM::get_A(int i, int j)
{
	return( this->A[i][j] );
}

double LogHMM::get_logP()
{
	return( this->logP );
}

void LogHMM::set_cutoff(int cutoff)
{
	this->cutoff = cutoff;
}

// Private ====================================================
// Methods ----------------------------------------------------
void LogHMM::forward()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// 	clock_t time = clock(), dtime;

	// Initialization
	for (int iN=0; iN<this->N; iN++)
	{
		this->logalpha[0][iN] = this->logproba[iN] + this->logdensities[iN][0];
		//FILE_LOG(logDEBUG4) << "logalpha[0]["<<iN<<"] = " << logalpha[0][iN];
	}
	// Induction
	for (int t=1; t<this->T; t++)
	{
		double temp = Max(&this->logalpha[t-1][0],this->N);
		for (int iN=0; iN<this->N; iN++)
		{
			double helpsum = 0.0;
			for (int jN=0; jN<this->N; jN++)
			{
				helpsum += exp( this->logalpha[t-1][jN]  + this->logA[jN][iN] - temp );
			}
			this->logalpha[t][iN] = temp + log(helpsum) + this->logdensities[iN][t];
			//FILE_LOG(logDEBUG4) << "logalpha["<<t<<"]["<<iN<<"] = " << logalpha[t][iN];
			// Security check for NANs
			if(isnan(this->logalpha[t][iN]))
			{
				//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				for (int jN=0; jN<this->N; jN++)
				{
					//FILE_LOG(logERROR) << "logalpha["<<t-1<<"]["<<jN<<"] = " << logalpha[t-1][jN];
					//FILE_LOG(logERROR) << "A["<<iN<<"]["<<jN<<"] = " << A[iN][jN];
				}
				//FILE_LOG(logERROR) << "temp = "<<temp << ", helpsum = "<<helpsum << ", log(helpsum) = "<<log(helpsum) << ", logdensities = "<<logdensities[iN][t];
				//FILE_LOG(logERROR) << "logalpha["<<t<<"]["<<iN<<"] = " << logalpha[t][iN];
				throw nan_detected;
			}
		}
	}

// 	dtime = clock() - time;
	//FILE_LOG(logDEBUG) << "forward(): " << dtime << " clicks";
}

void LogHMM::backward()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// 	clock_t time = clock(), dtime;

	// Initialization
	for (int iN=0; iN<this->N; iN++)
	{
		this->logbeta[T-1][iN] = 0.0; //=log(1)
		//FILE_LOG(logDEBUG4) << "logbeta["<<T-1<<"]["<<iN<<"] = " << logbeta[T-1][iN];
	}
	// Induction
	for (int t=this->T-2; t>=0; t--)
	{
		for (int iN=0; iN<this->N; iN++)
		{
			std::vector<double> tempvec(this->N);
			for(int jN=0; jN<this->N; jN++)
			{
				tempvec[jN] = this->logA[iN][jN] + this->logdensities[jN][t+1] + this->logbeta[t+1][jN];
			}
			double temp = *std::max_element(tempvec.begin(), tempvec.end());
			double helpsum = 0.0;
			for (int jN=0; jN<this->N; jN++)
			{
				helpsum += exp( tempvec[jN] - temp );
			}
			this->logbeta[t][iN] = log(helpsum) + temp;
			//FILE_LOG(logDEBUG4) << "logbeta["<<t<<"]["<<iN<<"] = " << logbeta[t][iN];
			// Security check for NANs
			if (isnan(this->logbeta[t][iN]))
			{
				//FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				for (int jN=0; jN<this->N; jN++)
				{
					//FILE_LOG(logERROR) << "logbeta["<<jN<<"]["<<t+1<<"] = " << scalebeta[t+1][jN];
					//FILE_LOG(logERROR) << "A["<<iN<<"]["<<jN<<"] = " << A[iN][jN];
					//FILE_LOG(logERROR) << "logdensities["<<jN<<"]["<<t+1<<"] = " << logdensities[jN][t+1];
				}
				//FILE_LOG(logERROR) << "temp = "<<temp << ", helpsum = "<<helpsum << ", log(helpsum) = "<<log(helpsum);
				//FILE_LOG(logERROR) << "logbeta["<<iN<<"]["<<t<<"] = " << logbeta[t][iN];
				throw nan_detected;
			}
		}
	}

// 	dtime = clock() - time;
	//FILE_LOG(logDEBUG) << "backward(): " << dtime << " clicks";
}

void LogHMM::calc_sumgamma()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// 	clock_t time = clock(), dtime;

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
			this->gamma[iN][t] = exp( this->logalpha[t][iN] + this->logbeta[t][iN] - this->logP );
			this->sumgamma[iN] += this->gamma[iN][t];
		}
	}
	// Subtract the last value because sumgamma goes only until T-1 and we computed until T to get also loggamma at T
	for (int iN=0; iN<this->N; iN++)
	{
		this->sumgamma[iN] -= this->gamma[iN][T-1];
	}

// 	dtime = clock() - time;
	//FILE_LOG(logDEBUG) << "calc_sumgamma(): " << dtime << " clicks";
}

void LogHMM::calc_sumxi()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// 	clock_t time = clock(), dtime;

	double logxi;
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
					logxi = this->logalpha[t][iN] + this->logA[iN][jN] + this->logdensities[jN][t+1] + this->logbeta[t+1][jN] - this->logP;
					this->sumxi[iN][jN] += exp( logxi );
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

// 	dtime = clock() - time;
	//FILE_LOG(logDEBUG) << "calc_sumxi(): " << dtime << " clicks";
}

void LogHMM::calc_loglikelihood()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// 	clock_t time = clock(), dtime;

	double helpsum = 0.0;
	std::vector<double> tempvec(this->N);
	
	for (int iN=0; iN<this->N; iN++)
	{
		tempvec[iN] = this->logalpha[this->T-1][iN];
	}
	double temp = *std::max_element(tempvec.begin(), tempvec.end());
	for (int iN=0; iN<this->N; iN++)
	{
		helpsum += exp( this->logalpha[this->T-1][iN] - temp );
	}
	this->logP = temp + log(helpsum);

// 	dtime = clock() - time;
	//FILE_LOG(logDEBUG) << "calc_loglikelihood(): " << dtime << " clicks";
}

void LogHMM::calc_densities()
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// 	clock_t time = clock(), dtime;
	#pragma omp parallel for
	for (int iN=0; iN<this->N; iN++)
	{
		//FILE_LOG(logDEBUG3) << "Calculating densities for state " << iN;
		this->densityFunctions[iN]->calc_logdensities(this->logdensities[iN]);
	}

// 	if (this->use_tdens)
// 	{
// 
// 		for (int t=0; t<this->T; t++)
// 		{
// 			for (int iN=0; iN<this->N; iN++)
// 			{
// 				this->tdensities[t][iN] = this->densities[iN][t];
// 			}
// 		}
// 
// 	}

// 	dtime = clock() - time;
	//FILE_LOG(logDEBUG) << "calc_densities(): " << dtime << " clicks";
}

void LogHMM::print_uni_iteration(int iteration)
{
	//FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->EMTime_real = difftime(time(NULL),this->EMStartTime_sec);
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
		snprintf(buffer, bs, "%10s%20s%20s%20s%*d", "0", "-inf", "-", "-", 15, this->EMTime_real);
	}
	else if (iteration == 1)
	{
		snprintf(buffer, bs, "%*d%*f%20s%*f%*d", 10, iteration, 20, this->logP, "inf", 20, this->sumdiff_posterior, 15, this->EMTime_real);
	}
	else
	{
		snprintf(buffer, bs, "%*d%*f%*f%*f%*d", 10, iteration, 20, this->logP, 20, this->dlogP, 20, this->sumdiff_posterior, 15, this->EMTime_real);
	}
	//FILE_LOG(logITERATION) << buffer;
	Rprintf("%s\n", buffer);

	// Flush Rprintf statements to R console
	R_FlushConsole();
}

void LogHMM::print_uni_params()
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


