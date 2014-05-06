#include "densities.h"

// ============================================================
// Normal (Gaussian) density
// ============================================================

// Constructor and Destructor ---------------------------------
Normal::Normal(int* observations, int T, double mean, double variance)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->obs = observations;
	this->T = T;
	this->mean = mean;
	this->variance = variance;
	this->sd = sqrt(variance);
}

Normal::~Normal()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

// Methods ----------------------------------------------------
void Normal::calc_densities(double* density)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int t=0; t<this->T; t++)
	{
		density[t] = dnorm(this->obs[t], this->mean, this->sd, 0);
	}
}

void Normal::calc_logdensities(double* logdensity)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int t=0; t<this->T; t++)
	{
		logdensity[t] = dnorm(this->obs[t], this->mean, this->sd, 1);
	}
}

void Normal::update(double* weights)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
// TODO
}

// Getter and Setter ------------------------------------------
DensityName Normal::get_name()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(NORMAL);
}

double Normal::get_mean()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->mean);
}

void Normal::set_mean(double mean)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->mean = mean;
}

double Normal::get_variance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->variance);
}

void Normal::set_variance(double variance)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->variance = variance;
	this->sd = sqrt(variance);
}

double Normal::get_stdev()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->sd);
}

void Normal::set_stdev(double stdev)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->sd = stdev;
	this->variance = stdev*stdev;
}


// ============================================================
// Negative Binomial density
// ============================================================

// Constructor and Destructor ---------------------------------
NegativeBinomial::NegativeBinomial(int* observations, int T, double mean, double variance)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->obs = observations;
	this->T = T;
	this->mean = mean;
	this->variance = variance;
	this->size = this->fsize(mean, variance);
	this->prob = this->fprob(mean, variance);
	this->max_obs = intMax(observations, T);
}

NegativeBinomial::~NegativeBinomial()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

// Methods ----------------------------------------------------
void NegativeBinomial::calc_densities(double* density)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int t=0; t<this->T; t++)
	{
		density[t] = dnbinom(this->obs[t], this->size, this->prob, 0);
	}
}

void NegativeBinomial::calc_logdensities(double* logdensity)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int t=0; t<this->T; t++)
	{
		logdensity[t] = dnbinom(this->obs[t], this->size, this->prob, 1);
	}
}

void NegativeBinomial::update(double* weights)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double eps = 1e-4, kmax;
	double numerator, denominator, rhere, dr, Fr, dFrdr, DigammaR, DigammaRplusDR, logp;
	// Update p
	numerator=denominator=0.0;
	clock_t time, dtime;
	time = clock();
	for (int t=0; t<this->T; t++)
	{
		numerator+=weights[t]*this->size;
		denominator+=weights[t]*(this->size+this->obs[t]);
	}
	this->prob = numerator/denominator; // Update of r is now done with updated p
	dtime = clock() - time;
	FILE_LOG(logDEBUG1) << "update prob: "<<dtime<< " clicks";
	// Update of r with Newton Method
	logp = log(this->prob);
	rhere = this->size;
	dr = 0.00001;
	kmax = 20;
	time = clock();
	// Select strategy for computing digammas
	if (this->max_obs <= this->T)
	{
		FILE_LOG(logDEBUG2) << "Precomputing digammas in " << __func__ << " for every obs[t], because max(obs)<=T";
		double DigammaRplusX[this->max_obs+1], DigammaRplusDRplusX[this->max_obs+1];
		for (int k=1; k<kmax; k++)
		{
			Fr=dFrdr=0.0;
			DigammaR = digamma(rhere);
			DigammaRplusDR = digamma(rhere + dr);
			// Precompute the digammas by iterating over all possible values of the observation vector
			for (int j=0; j<=this->max_obs; j++)
			{
				DigammaRplusX[j] = digamma(rhere+j);
				DigammaRplusDRplusX[j] = digamma(rhere+dr+j);
			}
			for(int t=0; t<this->T; t++)
			{
				if(this->obs[t]==0)
				{
					Fr+=weights[t]*logp;
					//dFrdr+=0;
				}
				if(this->obs[t]!=0)
				{
					Fr+=weights[t]*(logp-DigammaR+DigammaRplusX[(int)obs[t]]);
					dFrdr+=weights[t]/dr*(DigammaR-DigammaRplusDR+DigammaRplusDRplusX[(int)obs[t]]-DigammaRplusX[(int)obs[t]]);
				}
			}
			if(fabs(Fr)<eps)
{
				break;
			}
			if(Fr/dFrdr<rhere) rhere=rhere-Fr/dFrdr;
			if(Fr/dFrdr>rhere) rhere=rhere/2.0;
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing digammas in " << __func__ << " for every t, because max(obs)>T";
		double DigammaRplusX, DigammaRplusDRplusX;
		for (int k=1; k<kmax; k++)
		{
			Fr = dFrdr = 0.0;
			DigammaR = digamma(rhere);
			DigammaRplusDR = digamma(rhere + dr);
			for(int t=0; t<this->T; t++)
			{
				DigammaRplusX = digamma(rhere+this->obs[t]);
				DigammaRplusDRplusX = digamma(rhere+dr+this->obs[t]);
				if(this->obs[t]==0)
				{
					Fr+=weights[t]*logp;
				}
				if(this->obs[t]!=0)
				{
					Fr+=weights[t]*(logp-DigammaR+DigammaRplusX);
					dFrdr+=weights[t]/dr*(DigammaR-DigammaRplusDR+DigammaRplusDRplusX-DigammaRplusX);
				}
			}
			if(fabs(Fr)<eps)
			{
				break;
			}
			if(Fr/dFrdr<rhere) rhere=rhere-Fr/dFrdr;
			if(Fr/dFrdr>rhere) rhere=rhere/2.0;
		}
	}
	this->size = rhere;
	FILE_LOG(logDEBUG1) << "size = "<<this->size << ", prob = "<<this->prob;

	// Update mean and variance
	this->mean = this->fmean(this->size, this->prob);
	this->variance = this->fvariance(this->size, this->prob);

	dtime = clock() - time;
	FILE_LOG(logDEBUG1) << "update size: "<<dtime<< " clicks";

}

double NegativeBinomial::fsize(double mean, double variance)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return( mean*mean / (variance - mean) );
}

double NegativeBinomial::fprob(double mean, double variance)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return( mean / variance );
}

double NegativeBinomial::fmean(double size, double prob)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return( size/prob - size );
}

double NegativeBinomial::fvariance(double size, double prob)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return( (size - prob*size) / (prob*prob) );
}

// Getter and Setter ------------------------------------------
DensityName NegativeBinomial::get_name()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(NEGATIVE_BINOMIAL);
}

double NegativeBinomial::get_mean()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->mean);
}

void NegativeBinomial::set_mean(double mean)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->mean = mean;
	this->size = fsize(mean, this->variance);
	this->prob = fprob(mean, this->variance);
}

double NegativeBinomial::get_variance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->variance);
}

void NegativeBinomial::set_variance(double variance)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->variance = variance;
	this->size = fsize(mean, this->variance);
	this->prob = fprob(mean, this->variance);
}

double NegativeBinomial::get_size()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->size);
}

double NegativeBinomial::get_prob()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->prob);
}


// ============================================================
// Zero Inflation density
// ============================================================

// Constructor and Destructor ---------------------------------
ZeroInflation::ZeroInflation(int* observations, int T)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->obs = observations;
	this->T = T;
}

ZeroInflation::~ZeroInflation()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

// Methods ----------------------------------------------------
void ZeroInflation::calc_densities(double* density)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int t=0; t<this->T; t++)
	{
		if(this->obs[t]==0)
		{
			density[t] = 1.0;
		}
		if(this->obs[t]>0)
		{
			density[t] = 0.0;
		}
		FILE_LOG(logDEBUG4) << "density["<<t<<"] = " << density[t];
	}
}

void ZeroInflation::calc_logdensities(double* logdensity)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int t=0; t<this->T; t++)
	{
		if(this->obs[t]==0)
		{
			logdensity[t] = 0.0;
		};
		if(this->obs[t]>0)
		{
			logdensity[t] = -100; // -INFINITY gives nan's somewhere downstream
		}
		FILE_LOG(logDEBUG4) << "logdensity["<<t<<"] = " << logdensity[t];
	}
}

void ZeroInflation::update(double* weights)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

// Getter and Setter ------------------------------------------
DensityName ZeroInflation::get_name()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(ZERO_INFLATION);
}

double ZeroInflation::get_mean()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(0);
}

void ZeroInflation::set_mean(double mean)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

double ZeroInflation::get_variance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(0);
}

void ZeroInflation::set_variance(double variance)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}


