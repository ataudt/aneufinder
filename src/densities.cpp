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
NegativeBinomial::NegativeBinomial(int* observations, int T, double size, double prob)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->obs = observations;
	this->T = T;
	this->size = size;
	this->prob = prob;
	this->lxfactorials = NULL;
	// Precompute the lxfactorials that are used in computing the densities
	if (this->obs != NULL)
	{
		this->max_obs = intMax(observations, T);
		this->lxfactorials = (double*) calloc(max_obs+1, sizeof(double));
		this->lxfactorials[0] = 0.0;	// Not necessary, already 0 because of calloc
		this->lxfactorials[1] = 0.0;
		for (int j=2; j<=max_obs; j++)
		{
			this->lxfactorials[j] = this->lxfactorials[j-1] + log(j);
		}
	}
}

NegativeBinomial::~NegativeBinomial()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	if (this->lxfactorials != NULL)
	{
		free(this->lxfactorials);
	}
}

// Methods ----------------------------------------------------
void NegativeBinomial::calc_logdensities(double* logdens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->T)
	{
		FILE_LOG(logDEBUG2) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		double logdens_per_read [this->max_obs+1];
		for (int j=0; j<=this->max_obs; j++)
		{
			logdens_per_read[j] = lgamma(this->size + j) - lGammaR - lxfactorials[j] + this->size * logp + j * log1minusp;
		}
		for (int t=0; t<this->T; t++)
		{
			logdens[t] = logdens_per_read[(int) this->obs[t]];
			FILE_LOG(logDEBUG4) << "logdens["<<t<<"] = " << logdens[t];
			if (isnan(logdens[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "logdens["<<t<<"] = "<< logdens[t];
				throw nan_detected;
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGammaRplusX = lgamma(this->size + this->obs[t]);
			lxfactorial = this->lxfactorials[(int) this->obs[t]];
			logdens[t] = lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp;
			FILE_LOG(logDEBUG4) << "logdens["<<t<<"] = " << logdens[t];
			if (isnan(logdens[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "logdens["<<t<<"] = "<< logdens[t];
				throw nan_detected;
			}
		}
	}
} 

void NegativeBinomial::calc_densities(double* dens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->T)
	{
		FILE_LOG(logDEBUG2) << "Precomputing gammas in " << __func__ << " for every obs[t], because max(O)<=T";
		double dens_per_read [this->max_obs+1];
		for (int j=0; j<=this->max_obs; j++)
		{
			dens_per_read[j] = exp( lgamma(this->size + j) - lGammaR - lxfactorials[j] + this->size * logp + j * log1minusp );
		}
		for (int t=0; t<this->T; t++)
		{
			dens[t] = dens_per_read[(int) this->obs[t]];
			FILE_LOG(logDEBUG4) << "dens["<<t<<"] = " << dens[t];
			if (isnan(dens[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "dens["<<t<<"] = "<< dens[t];
				throw nan_detected;
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGammaRplusX = lgamma(this->size + this->obs[t]);
			lxfactorial = this->lxfactorials[(int) this->obs[t]];
			dens[t] = exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp );
			FILE_LOG(logDEBUG4) << "dens["<<t<<"] = " << dens[t];
			if (isnan(dens[t]))
			{
				FILE_LOG(logERROR) << __PRETTY_FUNCTION__;
				FILE_LOG(logERROR) << "dens["<<t<<"] = "<< dens[t];
				throw nan_detected;
			}
		}
	}
} 

void NegativeBinomial::update(double* weight)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double eps = 1e-4, kmax;
	double numerator, denominator, rhere, dr, Fr, dFrdr, DigammaR, DigammaRplusDR;
	// Update p
	numerator=denominator=0.0;
// 	clock_t time, dtime;
// 	time = clock();
	for (int t=0; t<this->T; t++)
	{
		numerator+=weight[t]*this->size;
		denominator+=weight[t]*(this->size+this->obs[t]);
	}
	this->prob = numerator/denominator; // Update of r is now done with updated p
	double logp = log(this->prob);
// 	dtime = clock() - time;
// 	FILE_LOG(logDEBUG1) << "updateP(): "<<dtime<< " clicks";
	// Update of r with Newton Method
	rhere = this->size;
	dr = 0.00001;
	kmax = 20;
// 	time = clock();
	// Select strategy for computing digammas
	if (this->max_obs <= this->T)
	{
		FILE_LOG(logDEBUG2) << "Precomputing digammas in " << __func__ << " for every obs[t], because max(O)<=T";
		double DigammaRplusX[this->max_obs+1], DigammaRplusDRplusX[this->max_obs+1];
		for (int k=1; k<kmax; k++)
		{
			Fr=dFrdr=0.0;
			DigammaR = digamma(rhere); // boost::math::digamma<>(rhere);
			DigammaRplusDR = digamma(rhere + dr); // boost::math::digamma<>(rhere+dr);
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
					Fr+=weight[t]*logp;
					//dFrdr+=0;
				}
				if(this->obs[t]!=0)
				{
					Fr+=weight[t]*(logp-DigammaR+DigammaRplusX[(int)obs[t]]);
					dFrdr+=weight[t]/dr*(DigammaR-DigammaRplusDR+DigammaRplusDRplusX[(int)obs[t]]-DigammaRplusX[(int)obs[t]]);
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
		FILE_LOG(logDEBUG2) << "Computing digammas in " << __func__ << " for every t, because max(O)>T";
		double DigammaRplusX, DigammaRplusDRplusX;
		for (int k=1; k<kmax; k++)
		{
			Fr = dFrdr = 0.0;
			DigammaR = digamma(rhere); // boost::math::digamma<>(rhere);
			DigammaRplusDR = digamma(rhere + dr); // boost::math::digamma<>(rhere+dr);
			for(int t=0; t<this->T; t++)
			{
				DigammaRplusX = digamma(rhere+this->obs[t]); //boost::math::digamma<>(rhere+this->obs[ti]);
				DigammaRplusDRplusX = digamma(rhere+dr+this->obs[t]); // boost::math::digamma<>(rhere+dr+this->obs[ti]);
				if(this->obs[t]==0)
				{
					Fr+=weight[t]*logp;
					//dFrdr+=0;
				}
				if(this->obs[t]!=0)
				{
					Fr+=weight[t]*(logp-DigammaR+DigammaRplusX);
					dFrdr+=weight[t]/dr*(DigammaR-DigammaRplusDR+DigammaRplusDRplusX-DigammaRplusX);
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
	FILE_LOG(logDEBUG1) << "r = "<<this->size << ", p = "<<this->prob;

// 	dtime = clock() - time;
// 	FILE_LOG(logDEBUG1) << "updateR(): "<<dtime<< " clicks";

}

// Getter and Setter ------------------------------------------
double NegativeBinomial::get_mean()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return this->size*(1-this->prob)/this->prob;
}

double NegativeBinomial::get_variance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return this->size*(1-this->prob)/this->prob/this->prob;
}

DensityName NegativeBinomial::get_name()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(NEGATIVE_BINOMIAL);
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


