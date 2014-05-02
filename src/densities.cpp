#include "densities.h"
#include "utility.h"

// ------------------------------------------------------------
// Zero inflated negative binomial
// ------------------------------------------------------------

ZiNB::ZiNB()
{
	this->lxfactorials = NULL;
}

ZiNB::ZiNB(int* observations, int T, double r, double p, double w)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->O = observations;
	this->T = T;
	this->p = p;
	this->r = r;
	this->w = w;
	this->lxfactorials = NULL;
	if (this->O != NULL)
	{
		this->max_O = intMax(observations, T);
		this->lxfactorials = (double*) calloc(max_O+1, sizeof(double));
		this->lxfactorials[0] = 0.0;	// Not necessary, already 0 because of calloc
		this->lxfactorials[1] = 0.0;
		for (int j=2; j<=max_O; j++)
		{
			this->lxfactorials[j] = this->lxfactorials[j-1] + log(j);
		}
	}
}

ZiNB::~ZiNB()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	if (this->O != NULL)
	{
		free(this->lxfactorials);
	}
}

void ZiNB::calc_logdensities(double* logdens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->p);
	double log1minusp = log(1-this->p);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->r);
	// Select strategy for computing gammas
	if (this->max_O <= this->T)
	{
		FILE_LOG(logDEBUG2) << "Precomputing gammas in " << __func__ << " for every O[t], because max(O)<=T";
		double lGammaRplusX[this->max_O+1];
		for (int j=0; j<=this->max_O; j++)
		{
			lGammaRplusX[j] = lgamma(this->r + j);
		}
		for (int t=0; t<this->T; t++)
		{
			lxfactorial = this->lxfactorials[(int) this->O[t]];
			if (O[t] == 0)
			{
				logdens[t] = log( this->w + (1-this->w) * exp( lGammaRplusX[(int) this->O[t]] - lGammaR - lxfactorial + this->r * logp + this->O[t] * log1minusp ) );
			}
			else
			{
				logdens[t] = log(1-this->w) + lGammaRplusX[(int) this->O[t]] - lGammaR - lxfactorial + this->r * logp + this->O[t] * log1minusp;
			}
			if (isnan(logdens[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "logdens["<<t<<"] = "<< logdens[t];
				exit(1);
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGammaRplusX = lgamma(this->r + this->O[t]);
			lxfactorial = this->lxfactorials[(int) this->O[t]];
			if (O[t] == 0)
			{
				logdens[t] = log( this->w + (1-this->w) * exp( lGammaRplusX - lGammaR - lxfactorial + this->r * logp + this->O[t] * log1minusp ) );
			}
			else
			{
				logdens[t] = log(1-this->w) + lGammaRplusX - lGammaR - lxfactorial + this->r * logp + this->O[t] * log1minusp;
			}
			if (isnan(logdens[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "logdens["<<t<<"] = "<< logdens[t];
				exit(1);
			}
		}
	}
}

void ZiNB::calc_densities(double* dens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->p);
	double log1minusp = log(1-this->p);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->r);
	// Select strategy for computing gammas
	if (this->max_O <= this->T)
	{
		FILE_LOG(logDEBUG2) << "Precomputing gammas in " << __func__ << " for every O[t], because max(O)<=T";
		double lGammaRplusX[this->max_O+1];
		for (int j=0; j<=this->max_O; j++)
		{
			lGammaRplusX[j] = lgamma(this->r + j);
		}
		for (int t=0; t<this->T; t++)
		{
			lxfactorial = this->lxfactorials[(int) this->O[t]];
			if (O[t] == 0)
			{
				dens[t] = this->w + (1-this->w) * exp( lGammaRplusX[(int) this->O[t]] - lGammaR - lxfactorial + this->r * logp + this->O[t] * log1minusp );
			}
			else
			{
				dens[t] = (1-this->w) * exp( lGammaRplusX[(int) this->O[t]] - lGammaR - lxfactorial + this->r * logp + this->O[t] * log1minusp );
			}
			if (isnan(dens[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "dens["<<t<<"] = "<< dens[t];
				exit(1);
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGammaRplusX = lgamma(this->r + this->O[t]);
			lxfactorial = this->lxfactorials[(int) this->O[t]];
			if (O[t] == 0)
			{
				dens[t] = this->w + (1-this->w) * exp( lGammaRplusX - lGammaR - lxfactorial + this->r * logp + this->O[t] * log1minusp );
			}
			else
			{
				dens[t] = (1-this->w) * exp( lGammaRplusX - lGammaR - lxfactorial + this->r * logp + this->O[t] * log1minusp );
			}
			if (isnan(dens[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "dens["<<t<<"] = "<< dens[t];
				exit(1);
			}
		}
	}
}

void ZiNB::calc_CDFs(double* CDF)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double log1minusp = log(1-this->p);
	double lGammaR = lgamma(this->r);
	double lppowerr = this->r * log(this->p);
	// No selection strategy here, because we must precompute CDF to deal with 1s by shifting them
// 	// Select strategy for computing gammas
// 	if (this->max_O <= this->T)
//	{
		FILE_LOG(logDEBUG2) << "Precomputing gammas in " << __func__ << " for every O[t], because max(O)<=T";
		double lGamma1plusRplusX[this->max_O+1], lGamma2plusX[this->max_O+1], lHyper[this->max_O+1], lppowert[this->max_O+1];
		double precomputed_CDF[this->max_O+1];
		for (int j=0; j<=this->max_O; j++)
		{
			lGamma1plusRplusX[j] = lgamma(1 + this->r + j);
			lGamma2plusX[j] = lgamma(2 + j);
			lHyper[j] = log(gsl_sf_hyperg_2F1(1, 1 + this->r + j, 2 + j, 1-this->p));
			lppowert[j] = (1+j) * log1minusp;
			precomputed_CDF[j] = 1 - exp( log(1-this->w) + lppowerr + lppowert[j] + lGamma1plusRplusX[j] + lHyper[j] - lGammaR - lGamma2plusX[j] );
			if (precomputed_CDF[j] == 1)
			{
				FILE_LOG(logDEBUG4) << "CDF = 1 for O[t] = "<<j<< ", shifting to value of O[t] = "<<j-1;
				precomputed_CDF[j] = precomputed_CDF[j-1]; 
			}
		}
		for (int t=0; t<this->T; t++)
		{
			CDF[t] = precomputed_CDF[(int)O[t]];
			if (isnan(CDF[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "CDF["<<t<<"] = "<< CDF[t];
				exit(1);
			}
		}
// 	}
// 	else
// 	{
// 		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
// 		double lGamma1plusRplusX, lGamma2plusX, lHyper, lppowert;
// 		for (int t=0; t<this->T; t++)
//		{
// 			lGamma1plusRplusX = lgamma(1 + this->r + this->O[t]);
// 			lGamma2plusX = lgamma(2 + this->O[t]);
// 			lHyper = log(gsl_sf_hyperg_2F1(1, 1 + this->r + this->O[t], 2 + this->O[t], 1-this->p));
// 			lppowert = (1+this->O[t]) * log1minusp;
// 			CDF[t] = 1 - exp( log(1-this->w) + lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX ); //TODO: Check formula for log
// 			if(CDF[t] == 0)
// 			{
// 				FILE_LOG(logWARNING) << "CDF["<<t<<"] = "<< CDF[t]; //TODO: Check if this works
// // 				cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);
// 			}
// 			if (isnan(CDF[t]))
// 			{
// 				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
// 				FILE_LOG(logWARNING) << "CDF["<<t<<"] = "<< CDF[t];
// 				exit(1);
// 			}
// 		}
// 	}
}

void ZiNB::calc_logCDFs(double* logCDF)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double log1minusp = log(1-this->p);
	double lGamma1plusRplusX, lHyper, lGammaR, lGamma2plusX, lppowert, lppowerr;
	lGammaR = lgamma(this->r);
	lppowerr = this->r * log(this->p);
	// Select strategy for computing gammas
	if (this->max_O <= this->T)
	{
		FILE_LOG(logDEBUG2) << "Precomputing gammas in " << __func__ << " for every O[t], because max(O)<=T";
		double lGamma1plusRplusX[this->max_O+1], lGamma2plusX[this->max_O+1], lHyper[this->max_O+1], lppowert[this->max_O+1];
		for (int j=0; j<=this->max_O; j++)
		{
			lGamma1plusRplusX[j] = lgamma(1 + this->r + j);
			lGamma2plusX[j] = lgamma(2 + j);
			lHyper[j] = log(gsl_sf_hyperg_2F1(1, 1 + this->r + j, 2 + j, 1-this->p));
			lppowert[j] = (1+j) * log1minusp;
		}
		for (int t=0; t<this->T; t++)
		{
			logCDF[t] = log(1 - exp( log(1-this->w) + lppowerr + lppowert[(int)O[t]] + lGamma1plusRplusX[(int)O[t]] + lHyper[(int)O[t]] - lGammaR - lGamma2plusX[(int)O[t]] ));
			if(logCDF[t] == 0)
			{
				FILE_LOG(logWARNING) << "logCDF["<<t<<"] = 0";
// 				cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);
			}
			if (isnan(logCDF[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "logCDF["<<t<<"] = "<< logCDF[t];
				exit(1);
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGamma1plusRplusX = lgamma(1 + this->r + this->O[t]);
			lGamma2plusX = lgamma(2 + this->O[t]);
			lHyper = log(gsl_sf_hyperg_2F1(1, 1 + this->r + this->O[t], 2 + this->O[t], 1-this->p));
			lppowert = (1+this->O[t]) * log1minusp;
			logCDF[t] = log(1 - exp( log(1-this->w) + lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX )); //TODO: Check formula for log
			if(logCDF[t] == 0)
			{
				FILE_LOG(logWARNING) << "logCDF["<<t<<"] = 0"; //TODO: Check if this works
// 				cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);
			}
			if (isnan(logCDF[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "logCDF["<<t<<"] = "<< logCDF[t];
				exit(1);
			}
		}
	}
}

double ZiNB::getMean()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return (1-this->w)*this->r*(1-this->p)/this->p;
}

double ZiNB::getVariance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return (1-this->w)*this->r*(1-this->p)/this->p/this->p; //TODO: Is this correct?
}

DensityName ZiNB::getType()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(ZINB);
}

void ZiNB::setR (double newR)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
  this->r = newR;
}

double ZiNB::getR()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
  return(this->r);
}

void ZiNB::setP (double newP)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
  this->p = newP;
}

double ZiNB::getP()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
  return(this->p);
}

void ZiNB::setW (double newW)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
    this->w = newW;
}

double ZiNB::getW()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->w);
}

void ZiNB::copy(Density* other)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	ZiNB* o = (ZiNB*)other;
	this->p = o->p;
	this->r = o->r;
	this->O = o->O;
	this->w = o->w;
}

double ZiNB::getLogDensityAt10Variance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->p);
	double log1minusp = log(1-this->p);
	double lGammaR,lGammaRplusX,lxfactorial;
	double logdens;
	// Calculate variance
	double mean = 0, variance = 0;
	for(int t=0; t<this->T; t++)
	{
		mean += O[t];
	}
	mean = mean / this->T;
	for(int t=0; t<this->T; t++)
	{
		variance += pow(O[t] - mean, 2);
	}
	variance = variance / this->T;
	// Calculate logdensity
	int x = 10*((int)variance+1);
	lGammaR=lgamma(this->r);
	lGammaRplusX = lgamma(this->r + x);
	lxfactorial = this->lxfactorials[x];
	if (x == 0)
	{
		logdens = log( this->w + (1-this->w) * exp( lGammaRplusX - lGammaR - lxfactorial + this->r * logp + x * log1minusp ) );
	}
	else
	{
		logdens = log(1-this->w) + lGammaRplusX - lGammaR - lxfactorial + this->r * logp + x * log1minusp;
	}
	if (isnan(logdens))
	{
		FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
		FILE_LOG(logWARNING) << "logdens = "<< logdens;
		exit(1);
	}
	
	return(logdens);
}


// ------------------------------------------------------------
// Negative binomial
// ------------------------------------------------------------

NegativeBinomial::NegativeBinomial()
{
	this->lxfactorials = NULL;
}

NegativeBinomial::NegativeBinomial(int* observations, int T, double r, double p)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->O = observations;
	this->T = T;
	this->r = r;
	this->p = p;
	this->lxfactorials = NULL;
	// Precompute the lxfactorials that are used in computing the densities
	if (this->O != NULL)
	{
		this->max_O = intMax(observations, T);
		this->lxfactorials = (double*) calloc(max_O+1, sizeof(double));
		this->lxfactorials[0] = 0.0;	// Not necessary, already 0 because of calloc
		this->lxfactorials[1] = 0.0;
		for (int j=2; j<=max_O; j++)
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

void NegativeBinomial::calc_logdensities(double* logdens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->p);
	double log1minusp = log(1-this->p);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->r);
	// Select strategy for computing gammas
	if (this->max_O <= this->T)
	{
		FILE_LOG(logDEBUG2) << "Precomputing gammas in " << __func__ << " for every O[t], because max(O)<=T";
		double logdens_per_read [this->max_O+1];
		for (int j=0; j<=this->max_O; j++)
		{
			logdens_per_read[j] = lgamma(this->r + j) - lGammaR - lxfactorials[j] + this->r * logp + j * log1minusp;
		}
		for (int t=0; t<this->T; t++)
		{
			logdens[t] = logdens_per_read[(int) this->O[t]];
			FILE_LOG(logDEBUG4) << "logdens["<<t<<"] = " << logdens[t];
			if (isnan(logdens[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "logdens["<<t<<"] = "<< logdens[t];
				exit(1);
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGammaRplusX = lgamma(this->r + this->O[t]);
			lxfactorial = this->lxfactorials[(int) this->O[t]];
			logdens[t] = lGammaRplusX - lGammaR - lxfactorial + this->r * logp + this->O[t] * log1minusp;
			FILE_LOG(logDEBUG4) << "logdens["<<t<<"] = " << logdens[t];
			if (isnan(logdens[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "logdens["<<t<<"] = "<< logdens[t];
				exit(1);
			}
		}
	}
} 

void NegativeBinomial::calc_densities(double* dens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->p);
	double log1minusp = log(1-this->p);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->r);
	// Select strategy for computing gammas
	if (this->max_O <= this->T)
	{
		FILE_LOG(logDEBUG2) << "Precomputing gammas in " << __func__ << " for every O[t], because max(O)<=T";
		double dens_per_read [this->max_O+1];
		for (int j=0; j<=this->max_O; j++)
		{
			dens_per_read[j] = exp( lgamma(this->r + j) - lGammaR - lxfactorials[j] + this->r * logp + j * log1minusp );
		}
		for (int t=0; t<this->T; t++)
		{
			dens[t] = dens_per_read[(int) this->O[t]];
			FILE_LOG(logDEBUG4) << "dens["<<t<<"] = " << dens[t];
			if (isnan(dens[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "dens["<<t<<"] = "<< dens[t];
				exit(1);
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
					lGammaRplusX = lgamma(this->r + this->O[t]);
			lxfactorial = this->lxfactorials[(int) this->O[t]];
			dens[t] = exp( lGammaRplusX - lGammaR - lxfactorial + this->r * logp + this->O[t] * log1minusp );
			FILE_LOG(logDEBUG4) << "dens["<<t<<"] = " << dens[t];
			if (isnan(dens[t]))
			{
							FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "dens["<<t<<"] = "<< dens[t];
				exit(1);
			}
		}
	}
} 

void NegativeBinomial::calc_CDFs(double* CDF)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double log1minusp = log(1-this->p);
	double lGammaR = lgamma(this->r);
	double lppowerr = this->r * log(this->p);
	// No selection strategy here, because we must precompute CDF to deal with 1s by shifting them
// 	// Select strategy for computing gammas
// 	if (this->max_O <= this->T)
// 	{
		FILE_LOG(logDEBUG2) << "Precomputing gammas in " << __func__ << " for every O[t], because max(O)<=T";
		double lGamma1plusRplusX[this->max_O+1], lGamma2plusX[this->max_O+1], lHyper[this->max_O+1], lppowert[this->max_O+1];
		double precomputed_CDF[this->max_O+1];
		for (int j=0; j<=this->max_O; j++)
		{
			lGamma1plusRplusX[j] = lgamma(1 + this->r + j);
			lGamma2plusX[j] = lgamma(2 + j);
			lHyper[j] = log(gsl_sf_hyperg_2F1(1, 1 + this->r + j, 2 + j, 1-this->p));
			lppowert[j] = (1+j) * log1minusp;
			precomputed_CDF[j] = 1 - exp( lppowerr + lppowert[j] + lGamma1plusRplusX[j] + lHyper[j] - lGammaR - lGamma2plusX[j] );
			if (precomputed_CDF[j] == 1)
			{
				FILE_LOG(logDEBUG4) << "CDF = 1 for O[t] = "<<j<< ", shifting to value of O[t] = "<<j-1;
				precomputed_CDF[j] = precomputed_CDF[j-1]; 
			}
		}
		for (int t=0; t<this->T; t++)
		{
			CDF[t] = precomputed_CDF[(int)O[t]];
			if (isnan(CDF[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "CDF["<<t<<"] = "<< CDF[t];
				exit(1);
			}
		}
// 	}
// 	else
// 	{
// 		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
// 		double lGamma1plusRplusX, lGamma2plusX, lHyper, lppowert;
// 		for (int t=0; t<this->T; t++)
// 		{
// 			lGamma1plusRplusX = lgamma(1 + this->r + this->O[t]);
// 			lGamma2plusX = lgamma(2 + this->O[t]);
// 			lHyper = log(gsl_sf_hyperg_2F1(1, 1 + this->r + this->O[t], 2 + this->O[t], 1-this->p));
// 			lppowert = (1+this->O[t]) * log1minusp;
// 			CDF[t] = 1 - exp( lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX ); //TODO: Check formula for log
// 			if(CDF[t] == 1)
// 			{
// 				FILE_LOG(logWARNING) << "CDF["<<t<<"] = "<< CDF[t]; //TODO: Check if this works
// // 				cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);
// 			}
// 			if (isnan(CDF[t]))
// 			{
// 				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
// 				FILE_LOG(logWARNING) << "CDF["<<t<<"] = "<< CDF[t];
// 				exit(1);
// 			}
// 		}
// 	}
}

void NegativeBinomial::calc_logCDFs(double* logCDF)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double log1minusp = log(1-this->p);
	double lGamma1plusRplusX, lHyper, lGammaR, lGamma2plusX, lppowert, lppowerr;
	lGammaR = lgamma(this->r);
	lppowerr = this->r * log(this->p);
	// Select strategy for computing gammas
	if (this->max_O <= this->T)
	{
		FILE_LOG(logDEBUG2) << "Precomputing gammas in " << __func__ << " for every O[t], because max(O)<=T";
		double lGamma1plusRplusX[this->max_O+1], lGamma2plusX[this->max_O+1], lHyper[this->max_O+1], lppowert[this->max_O+1];
		for (int j=0; j<=this->max_O; j++)
		{
			lGamma1plusRplusX[j] = lgamma(1 + this->r + j);
			lGamma2plusX[j] = lgamma(2 + j);
			lHyper[j] = log(gsl_sf_hyperg_2F1(1, 1 + this->r + j, 2 + j, 1-this->p));
			lppowert[j] = (1+j) * log1minusp;
		}
		for (int t=0; t<this->T; t++)
		{
			logCDF[t] = log(1 - exp( lppowerr + lppowert[(int)O[t]] + lGamma1plusRplusX[(int)O[t]] + lHyper[(int)O[t]] - lGammaR - lGamma2plusX[(int)O[t]] ));
			if(logCDF[t] == 0)
			{
				FILE_LOG(logWARNING) << "logCDF["<<t<<"] = 0";
// 				cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);
			}
			if (isnan(logCDF[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "logCDF["<<t<<"] = "<< logCDF[t];
				exit(1);
			}
		}
	}
	else
	{
		FILE_LOG(logDEBUG2) << "Computing gammas in " << __func__ << " for every t, because max(O)>T";
		for (int t=0; t<this->T; t++)
		{
			lGamma1plusRplusX = lgamma(1 + this->r + this->O[t]);
			lGamma2plusX = lgamma(2 + this->O[t]);
			lHyper = log(gsl_sf_hyperg_2F1(1, 1 + this->r + this->O[t], 2 + this->O[t], 1-this->p));
			lppowert = (1+this->O[t]) * log1minusp;
			logCDF[t] = log(1 - exp( lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX )); //TODO: Check formula for log
			if(logCDF[t] == 0)
			{
				FILE_LOG(logWARNING) << "logCDF["<<t<<"] = 1"; //TODO: Check if this works
// 				cout<<"CAUTION!!!! current ="<<current<<endl; current = this->cdf->at(i-1);
			}
			if (isnan(logCDF[t]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "logCDF["<<t<<"] = "<< logCDF[t];
				exit(1);
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
	clock_t time, dtime;
	time = clock();
	for (int t=0; t<this->T; t++)
	{
		numerator+=weight[t]*this->r;
		denominator+=weight[t]*(this->r+this->O[t]);
	}
	this->p = numerator/denominator; // Update of r is now done with updated p
	double logp = log(this->p);
	dtime = clock() - time;
	FILE_LOG(logDEBUG1) << "updateP(): "<<dtime<< " clicks";
	// Update of r with Newton Method
	rhere = this->r;
	dr = 0.00001;
	kmax = 20;
	time = clock();
	// Select strategy for computing digammas
	if (this->max_O <= this->T)
	{
		FILE_LOG(logDEBUG2) << "Precomputing digammas in " << __func__ << " for every O[t], because max(O)<=T";
		double DigammaRplusX[this->max_O+1], DigammaRplusDRplusX[this->max_O+1];
		for (int k=1; k<kmax; k++)
		{
			Fr=dFrdr=0.0;
			DigammaR = digamma(rhere); // boost::math::digamma<>(rhere);
			DigammaRplusDR = digamma(rhere + dr); // boost::math::digamma<>(rhere+dr);
			// Precompute the digammas by iterating over all possible values of the observation vector
			for (int j=0; j<=this->max_O; j++)
			{
				DigammaRplusX[j] = digamma(rhere+j);
				DigammaRplusDRplusX[j] = digamma(rhere+dr+j);
			}
			for(int t=0; t<this->T; t++)
			{
				if(this->O[t]==0)
				{
					Fr+=weight[t]*logp;
					//dFrdr+=0;
				}
				if(this->O[t]!=0)
				{
					Fr+=weight[t]*(logp-DigammaR+DigammaRplusX[(int)O[t]]);
					dFrdr+=weight[t]/dr*(DigammaR-DigammaRplusDR+DigammaRplusDRplusX[(int)O[t]]-DigammaRplusX[(int)O[t]]);
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
				DigammaRplusX = digamma(rhere+this->O[t]); //boost::math::digamma<>(rhere+this->O[ti]);
				DigammaRplusDRplusX = digamma(rhere+dr+this->O[t]); // boost::math::digamma<>(rhere+dr+this->O[ti]);
				if(this->O[t]==0)
				{
					Fr+=weight[t]*logp;
					//dFrdr+=0;
				}
				if(this->O[t]!=0)
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
	this->r = rhere;
	FILE_LOG(logDEBUG1) << "r = "<<this->r << ", p = "<<this->p;

	dtime = clock() - time;
	FILE_LOG(logDEBUG1) << "updateR(): "<<dtime<< " clicks";

}

double NegativeBinomial::getMean()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return this->r*(1-this->p)/this->p;
}

double NegativeBinomial::getVariance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return this->r*(1-this->p)/this->p/this->p;
}

DensityName NegativeBinomial::getType()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(NB);
}

void NegativeBinomial::setR (double newR)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->r = newR;
}

double NegativeBinomial::getR()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->r);
}

void NegativeBinomial::setP (double newP)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->p = newP;
}

double NegativeBinomial::getP()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(this->p);
}

void NegativeBinomial::copy(Density* other)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	NegativeBinomial* o = (NegativeBinomial*)other;
	this->p = o->p;
	this->r = o->r;
	this->O = o->O;
}

double NegativeBinomial::getLogDensityAt10Variance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logp = log(this->p);
	double log1minusp = log(1-this->p);
	double lGammaR,lGammaRplusX,lxfactorial;
	double logdens;
	// Calculate variance
	double mean = 0, variance = 0;
	for(int t=0; t<this->T; t++)
	{
		mean += O[t];
	}
	mean = mean / this->T;
	for(int t=0; t<this->T; t++)
	{
		variance += pow(O[t] - mean, 2);
	}
	variance = variance / this->T;
	// Calculate logdensity
	int x = 10*((int)variance+1);
	lGammaR=lgamma(this->r);
	lGammaRplusX = lgamma(this->r + x);
	lxfactorial = this->lxfactorials[x];
	logdens = lGammaRplusX - lGammaR - lxfactorial + this->r * logp + x * log1minusp;
	if (isnan(logdens))
	{
		FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
		FILE_LOG(logWARNING) << "logdens = "<< logdens;
		exit(1);
	}
	
	return(logdens);
}


// ------------------------------------------------------------
// Only zeros
// ------------------------------------------------------------

OnlyZeros::OnlyZeros() {}

OnlyZeros::OnlyZeros(int* observations, int T)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->O = observations;
	this->T = T;
}

OnlyZeros::~OnlyZeros()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

void OnlyZeros::calc_logdensities(double* logdens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int t=0; t<this->T; t++)
	{
		if(O[t]==0)
		{
			logdens[t] = 0.0;
		};
		if(O[t]>0)
		{
			logdens[t] = -100; // -INFINITY gives nan's somewhere downstream
		}
		FILE_LOG(logDEBUG4) << "logdens["<<t<<"] = " << logdens[t];
	}
}

void OnlyZeros::calc_densities(double* dens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int t=0; t<this->T; t++)
	{
		if(O[t]==0)
		{
			dens[t] = 1.0;
		}
		if(O[t]>0)
		{
			dens[t] = 0.0;
		}
		FILE_LOG(logDEBUG4) << "dens["<<t<<"] = " << dens[t];
	}
}

void OnlyZeros::copy(Density* other) {}

double OnlyZeros::getMean()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return 0;
}

double OnlyZeros::getVariance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return 0;
}

DensityName OnlyZeros::getType()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(ZI);
}

void OnlyZeros::update(double* weight)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

double OnlyZeros::getLogDensityAt10Variance()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double logdens;
	// Calculate variance
	double mean = 0, variance = 0;
	for(int t=0; t<this->T; t++)
	{
		mean += O[t];
	}
	mean = mean / this->T;
	for(int t=0; t<this->T; t++)
	{
		variance += pow(O[t] - mean, 2);
	}
	variance = variance / this->T;
	// Calculate logdensity
	int x = 10*((int)variance+1);
	if (x == 0)
	{
		logdens = 0;
	}
	else
	{
		logdens = -100;
	}
	
	return(logdens);
}

	
// ------------------------------------------------------------
// Multivariate Copula Approximation
// ------------------------------------------------------------

MVCopulaApproximation::MVCopulaApproximation(int** multiobservations, int T, vector<Density*> marginals, double* cor_matrix_inv, double cor_matrix_determinant)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->multiO = multiobservations;
	this->T = T;
	// these are the marginal distributions (we need their CDF function)
	this->marginals = marginals;
	this->Nmod = this->marginals.size();
	this->cor_matrix_inv = cor_matrix_inv;
	this->cor_matrix_determinant = cor_matrix_determinant;
}

MVCopulaApproximation::~MVCopulaApproximation()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	for (int imod; imod<this->Nmod; imod++)
	{
		delete this->marginals[imod];
	}
}

void MVCopulaApproximation::calc_logdensities(double* logdens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	// Calculate logdensities for marginals
	double** marginals_logdensities = allocDoubleMatrix(this->Nmod, this->T);
	double** marginals_CDFs = allocDoubleMatrix(this->Nmod, this->T);
	for (int imod=0; imod<this->Nmod; imod++)
	{
		FILE_LOG(logDEBUG2) << __func__ << ": calculating marginals for imod = " << imod;
		this->marginals[imod]->calc_logdensities(marginals_logdensities[imod]);
		this->marginals[imod]->calc_CDFs(marginals_CDFs[imod]);
	}
	// Calculate multivariate Copula approximation
	FILE_LOG(logDEBUG2) << __func__ << ": calculate Copula approximation";
	double sum, uniform, exponent, exponentTemp;
	double* z = (double*) calloc(this->Nmod, sizeof(double));
	for (int t=0; t<this->T; t++)
	{
		sum = 0.0;
		for (int imod=0; imod<this->Nmod; imod++)
		{
			sum += marginals_logdensities[imod][t];
			uniform = marginals_CDFs[imod][t];
			z[imod] = qnorm(uniform, 0, 1, 1, 0);
			if (isnan(z[imod]))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "uniform = "<< uniform;
				FILE_LOG(logWARNING) << "z[imod] = "<< z[imod];
			}
		}
		exponent = 0.0;
		for (int imod=0; imod<this->Nmod; imod++)
		{
			exponentTemp = 0.0;
			for(int jmod=0; jmod<Nmod; jmod++)
			{
				if(imod==jmod)
				{
					exponentTemp += z[jmod] * (this->cor_matrix_inv[imod * Nmod + jmod] - 1);
				}
				else
				{
					exponentTemp += z[jmod] * this->cor_matrix_inv[imod * Nmod + jmod];
				}
				if (isnan(exponentTemp))
				{
					FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
					FILE_LOG(logWARNING) << "exponentTemp = "<< exponentTemp;
					FILE_LOG(logWARNING) << "cor_matrix_inv = "<< cor_matrix_inv[imod * Nmod + jmod];
					FILE_LOG(logWARNING) << "z["<<jmod<<"] = "<< z[jmod];
					exit(1);
				}
			}
			exponent += exponentTemp * z[imod];
			if (isnan(exponent))
			{
				FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
				FILE_LOG(logWARNING) << "exponentTemp = "<< exponentTemp;
				FILE_LOG(logWARNING) << "z["<<imod<<"] = "<< z[imod];
				FILE_LOG(logWARNING) << "exponent = "<< exponent;
				exit(1);
			}
		}
		logdens[t] = log(1/sqrt(this->cor_matrix_determinant)) - 0.5 * exponent + sum;
		if (isnan(logdens[t]))
		{
			FILE_LOG(logWARNING) << __PRETTY_FUNCTION__;
			FILE_LOG(logWARNING) << "cor_matrix_determinant = " << this->cor_matrix_determinant;
			FILE_LOG(logWARNING) << "sum = " << sum;
			FILE_LOG(logWARNING) << "exponentTemp = " << exponentTemp;
			FILE_LOG(logWARNING) << "exponent = " << exponent;
			FILE_LOG(logWARNING) << "logdens["<<t<<"] = " << logdens[t];
			exit(1);
		}		
	}

	// Clean up
	freeDoubleMatrix(marginals_logdensities, this->Nmod);
	freeDoubleMatrix(marginals_CDFs, this->Nmod);
	free(z);
}

void MVCopulaApproximation::calc_densities(double* dens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->calc_logdensities(dens);
	for (int t=0; t<this->T; t++)
	{
		dens[t] = exp( dens[t] );
		if (dens[t] == 0)
		{
			dens[t] = 0.00000000001; // TODO: find a non arbitrary number
		}
	}
}

DensityName MVCopulaApproximation::getType()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(Other);
}

	
// ------------------------------------------------------------
// Multivariate Product of Bernoullis
// ------------------------------------------------------------

BernoulliProduct::BernoulliProduct(double** multiobservations, bool* binary_states, int T, int Nmod)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->multiO = multiobservations;
	this->binary_states = binary_states;
	this->T = T;
	this->Nmod = Nmod;
}

BernoulliProduct::~BernoulliProduct()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
}

void BernoulliProduct::calc_logdensities(double* logdens)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	double d, mult;
	double** tempPost = allocDoubleMatrix(this->Nmod, this->T);

	for (int t=0; t<this->T; t++)
	{
		d = 1.0;
		for (int imod=0; imod<this->Nmod; imod++)
		{
			//if state[iN] is such that modification imod is unmodified, multiProb[t][imod] is the univariate posterior of being unmodified. 
			//if state[iN] is such that modification imod is modified, multiProb[t][imod] is the univariate posterior of being modified
			if (binary_states[imod])
			{
				mult = 1-this->multiO[imod][t];
			}
			else
			{
				mult = this->multiO[imod][t];
			}
			if(mult>=1) mult=0.9999999999999;
			if(mult<=0) mult=0.0000000000001;
			d=d*mult;
		}
		logdens[t] = log(d);
	}
	freeDoubleMatrix(tempPost, this->Nmod);
}

DensityName BernoulliProduct::getType()
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	return(Other);
}

