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
	this->var = variance;
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
	return(this->var);
}

void Normal::set_variance(double variance)
{
	FILE_LOG(logDEBUG2) << __PRETTY_FUNCTION__;
	this->var = variance;
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
	this->var = stdev*stdev;
}

