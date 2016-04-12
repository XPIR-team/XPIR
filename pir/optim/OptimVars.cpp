/* Copyright (C) 2014 Carlos Aguilar Melchor, Joris Barrier, Marc-Olivier Killijian
 * This file is part of XPIR.
 *
 *  XPIR is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  XPIR is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XPIR.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "OptimVars.hpp"
/**
 *  Costs get from: https://aws.amazon.com/fr/ec2/pricing/ and we set Irland localisation.
 *  We take the first 4 cores configuration: https://aws.amazon.com/fr/ec2/instance-types/ named m1.xlarge
 **/
static const double kSecondsInHour = 3600;
static const double kOneHourComputingCost = 0.266; //Amortize cost for m1.xlarge with one years subscribtion
static const long double kOneByteUploadCost = 0.00000000005;
// Default constructor with fitness method to compute total costs : MAX
OptimVars::OptimVars():
  fitness(MAX)
{
  reset();
}


// More specific constructor
OptimVars::OptimVars(FitnessType fitness_):
  fitness(fitness_)
{
 reset();
}
OptimVars::OptimVars(FixedVars fixed_vars):
  fixedVars(fixed_vars)
{ OptimVars();}

OptimVars::OptimVars(FitnessType fitness_, FixedVars fixed_vars):
  fitness(fitness_),
  fixedVars(fixed_vars)
{ reset(); }

//Destructor
OptimVars::~OptimVars() {}

// Method to change the fitness method
void OptimVars::setType(FitnessType fitness_)
{
  fitness = fitness_;
}


// Getters and setters
double OptimVars::getGenQ() { return costs[GENQ]; }
void OptimVars::setGenQ(double GenQ) { costs[GENQ] = GenQ; }

double OptimVars::getSendQ() { return costs[SENDQ]; } 
void OptimVars::setSendQ(double SendQ) { costs[SENDQ] = SendQ; }

double OptimVars::getGenR() { return costs[GENR]; }
void OptimVars::setGenR(double GenR) { costs[GENR] = GenR; }

double OptimVars::getSendR() { return costs[SENDR]; }
void OptimVars::setSendR(double SendR) { costs[SENDR] = SendR; }

double OptimVars::getDecR() { return costs[DECR]; }
void OptimVars::setDecR(double DecR) { costs[DECR] = DecR; }

unsigned int OptimVars::getAlpha() { return alpha; }
void OptimVars::setAlpha(unsigned int alpha_) { alpha = alpha_; }

unsigned int OptimVars::getDim() { return d; }
void OptimVars::setDim(unsigned int d_) { d = d_; } 

FixedVars OptimVars::getFixedVars() { return fixedVars; }
void OptimVars::setFixedVars(FixedVars fixed_vars) { fixedVars = fixed_vars; }
// Getter for total value for the current fitness method
#include <iostream>
double OptimVars::getValue()
{
  switch(fitness)
  {
    case SUM :
      {
        return (getGenQ() + getSendQ()) + (getGenR()+ getSendR() + getDecR());
      }
    case MAX :
      {
        return std::max(getGenQ(), getSendQ()) + std::max(std::max(getGenR(), getSendR()), getDecR());
      }
    case CLOUD :
      {
        return getGenR() * (kOneHourComputingCost / kSecondsInHour) + (getSendR() * std::min(fixedVars.Tups, fixedVars.Tdoc)) * kOneByteUploadCost;
      }
    default :
      {
        return std::max(getGenQ(), getSendQ()) + std::max(std::max(getGenR(), getSendR()), getDecR());
      }
  }
}


// Access and comparison operators
double OptimVars::operator[](unsigned int i)
{
	assert(i < COST_NBR);
	return costs[i];
}
bool OptimVars::operator==(OptimVars &other) 
{
	return (getValue() ==  other.getValue()) ? true : false;
}
bool OptimVars::operator<(OptimVars &other) 
{
	return (getValue() < other.getValue()) ? true : false;
}
bool OptimVars::operator>(OptimVars &other) 
{
	return (getValue() > other.getValue()) ? true : false;
}
void OptimVars::operator=(OptimVars &other)
{
	for (unsigned i = 0 ; i  < COST_NBR ; i++) costs[i] = other.costs[i];
	alpha = other.alpha;
	d = other.d;
	crypto_params = std::string(other.crypto_params);
	fitness = other.fitness;
  fixedVars = other.fixedVars;
}


// Reset all costs to the maximum (to restart optimization)
void OptimVars::reset()
{
  alpha = d = 1;
  crypto_params = "No crypto params";

	for (unsigned i = 0 ; i < COST_NBR ; i++){
    costs[i] = std::numeric_limits<double>::max();
    limits[i] = std::numeric_limits<double>::min();
  }
}

