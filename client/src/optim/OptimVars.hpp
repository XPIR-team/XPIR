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

#ifndef DEF_OPTIVAR
#define DEF_OPTIVAR

#include <assert.h>
#include <limits>
#include <math.h>
#include <algorithm>
#include <string>

// Fitness method to compute total value
// SUM: sum all costs
// MAX: sum the max costs of the client and server pipelines
enum FitnessType {SUM, MAX, CLOUD};

// Struct with the variables fixed before calling the optimizer
// These are considered mandatory (alternatives will not be explored)
struct FixedVars
{
  // Number of files in the database
  uint64_t n;
  // Length in bits of the largest file in the database
  uint64_t l; 
  // Client and Server upload Throughputs
  double Tupc, Tups;
  // Client and Server download Throughputs
  double Tdoc, Tdos;  
  // Force crypto parameters manually
  std::string manual_crypto_params;
  // Security bits requested for the cryptosystem
  unsigned int k;
  // Minimum level of recursion tested
  unsigned int dMin;
  // Maximum level of recursion tested
  unsigned int dMax;
  // Maximum aggregation tested (0 = no limit, 1 = no aggregation)
  unsigned int alphaMax; 
  // Fitness method to compute total cost
  FitnessType fitness;
};

// Object with the variables TO BE fixed by the optimizer
// Also contains the costs associated to these variables:
// Query generation and transmission, reply generation, transmission and decryption
// Able to compute from these costs the total cost for a given fitness method
// Variables start uninitialized and costs to the maximum values of the corresponding types
// Each time a set of variables improves the total cost in a test, variables and costs are updated
class OptimVars
{
  private:
    // Enum type to make more readable the cost array indexes
    enum {GENQ, SENDQ, GENR, SENDR, DECR, COST_NBR}; 
    // Array containing the costs for the different operations of a PIR instance
    double costs[COST_NBR];
    double limits[COST_NBR];
    // Fitness method (redundant with the one of fixed vars but useful for
    // comparison operators to know what fitness method is used to compare)
    FitnessType fitness;
    FixedVars fixedVars;

  public:
    // Default constructor (fitness=MAX)
    OptimVars();
    // More specific constructor and destructor
    OptimVars(FitnessType fitness);
    OptimVars(FixedVars fixed_vars);
    OptimVars(FitnessType fitness_, FixedVars fixed_vars);
    ~OptimVars();
    // Method to change the fitness method
    void setType(FitnessType fitness_);

    // Getters and setters for individual costs
    double getGenQ();
    void setGenQ(double GenQ);
    double getSendQ();
    void setSendQ(double SendQ);
    double getGenR();
    void setGenR(double GenR);
    double getSendR();
    void setSendR(double SendR);
    double getDecR();
    void setDecR(double DecR);
    unsigned int getAlpha();
    void setAlpha(unsigned int alpha_);
    unsigned int getDim();
    void setDim(unsigned int d_);

    FixedVars getFixedVars();
    void setFixedVars(FixedVars fixed_vars);
    // Variables TO BE fixed by the optimizer (don't feel like doing getters and setters)
    unsigned int alpha;
    unsigned int d;
    std::string crypto_params;

    // Getter for total value for the current fitness method
    virtual double getValue();

    // Access and comparison operators
    double operator[](unsigned int i);
    bool operator<(OptimVars &other);
    bool operator>(OptimVars &other);
    bool operator==(OptimVars &other);
    void operator=(OptimVars &other);

    // Reset all costs to the maximum (to restart optimization)
    void reset();
};

#endif

