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

#include "PIRQueryGenerator.hpp"
#include "../crypto/NFLLWE.hpp"
/**
 *	Class constructor
 *	Params:
 *		- PIRParameters& pirParameters_ : PIRParameters reference shared with PIRClient.
 *		- crypto_ptr cryptoMethod_ 			: shared_pointer of Homomorphic crypto.
 **/
PIRQueryGenerator::PIRQueryGenerator(PIRParameters& pirParameters_,HomomorphicCrypto& cryptoMethod_) : 
	pirParams(pirParameters_),
	cryptoMethod(cryptoMethod_),
  queryBuffer("query_buffer"),
	mutex()
{}

/**
 * Generates asyncronously queries for each files.
 * Makes encrypted of 0 or 1.
 **/
void PIRQueryGenerator::generateQuery() 
{
  double start = omp_get_wtime();
	coord = new unsigned int[pirParams.d]();

	computeCoordinates();
	for (unsigned int j = 0 ; j < pirParams.d ; j++)
	{
		for (unsigned int i = 0 ; i < pirParams.n[j] ; i++) 
		{
			if (i == coord[j]) queryBuffer.push(cryptoMethod.encrypt(1, j + 1 ));
			else queryBuffer.push(cryptoMethod.encrypt(0, j + 1));
	  }
  std::cout << "PIRQueryGenerator: Generated a " << pirParams.n[j] << " element query" << std::endl;
  }
  double end = omp_get_wtime();
  
  std::cout << "PIRQueryGenerator: All the queries have been generated, total time is " << end - start << " seconds" << std::endl;
}

/**
 *	Compute coordinates of the chosen file.
 **/
void PIRQueryGenerator::computeCoordinates()
{
	uint64_t x = chosenElement;

	for (unsigned int i = 0 ; i < pirParams.d ; i++)
	{
		coord[i] = x % pirParams.n[i];
		x /= pirParams.n[i];
	}
}

/**
 * Starts computation in a new thread
 **/
void PIRQueryGenerator::startGenerateQuery()	
{
	queryThread = thread(&PIRQueryGenerator::generateQuery, this);
}

uint64_t PIRQueryGenerator::getChosenElement() 
{
	return chosenElement;
}

void PIRQueryGenerator::setChosenElement( uint64_t _chosenElement ) 
{
	chosenElement = _chosenElement;
}

void PIRQueryGenerator::setPIRParameters(PIRParameters& pirParams_)
{
  pirParams = pirParams_;
}

/**
 *	Join query thread if it's possible.
 **/
void PIRQueryGenerator::joinThread() 
{
	if(queryThread.joinable()) queryThread.join();
}

PIRQueryGenerator::~PIRQueryGenerator() 
{
	joinThread();
  cleanQueryBuffer();

	delete[] coord;
}

void PIRQueryGenerator::cleanQueryBuffer()
{
	while (!queryBuffer.empty())
		free(queryBuffer.pop_front());
}
