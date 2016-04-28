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

#ifndef DEF_PIRREPLYGENERATOR
#define DEF_PIRREPLYGENERATOR

#include <omp.h>
#include <string>

#include <boost/thread.hpp>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>

#include "pir/PIRParameters.hpp"
#include "pir/replyGenerator/GenericPIRReplyGenerator.hpp"
#include "crypto/PaillierAdapter.hpp"

using namespace std;

class PIRReplyGeneratorGMP : public GenericPIRReplyGenerator 
{
	private:
		mpz_t **datae;
		mpz_t **queriesBuf;

    bool firstTime;
    bool finished;
		PaillierAdapter* cryptoMethod;
		
		void computeMul (mpz_t query, mpz_t n, mpz_t res, int);
		void computeSum (mpz_t a, mpz_t b, int);
		
		void pushReply(mpz_t* replies, unsigned init_reply, unsigned replies_size);

		void generateReply(	mpz_t *queries, 
													mpz_t** data, int begin_data,
													int s,
													mpz_t* result); 
		void importData();
    void importFakeData(uint64_t plaintext_nbr);
    void clearFakeData(uint64_t plaintext_nbr);

    void generateReply();
    void cleanQueryBuffer();
    void freeResult();

	public:

		PIRReplyGeneratorGMP();
		PIRReplyGeneratorGMP(PIRParameters& param, DBHandler *db);
		~PIRReplyGeneratorGMP();

    imported_database_t generateReplyGeneric(bool keep_imported_data);
    void generateReplyGenericFromData(const imported_database_t database);
    double generateReplySimulation(const PIRParameters& pir_params, uint64_t plaintext_nbr);

		void initQueriesBuffer();
		void setCryptoMethod(CryptographicSystem* generator);
		void setPirParams(PIRParameters&);
		shared_queue< char*>& getRepliesBuffer();
		bool isFinished();
		void pushQuery(char* rawQuery, unsigned int size, int dim, int nbr);

		unsigned long computeReplySizeInChunks(unsigned long int);
};
#endif
