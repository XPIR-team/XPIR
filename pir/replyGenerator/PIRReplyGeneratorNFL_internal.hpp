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

#ifndef DEF_PIRREPLYGENERATORNFL
#define DEF_PIRREPLYGENERATORNFL

#include <unistd.h>
#include <omp.h>
#include "pir/replyGenerator/GenericPIRReplyGenerator.hpp"
#include "crypto/HomomorphicCrypto.hpp"
#include "crypto/NFLLWE.hpp"
#include "pir/dbhandlers/DBGenerator.hpp"

class PIRReplyGeneratorNFL_internal : public GenericPIRReplyGenerator
{
public:
    void setPirParams(PIRParameters& param);
private:
    bool lwe;
    uint64_t current_query_index;
    uint64_t current_dim_index;
#ifdef SHOUP
    lwe_query ***queriesBuf;
#else
    lwe_query **queriesBuf;
#endif
    //NFLLWE* cryptoMethod;
    void generateReply(	lwe_query *queries,
                       lwe_in_data* data,
                       unsigned int lvl,
                       lwe_cipher* result);
    char** exportResult(lwe_cipher** inter_reply, unsigned int absp_nbr, unsigned int max_chunk, unsigned int rec_lvl);
    lwe_cipher** importDatabase(char** database, unsigned int x_size, unsigned int y_size, unsigned int absp_nbr);
    lwe_in_data* fromResulttoInData(lwe_cipher** inter_reply, uint64_t reply_elt_nbr, unsigned int reply_rec_lvl);
    void importFakeData(uint64_t plaintext_nbr);
    void pushFakeQuery();
    void freeInputData();
    void freeFakeInputData();
protected:
    void freeResult();
    void freeQueries();
    void freeQueriesBuffer();
    void generateReply();
    imported_database_t generateReplyGeneric(bool keep_imported_data = false);
    void generateReplyGenericFromData(const imported_database_t database);
    size_t getTotalSystemMemory();
    lwe_in_data* input_data;
    LatticesBasedCryptosystem* cryptoMethod;
    uint64_t currentMaxNbPolys;
	
public:
    PIRReplyGeneratorNFL_internal();
    PIRReplyGeneratorNFL_internal(PIRParameters& param, DBHandler *db);
    ~PIRReplyGeneratorNFL_internal();
    void importDataNFL(uint64_t offset, uint64_t bytes_per_file);
    void initQueriesBuffer();
    void generateReplyExternal(imported_database_t* database);
    double generateReplySimulation(const PIRParameters& pir_params, uint64_t plaintext_nbr);
    double precomputationSimulation(const PIRParameters& pir_params, uint64_t plaintext_nbr);
    unsigned long computeReplySizeInChunks(unsigned long int);
    void pushQuery(char* rawQuery, unsigned int size, int dim, int nbr);
    void pushQuery(char* rawQuery);
    void setCryptoMethod(CryptographicSystem* cm);
};

#endif
