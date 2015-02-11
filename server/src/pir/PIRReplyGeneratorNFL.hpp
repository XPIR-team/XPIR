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
#include "GenericPIRReplyGenerator.hpp"
#include "../crypto/HomomorphicCrypto.hpp"
#include "../crypto/NFLLWE.hpp"
#include "../server/DBGenerator.hpp"

class PIRReplyGeneratorNFL : public GenericPIRReplyGenerator
{
private:
    bool lwe;
    uint64_t currentMaxNbPolys;
#ifdef SHOUP
    lwe_query ***queriesBuf;
#else
    lwe_query **queriesBuf;
#endif
    lwe_in_data* input_data;
    //NFLLWE* cryptoMethod;
    LatticesBasedCryptosystem* cryptoMethod;
    void generateReply(	lwe_query *queries,
                       lwe_in_data* data,
                       unsigned int lvl,
                       lwe_cipher* result);
    char** exportResult(lwe_cipher** inter_reply, unsigned int absp_nbr, unsigned int max_chunk, unsigned int rec_lvl);
    lwe_cipher** importDatabase(char** database, unsigned int x_size, unsigned int y_size, unsigned int absp_nbr);
    lwe_in_data* fromResulttoInData(lwe_cipher** inter_reply, uint64_t reply_elt_nbr, unsigned int reply_rec_lvl);
    void setPirParams(PIRParameters& param);
    void importFakeData(uint64_t plaintext_nbr);
    void pushFakeQuery();
    void generateReply();
    size_t getTotalSystemMemory();
    void freeInputData();
    void freeFakeInputData();
    void freeResult();
    void freeQuery();
public:
    PIRReplyGeneratorNFL();
    PIRReplyGeneratorNFL(vector <std::string>& database_, PIRParameters& param, DBHandler *db);
    ~PIRReplyGeneratorNFL();
    void importDataNFL(uint64_t offset, uint64_t bytes_per_file);
    void initQueriesBuffer();
    database_t generateReplyGeneric(bool keep_imported_data = false);
    void generateReplyGenericFromData(const database_t database);
    double generateReplySimulation(const PIRParameters& pir_params, uint64_t plaintext_nbr);
    double precomputationSimulation(const PIRParameters& pir_params, uint64_t plaintext_nbr);
    unsigned long computeReplySizeInChunks(unsigned long int);
    void pushQuery(char* rawQuery, unsigned int size, int dim, int nbr);
    void setCryptoMethod(CryptographicSystem* cm);
};

#endif
