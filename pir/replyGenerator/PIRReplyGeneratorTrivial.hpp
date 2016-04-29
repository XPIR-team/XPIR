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

#ifndef DEF_PIRREPLYGENERATORTRIVIAL
#define DEF_PIRREPLYGENERATORTRIVIAL

#include "pir/replyGenerator/GenericPIRReplyGenerator.hpp"
#include "crypto/HomomorphicCrypto.hpp"
#include "crypto/NoCryptography.hpp"

class PIRReplyGeneratorTrivial : public GenericPIRReplyGenerator
{
private:
//  boost::mutex size_mutex;
//  bool lwe;
//    uint64_t currentMaxNbPolys;
//    lwe_query **queriesBuf;
    uint64_t totalNbChunks;
    NoCryptography * cryptoMethod;
    char* input_data;
    bool firstTimeImport;
    //    LatticesBasedCryptosystem* cryptoMethod;
//    void generateReply(	lwe_query *queries,
//                       lwe_in_data* data,
//                       unsigned int lvl,
//                       lwe_cipher* result);
//    char** exportResult(lwe_cipher** inter_reply, unsigned int absp_nbr, unsigned int max_chunk, unsigned int rec_lvl);
//    lwe_cipher** importDatabase(char** database, unsigned int x_size, unsigned int y_size, unsigned int absp_nbr);
//    void importFakeData(uint64_t file_kb_size);
//    void pushFakeQuery();
    void generateReply();
    void freeResult();
//    void freeInputData();
//    void freeFakeInputData();
public:
//    PIRReplyGeneratorNFL_internal();
    PIRReplyGeneratorTrivial();
    PIRReplyGeneratorTrivial(PIRParameters& param, DBHandler *db);
    ~PIRReplyGeneratorTrivial();
    void importData();
    void initQueriesBuffer();
    imported_database_t generateReplyGeneric(bool keep_imported_data = false);
    void generateReplyGenericFromData(const imported_database_t database);
    double generateReplySimulation(const PIRParameters& pir_params, uint64_t plaintext_nbr);
    unsigned long computeReplySizeInChunks(unsigned long int);
    void pushQuery(char* rawQuery, unsigned int size, int dim, int nbr);
    void setCryptoMethod(CryptographicSystem* cm);
};

#endif
