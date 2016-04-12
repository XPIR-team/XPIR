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

#include "PIRReplyGeneratorTrivial.hpp"
#include "sys/time.h"

PIRReplyGeneratorTrivial::PIRReplyGeneratorTrivial():
  firstTimeImport(true)
{
}

PIRReplyGeneratorTrivial::PIRReplyGeneratorTrivial(PIRParameters& param, DBHandler *db):
  GenericPIRReplyGenerator(param, db),
  firstTimeImport(true)
{
}


//void PIRReplyGeneratorTrivial::importFakeData(uint64_t file_kb_size)
//{
//  uint64_t files_nbr = 1;
//
//  for (unsigned int i = 0 ; i < pirParam.d ; i++) files_nbr *= pirParam.n[i];
//
//  input_data = new lwe_in_data[files_nbr];
//  currentMaxNbPolys = ceil(double(file_kb_size * 1024) / double(cryptoMethod->getPublicParameters().getAbsorptionBitsize()));
//
//  for (unsigned int i = 0 ; i < files_nbr ; i++)
//  {
//    input_data[i].p = new poly64[currentMaxNbPolys];
//    for (unsigned int j = 0 ; j < currentMaxNbPolys ; j++)
//    {
//      input_data[i].p[j] = cryptoMethod->getnflInstance().allocBoundedRandomPoly(
//          cryptoMethod->getpolyDegree(), cryptoMethod->getmoduli()[0]); 
//      input_data[i].nbPolys = currentMaxNbPolys;
//    }
//  }
//}
//
//


void PIRReplyGeneratorTrivial::importData() 
{
	uint64_t maxFileBytesize = dbhandler->getmaxFileBytesize();
  unsigned int size = static_cast<double>(cryptoMethod->getPublicParameters().getAbsorptionBitsize()/GlobalConstant::kBitsPerByte);
  char* dataptr;

  totalNbChunks = ceil((double)maxFileBytesize*dbhandler->getNbStream()/double(size));

  input_data = (char *)calloc(size*totalNbChunks,1);
  dataptr = input_data;

	//pour tous les fichiers.
	for (unsigned int i = 0 ; i < dbhandler->getNbStream() ; i++)
	{
    ifstream* stream = dbhandler->openStream(i, 0);
	  dbhandler->readStream(stream, dataptr, maxFileBytesize);
		dbhandler->closeStream(stream);
    dataptr += maxFileBytesize;
	}

}


void PIRReplyGeneratorTrivial::generateReply()
{
	unsigned int size = static_cast<double>(cryptoMethod->getPublicParameters().getAbsorptionBitsize()/GlobalConstant::kBitsPerByte);
  char* buffer;

#ifdef PERF_TIMERS
  double vtstart = omp_get_wtime();
#endif

  // Init the reply array
  repliesArray = (char**) calloc(totalNbChunks, sizeof(char*));
  repliesAmount = totalNbChunks;
  for (unsigned int j = 0 ; j < totalNbChunks ; j++)
  {
    buffer = (char*) malloc(size*sizeof(char));
		memcpy(buffer, input_data + j*size, size);
    repliesArray[j] = buffer;
#ifdef PERF_TIMERS
    // Give some feedback if it takes too long
    double vtstop = omp_get_wtime();
    if (vtstop - vtstart > 1) 
    {
      vtstart = vtstop;
      std::cout <<"PIRReplyGeneratorTrivial: Reply chunk " << j+1 << "/" << totalNbChunks << " generated\r" << std::flush;
    }
#endif

  }
  // Always print feedback for last chunk
  std::cout <<"PIRReplyGenerator: Reply chunk " << totalNbChunks << "/" << totalNbChunks << " generated" << std::endl;

}


imported_database_t PIRReplyGeneratorTrivial::generateReplyGeneric(bool keep_imported_data){
  imported_database_t database_wrapper;
  
  if(firstTimeImport)
  {
	  importData();
    database_wrapper.imported_database_ptr = (void*) input_data;
    database_wrapper.polysPerElement = totalNbChunks;
    firstTimeImport = false;
  }

  boost::mutex::scoped_lock l(mutex);
  generateReply();

  if(keep_imported_data==false)
  {
    free(input_data);
  }

  return database_wrapper;
}


void PIRReplyGeneratorTrivial::generateReplyGenericFromData(const imported_database_t database)
{
  input_data = (char*) database.imported_database_ptr;
  totalNbChunks = database.polysPerElement;
	boost::mutex::scoped_lock l(mutex);
  generateReply();
}


double PIRReplyGeneratorTrivial::generateReplySimulation(const PIRParameters& pir_params, uint64_t plaintext_nbr)
{
  return (double) 0;
}

///**
// * Compute Reply Size une chunks. 
// * WARNING blocking function.
// **/
unsigned long PIRReplyGeneratorTrivial::computeReplySizeInChunks(unsigned long int maxFileBytesize)
{
	double ciphBytesize = 
    (double)cryptoMethod->getPublicParameters().getAbsorptionBitsize()
    /(double)GlobalConstant::kBitsPerByte;
  return ceil(double(maxFileBytesize) * dbhandler->getNbStream()
      / double(ciphBytesize)); 
}


///**
// *	Overloaded fonction from  GenericPIRReplyGenerator.
// *	Initalise queriesBuf.
// **/
//
void PIRReplyGeneratorTrivial::initQueriesBuffer() {}

void PIRReplyGeneratorTrivial::pushQuery(char* rawQuery, unsigned int size, int dim, int nbr)
{
  std::cout << "PIRReplyGeneratorTrivial: Ignoring query element" << std::endl;
}

PIRReplyGeneratorTrivial::~PIRReplyGeneratorTrivial()
{
  // Ensure that mutex is always unlocked upon destruction
  mutex.try_lock();
  mutex.unlock();
  freeResult();
}

void PIRReplyGeneratorTrivial::setCryptoMethod(CryptographicSystem* cm) {
	cryptoMethod = (NoCryptography *)cm;
}

void PIRReplyGeneratorTrivial::freeResult()
{
  if(repliesArray!=NULL)
  {
    for(unsigned i=0 ; i < repliesAmount; i++)
    {
      if(repliesArray[i]!=NULL) free(repliesArray[i]);
      repliesArray[i] = NULL;
    }
    free(repliesArray);
    repliesArray=NULL;
  }
}

