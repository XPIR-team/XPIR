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

#include "PIRReplyGeneratorGMP.hpp"

PIRReplyGeneratorGMP::PIRReplyGeneratorGMP()
{}


PIRReplyGeneratorGMP::PIRReplyGeneratorGMP( PIRParameters& param, DBHandler *db) :
  GenericPIRReplyGenerator(param,db),
	finished(false) 
{}

/**
 * Read raw data from files and make padding if necessary
 **/
void PIRReplyGeneratorGMP::importData() 
{
	unsigned int size = static_cast<double>(cryptoMethod->getPublicParameters().getAbsorptionBitsize()/GlobalConstant::kBitsPerByte), theoretic_nbr_elements = 1 ;
	char rawBits[size];
  uint64_t nbFiles = dbhandler->getNbStream();

	maxChunkSize = ceil((double)(dbhandler->getmaxFileBytesize())/(double)size);

  // Ugly ugly trick until this class is improved
  (*(PaillierPublicParameters *) &cryptoMethod->getPublicParameters()).getPubKey()->complete_key(pirParam.d+1);
  
	for (unsigned int i = 0 ; i < pirParam.d ; i++)
		theoretic_nbr_elements *= pirParam.n[i];

	datae = new mpz_t*[theoretic_nbr_elements];

	// For each real file
	for (unsigned int i = 0 ; i < nbFiles ; i++)
	{
		if (i % pirParam.alpha == 0) datae[i/pirParam.alpha] = new mpz_t[maxChunkSize*pirParam.alpha];

    ifstream* stream=dbhandler->openStream(i, 0);

    // For each chunk of size "size" of the file
		for (unsigned int j = 0 ; j < maxChunkSize ; j++ )
		{
			dbhandler->readStream(stream,rawBits, size);
			mpz_init(datae[i/pirParam.alpha][j + (i % pirParam.alpha) * maxChunkSize]);
			mpz_import(datae[i/pirParam.alpha][j + (i % pirParam.alpha) * maxChunkSize], size, 1, sizeof(char), 0, 0, rawBits);
		}
    // Pad if needed the last group of pirParam.alpha files
    if (i == nbFiles - 1)
    {
      for (unsigned int j = 0 ; j < maxChunkSize * (pirParam.alpha - 1 - (i % pirParam.alpha)); j++) 
      {
			mpz_init_set_ui(datae[i/pirParam.alpha][j + ((i+1) % pirParam.alpha) * maxChunkSize], 0);
      }
    }
		dbhandler->closeStream(stream);
	}

  std::cout << "PIRReplyGeneratorGMP: " << pirParam.alpha*theoretic_nbr_elements - nbFiles << "  non-aggregated padding files need to be added ..." << std::endl;

	/* make file padding if necessary **/
	for (uint64_t i = ceil((double)nbFiles/pirParam.alpha) ; i < theoretic_nbr_elements ; i++) 
	{
		datae[i] = new mpz_t[maxChunkSize*pirParam.alpha];
		for (unsigned int k = 0 ; k < maxChunkSize*pirParam.alpha ; k++)
		{
			mpz_init_set_ui(datae[i][k], 0);
		}
	}
  // From now on all data sizes take into account aggregation
	maxChunkSize *= pirParam.alpha; 
}

void PIRReplyGeneratorGMP::importFakeData(uint64_t plaintext_nbr)
{
  unsigned int file_nbr = 1;
  unsigned int abs_size = cryptoMethod->getPublicParameters().getAbsorptionBitsize()/GlobalConstant::kBitsPerByte;
  unsigned int ciph_size = cryptoMethod->getPublicParameters().getCiphertextBitsize()/GlobalConstant::kBitsPerByte;
  char *raw_data, *raw_data2;
  
  for (unsigned int i = 0 ; i < pirParam.d ; i++)
		file_nbr *= pirParam.n[i];	
  
  datae = new mpz_t*[file_nbr];
  raw_data = (char*) malloc(abs_size + 1);	
  raw_data2 = (char*) malloc(ciph_size + 1);	
  memset(raw_data, 0xaa, abs_size);
  memset(raw_data2, 0xaa, ciph_size);
  maxChunkSize = plaintext_nbr;

  for (unsigned int i = 0 ; i < file_nbr; i++)
  {
    datae[i] = new mpz_t[maxChunkSize];
    for (unsigned int j = 0 ; j < maxChunkSize ; j++)
    {
      mpz_init(datae[i][j]);
      mpz_import(datae[i][j], abs_size, 1, sizeof(char), 0, 0, raw_data);
    }
  }

	for (unsigned int i = 0 ; i < pirParam.n[0] ; i++)
  {
    mpz_import(queriesBuf[0][i], ciph_size, 1, sizeof(char), 0, 0, raw_data2);
  }
  free(raw_data); 
  free(raw_data2); 
}

void PIRReplyGeneratorGMP::clearFakeData(uint64_t plaintext_nbr)
{
  unsigned int file_nbr = 1;
  unsigned int abs_size = cryptoMethod->getPublicParameters().getAbsorptionBitsize()/GlobalConstant::kBitsPerByte;
  char* raw_data;
  
  for (unsigned int i = 0 ; i < pirParam.d ; i++)
		file_nbr *= pirParam.n[i];	
  
  maxChunkSize = plaintext_nbr;

  for (unsigned int i = 0 ; i < file_nbr; i++)
  {
    for (unsigned int j = 0 ; j < maxChunkSize ; j++)
    {
      mpz_clear(datae[i][j]);
    }
  }
  delete[] datae;
}

double PIRReplyGeneratorGMP::generateReplySimulation(const PIRParameters& pir_params, uint64_t plaintext_nbr)
{
  setPirParams((PIRParameters&) pir_params);
  initQueriesBuffer();
  
  importFakeData(plaintext_nbr);

  double start = omp_get_wtime();
  generateReply();
  double stop = omp_get_wtime() - start;

  cleanQueryBuffer();
  clearFakeData(plaintext_nbr);
  freeResult();
  return stop;
}

imported_database_t PIRReplyGeneratorGMP::generateReplyGeneric(bool keep_imported_data = false)
{
  imported_database_t database_wrapper;
	importData();
	boost::mutex::scoped_lock l(mutex);
	//mutex.lock();
  generateReply();

  if(keep_imported_data) 
  {
    database_wrapper.imported_database_ptr = (void*)datae;
    database_wrapper.polysPerElement = maxChunkSize;
  } 
  else 
  {
    unsigned long theoretic_nbr_elements = 1;
	  for (unsigned int i = 0 ; i < pirParam.d ; i++)
		  theoretic_nbr_elements *= pirParam.n[i];
    for (unsigned int i = 0 ; i < theoretic_nbr_elements ; i++ ) 
    {
      for (unsigned int j = 0; j < maxChunkSize; j++)
      {
        mpz_clear(datae[i][j]);
      }
      delete[] datae[i];
    }
    delete[] datae;
	
    for (unsigned int i = 0 ; i < pirParam.d ; i++)
	  {
		  for (unsigned int j = 0 ; j < pirParam.n[i] ; j++)
      {
        mpz_clear(queriesBuf[i][j]);
      }
      delete[] queriesBuf[i];
	  }
    delete[] queriesBuf;
  }

  return database_wrapper;
}

void PIRReplyGeneratorGMP::generateReplyGenericFromData(const imported_database_t database)
{
  datae = (mpz_t**) database.imported_database_ptr;
  maxChunkSize = database.polysPerElement;
	boost::mutex::scoped_lock l(mutex);
  //mutex.lock();
  generateReply();
}
/**
 * Compute Lipmaa PIR Scheme with imported data.
 **/
void PIRReplyGeneratorGMP::generateReply() 
{
	uint64_t pir_nbr;
	double start;

	mpz_t **data_z = datae;
	mpz_t **reply_vec;
	mpz_t *queries;
	omp_set_nested(0);

	start = omp_get_wtime();
#ifdef PERF_TIMERS
  double vtstart = start;
  bool wasVerbose = false;
#endif

	for (unsigned int i = 0 ; i < pirParam.d ; i++) // For each dimension
	{
    pir_nbr = 1;
   	for (unsigned int j = i+1  ; j < pirParam.d ; j++)
		  pir_nbr *= pirParam.n[j];
    
    if (pir_nbr!=1) std::cout << "PIRReplyGeneratorGMP: Generating " << pir_nbr << " replies in recursion level " << i+1 << std::endl; 

		reply_vec = new mpz_t*[pir_nbr]; 

		queries = queriesBuf[i];

		if(pir_nbr == 1) omp_set_nested(1);
//#pragma omp parallel for
		for (unsigned int j = 0 ; j < pir_nbr ; j++) // Do pir_nbr PIR iterations
		{
			reply_vec[j] = new mpz_t[maxChunkSize];

			//Do a PIR iteration between a sub-list of elements and a query, i + 1 being the dimension
			generateReply(queries , data_z, pirParam.n[i] * j ,  i,  reply_vec[j]);

#ifdef PERF_TIMERS
      // Give some feedback if it takes too long
      double vtstop = omp_get_wtime();
      if (vtstop - vtstart > 1) 
      {
        vtstart = vtstop;
        if (pir_nbr!=1) std::cout <<"PIRReplyGeneratorGMP: Reply " << j+1 << "/" << pir_nbr << " generated\r" << std::flush;
        wasVerbose = true;
      }
#endif

    }
    // Always print feedback for last reply 
    if (pir_nbr!=1) std::cout <<"PIRReplyGeneratorGMP: Reply " << pir_nbr << "/" << pir_nbr << " generated" << std::endl;

    // Delete intermediate data obtained on the recursions
    if(i!=0)
    {
      for (unsigned int j = 0 ; j < pirParam.n[i] ; j++ ) delete[] data_z[j];
		  delete[] data_z;
    } 

		data_z = reply_vec;
	}

	printf( "PIRReplyGeneratorGMP: Global reply generation took %f seconds\n", omp_get_wtime() - start); 

  repliesArray = (char**)calloc(maxChunkSize,sizeof(char*));
  repliesAmount = maxChunkSize;
	pushReply(data_z[0], 0, maxChunkSize);

	for (unsigned int i = 0 ; i < maxChunkSize ; i++)
		mpz_clear(data_z[0][i]);

	delete[] data_z[0];
	delete[] data_z;
}

/**
 * Compute single parallelizable PIR with almost fully homomorphic crypto system
 * Params :
 *  - mpz_t* queries :  queries array to compute the pir ;
 *  - mpz_t** data :  two dimensionnal array for raw data to treat ;
 *  - int begin_data : begin data index
 *  - int s : current dimension
 *  - mpz_t* result : result array, no need to init
 **/
	void 
PIRReplyGeneratorGMP::generateReply(mpz_t *queries, 
		mpz_t** data, int begin_data,    
		int dimension,
		mpz_t* result) 
{
  int affichages = 1;
  unsigned int data_size =  pirParam.n[dimension];
	mpz_t replyTmp;
  int init_s = (*((PaillierPublicParameters*) &(cryptoMethod->getPublicParameters()))).getPubKey()->getinit_s();

#ifdef PERF_TIMERS
  bool wasVerbose = false;
  double vtstart = omp_get_wtime();
#endif

  //#pragma omp parallel for private(replyTmp) firstprivate(s, begin_data)
	for (unsigned int chunk = 0 ; chunk < maxChunkSize ; chunk++)
	{
		mpz_inits(result[chunk], replyTmp, NULL);
		computeMul(queries[0], data[begin_data][chunk], result[chunk], dimension+init_s+1);
    if(dimension != 0) mpz_clear(data[begin_data][chunk]);

//if ((chunk*chunkCost) << 21 == 0 || chunk == maxChunkSize - 1) std::cout <<"PIRReplyGeneratorGMP: Dealing with chunk " << chunk+1 << "/" << maxChunkSize << std::endl;

    for (unsigned int file = 1, k = begin_data + 1 ; file < data_size ; file++, k++)
		{
			computeMul(queries[file], data[k][chunk], replyTmp, dimension+init_s+1);
			
			if(dimension != 0) mpz_clear(data[k][chunk]);
			
			//We add the filechunks of index chunk for the files of index file
			// eg :  file[1]->chunk[1] + file[2]->chunk[1] + ... + file[file]->chunk[chunk]
			computeSum(result[chunk], replyTmp, dimension+init_s+1);
		}
		mpz_clear(replyTmp);

#ifdef CRYPTO_DEBUG
    gmp_printf("PIRReplyGeneratorGMP: Reply chunk generated %Zd\n\n",result[chunk]); 
#endif

#ifdef PERF_TIMERS
    // Give some feedback if it takes too long
    double vtstop = omp_get_wtime();
    if (vtstop - vtstart > 1) 
    {
      vtstart = vtstop;
      if(maxChunkSize != 1) std::cout <<"PIRReplyGeneratorGMP: Dealt with chunk " << chunk+1 << "/" << maxChunkSize << "\r" << std::flush;
      wasVerbose = true;
    }
#endif

  }

#ifdef PERF_TIMERS
  if (wasVerbose) std::cout <<"                                                     \r" << std::flush;
#endif

}

/**
 * Performs homomorphic multiplication between query and n then put the result in res
 * Params :
 *  - mpz_t query : a pir query
 *  - mpz_t n : raw data
 *  - mpz_t res : result of operation
 *  - int s : dimension
 **/
void PIRReplyGeneratorGMP::computeMul(mpz_t query, mpz_t n, mpz_t res, int modulus_index) 
{	
	cryptoMethod->e_mul_const(res, query, n, modulus_index );
#ifdef CRYPTO_DEBUG
  gmp_printf("PIRReplyGeneratorGMP: Raising %Zd\n To the power %Zd\nResult is %Zd\nModulus index is %d\n\n",query, n, res,modulus_index); 
#endif
}

/**
 * Performs homomorphic addition betwee two encrypted data and put the result in a
 * Params :
 * 	- mpz_t a : first data to multiply, store result also
 * 	- mpz_t b : second data to multoply
 * 	- int s : dimension
 **/
void PIRReplyGeneratorGMP::computeSum(mpz_t a, mpz_t b, int modulus_index) 
{
	cryptoMethod->e_add(a, a, b, modulus_index);
}

/**
 * Inits queries buffer, used by PIRServer before reveiving client queries.
 **/
void PIRReplyGeneratorGMP::initQueriesBuffer()
{
	queriesBuf = new mpz_t*[pirParam.d];

	for (unsigned int i = 0 ; i < pirParam.d ; i++)
	{
		queriesBuf[i] = new mpz_t[pirParam.n[i]];

		for (unsigned int j = 0 ; j < pirParam.n[i] ; j++)
			mpz_init(queriesBuf[i][j]);
	}

}

void PIRReplyGeneratorGMP::setCryptoMethod(CryptographicSystem* crypto_method) 
{
	cryptoMethod = (PaillierAdapter*) crypto_method;
}

void PIRReplyGeneratorGMP::pushQuery(char* rawQuery, unsigned int size, int dim, int nbr)
{
	mpz_import(queriesBuf[dim][nbr], size, 1, sizeof(char), 0, 0, rawQuery);
#ifdef CRYPTO_DEBUG
  gmp_printf("Imported query element: %Zd\n", queriesBuf[dim][nbr]);
#endif 
}

void PIRReplyGeneratorGMP::pushReply(mpz_t* replies, unsigned init_index, unsigned replies_nbr)
{
	size_t n;
	unsigned int size = cryptoMethod->getPublicParameters().getCiphBitsizeFromRecLvl(pirParam.d)/GlobalConstant::kBitsPerByte;
	char*ct, *tmp;

	for (int i = init_index ; i < init_index + replies_nbr; i++)
	{
		tmp = (char*)mpz_export(NULL, &n, 1, sizeof(char) , 0, 0, replies[i]);

		if (n < size)
		{
			ct = new char[size]();
			memcpy(ct+sizeof(char)*(size - n), tmp, n);
			repliesArray[i] = ct;
			free(tmp);
		}
		else
    {
			repliesArray[i] = tmp;
    }
	}
}

bool PIRReplyGeneratorGMP::isFinished()
{
	return finished;
}

unsigned long int PIRReplyGeneratorGMP::computeReplySizeInChunks(unsigned long int maxFileBytesize)
{

	float res = ceil(static_cast<float>(maxFileBytesize*pirParam.alpha) / static_cast<float>(cryptoMethod->getPublicParameters().getAbsorptionBitsize()/8.0));
	return static_cast<unsigned long int>(res);
}

void PIRReplyGeneratorGMP::setPirParams(PIRParameters& _pirParam)
{
	pirParam = _pirParam;
}

void PIRReplyGeneratorGMP::cleanQueryBuffer()
{
	for (unsigned int i = 0 ; i < pirParam.d; i++)
	{
		for (unsigned int j = 0 ; j < pirParam.n[i]; j++)
		{
			mpz_clear(queriesBuf[i][j]);
		}
		delete[] queriesBuf[i];
	}
	delete[] queriesBuf;
}

PIRReplyGeneratorGMP::~PIRReplyGeneratorGMP()
{
}

void PIRReplyGeneratorGMP::freeResult()
{
  for(unsigned i=0 ; i < repliesAmount; i++)
  {
    if(repliesArray[i]!=NULL) delete[] repliesArray[i];
    repliesArray[i] = NULL;
  }
  free(repliesArray);
  repliesArray=NULL;
}

