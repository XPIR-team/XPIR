/* Copyright (C) 2014 Carlos Aguilar Melchor, Joris Barrier, Marc-Olivier Killijian
 * This file is part of XPIRe.
 *
 *  XPIRe is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  XPIRe is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XPIRe.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "PIRReplyGeneratorNFL.hpp"
#include "sys/time.h"

#define PERFTIMERS
//#define SIMULATE_PRE_NFL_DATA //Use this to simulate imported data is in NFL form
//#define TEST_NFL_PERF_ITERATIONS //Loop to simulate very large databases
#define MULTI_THREAD
//#define SNIFFER_CHUNK_BYTESIZE 12800
//#define SNIFFER_CHUNK_BYTESIZE 2560
//#define SNIFFER_CHUNK_BYTESIZE 1500
//#define SNIFFER_CHUNK_BYTESIZE 600
//#define SNIFFER //Use this to activate a sniffer like behavior

PIRReplyGeneratorNFL::PIRReplyGeneratorNFL():
  lwe(false)
{
}

/**
 *	Constructor of the class.
 *	Params :
 *	- vector <File*>& database : reference of a File pointer vector.
 *	- PIRParameters& param : reference to a PIRParameters object.
 **/
PIRReplyGeneratorNFL::PIRReplyGeneratorNFL(vector <string>& database_, PIRParameters& param, DBHandler* db):
  lwe(false),
  currentMaxNbPolys(1),
  GenericPIRReplyGenerator(database_, param,db)
{
  // cryptoMethod will be set later by setCryptoMethod
}


void PIRReplyGeneratorNFL::importFakeData(uint64_t plaintext_nbr)
{
  uint64_t files_nbr = 1;
  for (unsigned int i = 0 ; i < pirParam.d ; i++) files_nbr *= pirParam.n[i];
  
  uint64_t plain_bytesize = cryptoMethod->getnflInstance().getpolyDegree()*cryptoMethod->getnflInstance().getnbModuli()*8; 
  dbhandler = new DBGenerator(files_nbr, plaintext_nbr*plain_bytesize, true);
  
  currentMaxNbPolys = plaintext_nbr;
  
  importDataNFL(0, plaintext_nbr*plain_bytesize);
}


/**
 *	Convert raw data from file in usable NFL data.
 **/
void PIRReplyGeneratorNFL::importDataNFL(uint64_t offset, uint64_t bytes_per_file)
{
	uint64_t fileByteSize = min(bytes_per_file, dbhandler->getmaxFileBytesize()-offset);
  uint64_t theoretical_files_nbr = 1;
  uint64_t nbFiles = dbhandler->getNbStream();

	for (unsigned int i = 0 ; i < pirParam.d ; i++) theoretical_files_nbr *= pirParam.n[i];

	input_data = new lwe_in_data[theoretical_files_nbr];
	// Because DBHandler does the filepadding (all files are maxFileSize) we don't need to calloc, malloc is enough
	char *rawBits = (char*)malloc(fileByteSize*pirParam.alpha*sizeof(char));

	currentMaxNbPolys=0;

#ifdef PERFTIMERS
  double vtstart = omp_get_wtime();
  bool wasVerbose = false;
  uint64_t lastindex = 0;
#endif
  // For global time measurement
  double start = omp_get_wtime();double now,delta;

  int nbruns=ceil((double)nbFiles/pirParam.alpha);
#ifdef MULTI_THREAD
  #pragma omp parallel for 
#endif
  for (int i=0; i < nbruns; i++)
	{
	  dbhandler->readAggregatedStream(i, pirParam.alpha, offset, bytes_per_file, rawBits);
    
#ifdef SIMULATE_PRE_NFL_DATA
    uint64_t abssize = cryptoMethod->getPublicParameters().getAbsorptionBitsize(); 
    uint64_t polysize = cryptoMethod->getpolyDegree() * cryptoMethod->getnbModuli()*sizeof(uint64_t);
    uint64_t nbpolys = ceil((double)fileByteSize * pirParam.alpha * 8 / abssize);
    input_data[i].p = (poly64*) malloc(nbpolys*sizeof(poly64*));
    input_data[i].p[0] = (poly64) malloc(nbpolys*polysize);
    for (unsigned j = 0; j < nbpolys ; j++) 
    {
      input_data[i].p[j] = input_data[i].p[0]+j*polysize/8;
      memcpy(input_data[i].p[j], rawBits, min(fileByteSize, polysize));
    }
    input_data[i].nbPolys = nbpolys; 
#else
    input_data[i].p = cryptoMethod->deserializeDataNFL((unsigned char**)&rawBits, (uint64_t) 1, fileByteSize*pirParam.alpha*GlobalConstant::kBitsPerByte, input_data[i].nbPolys);
#endif
    if (input_data[i].nbPolys>currentMaxNbPolys) currentMaxNbPolys=input_data[i].nbPolys;

#ifdef PERFTIMERS
    // Give some feedback if it takes too long
    double vtstop = omp_get_wtime();
    if (vtstop - vtstart > 1) 
    {
      vtstart = vtstop;
      std::cout <<"PIRReplyGeneratorNFL: Element " << i+1 << "/" << nbFiles << " imported\r" << std::flush;
     wasVerbose = true;
     lastindex = i+1;
    }
#endif
  }
		
#ifdef PERFTIMERS
    // If feedback was given say we finished
    if (wasVerbose && lastindex != nbFiles)  std::cout <<"PIRReplyGeneratorNFL: Element " << nbFiles << "/" << nbFiles << " imported" << std::endl;
#endif
  
  /** FILE PADDING **/
	for (uint64_t i = ceil((double)nbFiles/pirParam.alpha) ; i < theoretical_files_nbr ; i++) 
	{
	 	input_data[i].p = (poly64 *) malloc(currentMaxNbPolys*sizeof(poly64));
	 
	 	input_data[i].p[0] = (poly64) calloc(cryptoMethod->getpolyDegree()*cryptoMethod->getnbModuli()*currentMaxNbPolys,sizeof(uint64_t));
	 	for (uint64_t j = 1 ; j < currentMaxNbPolys ; j++) input_data[i].p[j] = input_data[i].p[0]+cryptoMethod->getpolyDegree()*cryptoMethod->getnbModuli()*j;
	 
	 	input_data[i].nbPolys = currentMaxNbPolys;
	}

	free(rawBits);
	std::cout<<"PIRReplyGeneratorNFL: Finished importing the database in " << omp_get_wtime() - start << " seconds" << std::endl;
}

#ifdef SNIFFER
database_t PIRReplyGeneratorNFL::generateReplyGeneric(bool keep_imported_data)
{
  database_t database_wrapper;
  boost::mutex::scoped_lock l(mutex);
  const uint64_t chunkBytesize = SNIFFER_CHUNK_BYTESIZE;
  const uint64_t iterations = dbhandler->getmaxFileBytesize()/(chunkBytesize+1);
  std::ifstream *is = dbhandler->openStream(0,0);
  const uint64_t nbFiles = dbhandler->getNbStream();
  const unsigned int polysize = cryptoMethod->getpolyDegree()*cryptoMethod->getnbModuli();
  const unsigned int jumpcipher = 2*polysize / sizeof(uint64_t);
	currentMaxNbPolys=0;
	lwe_in_data *input = new lwe_in_data[iterations];
  lwe_cipher *resul = new lwe_cipher[iterations];
  uint64_t index;
  lwe_query **queries;
	char **rawBits = (char**) malloc(iterations*sizeof(char*));

  cryptoMethod->setandgetAbsBitPerCiphertext(1,1);
  queries = queriesBuf[0];
  const uint64_t  mask = (1<<((int)log2(nbFiles)))-1;
	for (uint64_t it = 0 ; it < iterations ; it++)
  {
	  resul[it].a = (uint64_t *) malloc(2*polysize*sizeof(uint64_t));
	  resul[it].b = (uint64_t *) resul[it].a + polysize;
    rawBits[it] = (char*)malloc(((chunkBytesize)*sizeof(char)+sizeof(int)));
  }
  double start = omp_get_wtime();
  #pragma omp parallel for firstprivate(input,queries,resul)
	for (uint64_t it = 0 ; it < iterations ; it++)
  {
    dbhandler->readStream(0, rawBits[it], chunkBytesize+sizeof(int));
    input[it].p = cryptoMethod->deserializeDataNFL((unsigned char**)&(rawBits[it]), (uint64_t) 1, chunkBytesize*GlobalConstant::kBitsPerByte, input[it].nbPolys);
    index = *(int *)(rawBits[it]+chunkBytesize) & mask;
		cryptoMethod->mul(resul[it], input[it], queries[0][index],queries[1][index], 0, 0);
  }
  double end = omp_get_wtime();
	std::cout<<"PIRReplyGeneratorNFL: Finished processing the sniffed data in " << end - start << " seconds" << std::endl;
	std::cout<<"PIRReplyGeneratorNFL: Processing throughput " << (double)chunkBytesize*8*iterations / ((end - start)*1000000000ULL) << " Gbps" << std::endl;
  return database_wrapper;   
}

#else
database_t PIRReplyGeneratorNFL::generateReplyGeneric(bool keep_imported_data)
{
  database_t database_wrapper;
  uint64_t usable_memory, database_size, max_memory_per_file, max_readable_size, nbr_of_iterations;
  double start, end;

  // Init database_wrapper to NULL values so that we are able to know if it has been initialized
  database_wrapper.database_ptr = NULL;
  database_wrapper.elements_nbr = 0;

  // Don't use more than half of the computer's memory 
  usable_memory = getTotalSystemMemory()/2;
  database_size = dbhandler->getmaxFileBytesize() * dbhandler->getNbStream();
#ifndef TEST_NFL_PERF_ITERATIONS
  // This is the maximum amount of data per file we can get in memory
  max_memory_per_file = usable_memory / dbhandler->getNbStream();
  // Given the expansion factor of importation we get the max we can read per file
  max_readable_size = max_memory_per_file / 4 ; 
  // Reduce it so that we have full absorption in all the ciphertexts sent
  max_readable_size = (max_readable_size * GlobalConstant::kBitsPerByte / cryptoMethod->getPublicParameters().getAbsorptionBitsize(0)) * cryptoMethod->getPublicParameters().getAbsorptionBitsize(0)/GlobalConstant::kBitsPerByte;
#else
  // For our tests we will need to have databases of an integer amount of gigabits
  max_readable_size = 1280000000UL/dbhandler->getNbStream();
#endif
  // If we reduced it too much set it at least to a ciphertext
  if (max_readable_size == 0) max_readable_size = cryptoMethod->getPublicParameters().getAbsorptionBitsize(0);
  // Ensure it is not larger than maxfilebytesize
  max_readable_size = min(max_readable_size, dbhandler->getmaxFileBytesize());
  // Given readable size we get how many iterations we need
  nbr_of_iterations = ceil((double)dbhandler->getmaxFileBytesize()/max_readable_size);

#ifndef TEST_NFL_PERF_ITERATIONS
  // If aggregation is used we cannot iterate
  if ((pirParam.alpha != 1 || pirParam.d > 1) && nbr_of_iterations > 1)
  {
    std::cout << "PIRReplyGeneratorNFL: Cannot handle aggregation or dimensions on databases requiring multiple iterations" << std::endl;
    std::cout << "PIRReplyGeneratorNFL: Handling the database on a single iteration, this can cause memory issues ..." << std::endl;
    nbr_of_iterations = 1;
    max_readable_size = dbhandler->getmaxFileBytesize();
  }
  // If we cannot read the whole database we cannot store it precomputed
  if (nbr_of_iterations > 1) keep_imported_data = false;
#endif  

  // If we need to do more than an iteration say it
  if (nbr_of_iterations > 1)
  {
    std::cout << "PIRReplyGeneratorNFL: Database is considered too large, processing it in " 
      << nbr_of_iterations << " iterations" << std::endl; 
  }

  start = omp_get_wtime();
// #pragma omp parallel for
  for (unsigned iteration = 0; iteration < nbr_of_iterations; iteration++)
  {
    if (nbr_of_iterations > 1) cout << "PIRReplyGeneratorNFL: Iteration " << iteration << endl; 

    repliesIndex = computeReplySizeInChunks(iteration*max_readable_size);
    // Import a chunk of max_readable_size bytes per file with an adapted offset
    importDataNFL(iteration*max_readable_size, max_readable_size);
    if(keep_imported_data && iteration == nbr_of_iterations - 1)  // && added for Perf test but is no harmful
    {
      database_wrapper.elements_nbr = currentMaxNbPolys;
    } 

    boost::mutex::scoped_lock l(mutex);
    generateReply();
    end = omp_get_wtime();

    if(keep_imported_data && iteration == nbr_of_iterations - 1)  // && added for Perf test but is no harmful
    {
      database_wrapper.database_ptr = (void*)input_data;
    } 
    else 
    {
      freeInputData();
    }
  }

	std::cout<<"PIRReplyGeneratorNFL: Total process time " << end - start << " seconds" << std::endl;
	std::cout<<"PIRReplyGeneratorNFL: DB processing throughput " << 8*database_size/(end - start) << "bps" << std::endl;
	std::cout<<"PIRReplyGeneratorNFL: Client cleartext reception throughput  " << 8*dbhandler->getmaxFileBytesize()/(end - start) << "bps" << std::endl;
  freeQuery();

  return database_wrapper;
}
#endif
// Function used to generate a PIR reply if:
// - database is small enough to be kept in memory
// - it has already been imported to it
void PIRReplyGeneratorNFL::generateReplyGenericFromData(const database_t database)
{
#ifndef TEST_NFL_PERF_ITERATIONS
  input_data = (lwe_in_data*) database.database_ptr;
  currentMaxNbPolys = database.elements_nbr;
	boost::mutex::scoped_lock l(mutex);
  double start = omp_get_wtime();
  generateReply();
#else
  uint64_t max_readable_size, database_size, nbr_of_iterations;

  database_size = dbhandler->getmaxFileBytesize() * dbhandler->getNbStream();
  max_readable_size = 1280000000UL/dbhandler->getNbStream();
  // Ensure it is not larger than maxfilebytesize
  max_readable_size = min(max_readable_size, dbhandler->getmaxFileBytesize());
  // Given readable size we get how many iterations we need
  nbr_of_iterations = ceil((double)dbhandler->getmaxFileBytesize()/max_readable_size);


  boost::mutex::scoped_lock l(mutex);
  double start = omp_get_wtime();
  for (unsigned iteration = 0; iteration < nbr_of_iterations; iteration++)
  {

    input_data = (lwe_in_data*) database.database_ptr;
    currentMaxNbPolys = database.elements_nbr;
    generateReply();
  }
  freeInputData();
#endif
  double end = omp_get_wtime();
	std::cout<<"PIRReplyGeneratorNFL: Total process time " << end - start << " seconds" << std::endl;
	std::cout<<"PIRReplyGeneratorNFL: DB processing throughput " << 8*dbhandler->getmaxFileBytesize()*dbhandler->getNbStream()/(end - start) << "bps" << std::endl;
	std::cout<<"PIRReplyGeneratorNFL: Client cleartext reception throughput  " << 8*dbhandler->getmaxFileBytesize()/(end - start) << "bps" << std::endl;
  freeQuery();
}


/**
 *	Prepare reply and start absoptions.
 **/
void PIRReplyGeneratorNFL::generateReply()
{
  lwe_in_data *in_data = input_data;
  lwe_cipher **inter_reply;
#ifdef SHOUP
  lwe_query **queries;
#else
  lwe_query *queries;
#endif  
  uint64_t old_reply_elt_nbr = 0;
  uint64_t reply_elt_nbr = 1;
  uint64_t old_poly_nbr = 1;
  
  // Allocate memory for the reply array
  repliesAmount = computeReplySizeInChunks(dbhandler->getmaxFileBytesize());
  repliesArray = (char**)calloc(repliesAmount,sizeof(char*)); 


  // Start global timers
  double start = omp_get_wtime();

#ifdef PERFTIMERS
  double vtstart = start;
  bool wasVerbose = false;
#endif

  for (unsigned int i = 0 ; i < pirParam.d ; i++) // For each recursion level
  {
    old_reply_elt_nbr = reply_elt_nbr;
    reply_elt_nbr = 1;
    for (unsigned int j = i + 1 ; j < pirParam.d ; j++ ) reply_elt_nbr *= pirParam.n[j];
#ifdef DEBUG
  cout << "PIRReplyGeneratorNFL:  currentMaxNbPolys = " << currentMaxNbPolys << endl;
#endif

 
    inter_reply = new lwe_cipher*[reply_elt_nbr]();
    queries = queriesBuf[i];
    
    for (uint64_t j = 0 ; j < reply_elt_nbr ; j++) // Boucle de reply_elt_nbr PIR
    {
      inter_reply[j] = new lwe_cipher[currentMaxNbPolys];
	  // Warning of the trick in case SHOUP is defined : we cast quesries to a (lwe_query*) and will have to uncast it
      generateReply((lwe_query*)queries , in_data + (pirParam.n[i] * j ),  i,  inter_reply[j]);
#ifdef DEBUG_WITH_FILE_OUTPUT
      if (i ==0 && j==1) {
        std::ofstream file(std::string("output_level_"+ std::to_string(i)).c_str(), std::ios::out| std::ios::binary);
        for (int k = 0 ; k < currentMaxNbPolys ; k++)
        {
	        file.write((char*)inter_reply[j][k].a,1024*2*8);
        }
        file.close();
      }
#endif

#ifdef PERFTIMERS
      // Give some feedback if it takes too long
      double vtstop = omp_get_wtime();
      if (vtstop - vtstart > 1) 
      {
        vtstart = vtstop;
        std::cout <<"PIRReplyGeneratorNFL: Reply " << j+1 << "/" << reply_elt_nbr << " generated\r" << std::flush;
        wasVerbose = true;
      }
#endif

    }

    /*****************/
    /*MEMORY CLEANING*/
    /*****************/
    if ( i > 0)
    {
#ifdef DEBUG
      cout << "PIRReplyGeneratorNFL: reply_elt_nbr_OLD: " << old_reply_elt_nbr << endl;
#endif
     // for (unsigned int j = 0 ; j < old_reply_elt_nbr ; j++) {
     //   free(in_data[j].p[0]);
     //   free(in_data[j].p);
     // }
     // delete[] in_data;
    }
  // When i i=> 2 clean old in_data.
    if (i < pirParam.d - 1) { 
      old_poly_nbr = currentMaxNbPolys;
      in_data = fromResulttoInData(inter_reply, reply_elt_nbr, i);
    }

    for (uint64_t j = 0 ; j < reply_elt_nbr ; j++) {
      for (uint64_t k = 0 ; (k < old_poly_nbr) && (i < pirParam.d - 1); k++) free(inter_reply[j][k].a);

      delete[] inter_reply[j];
    }
    delete[] inter_reply; // allocated with a 'new' above. 
  }

  // Compute execution time
  printf( "PIRReplyGeneratorNFL: Global reply generation took %f (omp)seconds\n", omp_get_wtime() - start);
}


double PIRReplyGeneratorNFL::generateReplySimulation(const PIRParameters& pir_params, uint64_t plaintext_nbr)
{
  setPirParams((PIRParameters&)pir_params);
  initQueriesBuffer();
  pushFakeQuery();
  
  importFakeData(plaintext_nbr);


  uint64_t repliesAmount = computeReplySizeInChunks(cryptoMethod->getPublicParameters().getCiphertextBitsize() / CHAR_BIT);
  repliesArray = (char**)calloc(repliesAmount,sizeof(char*)); 
	repliesIndex = 0;

  double start = omp_get_wtime();
  generateReply();
  double result = omp_get_wtime() - start;

  freeQuery();
  freeInputData();
  freeResult();
  delete dbhandler;
  return result;
}



double PIRReplyGeneratorNFL::precomputationSimulation(const PIRParameters& pir_params, uint64_t plaintext_nbr)
{
  NFLlib *nflptr = &(cryptoMethod->getnflInstance());
  setPirParams((PIRParameters&)pir_params);
  initQueriesBuffer();
  pushFakeQuery();
  importFakeData(plaintext_nbr);

  uint64_t files_nbr = 1;
  for (unsigned int i = 0 ; i < pir_params.d ; i++) files_nbr *= pir_params.n[i];

  double start = omp_get_wtime();
  for (unsigned int i = 0 ; i < files_nbr ; i++)
  {
    {
      poly64 *tmp;
      tmp= cryptoMethod->deserializeDataNFL((unsigned char**)(input_data[i].p), (uint64_t) plaintext_nbr, cryptoMethod->getPublicParameters().getCiphertextBitsize()/2 , input_data[i].nbPolys);
	    free(tmp[0]);	
    }
  }
  double result = omp_get_wtime() - start;
  std::cout << "PIRReplyGeneratorNFL: Deserialize took " << result << " (omp)seconds" << std::endl;
  freeQuery();
  freeInputData();
  freeResult();
  delete dbhandler;
  return result;
}



/**
 *	Multiply each query parts by each files and sum the result.
 *	Params :
 *	- lwe_queries* : the query ;
 *	- lwe_in_data* : data to be processed.
 *	- int begin_data : index where begins the data absorption
 *	- int lvl : recursion level ;
 *	- lwe_cipher* result : Array to store the result.
 **/
void PIRReplyGeneratorNFL::generateReply(	lwe_query *queries_,
lwe_in_data* data,
unsigned int lvl,
lwe_cipher* result)
{
#ifdef SHOUP
	lwe_query **queries=(lwe_query**)queries_;
#else
	lwe_query *queries=queries_;
#endif	
	unsigned int query_size = pirParam.n[lvl];

#ifdef PERFTIMERS
  bool wasVerbose = false;
  double vtstart = omp_get_wtime();
#endif

  // In order to parallelize we must ensure replies are somehow ordered 
	// (see comment at the end of PIRReplyExtraction)
	//#pragma omp parallel for firstprivate(result,data, lvl, queries)
  for (unsigned int offset = 0; offset < query_size; offset += 200)
  {
#ifdef MULTI_THREAD
   # pragma omp parallel for
#endif
    for (unsigned int current_poly=0 ; current_poly < currentMaxNbPolys ; current_poly++)
	  { 
		  posix_memalign((void**) &(result[current_poly].a), 32, 
          2*cryptoMethod->getpolyDegree()*cryptoMethod->getnbModuli()*sizeof(uint64_t));
		  memset(result[current_poly].a,0,
          2*cryptoMethod->getpolyDegree()*cryptoMethod->getnbModuli()*sizeof(uint64_t));
		  result[current_poly].b = (uint64_t *) result[current_poly].a +
        cryptoMethod->getpolyDegree()*cryptoMethod->getnbModuli();
      for (unsigned int query_index = 0, ggg=0; query_index < query_size && ggg < 200 ; 
          query_index++, ggg++)
		  {
#ifdef SHOUP
#ifdef CRYPTO_DEBUG  
		    if(current_poly==0)
        {
			    std::cout<<"Query poped.a  ";NFLTools::print_poly64hex(queries[0][query_index].a,4);	
			    if (lwe)	
          {
            std::cout<<"Query poped.b  ";NFLTools::print_poly64hex(queries[0][query_index].b,4);
          }
          std::cout<<"Query poped.a' ";NFLTools::print_poly64hex(queries[1][query_index].a,4);	
          if (lwe) 
          {
            std::cout<<"Query poped.b' ";NFLTools::print_poly64hex(queries[1][query_index].b,4);
          }
		    }
#endif
			  cryptoMethod->mulandadd(result[current_poly], data[query_index], queries[0][query_index],
          queries[1][query_index], current_poly, lvl);
#else
			  cryptoMethod->mulandadd(result[current_poly], data[query_index], queries[query_index], 
          current_poly, lvl);
#endif				
		  }
      if ( lvl == pirParam.d-1 && offset + 200 >= query_size) 
      {
        // Watchout lwe_cipher.a and .b need to be allocated contiguously
        repliesArray[repliesIndex+current_poly] = (char*)result[current_poly].a; 
      }

#ifdef PERFTIMERS
      // Give some feedback if it takes too long
      double vtstop = omp_get_wtime();
      if (vtstop - vtstart > 1) 
      {
        vtstart = vtstop;
        if(currentMaxNbPolys != 1) std::cout <<"PIRReplyGeneratorNFL: Dealt with chunk " << 
          current_poly+1 << "/" << currentMaxNbPolys << "\r" << std::flush;
        wasVerbose = true;
      }
#endif

    }
  }

#ifdef PERFTIMERS
  if (wasVerbose) std::cout <<"                                                     \r" << std::flush;
#endif

}


// New version using the multiple buffer serialize function
lwe_in_data* PIRReplyGeneratorNFL::fromResulttoInData(lwe_cipher** inter_reply, uint64_t reply_elt_nbr, unsigned int reply_rec_lvl)
{
    uint64_t in_data2b_bytes = cryptoMethod->getPublicParameters().getAbsorptionBitsize()/8;
    uint64_t in_data2b_nbr_polys = ceil((double(currentMaxNbPolys * cryptoMethod->getPublicParameters().getCiphertextBitsize())/8.)/double(in_data2b_bytes));

    lwe_in_data *in_data2b  = new lwe_in_data[reply_elt_nbr]();
    uint64_t **bufferOfBuffers = (uint64_t **) calloc(currentMaxNbPolys,sizeof(uint64_t*));

    //For each element in the reply
    for (uint64_t i = 0 ; i < reply_elt_nbr ; i++)
    {
      //Build the buffer of buffers
      for (uint64_t j = 0 ; j < currentMaxNbPolys ; j++)
      {
        bufferOfBuffers[j]=inter_reply[i][j].a;
      }

	    // Ciphertexts can be serialized in a single block as a,b are allocatted contiguously
	    in_data2b[i].p = cryptoMethod->deserializeDataNFL((unsigned char**)bufferOfBuffers, 
														currentMaxNbPolys, 
														cryptoMethod->getPublicParameters().getCiphertextBitsize(), 
														in_data2b[i].nbPolys);
      //delete[] inter_reply[i]; free in generateReplyGeneric
    }
    //delete[] inter_reply;
    free(bufferOfBuffers);

    currentMaxNbPolys = in_data2b_nbr_polys; 
    return in_data2b;
}

//// Original function
//lwe_in_data* PIRReplyGeneratorNFL::fromResulttoInData(lwe_cipher** inter_reply, uint64_t reply_elt_nbr, unsigned int reply_rec_lvl)
//{
//    uint64_t in_data2b_bytes = cryptoMethod->getPublicParameters().getAbsorptionBitsize()/8;
//    uint64_t in_data2b_polys_per_reply_poly = ceil((double)(cryptoMethod->getPublicParameters().getCiphertextBitsize()/8)/in_data2b_bytes);
//    uint64_t in_data2b_nbr_polys = currentMaxNbPolys * in_data2b_polys_per_reply_poly;
//
//    lwe_in_data *in_data2b  = new lwe_in_data[reply_elt_nbr]();
//    lwe_in_data tmp_in_data;
//
//    //For each element in the reply
//    for (uint64_t i = 0 ; i < reply_elt_nbr ; i++)
//    {
//    	in_data2b[i].p = (poly64 *) malloc(in_data2b_nbr_polys*sizeof(poly64));
//	    in_data2b[i].nbPolys = 0;
//
//     	// For each polynomial in a reply element
//     	for (uint64_t j = 0 ; j < currentMaxNbPolys ; j++)
//     	{
//	      // Ciphertexts can be serialized in a single block as a,b are allocatted contiguously
//	      tmp_in_data.p = cryptoMethod->deserializeDataNFL((unsigned char*)inter_reply[i][j].a, cryptoMethod->getPublicParameters().getCiphertextBitsize(), tmp_in_data.nbPolys);
//	      for (uint64_t k = 0 ; k < in_data2b_polys_per_reply_poly; k++)
//   	    {
//		      in_data2b[i].p[k + j * in_data2b_polys_per_reply_poly] = tmp_in_data.p[k];
//	      }
//	      in_data2b[i].nbPolys += in_data2b_polys_per_reply_poly;
//	    }
//      delete[] inter_reply[i];
//    }
//    delete[] inter_reply;
//
//    currentMaxNbPolys = in_data2b_nbr_polys; 
//    return in_data2b;
//}

/**
 * Compute Reply Size une chunks. 
 * WARNING blocking function.
 **/
unsigned long PIRReplyGeneratorNFL::computeReplySizeInChunks(unsigned long int maxFileBytesize)
{
  using namespace GlobalConstant;
  
  unsigned int out = ceil((double)maxFileBytesize*kBitsPerByte*pirParam.alpha/cryptoMethod->getPublicParameters().getAbsorptionBitsize(0));

  for (unsigned int i = 1; i < pirParam.d; i++) {
    out = ceil(out  * double(cryptoMethod->getPublicParameters().getCiphBitsizeFromRecLvl(i)/kBitsPerByte) / double(cryptoMethod->getPublicParameters().getAbsorptionBitsize(i) / kBitsPerByte));
      
  }

  return out;
}

/**
 *	Overloaded fonction from  GenericPIRReplyGenerator.
 *	Initalise queriesBuf.
 **/

void PIRReplyGeneratorNFL::initQueriesBuffer() {
	const unsigned int nbQueriesBuf=pirParam.d;;
#ifdef SHOUP
	queriesBuf = new lwe_query**[nbQueriesBuf]();
	for (unsigned int i = 0 ; i < nbQueriesBuf ; i++)
	{
		queriesBuf[i] = new lwe_query*[2];
		queriesBuf[i][0] = new lwe_query[pirParam.n[i]]();
		queriesBuf[i][1] = new lwe_query[pirParam.n[i]]();
	}
#else
	queriesBuf = new lwe_query*[nbQueriesBuf]();
	for (unsigned int i = 0 ; i < nbQueriesBuf ; i++)
	{
		queriesBuf[i] = new lwe_query[pirParam.n[i]]();
	}
#endif
#ifdef DEBUG
	std::cout<<"Created a queriesBuf for "<<nbQueriesBuf<<" queries"<<std::endl;	
#endif

}

void PIRReplyGeneratorNFL::pushFakeQuery()
{
  char* query_element = cryptoMethod->encrypt(0, 1); 

  for (unsigned int dim  = 0 ; dim < pirParam.d ; dim++) {
    for(unsigned int j = 0 ; j < pirParam.n[dim] ; j++) {
      pushQuery(query_element, cryptoMethod->getPublicParameters().getCiphertextBitsize()/8, dim, j); 
    }
  }

  free(query_element);
}

void PIRReplyGeneratorNFL::pushQuery(char* rawQuery, unsigned int size, int dim, int nbr)
{
  unsigned int polyDegree = cryptoMethod->getpolyDegree();
  unsigned int nbModuli = cryptoMethod->getnbModuli();
  // Trick, we get both a and b at the same time, b needs to be set afterwards
  uint64_t *a,*b;
  a = (poly64) calloc(size, 1);
  memcpy(a,rawQuery,size);
  if (lwe) b = a+nbModuli*polyDegree;
#ifdef CRYPTO_DEBUG
	std::cout<<"\nQuery received.a ";NFLTools::print_poly64(a,4);	
	if (lwe) {std::cout<<"Query received.b ";NFLTools::print_poly64hex(b,4);}
#endif
#ifdef SHOUP

  uint64_t *ap,*bp;
  ap = (poly64) calloc(size, 1);
  if (lwe) bp = ap+nbModuli*polyDegree;
  
  for (unsigned int cm = 0 ; cm < nbModuli ; cm++)
  { 
    for (unsigned i = 0 ; i < polyDegree ;i++) 
    {
  	  ap[i+cm*polyDegree] = ((uint128_t) a[i+cm*polyDegree] << 64) / cryptoMethod->getmoduli()[cm];
  	  if (lwe) bp[i+cm*polyDegree] = ((uint128_t) b[i+cm*polyDegree] << 64) / cryptoMethod->getmoduli()[cm];
    }  
  }
  queriesBuf[dim][0][nbr].a = a;
  queriesBuf[dim][0][nbr].b = b;
  queriesBuf[dim][1][nbr].a = ap;
  queriesBuf[dim][1][nbr].b = bp;

#ifdef CRYPTO_DEBUG
	std::cout << "Query NFL pushed.a' "; NFLTools::print_poly64hex(queriesBuf[dim][1][nbr].a,4);	
	if (lwe) { std::cout << "Query NFL pushed.b' "; NFLTools::print_poly64hex(queriesBuf[dim][1][nbr].b,4);}
#endif
#else
  queriesBuf[dim][nbr].a = a;
  queriesBuf[dim][nbr].b = b;
#endif
}

size_t PIRReplyGeneratorNFL::getTotalSystemMemory()
{
#ifdef OSX
  int m[2];
  m[0] = CTL_HW;
  m[1] = HW_MEMSIZE;

  int64_t size = 0;               
  size_t len = sizeof( size );

  sysctl( m, 2, &size, &len, NULL, 0 );
    
  return (size_t)size;

#else
  long pages = /*get_phys_pages();*/sysconf(_SC_PHYS_PAGES);
  long page_size = /*getpagesize();*/sysconf(_SC_PAGE_SIZE);
  return pages * page_size;
#endif
}

PIRReplyGeneratorNFL::~PIRReplyGeneratorNFL()
{
  freeResult();
}


void PIRReplyGeneratorNFL::setPirParams(PIRParameters& param)
{
  pirParam = param;
  cryptoMethod->setandgetAbsBitPerCiphertext(1, pirParam.n[0]);
}


void PIRReplyGeneratorNFL::setCryptoMethod(CryptographicSystem* cm)
{
  //cryptoMethod = (NFLLWE*) cm;
  cryptoMethod = (LatticesBasedCryptosystem*) cm;
  lwe = (cryptoMethod->toString() == "LWE") ? true : false; 
}

void PIRReplyGeneratorNFL::freeInputData()
{
  uint64_t theoretical_files_nbr = 1;
	for (unsigned int i = 0 ; i < pirParam.d ; i++) theoretical_files_nbr *= pirParam.n[i];

  for (unsigned int i = 0 ; i < theoretical_files_nbr ; i++){
#ifdef DEBUG
  printf( "PIRReplyGeneratorNFL: freeing input_data[%d]\n",i);
#endif
    free(input_data[i].p[0]);
    free(input_data[i].p);
  }
  delete[] input_data;

#ifdef DEBUG
  printf( "PIRReplyGeneratorNFL: input_data freed\n");
#endif
}

void PIRReplyGeneratorNFL::freeQuery()
{
  for (unsigned int i = 0; i < pirParam.d; i++)
  {
#ifdef SHOUP
      for (unsigned int j = 0 ; j < pirParam.n[i] ; j++) {
		  free(queriesBuf[i][0][j].a); //only free a because a and b and contingus, see pushQuery
		  free(queriesBuf[i][1][j].a); //only free a because a and b and contingus, see pushQuery
	  }
      delete[] queriesBuf[i][0]; //allocated in intQueriesBuf with new.
      delete[] queriesBuf[i][1]; //allocated in intQueriesBuf with new.
      delete[] queriesBuf[i];
#else
	  for (unsigned int j = 0 ; j < pirParam.n[i] ; j++) free(queriesBuf[i][j].a); //only free a because a and b and contingus, see pushQuery
      delete[] queriesBuf[i]; //allocated in intQueriesBuf with new.
#endif
  }
  delete[] queriesBuf;//allocated in intQueriesBuf with new.
#ifdef DEBUG
  printf( "queriesBuf freed\n");
#endif
  
}

void PIRReplyGeneratorNFL::freeResult()
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
