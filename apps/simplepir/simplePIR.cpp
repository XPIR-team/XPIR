#include "libpir.hpp"
#include "apps/server/DBDirectoryProcessor.hpp"

bool run(DBHandler *db, uint64_t chosen_element, PIRParameters params){


  /******************************************************************************
  * PIR and Crypto Setup (must be done by both the client and the server)      
  * In a real application the client and server must agree on the parameters
  * For example the client chooses and sends them to the server (or inversely)
  ******************************************************************************/
	
	HomomorphicCrypto *crypto = HomomorphicCryptoFactory::getCryptoMethod(params.crypto_params);

  // Absorption capacity of an LWE encryption scheme depends on the number of sums that are going
  // to be done in the PIR protocol, it must therefore be initialized
	// Warning here we suppose the biggest dimension is in d[0] 
  // otherwise absorbtion needs to be computed accordingly
  crypto->setandgetAbsBitPerCiphertext(params.n[0]);


  /******************************************************************************
  * Query generation phase (client-side)                                       
  ******************************************************************************/
	
  // Create the query generator object
	PIRQueryGenerator q_generator(params,*crypto);
  std::cout << "SimplePIR: Generating query ..." << std::endl;
  // Generate a query to get the FOURTH element in the database (indexes begin at 0)
  // Warning : if we had set params.alpha=2 elements would be aggregated 2 by 2 and 
  // generatequery would only accept as input 0 (the two first elements) or 1 (the other two)
	q_generator.generateQuery(chosen_element);
  std::cout << "SimplePIR: Query generated" << std::endl;
	

  /******************************************************************************
  * Reply generation phase (server-side)
  ******************************************************************************/
	
  // Create the reply generator object
  // We could have also defined PIRReplyGenerator *r_generator(params,*crypto,db);
  // But we prefer a pointer to show (below) how to use multiple generators for a given db
  PIRReplyGenerator *r_generator = new PIRReplyGenerator(params,*crypto,db);
  r_generator->setPirParams(params);

  // In a real application the client would pop the queries from q with popQuery and 
  // send them through the network and the server would receive and push them into s 
  // using pushQuery
  char* query_element;
  while (q_generator.popQuery(&query_element))
  {
    r_generator->pushQuery(query_element);
  }
 
	// Import database
  // This could have been done on the "Database setup" phase if:
  //  - the contents are static
  //  - AND the imported database fits in RAM
  //  - AND the server knows in advance the PIR and crypto parameters (e.g. chosen by him)
  std::cout << "SimplePIR: Importing database ..." << std::endl;
  // Warning aggregation is dealt with internally the bytes_per_db_element parameter here
  // is to be given WITHOUT multiplying it by params.alpha
	imported_database* imported_db = r_generator->importData(/* uint64_t offset*/ 0, /*uint64_t
    bytes_per_db_element */ db->getmaxFileBytesize());
  std::cout << "SimplePIR: Database imported" << std::endl;

	// Once the query is known and the database imported launch the reply generation
  std::cout << "SimplePIR: Generating reply ..." << std::endl;
	double start = omp_get_wtime();
	r_generator->generateReply(imported_db);
	double end = omp_get_wtime();
  std::cout << "SimplePIR: Reply generated in " << end-start << " seconds" << std::endl;

  /********************************************************************************
   * Advanced example: uncomment it to test
   * The object imported_db is separated from r_generator in purpose
   * Here is an example on how to use the same imported_db for multiple queries
   * DO NOT try to use the same reply generator more than once, this causes issues
   * ******************************************************************************/

  // Generate 3 replies from 3 queries
  for (int i = 0 ; i < 3 ; i++){

    // Pop (and drop for this simple example) the generated reply
    char* reply_element_tmp;
    while (r_generator->popReply(&reply_element_tmp)){
      free(reply_element_tmp);
    }

    // Free memory and the r_generator object
    r_generator->freeQueries();
    //r_generator->freeResult();
    //delete r_generator;
	
    // Create and initialize a new generator object
    //r_generator = new PIRReplyGenerator(params,*crypto,db);
    //r_generator->setPirParams(params);
	
    // Generate a new query
  	q_generator.generateQuery(chosen_element);

    // Push it to the reply generator
    while (q_generator.popQuery(&query_element))
    {
      r_generator->pushQuery(query_element);
    }

    // Generate again the reply
	  r_generator->generateReply(imported_db);
  }


  /******************************************************************************
  * Reply extraction phase (client-side)
  ******************************************************************************/
	
	PIRReplyExtraction *r_extractor=new PIRReplyExtraction(params,*crypto);

  // In a real application the server would pop the replies from s with popReply and 
  // send them through the network together with nbRepliesGenerated and aggregated_maxFileSize 
  // and the client would receive the replies and push them into r using pushEncryptedReply
  std::cout << "SimplePIR: "<< r_generator->getnbRepliesGenerated()<< " Replies generated " << std::endl;

  uint64_t clientside_maxFileBytesize = db->getmaxFileBytesize();
  char* reply_element;
  while (r_generator->popReply(&reply_element))
  {
    r_extractor->pushEncryptedReply(reply_element);
  }

  std::cout << "SimplePIR: Extracting reply ..." << std::endl;
	r_extractor->extractReply(clientside_maxFileBytesize);
  std::cout << "SimplePIR: Reply extracted" << std::endl;

  // In a real application instead of writing to a buffer we could write to an output file
  char *outptr, *result, *tmp;
  outptr = result = (char*)calloc(r_extractor->getnbPlaintextReplies(clientside_maxFileBytesize)*r_extractor->getPlaintextReplyBytesize(), sizeof(char));
  while (r_extractor->popPlaintextResult(&tmp)) 
  {
    memcpy(outptr, tmp, r_extractor->getPlaintextReplyBytesize()); 
    outptr+=r_extractor->getPlaintextReplyBytesize();
    free(tmp);
  }
  // Result is in ... result
  

  /******************************************************************************
  * Test correctness
  ******************************************************************************/
	char *db_element = (char*)calloc(clientside_maxFileBytesize*params.alpha, sizeof(char));
  bool fail = false;
	db->readAggregatedStream(chosen_element, params.alpha, 0, clientside_maxFileBytesize, db_element);
  if (memcmp(result, db_element, clientside_maxFileBytesize*params.alpha))
  {
    std::cout << "SimplePIR: Test failed, the retrieved element is not correct" << std::endl;
    fail = true;
  }
  else
  {
    std::cout << "SimplePIR: Test succeeded !!!!!!!!!!!!!!!!!!!!!!!!"  << std::endl<< std::endl;
    fail = false;
  }

  /******************************************************************************
  * Cleanup
  ******************************************************************************/
  
  delete imported_db;
  r_generator->freeQueries();
  delete r_generator;
  free(result);
  free(db_element);
  
	return fail;
	
}


int main(int argc, char * argv[]) {

  uint64_t database_size, nb_files, chosen_element, maxFileBytesize;
	PIRParameters params; 
  bool tests_failed = false;
  
/******************************************************************************
  * Database setup (server-side)                              
  ******************************************************************************/


  // To Create the database generator object
  // it can be a DBGenerator that simulate nb_files files of size streamBytesize
  // database_size = 1ULL<<25; nb_files = 4; maxFileBytesize = database_size/nb_files;
  // DBGenerator db(nb_files, maxFileBytesize, /*bool silent*/ false); 
  // 
  // OR it can be a DBDirectoryProcessor that reads a real file in the ./db directory 
  // and splits it into nb_files virtual files
  // nb_files = 4;
  // DBDirectoryProcessor db(nb_files);
  // database_size = db.getDBSizeinbits();maxFileBytesize = database_size/nb_files;
  //
  // OR it can be a DBDirectoryProcessor that reads the real files in the ./db directory
  // DBDirectoryProcessor db;
  // nb_files=db.getNbStream();database_size = db.getDBSizeinbits();
  // maxFileBytesize = database_size/nb_files;

  // Simple test
  std::cout << "======================================================================" << std::endl;
  std::cout << "Test 1/7: database_size = 1ULL<<30; nb_files = 20;" << std::endl;
  std::cout << "params.alpha = 1; params.d = 1; crypto_params = LWE:80:2048:120;" << std::endl; 
  std::cout << "======================================================================" << std::endl;
  database_size = 1ULL<<20; nb_files = 20; maxFileBytesize = database_size/nb_files;
  DBGenerator db(nb_files, maxFileBytesize, /*bool silent*/ false); 
  chosen_element = 3;
  params.alpha = 1; params.d = 1; params.n[0] = nb_files; 
  // The crypto parameters can be set to other values
  // You can get a list of all available cryptographic parameters with this function call 
  // HomomorphicCryptoFactory::printAllCryptoParams();
  params.crypto_params = "LWE:80:2048:120"; 
  tests_failed |= run(&db, chosen_element, params);
  
  // Test with aggregation
  // WARNING we must provide the representation of the database GIVEN recursion and aggregation
  // as here we have 100 elements and aggregate them in a unique group we have params.n[0]=1
  std::cout << "======================================================================" << std::endl;
  std::cout << "Test 2/7: database_size = 1ULL<<25; nb_files = 100;" << std::endl;
  std::cout << "params.alpha = 100; params.d = 1; crypto_params = LWE:80:2048:120;" << std::endl; 
  std::cout << "======================================================================" << std::endl;
  database_size = 1ULL<<25; nb_files = 100; maxFileBytesize = database_size/nb_files;
  DBGenerator db2(nb_files, maxFileBytesize, /*bool silent*/ false); 
  chosen_element = 0;
  params.alpha = 100; params.d = 1; params.n[0] = 1; 
  params.crypto_params = "LWE:80:2048:120";
  tests_failed |= run(&db2, chosen_element, params);
  
  // Test with recursion 2
  std::cout << "======================================================================" << std::endl;
  std::cout << "Test 3/7: database_size = 1ULL<<25; nb_files = 100;" << std::endl;
  std::cout << "params.alpha = 1; params.d = 2; crypto_params = LWE:80:2048:120;" << std::endl; 
  std::cout << "======================================================================" << std::endl;
  database_size = 1ULL<<25; nb_files = 100; maxFileBytesize = database_size/nb_files;
  DBGenerator db3(nb_files, maxFileBytesize, /*bool silent*/ false); 
  chosen_element = 3;
  params.alpha = 1; params.d = 2; params.n[0] = 50; params.n[1] = 2; 
  params.crypto_params = "LWE:80:2048:120";
  tests_failed |= run(&db3, chosen_element, params);
  
  // Test with recursion 2 and aggregation
  std::cout << "======================================================================" << std::endl;
  std::cout << "Test 4/7: database_size = 1ULL<<25; nb_files = 100;" << std::endl;
  std::cout << "params.alpha = 2; params.d = 2; crypto_params = LWE:80:2048:120;" << std::endl; 
  std::cout << "======================================================================" << std::endl;
  database_size = 1ULL<<25; nb_files = 100; maxFileBytesize = database_size/nb_files;
  DBGenerator db4(nb_files, maxFileBytesize, /*bool silent*/ false); 
  chosen_element = 3;
  params.alpha = 2; params.d = 2; params.n[0] = 25; params.n[1] = 2; 
  params.crypto_params = "LWE:80:2048:120";
  tests_failed |= run(&db4, chosen_element, params);
  
  // Test with recursion 3
  std::cout << "======================================================================" << std::endl;
  std::cout << "Test 5/7: database_size = 1ULL<<25; nb_files = 100;" << std::endl;
  std::cout << "params.alpha = 1; params.d = 3; crypto_params = LWE:80:2048:120;" << std::endl; 
  std::cout << "======================================================================" << std::endl;
  database_size = 1ULL<<25; nb_files = 100; maxFileBytesize = database_size/nb_files;
  DBGenerator db5(nb_files, maxFileBytesize, /*bool silent*/ false); 
  chosen_element = 3;
  params.alpha = 1; params.d = 3; params.n[0] = 5; params.n[1] = 5; params.n[2] = 4; 
  params.crypto_params = "LWE:80:2048:120";
  tests_failed |= run(&db5, chosen_element, params);
  
  // Test with a DBDirectoryProcessor splitting a big real file
  std::cout << "======================================================================" << std::endl;
  std::cout << "Test 6/7: DBDirectoryProcessor with split; database_size = 1ULL<<25; nb_files = 4;" << std::endl;
  std::cout << "params.alpha = 1; params.d = 1; crypto_params = LWE:80:2048:120;" << std::endl; 
  std::cout << "======================================================================" << std::endl;
  database_size = 1ULL<<25; nb_files = 4; maxFileBytesize = database_size/nb_files;
  DBDirectoryProcessor db6(/*split the first file in*/ nb_files /*files*/);
  if (db6.getErrorStatus()==true){
    std::cout << "SimplePIR : Error with db directory skipping test ..." << std::endl << std::endl;
  } else {
    chosen_element = 3;
    params.alpha = 1; params.d = 1; params.n[0] = nb_files; 
    params.crypto_params = "LWE:80:2048:120";
    tests_failed |= run(&db6, chosen_element, params);
  }
  
  // Test with a DBDirectoryProcessor reading real files
  std::cout << "======================================================================" << std::endl;
  std::cout << "Test 7/7: DBDirectoryProcessor without split;" << std::endl;
  std::cout << "params.alpha = 1; params.d = 1; crypto_params = LWE:80:2048:120;" << std::endl; 
  std::cout << "======================================================================" << std::endl;
  DBDirectoryProcessor db7;
  if (db6.getErrorStatus()==true){
    std::cout << "SimplePIR : Error with db directory skipping test ..." << std::endl << std::endl;
  } else {
    database_size = db7.getDBSizeBits()/8; nb_files = db7.getNbStream(); 
    maxFileBytesize = database_size/nb_files;
    chosen_element = 0;
    params.alpha = 1; params.d = 1; params.n[0] = nb_files; 
    params.crypto_params = "LWE:80:2048:120";
    tests_failed |= run(&db7, chosen_element, params);
  }

  if (tests_failed) 
  {
    std::cout << "WARNING : at least one tests failed" << std::endl;
    return 1;
  }
  else
  {
    std::cout << "All tests succeeded" << std::endl;
    return 0;
  }
}



