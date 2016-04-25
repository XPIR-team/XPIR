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

#include "PIRSession.hpp"
#include "pir/replyGenerator/PIROptimizer.hpp"

#define NDSS_UPLOAD_SPEED 100000000UL
PIRSession::pointer PIRSession::create(boost::asio::io_service& ios)
{
  return PIRSession::pointer(new PIRSession(ios));
}

PIRSession::PIRSession(boost::asio::io_service& ios) :
  sessionSocket(ios),
  handmadeExceptionRaised(false),
  finished(false),
  cryptoMethod(NULL),
  generator(NULL),
  no_pipeline_mode(false)
{
}

void PIRSession::setDBHandler(DBHandler *db)
{
  dbhandler = db;
}

/**
 * Do the PIR protocol
 **/
bool PIRSession::start(session_option_t session_option)
{
  
  uint64_t nbFiles=dbhandler->getNbStream();
  maxFileBytesize = dbhandler->getmaxFileBytesize();
  short exchange_method = (session_option.driven_mode) ? CLIENT_DRIVEN : SERVER_DRIVEN;
  bool rcv_paramsandkey = true;

  // Deal with Optimizer dry-run queries if requested
  if(!rcvIsClient())
  {
    std::cout << "PIRSession: Incoming dry-run optimizer request, dealing with it ..." << std::endl;
    PIROptimizer optimizer(dbhandler);
    optimizer.prepareOptimData();
    optimizer.controlAndCommand(sessionSocket);
    std::cout << "PIRSession: Session finished" << std::endl << std::endl;
    finished = true;
    return 1; // Did deal with a dry-run query
  }

  sendCatalog();
  sendPIRParamsExchangeMethod(exchange_method);

  if(session_option.driven_mode) {
    PIROptimizer optimizer(dbhandler);
    optimizer.prepareOptimData();
    optimizer.controlAndCommand(sessionSocket);
    rcvCryptoParams(rcv_paramsandkey);
    rcvPirParams();
  }
  else //not driven
  {
    std::vector<std::string> fields;
    boost::algorithm::split(fields, pirParam.crypto_params, boost::algorithm::is_any_of(GlobalConstant::kDelim));
    cryptoMethod = HomomorphicCryptoFactory_internal::getCrypto(fields.at(0));
    cryptoMethod->setNewParameters(pirParam.crypto_params);

    generator = PIRReplyGeneratorFactory::getPIRReplyGenerator(fields.at(0), pirParam,dbhandler);
    generator->setCryptoMethod(cryptoMethod);
    sendCryptoParams();
    rcvCryptoParams(!rcv_paramsandkey);
    sendPirParams();
  }

  // If one of the functions above generates an error, handmadeExceptionRaised is set 
  if (!handmadeExceptionRaised)
  {
    // This is just a download thread. Reply generation is unlocked (by a mutex)
    // when this thread finishes.
    startProcessQuery(); 

    // Import the database
    // Start reply generation when mutex unlocked
    // Start a thread which uploads the reply as it is generated
    startProcessResult(session_option);
  }

  // Wait for child threads
  if (upThread.joinable())  upThread.join();
  if (downThread.joinable()) downThread.join();

  std::cout << "PIRSession: Session finished" << std::endl << std::endl;
  // Note that we have finished so that the PIRServer can do garbage collection
  finished = true;
  if (generator != NULL) delete generator;
  return 0; // Did deal with a client query
}


/** 
 * Getter for other classes (such as PIRServer) 
 **/
tcp::socket& PIRSession::getSessionSocket()
{
  return sessionSocket;
}

/**
 * Send the catalog to the client as a string. 
 * Format : file_list_size \n filename1 \n filesize1 \n filename2 \n ... filesizeN \n
 *          if SEND_CATALOG is defined and catalog size is less than 1000, the full catalog is sent
 *			otherwise only the size of the catalog is sent
 **/
void PIRSession::sendCatalog() 
{
	string buf;
#ifdef SEND_CATALOG
    if(dbhandler->getNbStream()>1000) {
    	buf=dbhandler->getCatalog(false);
    } else {
      	buf=dbhandler->getCatalog(true);
    }
#else
    buf = dbhandler->getCatalog(false);
#endif
  	// Send the buffer
  	const boost::uint64_t size = buf.size();
	  
  	try
  	{
    	if (write(sessionSocket, boost::asio::buffer(&size, sizeof(size))) <= 0)
      	 	exitWithErrorMessage(__FUNCTION__,"Error sending catalog size");
    	if (write(sessionSocket, boost::asio::buffer(buf.c_str(), size)) < size) 
      	  	exitWithErrorMessage(__FUNCTION__,"Error sending catalog");
  	} catch (std::exception const& ex) {
    	exitWithErrorMessage(__FUNCTION__,"Error sending catalog: " + string(ex.what()));
  	}
  	writeWarningMessage(__FUNCTION__ , "done.");
}

void PIRSession::sendCryptoParams() 
{
  std::cout << "PIRSession: Mandatory crypto params sent to the client are " << pirParam.crypto_params << std::endl;

  try 
  {
    int crypto_params_size = pirParam.crypto_params.size();

    write(sessionSocket, boost::asio::buffer(&crypto_params_size, sizeof(crypto_params_size)));

    write(sessionSocket, boost::asio::buffer(pirParam.crypto_params));
  } 
  catch(std::exception const& ex) 
  {
    exitWithErrorMessage(__FUNCTION__, string(ex.what()));
  }
}

bool PIRSession::rcvIsClient() 
{
  try{
    int is_client;
    while ( read(sessionSocket, boost::asio::buffer(&(is_client), sizeof(int))) == 0)
      boost::this_thread::yield();
      //exitWithErrorMessage(__FUNCTION__, "No client or optim choice recieved");

    return is_client == 1; 
  }catch (std::exception const& ex)
  {
    exitWithErrorMessage(__FUNCTION__, string(ex.what()));
  }
}

/**
 * Receive client's cryptographic parameters as:
 * - An int describing a string size
 * - A string of the given size containing the public parameters
 * - A second int describing a byte size
 * - A buffer of bytes of the given size with key material
 **/
void PIRSession::rcvCryptoParams(bool paramsandkey)
{
  unsigned int size = 0;

  try{
    std::vector<std::string> fields;
    
    if(paramsandkey == true)
    { 
      // First get the int and allocate space for the string (plus the null caracter) 
      read(sessionSocket, boost::asio::buffer(&size, sizeof(int)));
      char params_buf[size + 1];	

      // Get the string in the buffer and add a null character
      read(sessionSocket, boost::asio::buffer(params_buf, size));
      params_buf[size] = '\0';

      cout << "PIRSession: Received crypto parameters " << params_buf << " processing them ..." << endl;
#ifdef DEBUG
      cout << "PIRSession: Parameter string size is " << size << endl;
#endif

      // Extract cryptosystem's name
      string crypto_system_desc(params_buf);
      boost::algorithm::split(fields, crypto_system_desc, boost::algorithm::is_any_of(":"));

      // Create cryptosystem using a factory and the extracted name
      cryptoMethod = HomomorphicCryptoFactory_internal::getCrypto(fields[0]); 
    
      // Set cryptosystem with received parameters and key material
      cryptoMethod->setNewParameters(params_buf);

      // Create the PIR reply generator object using a factory to have 
      // the correct object given the cryptosystem used (for optimization 
      // reply generation is cryptosystem dependent)
      generator = PIRReplyGeneratorFactory::getPIRReplyGenerator(fields[0], pirParam, dbhandler);
      if (generator == NULL) 
      {
        std::cout << "PIRSession: CRITICAL no reply generator found, exiting session" << std::endl;
        pthread_exit(this);
      }
      else
      {
        generator->setCryptoMethod(cryptoMethod);
      } 
    }
    
    // Use again size to describe the key size 
    if (read(sessionSocket, boost::asio::buffer(&size ,sizeof(size))) <= 0)
      exitWithErrorMessage(__FUNCTION__,"No key received, abort.");

#ifdef DEBUG
    cout << "PIRSession: Size of received key material is " << size << endl;
#endif
    // Get the key material only if there is some  
    if(size > 0)
    {
      // This time we don't use a string so no need for an extra character
      char buf[size];

      if (read(sessionSocket, boost::asio::buffer(buf, size)) < size)
        exitWithErrorMessage(__FUNCTION__,"No parameters received, abort.");
    
      cryptoMethod->getPublicParameters().setModulus(buf);
    }

  }catch(std::exception const& ex)
  {
    exitWithErrorMessage(__FUNCTION__, string(ex.what()));
  }

  std::cout << "PIRSession: Finished processing crypto parameters" << std::endl;
}

void PIRSession::sendPIRParamsExchangeMethod(short exchange_method)
{
  if (exchange_method == CLIENT_DRIVEN)
  {
    cout << "PIRSession: Notifying the client this is a client-driven session" << endl;
  }
  else
  {
    cout << "PIRSession: Notifying the client this is a server-driven session" << endl;
  }
  try
  {
    write(sessionSocket, boost::asio::buffer(&exchange_method, sizeof(exchange_method)));
  }catch(std::exception const& ex)
  {
    exitWithErrorMessage(__FUNCTION__, string(ex.what()));
  }
}

/**
 * Receive client's PIR parameters.
 **/
void PIRSession::rcvPirParams() 
{
  try{
    // First we get an int with the recursion level
    if ( read(sessionSocket, boost::asio::buffer(&(pirParam.d), sizeof(int))) <= 0) 
      exitWithErrorMessage(__FUNCTION__, "No pir param recieved");

    // Then we get an int with the aggregation
    if ( read(sessionSocket, boost::asio::buffer(&(pirParam.alpha), sizeof(int))) <= 0) 
      exitWithErrorMessage(__FUNCTION__, "No pir param recieved");

    // Finally for each level we get an int withe the corresponding dimension size
    for (unsigned int i = 0 ; i < pirParam.d ; i++)
    {
      if ( read(sessionSocket, boost::asio::buffer(	&(pirParam.n[i]), sizeof(int))) <= 0) 
        exitWithErrorMessage(__FUNCTION__, "No pir param recieved");
    }
    // The last dimension + 1 is set to 1 (used in some functions to compute the number of
    // elements after a PIR recursion level)
    pirParam.n[pirParam.d] = 1;

    // Pass the PIR parameters to the reply generator
    generator->setPirParams(pirParam);

    cout << "PIRSession: PIR params recieved from the client are d=" << pirParam.d << ", alpha=" << pirParam.alpha << ", data_layout=";
    for (unsigned int i = 0; i < pirParam.d; i++)
    {
      if (pirParam.d == 1) cout << "1x";
      cout << pirParam.n[i];
      if (i != pirParam.d-1) cout << "x";
    }
    cout << endl;
  }catch (std::exception const& ex)
  {
    exitWithErrorMessage(__FUNCTION__, string(ex.what()));
  }
}

void PIRSession::sendPirParams() 
{
  cout << "PIRSession: Mandatory PIR params sent to the client are d=" << pirParam.d << ", alpha=" << pirParam.alpha << ", data_layout=";
  for (unsigned int i = 0; i < pirParam.d; i++)
  {
    if (pirParam.d == 1) cout << "1x";
    cout << pirParam.n[i];
    if (i != pirParam.d-1) cout << "x";
  }
  cout << endl;
  try{
    // First we send an int with the recursion level
    if ( write(sessionSocket, boost::asio::buffer(&(pirParam.d), sizeof(int))) <= 0) 
      exitWithErrorMessage(__FUNCTION__, "No pir param sended");

    // Then we send an int with the aggregation 
    if ( write(sessionSocket, boost::asio::buffer(&(pirParam.alpha), sizeof(int))) <= 0) 
      exitWithErrorMessage(__FUNCTION__, "No pir param sended");

    // Finally for each level we send an int withe the corresponding dimension size
    for (unsigned int i = 0 ; i < pirParam.d ; i++)
    {
      if ( write(sessionSocket, boost::asio::buffer(	&(pirParam.n[i]), sizeof(int))) <= 0) 
        exitWithErrorMessage(__FUNCTION__, "No pir param sended");
    }
    // The last dimension + 1 is set to 1 (used in some functions to compute the number of
    // elements after a PIR recursion level)
    pirParam.n[pirParam.d] = 1;

    // Pass the PIR parameters to the reply generator
    generator->setPirParams(pirParam);

  }catch (std::exception const& ex)
  {
    exitWithErrorMessage(__FUNCTION__, string(ex.what()));
  }
}


/**
 * Start downloadworker in downThread.
 **/
void PIRSession::startProcessQuery ()
{
  if(no_pipeline_mode) {
    std::cout << "No pipeline in query processing." << std::endl;
    downloadWorker();
  } else {
    downThread = boost::thread(&PIRSession::downloadWorker, this);
  }
}


/**
 * Recieve queries n messages with n = nbr of files.
 **/
void blo(const boost::system::error_code& err) {
  std::cout <<"rec "<<omp_get_wtime()<<std::endl;
}
void PIRSession::downloadWorker()
{
  double start = omp_get_wtime();
  unsigned int msg_size = 0;

  // Allocate an array with d dimensions with pointers to arrays of n[i] lwe_query elements 
  generator->initQueriesBuffer();

#ifdef PERF_TIMERS
  double vtstart = omp_get_wtime();
  bool wasVerbose = false;
  unsigned int previous_elts = 0;
  unsigned int total_elts = 0;
  for (unsigned int k = 0 ;  k < pirParam.d ; k++) total_elts += pirParam.n[k];
#endif

  try{
    for (unsigned int j = 0 ; j < pirParam.d ; j++)
    {
      // Compute and allocate the size in bytes of a query ELEMENT of dimension j 
      msg_size = cryptoMethod->getPublicParameters().getQuerySizeFromRecLvl(j+1) / 8;
      boost::asio::socket_base::receive_buffer_size opt(65535);
      sessionSocket.set_option(opt);
  //    boost_buffer = new boost::asio::buffer(buf, msg_size); 
#ifdef DEBUG
      cout << "PIRSession: Size of the query element to be received is " << msg_size << endl;
      cout << "PIRSession: Number of query elements to be received is " << pirParam.n[j] << endl;
#endif

      // Iterate over all the elements of the query corresponding to the j-th dimension
      for (unsigned int i = 0; i < pirParam.n[j]; i++)
      {
      char *buf = (char *) malloc(msg_size*sizeof(char));
      auto boost_buffer = boost::asio::buffer(buf,msg_size); 
      if (i==0 && j == 0) cout << "PIRSession: Waiting for query elements ..." << endl;
        // Get a query element 
         //( async_read(sessionSocket, boost_buffer,boost::bind(&blo,boost::asio::placeholders::error)) );
        if (read(sessionSocket, boost_buffer) < msg_size )
          writeWarningMessage(__FUNCTION__, "Query element not entirely recieved");	
//  std::cout <<"PIRSession: " << total_elts << " query elements received in " << omp_get_wtime() - start << std::endl;

        // Allocate the memory for the element, copy it, and point to it with the query buffer  
        if (i==0 && j == 0) cout << "PIRSession: Starting query element reception"  << endl;

#ifdef PERF_TIMERS
        // Give some feedback if it takes too long
        double vtstop = omp_get_wtime();
        if (vtstop - vtstart > 1) 
        {
          vtstart = vtstop;
          previous_elts = 0;
          for (unsigned int k = 0 ;  k < j ; k++) previous_elts += pirParam.n[k];
          std::cout <<"PIRSession: Query element " << i+1+previous_elts << "/" << total_elts << " received\r" << std::flush;
          wasVerbose = true;
        }
#endif
        
        generator->pushQuery(buf, msg_size, j, i);
      }

    }
  }catch (std::exception const& ex)
  {
    exitWithErrorMessage(__FUNCTION__, string(ex.what()));
    return;
  }

#ifdef PERF_TIMERS
  std::cout <<"PIRSession: Query element " << total_elts << "/" << total_elts << " received" << std::endl;
  std::cout <<"PIRSession: " << total_elts << " query elements received in " << omp_get_wtime() - start << std::endl;
#endif
  
  // All the query elements received, unlock reply generation
  generator->mutex.unlock();

  // Output we are done
  writeWarningMessage(__FUNCTION__, "done.");
}


/**
 * Start uploadWorker in a thread.
 **/
void PIRSession::startProcessResult(session_option_t session_option)
{
  if (no_pipeline_mode) {
    std::cout << "No pipeline in query Generation." << std::endl;
  }
  else {
    upThread = boost::thread(&PIRSession::uploadWorker, this);
  }
  // Import and generate reply once unlocked by the query downloader thread
  // If we got a preimported database generate reply directly from it 
  if (session_option.got_preimported_database == true) 
  {
    std::cout << "PIRSession: Already got an imported database available, using it" << std::endl;
    generator->generateReplyGenericFromData(session_option.data);
  }
  else if (session_option.keep_database == true) {
    savedDatabase = generator->generateReplyGeneric(true);
  } 
  else if (session_option.keep_database == false) 
  {
    generator->generateReplyGeneric(false);
  }
  
  if (no_pipeline_mode) {
    uploadWorker();
  } 
}

void sleepForBytes(unsigned int bytes) {
#ifdef NDSS_UPLOAD_SPEED
	uint64_t seconds=(bytes*8)/NDSS_UPLOAD_SPEED;
	uint64_t nanoseconds=((((double)bytes*8.)/(double)NDSS_UPLOAD_SPEED)-(double)seconds)*1000000000UL;
	  struct timespec req={0},rem={0};
	  req.tv_sec=seconds;
	  req.tv_nsec=nanoseconds;
	  nanosleep(&req,&rem);
#endif
}

/**
 * Send PIR's result, asynchronously. 
 **/
void PIRSession::uploadWorker() 
{
  // Ciphertext byte size 
  unsigned int byteSize = cryptoMethod->getPublicParameters().getCiphBitsizeFromRecLvl(pirParam.d)/GlobalConstant::kBitsPerByte;
  uint64_t totalbytesent=0;
  // Number of ciphertexts in the reply
  unsigned long reply_nbr = generator->computeReplySizeInChunks(maxFileBytesize), i = 0;
#ifdef DEBUG
  cout << "PIRSession: Number of ciphertexts to send is " << reply_nbr << endl;
  cout << "PIRSession: maxFileBytesize is  " << maxFileBytesize << endl;
  cout << "PIRSession: Ciphertext bytesize is " << byteSize << endl;
#endif

  try
  {
    // Pointer for the ciphertexts to be sent
    char *ptr;

    // For each ciphertext in the reply
    for (unsigned i = 0 ; i < reply_nbr ; i++)
    {
      while (generator->repliesArray == NULL || generator->repliesArray[i] == NULL) 
      { 
        boost::this_thread::sleep(boost::posix_time::milliseconds(10));
      }
      ptr = generator->repliesArray[i];
      
      // Send it
      //int byteSent=sessionSocket.send(boost::asio::buffer(ptr, byteSize));
      if (write(sessionSocket,boost::asio::buffer(ptr, byteSize)) <= 0)
        exitWithErrorMessage(__FUNCTION__,"Error sending request" );
	  	totalbytesent+=byteSize;
#ifdef NDSS_UPLOAD_SPEED
      sleepForBytes(byteSize);
#endif
      // Free its memory
      free(ptr);
      generator->repliesArray[i]=NULL;
    }

    // When everythind is send close the socket
    sessionSocket.close();
  }
  // If there was a problem sending the reply
  catch (std::exception const& ex)
  {
#ifdef DEBUG
    std::cerr << "Number of chunks sent: " << i << "/" << reply_nbr << std::endl;
#endif
    exitWithErrorMessage(__FUNCTION__, string(ex.what()));
    return;
  }


  // Tell when we have finished
  writeWarningMessage(__FUNCTION__ , "done.");;
}


// Functions for displaying logs 
void PIRSession::writeErrorMessage(string funcName, string message)
{
  cerr << BOLD  << funcName << " : " << RESET_COLOR << RED << message << RESET_COLOR <<endl;
}
void PIRSession::writeWarningMessage(string funcName, string message) 
{
  cerr << BOLD  << funcName << " : " << RESET_COLOR << ORANGE << message << RESET_COLOR <<endl;
}

// For critical erros
void PIRSession::exitWithErrorMessage(string funcName, string message) 
{
  writeErrorMessage(funcName, message);

  // This is used in the main function (start()) to skip costly operations if an error occurred
  handmadeExceptionRaised = true;	
}


// Used by the PIRServer for garbage collection
bool PIRSession::isFinished()
{
  return finished;
}

// Destructor
PIRSession::~PIRSession()
{
  if (cryptoMethod != NULL) delete cryptoMethod;

}

PIRParameters PIRSession::getPIRParams()
{
  pirParam.crypto_params = cryptoMethod->getPublicParameters().getSerializedParams(false);
  return pirParam;
}

void PIRSession::setPIRParams(PIRParameters pir_parameters)
{
  pirParam = pir_parameters;
}

imported_database_t PIRSession::getSavedDatabase()  {
   return  savedDatabase;
}

void PIRSession::no_pipeline(bool b)
{
  no_pipeline_mode = b;
}
