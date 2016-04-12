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

#include "PIRClient.hpp"
//Use this to limit the upload speed to UPLOAD_LIMIT bits per second
//#define UPLOAD_LIMIT 100000000UL


PIRClientSimple::PIRClientSimple(boost::asio::io_service& ios, ClientParams params, FixedVars vars):
	socket_up(ios),
  replyWriter(pirParams, writeListeners, messageListeners),
	catalog(messageListeners),
  clientParams(params),
  fixedVars(vars),
  optimum(vars),
  no_pipeline_mode(false)
{
  replyWriter.setdontWrite(params.dontwrite);
      //string lwe_params(crypto_params_const);
  //cryptoMethod.getPublicParameters().computeNewParameters(lwe_params);
  //pirParams.d = kDimensionNumber;
}

void PIRClientSimple::joinAllThreads()
{
  replyExt->replyThread.join();
  replyWriter.join();
}

PIRClientSimple::~PIRClientSimple()
{
  joinAllThreads();
  delete cryptoMethod;
}

void PIRClientSimple::connect()
{	
	using namespace boost::asio::ip;

	try
	{
		socket_up.connect(tcp::endpoint(address::from_string(clientParams.server_ip), clientParams.port));
    
    // Tell the server we are a client
    int is_client = 1;
		write(socket_up, boost::asio::buffer(&is_client, sizeof(int)));

  }
	catch (const std::exception& ex)
	{
		exitWithErrorMessage(__FUNCTION__, ex.what());
	}
}
/**
 *  Downloads the file catalog from the server.
 *  Exception :
 *  - Stop the programme if any network trouble.
 **/
void PIRClientSimple::downloadCatalog()
{
	boost::uint64_t size = 1;
	try
	{
		read(socket_up, boost::asio::buffer(&size, sizeof(size)));
		std::cout<< "PIRClient: Catalog bytesize is "<<size<<std::endl;
		char* buf = new char[size + 1]();
	
		if (read(socket_up, boost::asio::buffer(buf, size)) < size )
      	  writeWarningMessage(__FUNCTION__, " Catalog has not been fully received.");
   
		catalog.makeMenu(buf);

		if (catalog.getFilesNum() == 0) 
		{
			writeWarningMessage(__FUNCTION__,"Empty catalog recieved ! Is there a db directory with files in the server/ directory ? Aborting ");
			exit(1);
		}
   
    delete[] buf;
    
    std::cout << "PIRClient: Catalog received with " << catalog.getFilesNum() << " files" << std::endl;
    
  }
	catch (const std::exception& ex)
	{
		exitWithErrorMessage(__FUNCTION__, ex.what());
	}
}

void PIRClientSimple::optimize()
{
  if (exchange_method == CLIENT_DRIVEN)
  {
    PIROptimizer optimizer;
    if (clientParams.verboseoptim == false) optimizer.silent = true;
    std::cout << "PIRClient: Starting optimization ..." << std::endl;
    OptimVars op = optimizer.optimizeFromServer(fixedVars, socket_up);
    optimum = op;
    std::cout << "PIRClient: Optimization returned alpha=" << optimum.alpha << " d=" << optimum.d << " crypto_params=" << optimum.crypto_params << std::endl;
    pirParams.d = optimum.d;
    pirParams.alpha = optimum.alpha;
    optimizer.getDimSize(optimum.getFixedVars().n, optimum.alpha, optimum.d, pirParams.n);
    pirParams.crypto_params = std::string(optimum.crypto_params);
  }
}

void PIRClientSimple::processPIRParams()
{
  if (exchange_method == CLIENT_DRIVEN)
  {
    sendPirParams();
  }
  else rcvPirParams();
}

/**
 *  Sends PIR parameters to the server.
 *  Exception :
 *  - Stop the programme if any network trouble.
 **/
void PIRClientSimple::sendPirParams()
{
  try 
	{
		write(socket_up, boost::asio::buffer(&pirParams.d, sizeof(int)));
		write(socket_up, boost::asio::buffer(&pirParams.alpha, sizeof(int)));

		for (unsigned int i = 0 ; i < pirParams.d ; i++)
		{
			write(socket_up, boost::asio::buffer(&pirParams.n[i], sizeof(int)));
		}
	}
	catch (const std::exception& ex)
	{
		exitWithErrorMessage(__FUNCTION__, ex.what());
	}
}

void PIRClientSimple::rcvPirParams()
{
  try 
	{
		read(socket_up, boost::asio::buffer(&pirParams.d, sizeof(int)));
		read(socket_up, boost::asio::buffer(&pirParams.alpha, sizeof(int)));

		for (unsigned int i = 0 ; i < pirParams.d ; i++)
		{
			read(socket_up, boost::asio::buffer(&pirParams.n[i], sizeof(int)));
		}

    std::cout << "PIRClient: Received PIR parameters alpha=" << pirParams.alpha << " d=" << pirParams.d << std::endl; 
	}
	catch (const std::exception& ex)
	{
		exitWithErrorMessage(__FUNCTION__, ex.what());
	}
}

void PIRClientSimple::rcvPIRParamsExchangeMethod()
{
  try
  {
    read(socket_up, boost::asio::buffer(&exchange_method, sizeof(short)));
    if (exchange_method == CLIENT_DRIVEN) 
    {
      std::cout << "PIRClient: Client-driven mode used" << std::endl;
    } 
    else 
    {
      std::cout << "PIRClient: Server-driven mode used" << std::endl;
    }
  }
	catch (const std::exception& ex)
  {
		exitWithErrorMessage(__FUNCTION__, ex.what());
  }
}

void PIRClientSimple::processCryptoParams() {
  bool send_paramsandkey = true;
  if (exchange_method == CLIENT_DRIVEN) 
  {
    cryptoMethod = HomomorphicCryptoFactory_internal::getCryptoMethod(optimum.crypto_params);
    replyWriter.setCryptoMethod(cryptoMethod);
    sendCryptoParams(send_paramsandkey);
  }
  else 
  {
    rcvCryptoParams(); 
    cryptoMethod = HomomorphicCryptoFactory_internal::getCryptoMethod(pirParams.crypto_params);
    replyWriter.setCryptoMethod(cryptoMethod);
    sendCryptoParams(!send_paramsandkey);
  }
}

/**
 *  Send cryptographics parameters (e.g : key) to the server.
 *  Exception :
 *  - Stop the programme if any network trouble.
 **/
void PIRClientSimple::sendCryptoParams(bool paramsandkey)
{
	unsigned int size; 
	char* key;

	try
	{
    if (paramsandkey == true)
    {
      bool shortversion = true;

      string crypto_params = cryptoMethod->getSerializedCryptoParams(!shortversion);
      size = crypto_params.length();

		  write(socket_up, boost::asio::buffer(&size, sizeof(size)));
      write(socket_up, boost::asio::buffer(crypto_params));
    }

    size = ceil((double) cryptoMethod->getPublicParameters().getSerializedModulusBitsize() / GlobalConstant::kBitsPerByte);
		write(socket_up, boost::asio::buffer(&size, sizeof(size)));

		key = cryptoMethod->getPublicParameters().getByteModulus();
		write(socket_up, boost::asio::buffer(key, size));//Send key
		delete[] key;
	}
	catch (const std::exception& ex) 
	{
		exitWithErrorMessage(__FUNCTION__, ex.what());
	}
}

void PIRClientSimple::rcvCryptoParams() {

  try
  {
    int crypto_params_size = 0; 

    read(socket_up, boost::asio::buffer(&crypto_params_size, sizeof(crypto_params_size)));

    char crypto_params[crypto_params_size + 1];

    int read_size = read(socket_up, boost::asio::buffer(crypto_params,crypto_params_size));
    crypto_params[read_size] = '\0';

    pirParams.crypto_params = string(crypto_params);

    std::cout << "PIRClient: Received crypto parameters " << crypto_params << std::endl; 
  }
	catch (const std::exception& ex) 
	{
		exitWithErrorMessage(__FUNCTION__, ex.what());
	}

}

/**
 *  Launches the choose File action to the recorded view. If autochoice is set, 
 *  the first file is chosen.
 **/
void PIRClientSimple::chooseFile() 
{
	chosenElement = 0;

  if(!clientParams.autochoice)
	{
		CatalogEvent catalogEvent(catalog.getFileList());
		menuListeners(catalogEvent); // get user input and set chosenElement.
	}

	std::ostringstream oss;
	oss << "PIRClient: Retrieved file is nbr " + (chosenElement+1);
	oss << + ", \"" + catalog.getFileName(chosenElement) +  "\"";

	MessageEvent messageEvent(oss.str());
	messageListeners(messageEvent);
}

/**
 *  Set the query generator and launch parallely query generation and query upload.
 **/
void PIRClientSimple::startProcessQuery()
{
	PIRQueryGenerator_internal queryGen(pirParams, *cryptoMethod);

	queryGen.setChosenElement(chosenElement);


  if(no_pipeline_mode) {
    std::cout << "PIRClient: Calling query generation and upload without pipelining" << std::endl;
    queryGen.generateQuery(); 
  } else {
    queryGen.startGenerateQuery();
  }

	uploadWorker(queryGen);
}
void sleepForBytes(unsigned int bytes) {
#ifdef UPLOAD_LIMIT
	uint64_t seconds=(bytes*8)/UPLOAD_LIMIT;
	uint64_t nanoseconds=((((double)bytes*8.)/(double)UPLOAD_LIMIT)-(double)seconds)*1000000000UL;
	  struct timespec req={0},rem={0};
	  req.tv_sec=seconds;
	  req.tv_nsec=nanoseconds;
    double st=omp_get_wtime();
	  nanosleep(&req,&rem);
#endif
}
/**
 *  Upload query to the server et delete its parts from the shared query queue.
 *  Exception :
 *  - Stop the programme if any network trouble.
 **/
void PIRClientSimple::uploadWorker(PIRQueryGenerator_internal& queryGen)
{
	unsigned int size = 0;
	char *tmp;
	try 
	{
		for (unsigned int j = 1 ; j <= pirParams.d ; j++)
		{
			size = cryptoMethod->getPublicParameters().getQuerySizeFromRecLvl(j) / GlobalConstant::kBitsPerByte;

			for (unsigned int i = 0 ; i < pirParams.n[j - 1] ; i++)
			{
				tmp = queryGen.queryBuffer.pop_front();
				write(socket_up, boost::asio::buffer(tmp, size));
				free(tmp);
#ifdef UPLOAD_LIMIT
        sleepForBytes(size);
#endif
			}
		}
	} 
	catch (const std::exception& ex) 
	{
		exitWithErrorMessage(__FUNCTION__, ex.what());
	}

	
}

/**
 *  Sets reply extractpr and lauches parallely reply download and reply extraction.
 **/
void PIRClientSimple::startProcessResult()
{
	replyExt = new PIRReplyExtraction_internal(pirParams,*cryptoMethod);
	replyExt->setFileParam(catalog.getFileName(chosenElement), catalog.getFileSize(chosenElement)); 

  if(no_pipeline_mode)
  {
    std::cout << "PIRClient: Calling reply download, extraction and writting without pipelining" << std::endl;
    downloadWorker(*replyExt);
    replyExt->extractReply(catalog.getMaxFileSize()*pirParams.alpha, replyWriter.getClearDataQueue());
    // Tell reply writer we finished the extraction
    replyWriter.writeAggregatedFileSecurely(chosenElement, catalog);
  }
  else 
  {
	replyExt->startExtractReply(catalog.getMaxFileSize()*pirParams.alpha, replyWriter.getClearDataQueue());
  
  replyWriter.startFileWritting(chosenElement, catalog);
	downloadWorker(*replyExt);
  }
}

/**
 *  Download reply from the server and stores it chunks in shared replies queue.
 *  Exception :
 *  - Stop the program if any network trouble.
 **/
void PIRClientSimple::downloadWorker(PIRReplyExtraction_internal& turlututu)
{
  using  namespace GlobalConstant;
  unsigned int i,  ciph_siz = cryptoMethod->getPublicParameters().getCiphBitsizeFromRecLvl(pirParams.d)/kBitsPerByte;

  double paquet_nbr = ceil(static_cast<double>(catalog.getMaxFileSize()*pirParams.alpha) / double(cryptoMethod->getPublicParameters().getAbsorptionBitsize(0)/kBitsPerByte)); 

#ifdef DEBUG
  cout << "PIRClient: First layer nbr of ciphertexts is " << paquet_nbr << endl;
#endif

    for (unsigned int i = 1 ; i < pirParams.d; i++) paquet_nbr =  ceil(paquet_nbr * double(cryptoMethod->getPublicParameters().getCiphBitsizeFromRecLvl(i)/kBitsPerByte) / double(cryptoMethod->getPublicParameters().getAbsorptionBitsize(i) / kBitsPerByte));

	char* buf;

#ifdef DEBUG
  cout << "PIRClient: getMaxFileSize=" << catalog.getMaxFileSize() << " alpha=" << pirParams.alpha << " getAbsorptionBitsize=" << double(cryptoMethod->getPublicParameters().getAbsorptionBitsize(0))/*kBitsPerByte);*/ << endl;
  cout << "PIRClient: Last layer nbr of ciphertexts is " << paquet_nbr << endl;
  cout << "PIRClient: File size of the chosen element is " << catalog.getFileSize(chosenElement) << endl;
	cout << "PIRClient: Ciphertext bytesize is " << ciph_siz << endl;
#endif

	try
	{
		for (i = 0 ; i < paquet_nbr ; i++)
		{
			buf = (char*) malloc(ciph_siz * sizeof(char));
			read(socket_up,boost::asio::buffer(buf, ciph_siz));
			replyExt->repliesBuffer.push(buf);
		}
	}
	catch (const std::exception& ex) 
	{
#ifdef DEBUG
    cerr << "PIRReplyWriter: Number of downloaded chunks: " << i << "/" << paquet_nbr << endl;
#endif
		exitWithErrorMessage(__FUNCTION__, ex.what());
    return;
	}
  writeWarningMessage(__FUNCTION__, "done.");
}

/**
 *  Sets the chosen element.
 **/
void PIRClientSimple::setChosenElement(uint64_t choice)
{
	chosenElement = choice;
}

//void PIRClientSimple::setPIRParameters(PIRParameters& pirParameters)
//{
//  pirParams.alpha = pirParameters.alpha;
//  pirParams.d = pirParameters.d;
//  memcpy(pirParams.n, pirParameters.n, MAX_REC_LVL);
//  pirParams.r;
//}

/**
 *	Add a new message listener.
 *	Param :
 *		- messageListener::slot_function_type subscriber, a message listener to add.
 **/
boost::signals2::connection PIRClientSimple::addMessageListener(messageListener::slot_function_type subscriber)
{
	return messageListeners.connect(subscriber);
}

/**
 *	Add a new menu listener.
 *	Param :
 *		- menuListener::slot_function_type subscriber, a menu listener to add.
 **/
boost::signals2::connection  PIRClientSimple::addMenuListener(menuListener::slot_function_type subscriber)
{
	return menuListeners.connect(subscriber);
}

/**
 *	Add a new wite listener
 *	Param :
 *		- writeListener::slot_function_type subscriber, a write listener to add.
 **/
boost::signals2::connection PIRClientSimple::addWriteListener(writeListener::slot_function_type subscriber)
{
	return writeListeners.connect(subscriber);
}

void PIRClientSimple::exitWithErrorMessage(string s1, string s2)
{
	writeErrorMessage(s1 ,s2);
	socket_up.close();
	exit(errno);
}

/**
 *	Send ERROR message to listeners
 *	Params :
 *		- string funcName :  usually the function name who throw the error message ;
 *		- string message  :  message to show.
 **/
void PIRClientSimple::writeErrorMessage(string funcName, string message)
{
	MessageEvent event(ERROR, message, funcName);
	messageListeners(event);	
}

/**
 *	Send WARNING message to listeners
 *	Params :
 *		- string funcName :  usually the function name who throw the warning message ;
 *		- string message  :  message to show.
 **/
void PIRClientSimple::writeWarningMessage(string funcName, string message) 
{
	MessageEvent event(WARNING, message, funcName);
	messageListeners(event);
}

void PIRClientSimple::no_pipeline(bool b) 
{
  no_pipeline_mode = b;
}
