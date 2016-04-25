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

#include "PIROptimizer.hpp"

static const unsigned int kPrecision = 5;

using boost::asio::ip::tcp;

PIROptimizer::PIROptimizer(DBHandler *db) :
  pir("lipmaa"), 
  filePath("optim/preCompute.abs")
{
  fileCount = db->getNbStream();
  maxFileBytesize = db->getmaxFileBytesize();
}

PIROptimizer::~PIROptimizer()
{}

void PIROptimizer::speedTest(boost::asio::ip::tcp::socket& s)
{
  try
  {
    unsigned int loop = 1;
    char* msg = (char *) malloc(MEGA_BYTE); 

    boost::system::error_code error;

    double start = omp_get_wtime();

    for (unsigned int i = 0 ; i < loop ; i++) boost::asio::read(s, boost::asio::buffer(msg, MEGA_BYTE));

    cout << "PIROptimizer: Speed test gives upload=" <<  double(loop / (omp_get_wtime() - start) ) << " MB/s";

    start = omp_get_wtime();

    for (unsigned int i = 0 ; i < loop ; i++) write(s, boost::asio::buffer(msg, MEGA_BYTE));

    cout << ", download=" <<  double(loop / (omp_get_wtime() - start) ) << " MB/s" << endl;
    free(msg);
  }
  catch(std::exception const& ex)
  {
    cout << "Error : " << ex.what() << endl;
  }
}

unsigned int PIROptimizer::sendAbsAndPrecomputeCaches(boost::asio::ip::tcp::socket& s)
{

  try
  {
    unsigned int crypto_name_size = 0;
    read(s, boost::asio::buffer(&crypto_name_size, sizeof(int)));

    char crypto_name[crypto_name_size + 1];

    crypto_name_size = read(s, boost::asio::buffer(crypto_name, crypto_name_size));
    crypto_name[crypto_name_size] = '\0';
    cout << "PIROptimizer: Sending absorption and precompute costs for " << crypto_name << endl;

    std::string file_path(OptimService::folderName + OptimService::fileName + crypto_name + OptimService::absFileExtension); 
    std::string file_content;
    
    if (OptimService::readEntireFile(file_content, file_path) < 0) { cout << "PIROptimizer: Error when reading file : " << file_path << ", abort." << endl; return 0;}
    
    int file_content_size = file_content.size();
 
#ifdef DEBUG
    std::cout << "PIROptimizer: Absorption and precompute costs file size " <<  file_content_size << std::endl;
#endif

    write(s, boost::asio::buffer(&file_content_size, sizeof(file_content_size)));

    write(s, boost::asio::buffer(file_content));

  }
  catch(std::exception const& ex)
  {
    cout << "Error : " << ex.what() << endl;
  }
  return 0;
}


/**
 * Build cache if necessary.
 **/
void PIROptimizer::prepareOptimData()
{
  for (auto crypto_name : HomomorphicCryptoFactory_internal::crypto_method_name_vec)
  {
    std::string file_path(OptimService::folderName + OptimService::fileName + crypto_name + OptimService::absFileExtension); 
    if (OptimService::fileOutdated(crypto_name, OptimService::absFileExtension))
    {
      std::cout << "PIROptimizer: Absorption and precompute performance cache is outdated, regenerating it" << std::endl; 
      std::string optim_data2write = computeOptimData(crypto_name);

      if(OptimService::writeOptimDataBuffer(optim_data2write, file_path)) {std::cout << "PIROptimizer: Error when writing optimization data, aborting." << std::endl; exit(1);}
      std::cout << "PIROptimizer: Finished generating the absorption and precompute performance cache" << std::endl; 
    }
  }
}

std::string PIROptimizer::computeOptimData(const std::string& crypto_name)
{
  double abs1plaintext_time, precompute1plaintext_time;

  GenericPIRReplyGenerator* generator_ptr = PIRReplyGeneratorFactory::getPIRReplyGenerator(crypto_name);
  HomomorphicCrypto* crypto_ptr = HomomorphicCryptoFactory_internal::getCrypto(crypto_name);

  std::set<std::string> crypto_params_set;
  string optim_data2write;

  generator_ptr->setCryptoMethod(crypto_ptr);
  crypto_ptr->getAllCryptoParams(crypto_params_set);

  for (auto crypto_param : crypto_params_set)
  {
      std::cout << "PIROptimizer: Generating cache for " << crypto_param << std::endl; 

    crypto_ptr->getPublicParameters().computeNewParameters(crypto_param);    
    crypto_ptr->getPublicParameters().setMockedPubKey();

    abs1plaintext_time = getAbs1PlaintextTime(crypto_ptr, generator_ptr);
    precompute1plaintext_time = getPrecompute1PlaintextTime(crypto_ptr, generator_ptr);

    std::ostringstream out;
    out << std::setprecision(kPrecision) << abs1plaintext_time << " " << precompute1plaintext_time;
    optim_data2write += crypto_param + " " + out.str() + "\n";
  }

  delete generator_ptr;
  delete crypto_ptr;

  return optim_data2write;
}


double PIROptimizer::getAbs1PlaintextTime(HomomorphicCrypto* crypto_ptr, GenericPIRReplyGenerator* generator)
{
  double result;
  uint64_t plaintext_nbr; 
  PIRParameters pir_params;
  pir_params.d = 1;
  pir_params.alpha = 1;
  pir_params.n[0] = 4;
  
  crypto_ptr->setandgetAbsBitPerCiphertext(pir_params.n[0]); // Set best absorption possible 
  plaintext_nbr = 4;

  do
  {
    generator->mutex.try_lock();
    generator->mutex.unlock();
    result = generator->generateReplySimulation(pir_params, plaintext_nbr);
    plaintext_nbr *= 2;
  }
  while(result < 0.5 && plaintext_nbr*pir_params.n[0]*crypto_ptr->getPublicParameters().getCiphertextBitsize() < (1UL<<30));

  plaintext_nbr /= 2;

  double plaintexts_in_database = (double)pir_params.n[0] * plaintext_nbr;

  return result / plaintexts_in_database;
}

double PIROptimizer::getPrecompute1PlaintextTime(HomomorphicCrypto* crypto_ptr, GenericPIRReplyGenerator* generator)
{
  double result;
  uint64_t plaintext_nbr; 
  PIRParameters pir_params;
  pir_params.d = 1;
  pir_params.alpha = 1;
  pir_params.n[0] = 4;
  
  crypto_ptr->setandgetAbsBitPerCiphertext(pir_params.n[0]); // Set best absorption possible 
  plaintext_nbr = 128;

  do
  {
    generator->mutex.try_lock();
    generator->mutex.unlock();
    result = generator->precomputationSimulation(pir_params, plaintext_nbr);
    plaintext_nbr *= 2;
  }
  while(result < 0.5 && result != 0.0 && plaintext_nbr*pir_params.n[0]*crypto_ptr->getPublicParameters().getCiphertextBitsize() < (1UL<<30));

  plaintext_nbr /= 2;

  double plaintexts_in_database = (double)pir_params.n[0] * plaintext_nbr;

  return result / plaintexts_in_database;
}

void PIROptimizer::serve()
{
  tcp::acceptor acceptor(io_service, tcp::endpoint(tcp::v4(), COMMAND_AND_CONTROL_PORT));
  acceptor.listen(5);

  cout << "Waiting client..." << endl;
  while(true)
  {
    boost::asio::ip::tcp::socket socket(acceptor.get_io_service());
    acceptor.accept(socket);
   // boost::thread t(boost::bind(&PIROptimizer::controlAndCommand, this, socket));
    controlAndCommand(socket);
  }
}

void PIROptimizer::sendDatabaseInfos(boost::asio::ip::tcp::socket&s)
{
  cout << "PIROptimizer: Sending database infos file_count=" << fileCount << ", max_bytesize=" << maxFileBytesize << endl;
  try{
   
    /*Send number of files*/
    boost::asio::write(s, boost::asio::buffer(&fileCount, sizeof(fileCount)));

    /*Send max file size*/
    boost::asio::write(s, boost::asio::buffer(&maxFileBytesize, sizeof(maxFileBytesize)));
  }catch(std::exception const&  ex)
  {
    cerr << "Error when sending database informations : " << ex.what() << std::endl;
  }
}

void PIROptimizer::controlAndCommand(boost::asio::ip::tcp::socket& s)
{
  size_t cmd = NOP;
  try{
    while(cmd != EXIT)
    {
      boost::asio::read(s, boost::asio::buffer(&cmd, sizeof(size_t)));

      switch (cmd)
      {
        case ABS :
          {
            sendAbsAndPrecomputeCaches(s);
            break;
          }
        case SPEED :
          {
            speedTest(s);
            break;
          }
        case EXIT :
          {
            cout << "PIROptimizer: EXIT command recieved" << endl;
            break;
          }
        case DATA :
          {
            sendDatabaseInfos(s);
            break;
          }
        default :
          {
            cout << "PIROptimizer: Looping 30s as requested by client ..." << endl;
            sleep(30);
           	break;
          }
      }
    }
  }catch(std::exception const&  ex)
  {
    cerr << "Client quit without send EXIT code: "  << ex.what() << endl;
  }
}

void PIROptimizer::optimize()
{
  serve();
}
