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

#include "main.hpp"

// Client constants
static const std::string DEFAULT_IP("127.0.0.1");
static const int DEFAULT_PORT = 1234;
static const bool DEFAULT_AUTOCHOICE = false;

// Optimizer constants
static const int DEFAULT_SECURITY = 80;
static const unsigned long DEFAULT_N = 1000;
static const unsigned int DEFAULT_L = 24000000;
static const unsigned int DEFAULT_TUPC = 100000000;
static const unsigned int DEFAULT_TDOC = 100000000;
static const unsigned int DEFAULT_K = 80;
static const unsigned int DEFAULT_ALPHAMAX = 0; // 0 means no limit
static const unsigned int DEFAULT_DMAX = 4;
static const unsigned int DEFAULT_DMIN = 1; 
static const FitnessType DEFAULT_FITNESSMETHOD = MAX;
static const std::string DEFAULT_CRYPTO_PARAMS("LWE:80:1024:60:22");


//Global vars
string fileToOptimizeFrom;
bool no_pipeline = false;

void sighandler(int sig_num)
{
  FixedVars fixedVars;
  cerr.flush();
  if(sig_num == SIGPIPE)
    cerr << "Broken pipe detected";

  cerr << endl << "Exiting client..." << endl << endl;
  exit(EXIT_SUCCESS);       
}


void defineClientOptions(ClientParams* paramsPtr, po::options_description* odptr){
  odptr->add_options()
    ("help,h", 
     "help message") 
    
    ("serverip,ip", 
     po::value<string>(&paramsPtr->server_ip)->default_value(DEFAULT_IP), 
     "PIR server IP" )

    ("port,p", 
     po::value<int>(&paramsPtr->port)->default_value(DEFAULT_PORT), 
     "PIR server port")
    
    ("autochoice,c", 
     "Auto-choose the first file")
    
    ("dry-run", 
     "Enable dry-run mode")
    
    ("verbose-optim", 
     "Ask the optimizer to be more verbose")
    
    ("dont-write", 
     "Don't write result to a file")
    
    ("file,f", 
     po::value<std::string>(&fileToOptimizeFrom), 
     "Use a config file to test different optimizations in dry-run mode (see sample.conf)"); 
}


void defineOptimizerOptions(FixedVars* varsPtr, po::options_description* odptr){
  odptr->add_options()
    ("file-nbr,n", 
     po::value<uint64_t>(), 
     "Used in dry-run mode only: Number of database elements" )

    ("file-size,l", 
     po::value<uint64_t>(), 
     "Used in dry-run mode only: Database element size in bits")
    
    ("upload,up", 
     po::value<double>(),  
     "Force client upload speed in bits/s (bandwith test will be skipped)")
    
    ("download,down", 
     po::value<double>(), 
     "Force client download speed in bits/s (bandwidth test will be skipped)")
    
    ("crypto-params,r", 
     po::value<string>(&varsPtr->manual_crypto_params), 
     "Set cryptographic parameteres manually")
    
    ("security,k", 
     po::value<unsigned int>(&varsPtr->k)->default_value(DEFAULT_SECURITY), 
     "Security bits wanted")
    
    ("dmin", 
     po::value<unsigned int>(&varsPtr->dMin)->default_value(DEFAULT_DMIN), 
     "Min dimension value to test")
    
    ("dmax", 
     po::value<unsigned int>(&varsPtr->dMax)->default_value(DEFAULT_DMAX), 
     "Max dimension value to test")
    
    ("alphaMax,a", 
     po::value<unsigned int>(&varsPtr->alphaMax)->default_value(DEFAULT_ALPHAMAX), 
     "Max aggregation value to test (1 = no aggregation, 0 = no limit)")
    
    ("fitness,x",
     po::value<int>((int*)&varsPtr->fitness)->default_value((int)DEFAULT_FITNESSMETHOD),
     "Set fitness method to: \n0=SUM Sum of the times on each task\n1=MAX Max of server times + Max of client times\n2=CLOUD Dollars in a cloud model (see sourcecode)");    
}


void defineHiddenOptions(FixedVars* varsPtr, po::options_description* odptr){
  odptr->add_options()
    ("no-pipeline", "No pipeline mode\n")

    ("reclvl", 
     po::value<unsigned int>(), 
     "Number of dimension used for database representation"); 
}


void processOptions(FixedVars* varsPtr, ClientParams* paramsPtr, po::variables_map vm){

  // Client options
  if(vm.count("serverip")) 
  {
    std::cout << "CLI: Server ip set to " <<  paramsPtr->server_ip << std::endl;
  }
  
  if(vm.count("port")) 
  {
    std::cout << "CLI: Upload port set to " << paramsPtr->port << std::endl;
  }
  
  if(vm.count("autochoice")) 
  {
    std::cout << "CLI: Auto-choice activated" << std::endl;
    paramsPtr->autochoice = true;
  }
  
  if(vm.count("dry-run")) 
  {
    std::cout << "CLI: Dry-run mode activated" << std::endl;
    paramsPtr->dryrunmode = true;
    varsPtr->n = DEFAULT_N;
    varsPtr->l = DEFAULT_L;
    varsPtr->Tupc = DEFAULT_TUPC;
    varsPtr->Tdos = DEFAULT_TUPC;
    varsPtr->Tdoc = DEFAULT_TDOC;
    varsPtr->Tups = DEFAULT_TDOC;
    std::cout << "CLI: Setting default values for dry-run mode (a thousand mp3 files, 100MBit/s connection)" <<std::endl;
  }

  if(vm.count("verbose-optim")) 
  {
    std::cout << "CLI: Will ask the optimizer to be more verbose" << std::endl;
    paramsPtr->verboseoptim = true;
  }


  if(vm.count("dont-write")) 
  {
    std::cout << "CLI: Will ask PIRReplyWriter not to write to a file" << std::endl;
    paramsPtr->dontwrite = true;
  }
  else paramsPtr->dontwrite = false;

  // Optimizer options
  if(vm.count("file-nbr")) 
  {
    if(vm.count("dry-run"))
    {
      varsPtr->n = vm["file-nbr"].as<uint64_t>();
      std::cout << "CLI: Changing number of database elements to "\
        << varsPtr->n << std::endl;
    } 
    else 
    {
      std::cout << "CLI: Option -n,--file-nbr ignored as we are not in dry-run mode"\
        << std::endl; 
    }
  }

  if(vm.count("file-size")) {
    if(vm.count("dry-run"))
    {
      varsPtr->l = vm["file-size"].as<uint64_t>();
      std::cout << "CLI: Changing database element size to " << varsPtr->l <<std::endl;
    }
    else
    {
      std::cout << "CLI: Option -l,--file-size ignored as we are not in dry-run mode"
        << std::endl; 
    }
  }

  if (vm.count("upload"))
  {
    varsPtr->Tupc = vm["upload"].as<double>();
    varsPtr->Tdos = varsPtr->Tupc;
    cout << "CLI: Upload speed forced to "<< varsPtr->Tupc << endl;
  }
  
  if (vm.count("download"))
  {
    varsPtr->Tdoc = vm["download"].as<double>();
    varsPtr->Tups = varsPtr->Tdoc;
    cout << "CLI: Download speed forced to "<< varsPtr->Tdoc << endl;
  }
  
  if(vm.count("crypto-params")) 
  {
    std::cout << "CLI: Crypto parameters set to " << varsPtr->manual_crypto_params <<  std::endl;
    std::vector<std::string> fields;
    boost::algorithm::split(fields, varsPtr->manual_crypto_params, boost::algorithm::is_any_of(":"));
    if (fields.size()>1 && atoi(fields[1].c_str()) > 0)
    {
      varsPtr->k = atoi(fields[1].c_str());
    }
    if (fields.size()>4)
    { 
      std::cout << "CLI: WARNING Absorption size will be overriden by the optimizer" << varsPtr->manual_crypto_params <<  std::endl;
      
      varsPtr->k = 0;
    } 
  }

  if(vm.count("security")) 
  {
    std::cout << "CLI: Security set to " << varsPtr->k <<  std::endl;
  }

  if (vm.count("dmin")) 
  {
    cout << "CLI: Minimum recursion level set to "<< varsPtr->dMin << endl;
  }
  
  if (vm.count("dmax")) 
  {
    cout << "CLI: Maximum recursion level set to "<< varsPtr->dMax << endl;
  }
  
  if(vm.count("alphaMax")) 
  {
    varsPtr->alphaMax = vm["alphaMax"].as<unsigned int>();
    cout << "CLI: Max aggregation set to "<< varsPtr->alphaMax << endl;
  }
  
  if(vm.count("fitness")) 
  {
    cout << "CLI: Fitness method set to "<< varsPtr->fitness << endl;
  }   
  
  // Hidden options
  if (vm.count("reclvl")) 
  {
    varsPtr->dMax = vm["reclvl"].as<unsigned int>(); 
    varsPtr->dMin = vm["reclvl"].as<unsigned int>(); 
    cout << "CLI: Recursion level forced to "<< varsPtr->dMax << endl;
  }
  
  if (vm.count("no-pipeline"))
  {
    std::cout << "CLI: WARNING no pipeline mode activated" << std::endl;
    no_pipeline = true;
  }
}



int main(int argc, char** argv) 
{
  boost::asio::io_service ios;
  
  // Vars for the optimizer (and pot. client)
  FixedVars fixedVars = {}; // Inits to default value all fields
  
  // Vars for the client only
  ClientParams clientParams = {}; // Same here

  // Add and define options
  po::options_description od("Client options");
  po::options_description optimizeropts("Optimizer options");
  po::options_description hidden("Hidden options");
  defineOptimizerOptions(&fixedVars, &optimizeropts);
  defineClientOptions(&clientParams, &od);
  defineHiddenOptions(&fixedVars, &hidden);

  // Set which options are visible and which not
  po::options_description visible;
  visible.add(od).add(optimizeropts);
  po::options_description all;
  all.add(visible).add(hidden);
  po::variables_map vm;

  // Parse options from command line
  try {
    po::store(po::parse_command_line(argc, argv, all), vm);
  }
  catch (const std::exception& ex) 
  {
    std::cout << "CLI: Error checking program options: " << ex.what() << std::endl;
    std::cout << visible << std::endl;
    return 1;
  }
  po::notify(vm);

  // Show usage help if requested
  if(vm.count("help")) {
    std::cout << visible << endl;
    return 0;
  }

  // Set variables according to options
  processOptions(&fixedVars, &clientParams, vm);
 
  
  // If we are on dry-run mode run only the optimizer not the client
  if(vm.count("dry-run"))
  {
    std::cout << "CLI: Dry-run mode activated" << std::endl;

    // The optimizer connects to the same PIR server but on a higher port for 
    // optimization related exchanges
    PIROptimizer optimizer(clientParams.server_ip, clientParams.port, fixedVars.fitness);
    
    if (vm.count("file"))
    {
      int experience_nbr = OptimService::getNumberOfExperiences(fileToOptimizeFrom);
      if(experience_nbr == -1)
      {
        cout << "CLI: Unable to open : " << fileToOptimizeFrom << " aborting" << endl;
        cout << "CLI: Try exp/sample.conf to test a predefined set of configurations" << endl;
        return 1;
      }
      for (int exp_i = 0 ; exp_i <= experience_nbr ; exp_i++)
      {
        OptimService::readTestValues(exp_i, fixedVars, fileToOptimizeFrom);
        OptimService::writeHeadFile(exp_i, fixedVars);
        optimizer.optimize(fixedVars, exp_i);
      }
    }
    else
    {
      cout << "CLI: Default example : A thousand mp3 files, ADSL, no aggregation, k=80" << endl;
      cout << "CLI: n : " << fixedVars.n << " l : " << fixedVars.l;
      cout << "CLI: Tupc : " << fixedVars.Tupc << " Tdoc : " << fixedVars.Tdoc << endl;

      OptimService::writeHeadFile(0, fixedVars);
      optimizer.optimize(fixedVars, 0);
    }
    return 0;
  }

  // If we are not in dry-run mode create client and controller
  PIRClientSimple client(ios, clientParams, fixedVars);
  PIRController controller(client);
 
  // Set no_pipeline if needed
  if(no_pipeline) client.no_pipeline(true);

  // Start by connecting to the PIR server
  client.connect();

  // Downloads the file catalog
  client.downloadCatalog();

  // Get from the server whether we are in client-driven mode or not
  client.rcvPIRParamsExchangeMethod();

  // Use the optimizer to choose best parameters 
  // (returns immediately in server-driven mode)
  client.optimize();  

  // Send PIR and cryptographic parameters to the server in client-driven mode
  // and receive and process them in server-driven mode
  client.processCryptoParams();
  client.processPIRParams();

  /*User chooses the file here.*/
  client.chooseFile();
  
  double start = omp_get_wtime();
  
  
  /* Asynchronously generate and send the request
     separately in two threads*/
  client.startProcessQuery();
  /* Receive asynchronously the response from the server and
     asynchronously writes it */
  client.startProcessResult();

  client.joinAllThreads();
  double end = omp_get_wtime();
  cout << "CLI: Query RTT was " << end-start << " seconds" << endl;
  cout << "CLI: Exiting..." << endl;

  return EXIT_SUCCESS;
}

