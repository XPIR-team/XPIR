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
// Set this to true if the server database is stored in precomputed (e.g. NTT) form
// TODO this should be an option given on the command line
#define INITIAL_PRECOMPUTATION_DONE false

static const unsigned int kPrecision = 5;
const unsigned int PIROptimizer::kalphaBound = 100;

/**
 *	Class constructor.
 *	Params :
 *		- pirClient_ptr client : PIRClient shared_ptr ;
 *		- int security_bits		 : security bits wanted ;
 *		- string server_ip		 : IPv4 server adresse	;
 *		- io_service& ios			 : boost io_service reference, used for boost socket ;
 **/
PIROptimizer::PIROptimizer(string server_ip_, int port_, FitnessType fitnessMethod_) :
  serv_ip(server_ip_),
  server_optimport(port_),
  s(ios),
  // Number of clients the servers would like to handle in parallel : Fixed for the moment (the server will tell in future versions)
  nbc(1),
  optimum(fitnessMethod_),
  silent(false)
{  
  // Max latency for the socket when connecting to server
  struct timeval tv;
  tv.tv_sec  = 5; 
  tv.tv_usec = 0;         
  setsockopt(s.native(), SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));
  setsockopt(s.native(), SOL_SOCKET, SO_SNDTIMEO, &tv, sizeof(tv));
}

PIROptimizer::PIROptimizer():
  s(ios),
  nbc(1),
  silent(false)
{}


OptimVars PIROptimizer::optimizeFromServer(FixedVars& partial_fixed_vars, boost::asio::ip::tcp::socket& socket)
{
  fixedVars = partial_fixed_vars;
  getFixedVarsFromServer(fixedVars, socket);
  getNetworkInfos(fixedVars, socket);
  optimize(socket, fixedVars, 0);
  sendExitCommand(socket);
  return global_optimum;
}
void PIROptimizer::optimize(FixedVars& fixed_vars, unsigned int exp_nbr)
{

	try
	{
		s.connect(tcp::endpoint(boost::asio::ip::address::from_string(serv_ip), server_optimport));

    // Tell server we are an optimizer not a client
    int is_client = 0;
		write(s, boost::asio::buffer(&is_client, sizeof(int)));
	}
	catch (const std::exception& ex)
	{
    std::cout << "Optimizer: Unable to connect to "<< serv_ip << " port " << server_optimport <<std::endl;
	}
  optimize(s, fixed_vars, exp_nbr);
}

/**
 * Main function of PIROptimizer.
 * Loop on all crypto modules to finc possible parameters
 * Initialize caches with computing costs
 * Iterate on all parameters except aggregation
 * Use findAlphaDicho to explore best aggregation values through dichotomy 
 **/
void PIROptimizer::optimize(boost::asio::ip::tcp::socket& s, FixedVars& fixed_vars, unsigned int exp_nbr)
{
  std::regex *cp_rex;
  fixedVars = fixed_vars;
  // Do we optimize the sum of all the times or we consider the client and the server as pipelines (sum of max) ?
  fitnessMethod = fixedVars.fitness;
  std::cout << "Optimizer: Using fitness method " << fitnessMethod << " (use --help for code translation)" << std::endl;

  set<std::string> crypto_params_desc_set;
  vector<HomomorphicCrypto*> crypto_systems_vec;

  // Get all the crypto modules
  HomomorphicCryptoFactory_internal::getAllCryptoSystems(crypto_systems_vec);

  if(fixedVars.manual_crypto_params == "")
  {
    noCryptoSystemOptim(exp_nbr);
  }
  else
  {
    std::cout << "Optimizer: Crypto parameters manually forced to " << fixedVars.manual_crypto_params << std::endl;

    // Get a regular expression from the requested params
    cp_rex = new std::regex(fixedVars.manual_crypto_params); 

    // Only take into account NoCryptography as a possibility if requested
    if(std::regex_match("NoCryptography", *cp_rex))
    {
      noCryptoSystemOptim(exp_nbr);
    }
  }
//    else
//    {
//      std::vector<std::string> fields;
//      boost::algorithm::split(fields, fixedVars.manual_crypto_params, boost::algorithm::is_any_of(":"));
//      HomomorphicCrypto* h;
//
//      h = HomomorphicCryptoFactory_internal::getCrypto(fields[0]);
//      crypto_systems_vec.push_back(h);
//    }
//  }


  // Iterate on all crypto modules
  for (auto crypto_system_ptr : crypto_systems_vec)
  {
    crypto = crypto_system_ptr;
    //Get i-th crypto module's public parameter
    currentPublicParams = &crypto->getPublicParameters();
  
    // Get all the proposed parameter sets for the requested security on this crypto module
    crypto->getCryptoParams(fixedVars.k, crypto_params_desc_set); 
    
    // If special parameters requested, remove all those that don't match the regexp
    if(fixedVars.manual_crypto_params != "")
    {
		const   set<std::string> crypto_params_desc_set_copy=crypto_params_desc_set;
		
      for (auto cp_ptr : crypto_params_desc_set_copy)
      {
        if (!std::regex_match(cp_ptr, *cp_rex)) crypto_params_desc_set.erase(cp_ptr);
        //HomomorphicCrypto* h;

        //h = HomomorphicCryptoFactory_internal::getCrypto(fields[0]);
        //crypto_systems_vec.push_back(h);
      }
      if (crypto_params_desc_set.size() == 0) continue;
      //crypto_params_desc_set.insert(fixedVars.manual_crypto_params);
    }
    // Add to cost dictionnaries (abs_cache and precompute_cache) estimates from server 
    // or pre-established values if no server available
		getAbsAndPreCachesFromServer(s, crypto);
    // Add to cost dictionnaries (encryptcache, decryptcache) cost estimates (pot. precomputed)
    getPreComputedData(crypto);

    // Iterate on all the parameters for this crypto module
    for (std::string crypto_params_desc : crypto_params_desc_set)
    {
      // Set the current public parameters to the parameter set values
      currentPublicParams->computeNewParameters(crypto_params_desc);
      // Iterate on the dimensions
      for (unsigned int d = fixedVars.dMin ; d  <= fixedVars.dMax ; d++) 
      {
        // Define a OptimVars object (with the variables that we can choose for optimization)
        OptimVars vars(fixedVars.fitness, fixedVars);
        vars.crypto_params = crypto_params_desc;
        // Get best alpha for this iteration given the fitness method and set the variables in vars accordingly
        // - First do a test with no aggreggation
        findAlpha(d, 1, 1, vars, exp_nbr);
        // - Then search with dichotomy 
        findAlphaDicho(0, getMaxAlpha(), d, exp_nbr, vars);
      }
    }
    // Show and write the optimized values.
    processResults(exp_nbr);
    
    // Save best optimum inter cryptosystems
    std::cout << "Optimizer: Comparing (" << optimum.crypto_params << ", " << optimum.getValue() << ") with (" << global_optimum.crypto_params << ", " << global_optimum.getValue() << ")" << std::endl;
    if(optimum < global_optimum)
    {
      global_optimum = optimum;
    }
    std::cout << "Optimizer: Winner is (" << global_optimum.crypto_params << ", " << global_optimum.getValue() << ")" << std::endl;

    // Clean
    optimum.reset();
    crypto_params_desc_set.clear();
    delete crypto;
  }
  disconnect();
  
  // If nothing was found say it and exit
  if (global_optimum.crypto_params == "No crypto params")
  {
    std::cout << "Optimizer: No valid crypto parameters matching " << fixedVars.manual_crypto_params << " found, exiting ..." << std::endl;
    exit(0);
  }
}

/**
 *  Set optimum with values giving by a trivial PIR. (i.e download the entiere database) 
 **/
void PIROptimizer::noCryptoSystemOptim(unsigned int exp_nbr)
{
  optimum.crypto_params = "NoCryptography";
  optimum.setFixedVars(fixedVars);
  optimum.setGenQ(0);
  optimum.setSendQ(0);
  optimum.setGenR(0);
  optimum.setSendR((fixedVars.l * fixedVars.n) / fixedVars.Tdoc);
  optimum.setDecR(0);
  optimum.setDim(1);
  optimum.setAlpha(fixedVars.n);
  if (silent == false) showBestResults(exp_nbr);
  
  // Save this as the best result inter cryptosystems
  global_optimum = optimum;

  optimum.reset();
}


/**
 *  Find Iteratively the best alpha and set the optimal variables in var
 **/
void PIROptimizer::findAlpha(unsigned int d, OptimVars& vars, unsigned int exp_nbr)
{
  // alphaMax == 0 is a convention to note that no limit should be done in aggregation
  unsigned int alpha_max = (fixedVars.alphaMax == 0) ? fixedVars.n : fixedVars.alphaMax;
  //Really get best alpha between 0 and alpha_max
  findAlpha(d, 0, alpha_max, vars, exp_nbr);
}

void PIROptimizer::findAlpha(unsigned int d, unsigned int alpha_min, unsigned int alpha_max, OptimVars& vars, unsigned int exp_nbr)
{
  // Number of bits that can be absorbed by a ciphertext for our parameters
  long abs_bit;
  // If aggregation is done, try to use all plaintext space
  unsigned int minimum_reasonable_alpha = getMinAlpha(); 
  // Iterate on possible aggregation values
  for (unsigned int current_alpha = alpha_min ; current_alpha <= alpha_max ; current_alpha+=minimum_reasonable_alpha)
  {
    // In first iteration current_alpha=0 is treated as 1 (i.e. no aggregation)
    computeOptimVars(current_alpha, d, vars); 
    
    // If no bits can be absorbed for this choice ignore it
    abs_bit = crypto->getPublicParameters().getAbsorptionBitsize();
    if(abs_bit <= 0)
    {
      cout << "PIROptimizer: Unusable cryptoparams, ignoring them" << endl; 
      break;
    }
    // Write test result to output file
    OptimService::writeTestCurrentResult(1, alpha_max, current_alpha, 1, alpha_max, d, exp_nbr, vars);
  }
}

void PIROptimizer::findAlphaDicho(unsigned int inf, unsigned int sup, unsigned int d, unsigned int exp_nbr, OptimVars& vars)
{
  unsigned int min_alpha = getMinAlpha();

  //When the space to explore is relatively little we use an iterative method to steer clear of local optima.
  if ((sup - inf) <= 100)
  {
    unsigned int max_alpha_bound = 0;
		if (fixedVars.n < (sup * min_alpha))
			max_alpha_bound = fixedVars.n;  

		else 
			max_alpha_bound = min(fixedVars.alphaMax, max(kalphaBound, sup * min_alpha)); 

		findAlpha(d, inf * min_alpha, max_alpha_bound, vars, exp_nbr);
		return;
  }

  //Next value to test.
  unsigned int alpha_bound = inf + (sup - inf) / 2.0;
  //Return the state of the slope (up or down).
  int val = slop(alpha_bound, inf, sup, d, vars, exp_nbr);
  //Write test restult to output file.
  OptimService::writeTestCurrentResult(1, sup*getMinAlpha(), alpha_bound*getMinAlpha(), inf*getMinAlpha(), sup*getMinAlpha(), d, exp_nbr, vars);

  if(val == -1) {  
    findAlphaDicho(alpha_bound, sup, d, exp_nbr, vars);
    return;
  }
  findAlphaDicho(inf, alpha_bound, d, exp_nbr, vars);
}


int PIROptimizer::slop(unsigned int alphaMul, unsigned int inf, unsigned int sup,  unsigned int d, OptimVars& vars, unsigned int exp_nbr)
{
  // Try ten uniformly spaced values ahead to try to find a better result
  unsigned int deltaright = (sup-1-alphaMul)/10;
  unsigned int deltaleft = (alphaMul - inf-1)/10;
  unsigned int minalpha = getMinAlpha();
  double vright, vleft, vtmp;
  // Ten increasing tests take the best
  // ten decreasings tests take the best
  // compare the two bests to determine the slope
  // each side must have alphaMul+-1 at least
  vright = computeOptimVars((alphaMul + 1) * minalpha, d, vars); 
  for (unsigned int i = 1 ; i <= 10 ; i++)
  {
    vtmp = computeOptimVars((alphaMul + 1 + i * deltaright) * minalpha, d, vars);
    if (vtmp <= vright) vright = vtmp; 
  }
  vleft = computeOptimVars((alphaMul - 1) * minalpha, d, vars); 
  for (unsigned int i = 1 ; i <= 10 ; i++)
  {
    vtmp = computeOptimVars((alphaMul - 1 - i * deltaleft) * minalpha, d, vars);
    if (vtmp <= vleft) vleft = vtmp; 
  }

  return (vright < vleft - 10e-8) ? -1 : 1 ;
}


/**
 * Compute dimension sizes for a d-dimensional hypercube of n elements.
 * Params:
 * 	- unsigned int n: the number of elements
 * 	- unsigned int alpha: the aggregation value
 * 	- unsigned int d: the dimension
 * 	- unsigned int *dn: computed dimension sizes (output) 
 **/
void PIROptimizer::getDimSize(unsigned int n, unsigned int alpha, unsigned int d, unsigned int *dn)
{

  unsigned int prod = 1, j = 0;

  // Elements in the database after the aggregation 
  unsigned int new_n = ceil(static_cast<double>(n)/static_cast<double>(alpha)); //PAtch tres sale à reprendre
  // Dimension sizes. Product must be > n/alpha
  unsigned int factors[d];

  // Lower bound on the dimensions sizes needed. Correct only if n/alpha is a d-power.
  for (unsigned int i = 0 ; i < d ; i++) factors[i] = floor(pow(new_n,1./d));

  // If n/alpha is not a d-power 
  if (static_cast<double>(factors[0]) != pow(new_n, static_cast<double>(1.0 /  static_cast<double>(d))) )
  {
    // Increment each factor by one until the product of the dimension sizes is above n/alpha
    while (prod < new_n && j < d)
    {
      prod = 1;
      factors[j++]++;

      for (unsigned int i = 0 ; i < d ; i++)
        prod *= factors[i];
    }
  }

  // Copy the result to the output
  memcpy(dn, factors, sizeof(factors));
}


/**
 *	Compute time needed to send the request.
 * 	Params:
 * 		- double Tupc : client upload throughput
 * 		- double Tdos : server download throughput
 * 		- unsigned int nbc : number of clients that share server throughput (currently always 1);
 * 		- unsigned int d   : number of dimensions used for recursion
 *
 * 	Return:
 * 		- double : time in seconds to send the PIR query.
 **/
double PIROptimizer::getSendCost(double Tupc, double Tdos, unsigned int nbc, unsigned int d, unsigned int* dn)
{
  // Available throughput for the transfer
  double min_throughput = min(Tupc, Tdos/nbc );
  // Time needed to send the query
  double SenQ = 0.0;

  // Sum the times needed for the subqueries of each dimension in the recursion 
  for (unsigned int i = 0 ; i < d ; i++)
  {
    // For a given dimension time = size / throughput and size = (nb_elements * size_of_an_element)
    SenQ +=	(dn[i] * Tcq(i)) / min_throughput;
  }
  //std::cout << "dn[0]="<<dn[0]<<" d=" << d << " Tcq(0)=" << Tcq(0) << " Tupc=" << Tupc << " Tdos=" << Tdos << " SenQ=" << SenQ << std::endl;

  return SenQ;
}


/**
 *	Compute reply generation cost.
 *	Params:
 * 		- unsigned int d   : number of dimension used for recursion
 * 		- unsigned int Tmi : plaintext size on the i-th level of recursion
 *
 * 	Return:
 * 		- double : time in second to generate a reply.
 **/
double PIROptimizer::getReplyGenCost(unsigned int alpha, unsigned int d, unsigned int *dn)
{
  long double prod, GenR = 0;
  uint64_t plaintext_nbr_post_ntt, plaintext_nbr_pre_ntt;
  bool initialprecomputationdone = INITIAL_PRECOMPUTATION_DONE;
  
  // In order to compute absorption size we must define it first
  crypto->setandgetAbsBitPerCiphertext(dn[0]);

  for (unsigned int i = 0 ; i < d ; i++)
  {
    // Compute how many elements has the intermediate database
    prod = 1;
    for (unsigned int j = i ; j < d ; j++)
    { 
      prod *= dn[j];
    }

    // Compute how large they are
    plaintext_nbr_post_ntt = ceil((double) eltSize(alpha,i) / (double) crypto->getPublicParameters().getAbsorptionBitsize(i)); 
    plaintext_nbr_pre_ntt = ceil((double) eltSize(alpha,i) / (double) crypto->getPublicParameters().getCiphertextBitsize()/2); 

    
    // Consider absorption cost
    GenR += prod * plaintext_nbr_post_ntt * getAbs1PlaintextTime(crypto);
    
    // And potentially precomputation cost
    if (i > 0 || initialprecomputationdone == false)
    {
      GenR += prod * plaintext_nbr_pre_ntt * getPrecompute1PlaintextTime(crypto);
    }
  }
  return GenR;
}

/**
 *	Compute the size of a element at dimension "i".
 *	Params :
 *		- unsigned int i   : dimension ;
 *		- unsigned int Tmi : encrypted size.
 *
 *	Return :
 *		- double : the size of an element at dimension "i".
 **/
double PIROptimizer::eltSize(unsigned int alpha, unsigned int i)
{
  if (i == 0)
  {
    return static_cast<double>(fixedVars.l * alpha); 
  }
  return ceil( static_cast<double>(eltSize(alpha, i - 1) /  Tp(i-1) )) * Tcr(i-1); 
}

double PIROptimizer::Tp(unsigned int i)
{
  return currentPublicParams->getAbsorptionBitsize(i);
}

/**
 *	Compute decryption cost.
 *	Params :
 *		- unsigned int d   : dimension ;
 *		- unsigned int Tmi : encrypted size.
 *
 *	Return :
 *		- double : time in seconds to decrypt something with given paramtes.
 **/
double PIROptimizer::getDecryptCost(unsigned int alpha, unsigned int d)
{
  double DecR = 0;

  for (unsigned int i = 0 ; i < d ; i++)
  {
    // Yet another strange formula
    DecR += ceill((eltSize(alpha, i + 1) / Tcr(i))) * getDecCost();
  }

  return DecR;
}

/**
 * Return the size of reply for dimension "i"
 **/
inline double PIROptimizer::Tcr(unsigned i)
{
  return currentPublicParams->getCiphBitsizeFromRecLvl(i+1);
}

/**
 * Return the size of a query for dimension "i"
 **/
inline double PIROptimizer::Tcq(unsigned i)
{
  return Tcr(i);
}

/**
 * Return the time in secondes to decrypt data  at dimension "i"
 **/
double PIROptimizer::getDecCost(unsigned int d, crypto_ptr crypto)
{
  PIRParameters pirParams; //Works for d = 1 
  pirParams.d = 1;
  pirParams.alpha = 1;
  unsigned int chunks = 1;
  double start, stop, elapsed_time = 0;

  // Needed for calls to getAbsorptionBitsize
  crypto->setandgetAbsBitPerCiphertext(1); // Set internally best absorption possible 
  
  // Use a special function to generate a ciphertext to be sure its decryption cost is average
  char* encrypted_data = (char*) crypto->encrypt_perftest();
  PIRReplyExtraction_internal replyExt(pirParams,*crypto); 

  shared_queue<char*> clearChunks("clear_chunks");

  do{
    chunks *= 2;
    for(unsigned int i = 0 ; i < chunks ; i++) //Fill, reply buffer with copy of "encrypted_data".
    {
      char* encrypted_data_copy = (char*) malloc((crypto->getPublicParameters().getCiphertextBitsize()/8) * sizeof(char));
      memcpy(encrypted_data_copy, encrypted_data, crypto->getPublicParameters().getCiphertextBitsize()/8);
      replyExt.repliesBuffer.push(encrypted_data_copy);
    }


    start = omp_get_wtime();
    replyExt.extractReply((currentPublicParams->getAbsorptionBitsize()/8) * chunks, &clearChunks); //Do a PIR data extraction.
    stop  = omp_get_wtime();

    while(!clearChunks.empty())
      free(clearChunks.pop_front()); //free clear data.

  }while((elapsed_time = (stop - start)) < 0.20);

  double result = elapsed_time / chunks;

  return result;
}

double PIROptimizer::getDecCost()
{
  bool shortversion=true;
  return decrypt_cache[currentPublicParams->getSerializedParams(shortversion)];
}

/**
 *	Compute SendR.
 *	Params :
 *		- double Tups : server upload time/bit ;
 *		- double Tdoc : client download time/bit;
 *		- unsigned int nbc : number of client ;
 *		- unsigned int Tmi : encrypted size.
 *	Return :
 *	  - double : Duration in seconds.
 **/
double PIROptimizer::getReceiveCost(unsigned int alpha, double Tups, double Tdoc, unsigned int nbc, unsigned int d)
{
  // This is a bad clone of eltSize. To be removed.
  //double SendR = eltSize(alpha, 0);
  //
  //for (unsigned int i = 0 ; i < d; i++)
  //{
  //  SendR = ceil(SendR / Tp(i)) * Tcr(i);
  //}
  //SendR /= min(Tups/nbc, Tdoc); 
  //
  //return SendR;

  // eltSize computes how db elements grow with recursion. 
  // Reply size is what would be the element size for the d-th recursion (rec levels start at 0)
  return eltSize(alpha,d)/min(Tups/nbc, Tdoc);
}

/**
 *	Get absorption time for 1 bit.
 *	Param :
 *		- unsigned int d : dimension.
 **/
void PIROptimizer::getAbsAndPreCachesFromServer(boost::asio::ip::tcp::socket& s, crypto_ptr crypto_system)
{
  size_t cmd = ABS;
  try
  {

    write(s, boost::asio::buffer(&cmd, sizeof(size_t)));

    std::string crypto_name = crypto_system->toString();
    unsigned int crypto_name_size = crypto_name.size();
    std::cout << "Optimizer: Requesting absorption and precomputation costs file for " << crypto_name << std::endl;

    write(s, boost::asio::buffer(&crypto_name_size, sizeof(int)));
    write(s, boost::asio::buffer(crypto_name));

    int file_content_size = 0;
    read(s, boost::asio::buffer(&file_content_size, sizeof(int)));

#ifdef DEBUG
    std::cout << "Optimizer: File contents size " << file_content_size << std::endl;
#endif 

    char file_content[file_content_size + 1];
    file_content_size = read(s, boost::asio::buffer(file_content, file_content_size));
    file_content[file_content_size] = '\0';

    makeAbsAndPrecomputeCaches(file_content);
  }
  catch(std::exception const& ex)
  {
    cout << "Optimizer: No server, using reference values for " << crypto_system->toString() << endl;
    abs_cache.clear(); 
    precompute_cache.clear(); 

    set<string> crypto_params_set;
    crypto_system->getAllCryptoParams(crypto_params_set);

    for (auto crypto_param : crypto_params_set)
    {
      abs_cache[crypto_param] = crypto_system->estimateAbsTime(crypto_param);
      precompute_cache[crypto_param] = crypto_system->estimatePrecomputeTime(crypto_param);
    }
  }
}

void PIROptimizer::makeAbsAndPrecomputeCaches(char* serialized_cache)
{
  abs_cache.clear();
  precompute_cache.clear();
  std::vector<string> lines;
  boost::algorithm::split(lines, serialized_cache, boost::algorithm::is_any_of("\n"));

  for (auto line : lines)
  {
    std::cout << line << std::endl;
    vector<string> fields;
    boost::algorithm::split(fields, line, boost::algorithm::is_any_of(" "));

    if (fields.size() == 3)
    {
      abs_cache[fields.at(0)] = atof(fields.at(1).c_str());
      precompute_cache[fields.at(0)] = atof(fields.at(2).c_str());
    }
    else if (line.find_first_not_of(' ') != std::string::npos)// if it is not a blank line
    {
      std::cout << "PIROptimizer: Absorption and precompute cache corrupted" << std::endl;
      std::cout << "PIROptimizer: Remove files in exp/*.abs in the server before proceeding" << std::endl;
      exit(1);
    }
  }
}

/**
 *	Compute the time for generate the complete query.
 *	Params :
 *		 - unsigned int d : dimension ;
 *		 - crypto_ptr	crypto : HomomorphicCrypto shared_ptr ;
 *
 *	Return :
 *		- double : Complete query duration.
 **/
double PIROptimizer::getGenQueryCost(unsigned int d, unsigned int *dn)
{
  double GenQ = 0.0;
  bool shortversion = true;
  string current_crypto_params = currentPublicParams->getSerializedParams(shortversion);

  for (unsigned int i = 0 ; i < d ; i++)
  {
    GenQ += dn[i]  * encrypt_cache[current_crypto_params];
  }

  return GenQ;
}

/**
 *	Get time to encrypt a query for a given dimension.
 * 	Params :
 * 		- unsigned int d : dimension ;
 *		- crypto_ptr	crypto : HomomorphicCrypto shared_ptr ;
 *
 *	Return :
 *		- double : Duration in seconds.
 **/
double PIROptimizer::getQueryElemGenCost(unsigned int d, crypto_ptr crypto)
{
  double start, stop; 
  double elapsed_time = 0;
  unsigned int query_elts_nbr = 0;
  PIRParameters pirParams;

  pirParams.d = d;
  pirParams.alpha = 1;

  for (unsigned int i = 0 ; i < d; i++)
    pirParams.n[i] = 1;

  // Needed for calls to getAbsorptionBitsize
  crypto->setandgetAbsBitPerCiphertext(1); // Set internally best absorption possible 
  
  PIRQueryGenerator_internal queryGenerator(pirParams, *crypto); 
  queryGenerator.setChosenElement(1);

  do{
    for (unsigned int i = 0 ; i < d; i++) pirParams.n[i] *= 2;

    start = omp_get_wtime();
    queryGenerator.generateQuery();
    stop = omp_get_wtime();

    queryGenerator.cleanQueryBuffer();

  }while((elapsed_time = (stop - start)) <= 0.5);

  for(unsigned int i = 0 ; i < d ; i++) query_elts_nbr += pirParams.n[i];

  double result = elapsed_time / query_elts_nbr;
  return result;
}

const PIRParameters& PIROptimizer::getParameters()
{
  return pirParameters;
}

double PIROptimizer::getAbs1PlaintextTime(HomomorphicCrypto* crypto_system)
{
  bool shortversion = true; 
  return abs_cache[currentPublicParams->getSerializedParams(shortversion)];
}


double PIROptimizer::getPrecompute1PlaintextTime(HomomorphicCrypto* crypto_system)
{
  bool shortversion = true; 
  return precompute_cache[currentPublicParams->getSerializedParams(shortversion)];
}


/**
 *	Return natural aggregation value.
 **/
unsigned int PIROptimizer::getMinAlpha()
{
  double r = ceil(Tp(0) / static_cast<double>(fixedVars.l));
  r = (r > fixedVars.n) ? fixedVars.n : r; 
  unsigned ret = (r < 1.0) ? 1 : unsigned(r);

	/*No limit for agregation*/
	if (fixedVars.alphaMax == 0)
		return ret;

	else
		return min(ret, fixedVars.alphaMax);
}

unsigned int PIROptimizer::getMaxAlpha()
{
  return ((fixedVars.alphaMax == 0) || (fixedVars.n < fixedVars.alphaMax)) ? fixedVars.n : fixedVars.alphaMax;
}


/**
 * Given the fixed vars and the choices done in the optimize loop, estimate costs.
 * If total cost is better than the previus optimum, replace it.  
 * Return total cost of the PIR retrieval (given the fitness method in vars).
 **/
double PIROptimizer::computeOptimVars(unsigned int alpha, unsigned int d, OptimVars& vars)
{
  if (alpha == 0) alpha = 1;

  unsigned int dn[fixedVars.dMax];

  // Compute dimension sizes given the number of files, aggregation, and recursion levels
  getDimSize(fixedVars.n, alpha, d, dn);
  //Save alpha and d values.
  vars.setAlpha(alpha);
  vars.setDim(d);
  
  crypto->setandgetAbsBitPerCiphertext(dn[0]);


  // Get costs for query generation, emission, reply generation, reception, and decryption
  vars.setGenQ(getGenQueryCost(d, dn));
  vars.setSendQ(getSendCost(fixedVars.Tupc, fixedVars.Tdos, nbc, d, dn));
  vars.setGenR(getReplyGenCost(alpha, d, dn)); 
  vars.setSendR(getReceiveCost(alpha, fixedVars.Tups, fixedVars.Tdoc, nbc, d)); 
  vars.setDecR(getDecryptCost(alpha, d)); 

  // Decide whether this test is better than the optimum. If so replace it.
  analyseResult(alpha, d, vars); 

  // Return total cost of the PIR retrieval (given the fitness method in vars)
  return vars.getValue();
}


void PIROptimizer::processResults(unsigned int i)
{
  if (silent == false) showBestResults(i);
  writeTestBestResult(i);
}

/**
 *	Test current result to get the optimum.
 **/
void PIROptimizer::analyseResult(unsigned int alpha, unsigned int d, OptimVars& vars)
{
  bool shortversion = true;

  // < operator is defined in OptimVars object, compares total cost value
  if(vars < optimum)
  {
    // If vars is better than the current optimum redefine it
    optimum = vars;

    // Get a final version of the parameters with the absorption size
    optimum.crypto_params = crypto->getSerializedCryptoParams(!shortversion);
#ifdef DEBUG
    cout << "PIROptimizer: New optimum cryptoparams " << optimum.crypto_params << endl;
#endif
  }
}

void PIROptimizer::getFixedVarsFromServer(FixedVars& fixed_vars, boost::asio::ip::tcp::socket& s)
{
  try{
    /*Send command*/
    size_t cmd = DATA;
    write(s, boost::asio::buffer(&cmd, sizeof(cmd)));
    /*Get number of files*/
    read(s, boost::asio::buffer(&fixed_vars.n, sizeof(fixed_vars.n)));
    cout << "Optimizer: Number of files is " << fixed_vars.n << endl;

    /*Get max file size*/
    read(s, boost::asio::buffer(&fixed_vars.l, sizeof(fixed_vars.l)));
    cout << "Optimizer: File bytesize is " << fixed_vars.l << endl;
    // fixed_vars.l should be in bits
    fixed_vars.l*=8;
  }catch(const std::exception& ex)
  {
    std::cerr << "Catched exception in " << __FUNCTION__ << ": " << ex.what() << std::endl;
  }
}

/**
 *	Print the best result.
 **/
void PIROptimizer::showBestResults(unsigned int i)
{
  writeBigMessage(" RESULTS exp " + to_string(i+1) + " ");
  cout << "\tparams : " << optimum.crypto_params  << endl;
  cout << "\td      : "	<< optimum.getDim()	      << endl;
  cout << "\talpha  : "	<< optimum.getAlpha()     << endl;
  cout << "\tGenQ   : "	<< optimum.getGenQ()      << " sec" <<  endl;
  cout << "\tSenQ   : " << optimum.getSendQ()     << " sec" <<  endl;
  cout << "\tGenR   : "	<< optimum.getGenR()      << " sec" <<  endl;
  cout << "\tSendR  : " << optimum.getSendR()	    << " sec" <<  endl;
  cout << "\tDecR   : "	<< optimum.getDecR() 	    << " sec" <<  endl;
  cout << "\tTotal Time : "<< optimum.getValue()  << " sec" <<  endl;
  cout << "\t--Fixed Vars-- "                     << endl;
  cout << "\tl      : "  << fixedVars.l << " bits"<< endl;
  cout << "\tn      : "  << fixedVars.n           << endl;
  cout << "\tTupc/dos : "<< fixedVars.Tupc << " b/s" << endl;
  cout << "\tTdoc/ups : "<< fixedVars.Tdoc << " b/s" << endl;

  writeBigMessage("##############");
}

void PIROptimizer::getPreComputedData(HomomorphicCrypto* crypto)
{
  std::string	fdec_path(OptimService::folderName + OptimService::fileName + crypto->toString() 
      + OptimService::decFileExtension);
  std::string	fenc_path(OptimService::folderName + OptimService::fileName + crypto->toString() 
      + OptimService::encFileExtension);

  decrypt_cache.clear();
  encrypt_cache.clear();

  set<string> crypto_params_set;
  crypto->getAllCryptoParams(crypto_params_set);

  //if (OptimService::verifyOptimData(crypto_params_set, fdec_path, fenc_path))
  if (OptimService::fileOutdated(crypto->toString(), OptimService::encFileExtension) || 
      OptimService::fileOutdated(crypto->toString(), OptimService::decFileExtension) )
  {
    cout << "PIROptimizer:: Computing cache data for " <<  crypto->toString() << ", this can take a while..." << endl;
    computeOptimData(crypto_params_set);
    cout << "PIROptimizer: Finished computing cache data for "<< crypto->toString() << endl;
  }
  else
  {
    OptimService::readOptimData(decrypt_cache, fdec_path);
    OptimService::readOptimData(encrypt_cache, fenc_path);
  }
}

void PIROptimizer::computeOptimData(set<string> &crypto_params_set )
{
  unsigned int i = 1;
  unsigned int crypto_params_nbr = crypto_params_set.size();
  string optim_data2write;

  // Generate the encryption performance cache
  std::string file_path(OptimService::folderName + OptimService::fileName + crypto->toString());
  for (auto crypto_param : crypto_params_set)
  {
    cout << "PIROptimizer: Encrypt cache generation for " << crypto_param << " " << i++ 
      << "/" << crypto_params_nbr << "." << flush;
    crypto->setNewParameters(crypto_param);
    cout << "." << endl;

    encrypt_cache[crypto_param] = getQueryElemGenCost(1, crypto);
    
    std::ostringstream out;
    out << std::setprecision(kPrecision) << encrypt_cache[crypto_param];
    optim_data2write += crypto_param + " " + out.str() + "\n";
  }	
  if(OptimService::writeOptimDataBuffer(optim_data2write, file_path+OptimService::encFileExtension)) 
  {
    std::cout << "PIROptimizer: Error when writing optimization data, aborting." << std::endl;
    exit(1);
  }
  optim_data2write = "";
  i=1;

  // Generate the decryption performance cache
  for (auto crypto_param : crypto_params_set)
  {
    cout << "PIROptimizer: Decrypt cache generation for " << crypto_param << " " << i++ 
      << "/" << crypto_params_nbr << "." << flush;
    crypto->setNewParameters(crypto_param);
    cout << "." << endl;

    decrypt_cache[crypto_param] = getDecCost(1, crypto);
    
    std::ostringstream out;
    out << std::setprecision(kPrecision) << decrypt_cache[crypto_param];
    optim_data2write += crypto_param + " " + out.str() + "\n";
  }	
  if(OptimService::writeOptimDataBuffer(optim_data2write, file_path+OptimService::decFileExtension)) 
  {
    std::cout << "PIROptimizer: Error when writing optimization data, aborting." << std::endl;
    exit(1);
  }

}

/**
 *	Writes values of the best result.
 **/
void PIROptimizer::writeTestBestResult(unsigned int i)
{
  OptimService::writeMessage(i, "### BEST RESULT ###");
  OptimService::writeFootFile(i);
  OptimService::writeConfigFile(getMinAlpha(),  optimum.getAlpha(), optimum.getDim(), i);
}

/**
 *  Disconect to the sever's optimizer.
 **/
void PIROptimizer::disconnect()
{
  try
  {
    size_t cmd = EXIT;
    s.send(boost::asio::buffer(&cmd, sizeof(size_t)));// exit function
    s.close();
  }catch (const std::exception& ex)
  {}
}

void PIROptimizer::connect()
{
  try
  {
    tcp::endpoint endpoint(boost::asio::ip::address::from_string(serv_ip), server_optimport);
    s.connect(endpoint);
  }catch(const std::exception& ex )
  {

  }

}

/**
 *	Print "big" message on the terminal.
 **/
void PIROptimizer::writeBigMessage(string msg)
{
  std::transform(msg.begin(), msg.end(), msg.begin(), ::toupper);
  cout << "########" << msg << "########" << endl;
}

PIROptimizer::~PIROptimizer()
{
}

/**
 * Do a network speed test with the server.
 * Params :
 * 	- int port : server port
 **/
void PIROptimizer::getNetworkInfos(FixedVars& fixedVars, boost::asio::ip::tcp::socket& s)
{
  double start = 0, end=0, loop = 1;
  char *msg = (char *) malloc(MEGA_BYTE+4);
  unsigned int read_size = 0;
  memset(msg, 1, MEGA_BYTE);

  size_t cmd = SPEED;

  /*Send command*/
  write(s, boost::asio::buffer(&cmd, sizeof(size_t)));

  // Do upload test
  start = omp_get_wtime();
  for (int i = 0 ; i < loop ; i++) write(s, boost::asio::buffer(msg, MEGA_BYTE));
  end = omp_get_wtime();
  
  // Ignore it if value was forced
  if (fixedVars.Tupc != 0)
  {
    cout << "Optimizer: Upload speed value forced, dropping network test result" << std::endl;
  }
  else
  {
    fixedVars.Tupc = (loop /(end - start)) * 8 * MEGA_BYTE;
    fixedVars.Tdos = fixedVars.Tupc;  

    cout << "Optimizer: Upload speed test gives " << fixedVars.Tupc << " bits/s" << endl;
  }

  // Do download tests
  start = omp_get_wtime();
  for (int i = 0 ; i < loop ; i++) read(s, boost::asio::buffer(msg, MEGA_BYTE));
  end = omp_get_wtime();

  if (fixedVars.Tdoc != 0)
  {
    cout << "Optimizer: Download speed value forced, dropping network test result" << endl;
  }
  else
  {
    fixedVars.Tdoc =  loop / (end - start) * 8 * MEGA_BYTE;
    fixedVars.Tups = fixedVars.Tdoc;

    cout << "Optimizer: Download speed test gives " << fixedVars.Tdoc << " bits/s" << endl;
  }
  free(msg);
}

void PIROptimizer::sendExitCommand(boost::asio::ip::tcp::socket& s)
{
  try{
    size_t cmd = EXIT;
    write(s, boost::asio::buffer(&cmd, sizeof(cmd)));
  }catch(std::exception const& ex)
  {
    cerr << "Error when sending EXIT command to the server: " << ex.what() << endl;
  }
}
