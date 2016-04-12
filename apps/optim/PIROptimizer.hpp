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

#ifndef DEF_PIROPTIMIZER
#define DEF_PIROPTIMIZER

#include <math.h>
#include <map>
#include <limits>
#include <iomanip>
#include <string>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <regex>
#include <boost/asio.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

#include "pir/optim/OptimVars.hpp"
#include "pir/optim/OptimService.hpp"
#include "pir/PIRParameters.hpp"
#include "pir/queryGen/PIRQueryGenerator_internal.hpp"
#include "pir/replyExtraction/PIRReplyExtraction_internal.hpp"
#include "crypto/HomomorphicCryptoFactory_internal.hpp"
#include "crypto/NFLLWEPublicParameters.hpp"
#include "pir/ClientDefines.hpp"
#include "pir/PIROptimizerCommand.hpp"

#ifndef MODULBITS
#define MODULBITS 1024
#endif

#define MAX_D 10

#define MEGA_BYTE 1000000

using namespace std;
using boost::asio::ip::tcp;

typedef HomomorphicCrypto* crypto_ptr;

class PIROptimizer {

	protected:
    FitnessType fitnessMethod;
		PIRParameters pirParameters;
		HomomorphicCrypto* crypto;
    AbstractPublicParameters* currentPublicParams;
		string serv_ip;
    int server_optimport;
    boost::asio::io_service ios;
		boost::asio::ip::tcp::socket s;
		ofstream file;
    map<string, double> encrypt_cache;
    map<string, double> decrypt_cache;
    map<string, double> abs_cache;
    map<string, double> precompute_cache;

		/** Fixed Parameters **/
    FixedVars fixedVars;
		unsigned int nbc; /** number of client **/
		unsigned int dn[MAX_D];
    unsigned int security_bits;
    double ax[MAX_D], ayDec[MAX_D], ayEnc[MAX_D], axAbs[MAX_D], ayAbs[MAX_D];

		/** Free PIR Parameters **/
		string pir;
		/**Outputs**/
		OptimVars optimum, global_optimum;
		double Tcr(unsigned i);
		double Tcq(unsigned i);
    double Tp(unsigned i);
		
		void getAbsAndPreCachesFromServer(boost::asio::ip::tcp::socket& s, crypto_ptr crypto_system);
    void makeAbsAndPrecomputeCaches(char* serialized_cache);
		double getAbs1PlaintextTime(HomomorphicCrypto* crypto_system);
		double getPrecompute1PlaintextTime(HomomorphicCrypto* crypto_system);
		double getQueryElemGenCost(unsigned int d, crypto_ptr crypto);
		double eltSize(unsigned int alpha, unsigned int s);
		double getDecCost(unsigned int i, crypto_ptr);
		double getDecCost();
		
		void analyseResult(unsigned int alpha, unsigned int d, OptimVars& vars);
		void showBestResults(unsigned int i);

		void writeTestBestResult(unsigned int i);

		void getPreComputedData(HomomorphicCrypto* crypto);
		void computeOptimData(set<string>& crypto_params_set);

    void sendExitCommand(boost::asio::ip::tcp::socket& s);
		void disconnect();
		void connect(); //const throw(std::exception);


	public:
    bool silent;
		PIROptimizer(string server_ip, int port, FitnessType fitnessMethod);
		PIROptimizer();
	  virtual	~PIROptimizer();
    void optimize(boost::asio::ip::tcp::socket& s, FixedVars& fixed_vars, unsigned int exp_nbr);
    static const unsigned int kalphaBound;
    void optimize(FixedVars& fixed_vars, unsigned int exp_nbr);
    OptimVars optimizeFromServer(FixedVars& partial_fixed_vars,boost::asio::ip::tcp::socket& socket);
    void optimizeFromConfiguredCryptosystem(FixedVars& fixed_vars, HomomorphicCrypto* crypto_system);
    void noCryptoSystemOptim(unsigned int exp_nbr);
    void findAlpha(unsigned int d, OptimVars& vars, unsigned int exp_nbr);
    void findAlpha(unsigned int d, unsigned int alpha_min, unsigned int alpha_max, OptimVars& vars, unsigned int exp_nbr);

    void findAlphaDicho(unsigned int inf, unsigned int sup, unsigned int d, unsigned int exp_nbr, OptimVars& vars);
    int slop(unsigned int alphaMul, unsigned int inf, unsigned int sup, unsigned int d, OptimVars& vars, unsigned int exp_nbr);

		const PIRParameters& getParameters();

		static void getDimSize(unsigned int n, unsigned int alpha, unsigned int d, unsigned int *dn);
    unsigned int getMaxAlpha();
    unsigned int getMinAlpha();

    double getGenQueryCost(unsigned int  d, unsigned int *dn);
		double getSendCost(double Tupc, double Tdoc, unsigned int nbc, unsigned int d, unsigned int *dn);
		double getReplyGenCost(unsigned int alpha, unsigned int d, unsigned int *dn);
		double getReceiveCost(unsigned int alpha, double Tups, double Tdoc, unsigned int nbc, unsigned int d);
		double getDecryptCost(unsigned int alpha, unsigned int d);
    void getFixedVarsFromServer(FixedVars& fixed_vars0,boost::asio::ip::tcp::socket& socket);
    void getNetworkInfos(FixedVars& fixed_vars, boost::asio::ip::tcp::socket& s);

		static void writeBigMessage(string msg);
    void processResults(unsigned int i);
    double computeOptimVars(unsigned int alphalt, unsigned int d, OptimVars& vars);

};

#endif
