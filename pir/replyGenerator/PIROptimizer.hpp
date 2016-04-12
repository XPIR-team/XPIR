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

#ifndef DEF_PIROPTIM
#define DEF_PITOPTIM

#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include "crypto/HomomorphicCryptoFactory_internal.hpp"
#include "pir/optim/OptimService.hpp"
#include "apps/server/PIRServer.hpp"
#include "pir/PIROptimizerCommand.hpp"

#define COMMAND_AND_CONTROL_PORT 1235

#define MEGA_BYTE 1000000
#define MAX_D 10
using boost::asio::ip::tcp;

typedef boost::shared_ptr<boost::asio::ip::tcp::socket> socket_ptr;

class PIROptimizer
{
	private:
    std::vector<std::string> abs1bitCacheLines;
    boost::asio::io_service io_service;
		PIRParameters pirParameters;
    uint64_t fileCount;
    uint64_t maxFileBytesize;
		string pir;
		string filePath;
		void serve();

    double getAbs1PlaintextTime(HomomorphicCrypto* crypto_ptr, GenericPIRReplyGenerator* generator_ptr);
    double getPrecompute1PlaintextTime(HomomorphicCrypto* crypto_ptr, GenericPIRReplyGenerator* generator_ptr);
    void sendDatabaseInfos(boost::asio::ip::tcp::socket& s);
	public:
		PIROptimizer(DBHandler *db);
		~PIROptimizer();

		void speedTest(boost::asio::ip::tcp::socket& s);
    void controlAndCommand(boost::asio::ip::tcp::socket& s);
		unsigned int sendAbsAndPrecomputeCaches(boost::asio::ip::tcp::socket& s);
		double launcheAbs1BitBench();
    void processCommand(socket_ptr s, size_t& cmd);
		void optimize();

		void prepareOptimData();
    std::string computeOptimData(const std::string& crypto_name);
};
#endif
