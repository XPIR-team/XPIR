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

#ifndef DEF_PIRSESSION
#define DEF_PIRSESSION
#include <vector>
#include <iostream>
#include <dirent.h>
#include <boost/asio.hpp>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

#include "pir/dbhandlers/DBHandler.hpp"
#include "pir/PIRParameters.hpp"
#include "pir/replyGenerator/PIRReplyGeneratorFactory.hpp"
#include "crypto/HomomorphicCryptoFactory_internal.hpp"

#include "pir/PIRParametersExchangeMethods.hpp"
#include "pir/GlobalConstant.hpp"
#include "pir/shared_queue.hpp"

#define DEFAULT_DIR_NAME "db/"



#define MAX_ATTEMPT 10
#define DEFAULT_PORT 1234

#define RED "\033[31m"
#define BOLD "\033[1;30m"
#define ORANGE "\033[33m"
#define RESET_COLOR "\033[m"

using boost::asio::ip::tcp;

struct session_option_t {
  bool driven_mode;
  bool keep_database;
  bool got_preimported_database;
  imported_database_t data;
};

class PIRSession : public boost::enable_shared_from_this<PIRSession>
{
	private:
		DBHandler *dbhandler;
		tcp::socket sessionSocket;
		boost::thread upThread, downThread;
		bool handmadeExceptionRaised;
    	bool finished;

		PIRParameters 		 pirParam;
		CryptographicSystem* cryptoMethod;
		GenericPIRReplyGenerator*  generator;

		uint64_t maxFileBytesize; // La taille totale de la base de donnée en octets

		void serveForever();
		void rcvCryptoParams(bool paramsandkey);
    void sendCryptoParams(); 
    void sendPIRParamsExchangeMethod(short exchange_method);
		void rcvPirParams();
		void sendPirParams();
		void sendCatalog();
		void recieveKey();
		void startProcessQuery ();
		void startProcessResult(session_option_t session_option);
    bool rcvIsClient();

		void downloadWorker();
		void uploadWorker();
    imported_database_t savedDatabase;
     bool no_pipeline_mode;

	public:
    	typedef boost::shared_ptr<PIRSession> pointer;
    	static pointer create(boost::asio::io_service& ios);
		PIRSession(boost::asio::io_service& ios);
		~PIRSession();

    void setDBHandler(DBHandler *db);
		void writeErrorMessage(string, string);
		void writeWarningMessage(string, string);
		void exitWithErrorMessage(string, string);
		bool isFinished();
		tcp::socket& getSessionSocket();
		bool start(session_option_t session_option);
    PIRParameters getPIRParams();
    void setPIRParams(PIRParameters pir_parameters);
    imported_database_t getSavedDatabase();
    void no_pipeline(bool b);
};
#endif
