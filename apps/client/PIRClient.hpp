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

#ifndef DEF_PIRCLIENT
#define DEF_PIRCLIENT
//#define DEBUG
#include <string>
#include <iostream>

#include <boost/asio.hpp>
#include <boost/thread.hpp>
#include <boost/signals2.hpp>

#include "DESC.hpp"

#include "crypto/NFLLWE.hpp"

#include "pir/queryGen/PIRQueryGenerator_internal.hpp"
#include "pir/PIRParameters.hpp"
#include "pir/replyExtraction/PIRReplyExtraction_internal.hpp"
#include "pir/replyExtraction/PIRReplyWriter.hpp"

#include "apps/optim/PIROptimizer.hpp"

#include "pir/events/WriteEvent.hpp"
#include "pir/events/CatalogEvent.hpp"
#include "pir/events/MessageEvent.hpp"

#include "pir/PIRParametersExchangeMethods.hpp"
#include "pir/GlobalConstant.hpp"

typedef boost::signals2::signal<void (MessageEvent&)>   messageListener;
typedef boost::signals2::signal<void (CatalogEvent&)>   menuListener;
typedef boost::signals2::signal<void (WriteEvent&)>     writeListener;


struct ClientParams 
{
  // Networking parameters
  string server_ip;
  int port; 

  // Interface params
  bool autochoice;
  bool verboseoptim;
  bool dontwrite;

  // Dry-run mode
  bool dryrunmode;
};

class PIRClientSimple
{
  private:
    PIRReplyExtraction_internal *replyExt;
    ClientParams clientParams;
    FixedVars fixedVars;
    OptimVars optimum;
    HomomorphicCrypto* cryptoMethod;
    messageListener messageListeners; 
    menuListener 		menuListeners;
    writeListener		writeListeners;
    boost::asio::ip::tcp::socket socket_up;
  
    PIRReplyWriter replyWriter;
    PIRParameters pirParams;
    DESC catalog;
    short exchange_method;

    void uploadWorker(PIRQueryGenerator_internal& queryGen);
    void downloadWorker(PIRReplyExtraction_internal& replyExt);
    void sendPirParams();
    void rcvPirParams();

    void exitWithErrorMessage(string s1, string s2);
    void writeErrorMessage(string funcName, string message);
    void writeWarningMessage(string funcName, string message);

    bool no_pipeline_mode;
    uint64_t chosenElement;

  public:
    PIRClientSimple( boost::asio::io_service& ios, ClientParams params, FixedVars vars);
    ~PIRClientSimple();

    void connect();
    void chooseFile();
    void rcvPIRParamsExchangeMethod();
    void processPIRParams();
    void downloadCatalog();
    void optimize();
    void sendCryptoParams(bool paramsandkey);
    void rcvCryptoParams();
    void processCryptoParams(); 
    void startProcessQuery();
    void startProcessResult();

    void setChosenElement(uint64_t choice);
    void joinAllThreads();
    void no_pipeline(bool b);

    /*Add Observers*/
    boost::signals2::connection addMessageListener(messageListener::slot_function_type subscriber);
    boost::signals2::connection addMenuListener(menuListener::slot_function_type subscriber);
    boost::signals2::connection addWriteListener(writeListener::slot_function_type subscriber);
};

#endif
