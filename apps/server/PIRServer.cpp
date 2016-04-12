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

#include "PIRServer.hpp"
#include "PIRSession.hpp"

/**
 *	Class constructor:
 *	Params :
 *		- tcp::endpoint& endpoint_down : download enpoint ;
 *		- tcp::endpoint& endpoint_up   : upload endpoit.
 **/
PIRServer::PIRServer(boost::asio::io_service &ios, unsigned int port, uint64_t split_value, bool usedbgenerator, uint64_t dbgenerator_n, uint64_t dbgenerator_l) : 
	acceptor(ios, tcp::endpoint(tcp::v4(), port)),
  drivenMode(true),
  pirParamsFilePath("exp/PIRParams.cfg"),
  no_pipeline_mode(false)
{
	acceptor.set_option(tcp::acceptor::reuse_address(true));
  session_option.driven_mode = true;
  session_option.keep_database = false;
  session_option.got_preimported_database = false;
  if (usedbgenerator)
  {
    std::cout << "PIRServer: Launching DBGenerator" << std::endl;
    dbhandler = new DBGenerator(dbgenerator_n, dbgenerator_l, false); 
  }
  else if (split_value == 1)
  {
    std::cout << "PIRServer: Launching DBDirectoryProcessor without file splitting" << std::endl;
    dbhandler = new DBDirectoryProcessor();
  }
  else
  {
    std::cout << "PIRServer: Launching DBDirectoryProcessor with split_value=" << split_value 
      << std::endl;
    dbhandler = new DBDirectoryProcessor(split_value);
  }
}


/**
 * Wait new client and fork when new connexion.
 **/
void PIRServer::serve()
{
  PIRSession::pointer new_session =  PIRSession::create(acceptor.get_io_service());
  new_session->no_pipeline(no_pipeline_mode);
  new_session->setDBHandler(dbhandler);

  if (!drivenMode) new_session->setPIRParams(savedPIRParams);

	acceptor.async_accept(new_session->getSessionSocket(), 
    boost::bind( &PIRServer::handleAccept, this, new_session, boost::asio::placeholders::error ));
}

void PIRServer::processDrivenSession(const std::string & file_path)
{
	std::cout <<"PIRServer: Creating initial PIRSession" << std::endl;
  std::cout <<"PIRServer: Please launch a client to configure the server"<<std::endl;

  session_option.driven_mode = true;
  session_option.keep_database = true;

  // Loop while we get dry-run queries
  PIRSession* session_ptr = new PIRSession(acceptor.get_io_service());
  do 
  {
    // Trick to have only one PIRSession at a time ... must be a better way :)
    delete session_ptr;
    session_ptr = new PIRSession(acceptor.get_io_service());
    session_ptr->setDBHandler(dbhandler);
    session_ptr->no_pipeline(no_pipeline_mode);
    acceptor.accept(session_ptr->getSessionSocket());
  }while(session_ptr->start(session_option));
 
  savedPIRParams = session_ptr->getPIRParams();
  session_option.data = session_ptr->getSavedDatabase();

  session_option.driven_mode = false;
  session_option.keep_database = false;
  if (session_option.data.imported_database_ptr != NULL) session_option.got_preimported_database = true;

  if(ServerService::writePIRParameters(savedPIRParams,  file_path)) {
    cout << "Unable to write PIR parameters configuration data, maybe target repertory doesn't exist ?" << std::endl;
    return;
  }
	std::cout << "PIRServer: Configuration session successfuly ended" << endl << endl;
  drivenMode = false;
}

void PIRServer::handleAccept(PIRSession::pointer pir_session, const boost::system::error_code& error)
{
	std::cout<<"PIRServer: Creating new PIRSession"<<std::endl;
  if(!error)
  {
    boost::thread t(boost::bind(&PIRSession::start, pir_session, session_option));
  }

  serve();
}

void PIRServer::readPIRParamsFromFile(const std::string& file_path)
{
  if (ServerService::readPIRParameters(savedPIRParams, file_path)) {
        cout << "PIR parameters configuration file not found, try to launch the server in driven mode." << std::endl;
        exit(1);
    }
  drivenMode = false;
}

void PIRServer::no_pipeline(bool b)
{
  no_pipeline_mode = b;
}

PIRServer::~PIRServer()
{
  delete dbhandler; 
}
