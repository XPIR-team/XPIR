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

const std::string kDefault_file_path = "exp/PIRParams.cfg";

int main(int argc, char* argv[])
{
  std::cout << std::endl << "########################################################" << std::endl;
  std::cout << "# XPIR ... Private Information Retrieval for everyone #" << std::endl;
  std::cout << "########################################################" << std::endl << std::endl;

  int port = 1234;
  bool usedbgenerator = false;
  uint64_t dbgenerator_n = 10, dbgenerator_l = 12800000;
  int driven = -1;
  uint64_t split_value = 1; 
  std::string pir_params_file_path(kDefault_file_path);

	po::options_description od("Available options");	
	od.add_options()
		("help,h", "help message")
    ("driven,z", po::value<std::string>()->implicit_value(pir_params_file_path),"Server-driven mode : configure crypto and PIR parameters using the first client, make them mandatory for other clients and save them into a persistent file.")
		("load_file,L", po::value<std::string>(), "Load PIR parameters from file.")
    ("split_file,s", po::value<uint64_t>(&split_value)->default_value(1), "Only use first file in db directory and split it in arg elements.")
		("port,p", po::value<int>(&port)->default_value(1234), "Port used, default")
		("db-generator-files,n", po::value<uint64_t>(&dbgenerator_n)->default_value(10), "Number of files for the DB generator")
		("db-generator-filesize,l", po::value<uint64_t>(&dbgenerator_l)->default_value(12800000), "Size of file for the DB generator")
		("db-generator", "Generate the database instead of reading it from a directory")
    ("no-pipeline", "No pipeline mode");

	po::variables_map vm;
	try 
	{
		po::store(po::parse_command_line(argc, argv, od), vm);
	}
	catch (const std::exception& ex) 
	{
		std::cout << "Error checking program options: " << ex.what() << std::endl;
		std::cout << od << std::endl;
		return 1;
	}

	if (!vm.size()) 
	{
		std::cout << "No option specified : " << std::endl;
		std::cout << od << std::endl;
	}
	if (vm.count("help")) 
	{
		std::cout << od << std::endl;
		return 0;
	}
  if (vm.count("driven"))
  {
    driven = 1;
    pir_params_file_path =  vm["driven"].as<std::string>();
  }
  if (vm.count("load_file"))
  {
    pir_params_file_path =  vm["load_file"].as<std::string>();
    driven = 0;
  }
  if (vm.count("split_file"))
  {
    split_value = vm["split_file"].as<uint64_t>();
  }
  if (vm.count("db-generator-files"))
  {
    dbgenerator_n = vm["db-generator-files"].as<uint64_t>();
  }
  if (vm.count("db-generator-filesize"))
  {
    dbgenerator_l = vm["db-generator-filesize"].as<uint64_t>();
  }
  if (vm.count("db-generator"))
  {
    std::cout << "CLI: DBGenerator requested with params n=" << dbgenerator_n
     << " l=" << dbgenerator_l << std::endl;
    usedbgenerator = true;
  }
    
  try
  {
    boost::asio::io_service io_service;
    boost::asio::signal_set signals(io_service, SIGINT, SIGTERM);
    signals.async_wait(
    boost::bind(&boost::asio::io_service::stop, &io_service));

    PIRServer server(io_service, port, split_value, usedbgenerator, dbgenerator_n, dbgenerator_l);

    if (vm.count("no-pipeline"))
    {
    std::cout << "Main: WARNING no pipeline mode activated" << std::endl;
    server.no_pipeline(true);
    }

    if (driven == 1) server.processDrivenSession(pir_params_file_path);

    else if (driven == 0) server.readPIRParamsFromFile(pir_params_file_path);

    server.serve();
    io_service.run();

    cout << "Main: PIRServer has stopped" << endl;
  }
  catch (std::exception& e)
  {
    std::cerr << "Exception catched in main.cpp " << e.what() << std::endl;
  }
}
