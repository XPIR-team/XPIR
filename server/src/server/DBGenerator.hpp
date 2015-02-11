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

#ifndef DEF_DBGEN
#define DEF_DBGEN

#include "DBHandler.hpp"
#include <cstring>
#include <random>
#include <omp.h>

class DBGenerator : public DBHandler
{
private:
  std::vector <std::string> file_list; // the output file list
public:
  DBGenerator(uint64_t nbStreams, uint64_t streamBytesize, bool silent); 
  ~DBGenerator();
  
  std::string getCatalog(const bool typeOfCatalog);
  
  uint64_t getNbStream();
  uint64_t getmaxFileBytesize();
  
  std::ifstream* openStream(uint64_t streamNb, uint64_t requested_offset);
  uint64_t readStream(std::ifstream* s, char * buf, uint64_t size);
  void readAggregatedStream(uint64_t streamNb, uint64_t alpha, uint64_t offset, uint64_t bytes_per_file, char* rawBits);
  void closeStream(std::ifstream* s);
  
private:
  std::mt19937_64 random_engine; // Fixed seed of 0
  std::uniform_int_distribution<> random_distribution;
  
 	uint64_t maxFileBytesize;
  uint64_t nbFiles;  
};

#endif //DEF_CATALOGMAKER
