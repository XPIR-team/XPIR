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

#ifndef DEF_DBHANDLER
#define DEF_DBHANDLER

#include <vector>
#include <string>
#include <iostream>
#include <string.h>

class DBHandler
{
private:
  std::vector <std::string> file_list; // the output file list
public:
  
  virtual std::string getCatalog(const bool typeOfCatalog)=0;
  
  virtual uint64_t getNbStream()=0;
  virtual uint64_t getmaxFileBytesize()=0;
  
  virtual bool openStream(uint64_t streamNb, uint64_t requested_offset)=0;
  virtual uint64_t readStream(uint64_t streamNb, char * buf, uint64_t size)=0;
  virtual void readAggregatedStream(uint64_t streamNb, uint64_t alpha, uint64_t offset, uint64_t bytes_per_file, char* rawBits);
  virtual void closeStream(uint64_t streamNb)=0;
  virtual ~DBHandler(){};
  
  
private:
  	uint64_t maxFileBytesize;
};

#endif 
