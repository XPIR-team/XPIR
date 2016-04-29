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

#include "DBGenerator.hpp"

/************************************************/
/* Default constructor : no splitting           */
/*  -> 1 input file -> 1 output stream          */
/************************************************/
DBGenerator::DBGenerator(uint64_t nbStreams, uint64_t streamBytesize, bool silent):
random_engine(0), // Fixed seed of 0
random_distribution()
{
	maxFileBytesize=streamBytesize;
	nbFiles = 0;
#ifdef SEND_CATALOG
  for (unsigned int i = 0 ; i < nbStreams ; i++) {
	std::string fileName= std::to_string(i);
    file_list.push_back( fileName );
	}
#endif
  nbFiles = nbStreams;
  if(!silent)
  {
	  std::cout << "DBGenerator: The size of the database is " << maxFileBytesize*nbFiles << " bytes" << std::endl;
	  std::cout << "DBGenerator: The number of elements in the catalog is " << nbFiles << std::endl;
  }
}

DBGenerator::~DBGenerator() {}

std::string DBGenerator::getCatalog(const bool typeOfCatalog) {
	std::string buf;
	if(typeOfCatalog) {
		// Start with the number of elements in the catalog
		buf = std::to_string((unsigned int)0)+ "\n";	
		buf += std::to_string(getNbStream())+ "\n";		
		// Then for each file contactenate (with newlines) filename and filesize
		for (auto f : file_list)
		{
			buf += f + "\n" + std::to_string(maxFileBytesize) + "\n";
		}
		return buf;
	} 
	// else we want a compact representation, i.e. nbFiles / fileSize
	buf = std::to_string((unsigned int)1)+ "\n";	
	buf += std::to_string(getNbStream())+ "\n";
	buf += std::to_string(maxFileBytesize)+ "\n";
	return buf;
}

//uint64_t DBGenerator::getDBSizeBits() {
//	return maxFileBytesize*nbFiles*8;
//}
uint64_t DBGenerator::getNbStream() {
	return  nbFiles;
}
uint64_t DBGenerator::getmaxFileBytesize() {
	return maxFileBytesize;
}

std::ifstream* DBGenerator::openStream(uint64_t streamNb, uint64_t requested_offset) {
	return NULL;
}

uint64_t DBGenerator::readStream(std::ifstream* s, char * buf, uint64_t size) {
  //for (unsigned char i = 0xaa, j = 0; j < size; i++, j++)
  //{
  //  buf[j] = i;
	//}
//#define NDSS_SNIFFER
#ifdef NDSS_SNIFFER
  for (int i = 0; i < size/4; i++) {
    buf[i<<2] = random_distribution(random_engine);
  }
#else
  char ccc=0xaa;
  memset(buf, ccc++, size);
#endif
  return size;
}

void DBGenerator::closeStream(std::ifstream* s) {}

void DBGenerator::readAggregatedStream(uint64_t streamNb, uint64_t alpha, uint64_t offset, uint64_t bytes_per_file, char* rawBits){
  readStream(NULL, NULL, 0);
	uint64_t fileByteSize = std::min(bytes_per_file, maxFileBytesize-offset);
  uint64_t startStream = streamNb*alpha;
  uint64_t endStream = std::min(streamNb*alpha + alpha - 1, getNbStream() - 1);
  uint64_t paddingStreams = std::max((long long)((streamNb)*alpha+alpha) - (long long)getNbStream(), (long long)0);
  memset(rawBits, 0xaa, fileByteSize*(endStream-startStream + 1));
  
  if(paddingStreams !=0)
  {
	  bzero(rawBits + (endStream % alpha) * fileByteSize, fileByteSize*paddingStreams);
  }
}
