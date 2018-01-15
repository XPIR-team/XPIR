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

#include "DBDirectoryProcessor.hpp"

/************************************************/
/* Default constructor : no splitting           */
/*  -> 1 input file -> 1 output stream          */
/************************************************/
DBDirectoryProcessor::DBDirectoryProcessor() : filesSplitting(false) {
	// TODO(feature) deal with sub-directories in the database !

	directory=std::string(DEFAULT_DIR_NAME);
	maxFileBytesize=0;

	// Then create the catalog and get the filenumbers and size
	DIR *dir = opendir (directory.c_str());
	struct dirent *ent = nullptr;

	// If there is no error when opening the directory
	if (dir != NULL) 
	{
		uint32_t i=0;
		// For each entry
		while ( (ent = readdir (dir)) != NULL) 
		{
			// Ignore files . and ..
			if (strcmp(ent->d_name, ".") > 0 && strcmp(ent->d_name, ".."))
			{
				// Count processed files (one out of 2**7?)
				if ((i << 25)==0) std::cout << "DBDirectoryProcessor: " << i+1 << " entries processed\r" << std::flush;i++;
				// Add File object on the file list
				std::string fileName= std::string( ent->d_name );
				file_list.push_back( fileName );
				uint64_t fileSize = getFileSize(directory + fileName);
				if (fileSize > maxFileBytesize)
					maxFileBytesize = fileSize;		
			}
		}
		std::cout << "DBDirectoryProcessor: " << i << " entries processed" << std::endl;
    if (i==0) {
      std::cout <<"DBDirectoryProcessor: No entries in the database" << std::endl;
      error = true;
    }
		closedir (dir);
	}
	else // If there was a problem opening the directory
	{
		std::cout << "DBDirectoryProcessor: Error opening database directory" << std::endl;
    error = true;
	}

	std::cout << "DBDirectoryProcessor: The size of the database is " << maxFileBytesize*file_list.size() << " bytes" << std::endl;
	std::cout << "DBDirectoryProcessor: The number of elements in the catalog is " << file_list.size() << std::endl;
}

// This constructor is called when we need File-splitting
DBDirectoryProcessor::DBDirectoryProcessor(uint64_t nbStreams) : filesSplitting(true) {

	directory=std::string(DEFAULT_DIR_NAME);
	maxFileBytesize=0;

	// Then create the catalog and get the filenumbers and size
	DIR *dir = opendir (directory.c_str());
	struct dirent *ent = nullptr;

	// If there is no error when opening the directory
	if (dir != NULL) 
	{
		ent = readdir (dir);
		// WARNING: In case of file-splitting, we deal only with the first file
		// On some filesystems, the dir contains also special files such as "." and "..", skip them
		while (ent->d_name == NULL || ent->d_type != DT_REG) {
			ent = readdir (dir);
		}

		// Add File object on the file list
		std::string fileName=directory + std::string( ent->d_name );
		realFileName=fileName;
		uint64_t realFileSize = getFileSize(realFileName);				
		maxFileBytesize = realFileSize/nbStreams;

		if(maxFileBytesize==0) {
			std::cout << "DBDirectoryProcessor: ERROR cannot split a file en less than one byte elements!" << std::endl;
			std::cout << "DBDirectoryProcessor: file " << realFileName << " is only "<< realFileSize << " long" << std::endl;
			error = true;
		}

		closedir (dir);
		for(int i=0;i<nbStreams;i++) {
			file_list.push_back( std::to_string(i) );
		}
	}
	else // If there was a problem opening the directory
	{
		std::cout << "DBDirectoryProcessor: Error when opening directory " <<directory<< std::endl;
	  error = true;
  }

#ifdef DEBUG
	std::cout << "maxFileBytesize." <<maxFileBytesize<< std::endl;
	std::cout << "file_list.size()." <<file_list.size()<< std::endl;

#endif
	std::cout << "DBDirectoryProcessor: The size of the database is " << maxFileBytesize*file_list.size() << " bytes" << std::endl;
	std::cout << "DBDirectoryProcessor: The number of elements in the catalog is " << file_list.size() << std::endl;
}

DBDirectoryProcessor::~DBDirectoryProcessor() {
    for(auto it : fdPool) {
        delete it.second;
    }
}

std::string DBDirectoryProcessor::getCatalog(const bool typeOfCatalog) {
	std::string buf;
	directory=std::string(DEFAULT_DIR_NAME);
	if(typeOfCatalog) {
		// Start with the number of elements in the catalog
		buf = std::to_string((unsigned int)0)+ "\n";	
		buf += std::to_string(getNbStream())+ "\n";		
		// Then for each file contactenate (with newlines) filename and filesize
		for (auto f : file_list)
		{
			if(!filesSplitting) {
				buf += f + "\n" + std::to_string(getFileSize(directory+f)) + "\n";
			} else {
				buf += f + "\n" + std::to_string(getmaxFileBytesize()) + "\n";
			}
		}
		return buf;
	} 
	// else we want a compact representation, i.e. nbFiles / fileSize
	buf = std::to_string((unsigned int)1)+ "\n";	
	buf += std::to_string(getNbStream())+ "\n";
	buf += std::to_string(maxFileBytesize)+ "\n";
	return buf;
}

uint64_t DBDirectoryProcessor::getDBSizeBits() {
	return maxFileBytesize*file_list.size()*8;
}
uint64_t DBDirectoryProcessor::getNbStream() {
	return  file_list.size();
}
uint64_t DBDirectoryProcessor::getmaxFileBytesize() {
	return maxFileBytesize;
}
bool DBDirectoryProcessor::getErrorStatus() {
	return error;
}

bool DBDirectoryProcessor::openStream(uint64_t streamNb, uint64_t requested_offset) {
    if(fdPool.count(streamNb)) {
        return false;
    }
	std::string local_directory(DEFAULT_DIR_NAME);

	std::ifstream* is = new std::ifstream();
	// When there is no splitting, each ifstream is associated with a real file 
	// (at least when no aggregation is done which is the case for now)
	if(!filesSplitting) {
		is->open( local_directory + file_list[streamNb], std::ios::binary );
		is->seekg(requested_offset);
	} else {
		// But when we are doing file splitting, we just need to position the ifstream at the correct position
		uint64_t splitting_offset=streamNb*getmaxFileBytesize();
		is->open( realFileName, std::ios::binary );
		is->seekg(splitting_offset + requested_offset);
	}
    fdPool.insert( std::pair<uint64_t, std::ifstream*>(streamNb, is));
	return true;
}

uint64_t DBDirectoryProcessor::readStream(uint64_t streamNb, char * buf, uint64_t size) {
    std::ifstream *s = fdPool[streamNb];
	uint64_t sizeRead=0;
	//std::cout << "sizeRead = "<<sizeRead<<" size = "<<size<<std::endl;
	while(sizeRead<size) {
		uint64_t readThisrun=s->readsome(buf+sizeRead,size-sizeRead);
		sizeRead+=readThisrun;
		// Check if we need to pad
		if(readThisrun==0 && sizeRead<size) {
			//		std::cout << "padding = "<<size-sizeRead<<std::endl;
			bzero(buf+sizeRead,size-sizeRead);
			sizeRead=size;
		}
	}
	return size;
}

void DBDirectoryProcessor::closeStream(uint64_t streamNb) {
    if(!fdPool.count(streamNb)) {
        return;
    }
    std::map<uint64_t, std::ifstream*>::iterator it = fdPool.find(streamNb);
    it->second->close();
    fdPool.erase(it);
}

std::streampos DBDirectoryProcessor::getFileSize( std::string filePath ){
	std::streampos fsize = 0;
	std::ifstream file( filePath.c_str(), std::ios::binary );
	fsize = file.tellg();
	file.seekg( 0, std::ios::end );
	fsize = file.tellg() - fsize;
	file.close();
	return fsize;
}
