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

#include "DESC.hpp"
//#define NO_CATALOG //Use this option for performance tests (only with auto-choice)

/**
 *	Class constructor.
 *	Param :
 *		-	messageListener& messageListeners_ : message listerner reference to display error message.
 **/
DESC::DESC(messageListener& messageListeners_):
	messageListeners(messageListeners_),
	maxFileSize(0)
{
    ; // This semicolon is necessary under some uncertain circumstances .... e.g. OSX/gcc4.8.1
}

/**
 *	Generate a menu from a char* buffer.
 *	Param :
 *		- char* receiveBuffer : buffer to parse.
 **/
void DESC::makeMenu(char* receivedBuffer)
{
	maxFileSize = 0;
	uint64_t fileCount = 0;
	char fileName[FILENAME_MAX_BYTE_SIZE + 1];
	char tmpSize[FILENAME_MAX_BYTE_SIZE + 1];
	stringstream ss (receivedBuffer);
	unsigned int catalogType;
#ifdef DEBUG
  std::cout << "DESC: catalog received : ..."<<std::endl<<receivedBuffer<<std::endl;
#endif
	if (ss.good ())
	{
		//The first line of the catalog is the catalog type
		ss >> catalogType;
		if(catalogType==0) {
#ifdef DEBUG
			std::cout << "DESC: catalog type 0 : "<<catalogType<<std::endl;
#endif		
			ss >> fileCount;
			//The rest of catalog alternates filename lines and size lines.
			ss.getline (fileName, FILENAME_MAX_BYTE_SIZE);

      nbFiles = 0;
			for (uint64_t i = 0; i < fileCount; i++) 
			{
				ss.getline (fileName, FILENAME_MAX_BYTE_SIZE);
				ss.getline (tmpSize, FILENAME_MAX_BYTE_SIZE);

				if (ss.gcount () > 0) 
				{
					fileList.push_back(fileName);
					fileSize.push_back(atoi(tmpSize));
          nbFiles++;
#ifdef DEBUG
          std::cout<<"DESC: "<<i<<"-"<<fileName<<"-"<<tmpSize<<std::endl;
#endif
          if(maxFileSize < fileSize.back())
						maxFileSize = fileSize.back();
				}
			}
		} else if(catalogType==1) {
#ifdef DEBUG
			std::cout << "DESC: catalog type 0 : "<<catalogType<<std::endl;
#endif		
			
			ss >> fileCount;
			ss.getline (tmpSize, FILENAME_MAX_BYTE_SIZE);
			ss.getline (tmpSize, FILENAME_MAX_BYTE_SIZE);
			maxFileSize=atoi(tmpSize);
#ifndef NO_CATALOG
      for (uint64_t i = 0; i < fileCount; i++) 
			{
					fileList.push_back(std::to_string(i));
					fileSize.push_back(maxFileSize);
			}
#endif
      nbFiles = fileCount;
		} else {
			MessageEvent event(ERROR,"Catalog type unknown "+catalogType, __FUNCTION__);
			messageListeners(event);
		}
	} else {
			MessageEvent event(ERROR,"Impossible to read catalog from the server.", __FUNCTION__);
			messageListeners(event);
	}
}

/**
 *	Get file name from a given index.
 *	Param :
 *		- int index : file index.
 *
 *	Return :
 *		- std::string : return fileName if exist and "None" else.
 **/
string DESC::getFileName(uint64_t index) 
{
#ifndef NO_CATALOG
	return (index < nbFiles) ?  fileList[index] : "None";
#else
	return (index < nbFiles) ?  std::to_string(index) : "None";
#endif
}

/**
 *	Get file size frome a given index.
 *	Param :
 *		- int index : file index.
 *
 *	Return :
 *		- int : return fileSize of exist and 0 else.
 *
 **/
uint64_t DESC::getFileSize(uint64_t index) 
{
#ifndef NO_CATALOG
	return (index < fileSize.size()) ? fileSize[index] : -1;
#else
	return (index < nbFiles) ? maxFileSize : -1;
#endif
}

/**
 *	Checks whether the provided file number exists. 
 *	Param : 
 *		- int dex : file index.
 *
 *	Return :
 *		- bool : true if file exist, false else.
 **/
bool DESC::file_exists (uint64_t index) 
{
	return (index < nbFiles) ? true : false; 
}

/**
 *	Get number of files.
 *	Return :
 *		- unsigned int : number of files.
 **/
uint64_t DESC::getFilesNum() 
{
	return nbFiles;
}

/**
 *	Get the biggest file.
 *	Return :
 *		- int : biggest file size.
 *		
 **/
uint64_t DESC::getMaxFileSize()
{
	return maxFileSize;
}

/**
 *	Get file list.
 *	Return :
 *		- const std::vector<string>& : constant reference to the fileList attribute. 
 **/
const vector<string>& DESC::getFileList()
{
	return fileList;
}

DESC::~DESC() 
{
	fileList.clear();
	fileSize.clear();
}
