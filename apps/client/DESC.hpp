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

#ifndef DEF_DESC
#define DEF_DESC

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <boost/signals2.hpp>

#include "pir/events/MessageEvent.hpp"

using namespace std;
#define FILENAME_MAX_BYTE_SIZE 1023	

typedef boost::signals2::signal<void (MessageEvent&)>   messageListener;
class DESC
{
    public :
				DESC (messageListener& messageListeners);
        void makeMenu(char* receivedBuffer);
        bool file_exists (uint64_t ID);
       
				string 			 getFileName(uint64_t index);
        uint64_t				 getFileSize(uint64_t index);
				const        vector<string>& getFileList();
        uint64_t getFilesNum();
				uint64_t 				 getMaxFileSize();
   
				~DESC();
		private:
				messageListener& messageListeners;
        vector<string> fileList;
        vector<uint64_t> fileSize;
				uint64_t maxFileSize;
        uint64_t nbFiles;
};

#endif
