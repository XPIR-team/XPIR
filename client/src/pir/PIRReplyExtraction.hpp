/* Copyright (C) 2014 Carlos Aguilar Melchor, Joris Barrier, Marc-Olivier Killijian
 * This file is part of XPIRe.
 *
 *  XPIRe is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  XPIRe is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XPIRe.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DEF_PIRRESULT
#define DEF_PIRRESULT

#include <iostream>
#include <boost/thread.hpp>
#include <boost/signals2.hpp>
#include <omp.h>

#include "PIRParameters.hpp"
#include "../crypto/PaillierAdapter.hpp"
#include "../events/WriteEvent.hpp"
#include "../events/MessageEvent.hpp"
#include "../common/ClientDefines.hpp"
#include "../common/shared_queue.hpp"
#include "../common/GlobalConstant.hpp"


class PIRReplyExtraction 
{
	private :

		string 				 filename, filePath;
		PIRParameters& pirParams;
		HomomorphicCrypto& cryptoMethod;
		
		int   fileSize,keySize;
		
		void writeFile();
	
	public :
		boost::thread  replyThread;
		shared_queue<char*> repliesBuffer;
		
    PIRReplyExtraction(HomomorphicCrypto& cryptoMethod_, PIRParameters& pirParams); 
		~PIRReplyExtraction();
	
		void startExtractReply(int aggregated_maxFileSize, shared_queue<char*>* clearChunks);
		void extractReply(int aggregated_maxFileSize, shared_queue<char*>* clearChunks);

		void setFileParam(string filename,int fileSize);
		
		void pushCipheredChunk(char* rawBytes);
		char* getCipheredChunk();
		int getChosenFileSize();
		
		void joinThread();
};

#endif
