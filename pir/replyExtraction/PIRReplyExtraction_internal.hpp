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

#ifndef DEF_PIRRESULT
#define DEF_PIRRESULT

#include <iostream>
#include <boost/thread.hpp>
#include <boost/signals2.hpp>
#include <omp.h>

#include "pir/PIRParameters.hpp"
#include "crypto/PaillierAdapter.hpp"
#include "pir/events/WriteEvent.hpp"
#include "pir/events/MessageEvent.hpp"
#include "pir/ClientDefines.hpp"
#include "pir/shared_queue.hpp"
#include "pir/GlobalConstant.hpp"


class PIRReplyExtraction_internal 
{
  protected:
		PIRParameters& pirParams;
		HomomorphicCrypto& cryptoMethod;

	private :

		string 				 filename, filePath;
		
		int   fileSize,keySize;
		
		void writeFile();
	
	public :
		boost::thread  replyThread;
		shared_queue<char*> repliesBuffer;
		
    	PIRReplyExtraction_internal(PIRParameters& pirParameters, HomomorphicCrypto& cryptoMethod_); 
		~PIRReplyExtraction_internal();
	
		void startExtractReply(int aggregated_maxFileSize, shared_queue<char*>* clearChunks);
		void extractReply(int aggregated_maxFileSize, shared_queue<char*>* clearChunks);

		void setFileParam(string filename,int fileSize);
		int getChosenFileSize();
		
		void joinThread();
};

#endif
