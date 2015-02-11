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

#ifndef DEF_PIRQUERY
#define DEF_PIRQUERY

#include <boost/thread.hpp>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>

#include "../common/shared_queue.hpp"
#include "../client/DESC.hpp"
#include "../common/ClientDefines.hpp"
#include "../pir/PIRReplyExtraction.hpp"
#include "../crypto/PaillierAdapter.hpp"


using namespace std;
using namespace boost;

class PIRQueryGenerator
{
	private:
		unsigned int* coord;
		uint64_t chosenElement;
		PIRParameters& pirParams;
	  HomomorphicCrypto& cryptoMethod;
		void computeCoordinates(void);

	public:
		boost::mutex mutex;
		thread queryThread;
		shared_queue <char*> queryBuffer;
	
		PIRQueryGenerator(PIRParameters& pirParameters, HomomorphicCrypto& cryptoMethod_);
		~PIRQueryGenerator();

		void startGenerateQuery();
		//void generateQuery(PIRParameters& PIRParameters, HomomorphicCrypto& cryptoMethod);
		void generateQuery();
		
		/******Getters*****/
		uint64_t getChosenElement();
		void setChosenElement(uint64_t chosenElement);
    void setPIRParameters(PIRParameters& pirParams_);
		void joinThread();
    void cleanQueryBuffer();
};
#endif
