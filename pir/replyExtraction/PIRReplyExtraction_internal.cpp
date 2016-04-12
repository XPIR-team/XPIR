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

#include "PIRReplyExtraction_internal.hpp"
#include "apps/client/PIRViewCLI.hpp"
//#define DEBUG

/**
 *	Class constructor.
 *	Params :
 *		- crypto_ptr crypto 	 : HomomorphicCrypto shared_ptr ;
 *		- PIRParameters& param : PIRParameters reference ;
 *		- signal1<void, WriteEvent> 	: writeEvent listeners reference shared with PIRClient ;
 *		- ignal1<void, MessageEvent&> : messageEvent listeners reference shared with PIRClient.
 **/
PIRReplyExtraction_internal::PIRReplyExtraction_internal(PIRParameters& param_, HomomorphicCrypto& cryptoMethod_) :
  filePath("reception"),
  pirParams(param_),
  cryptoMethod(cryptoMethod_),
  repliesBuffer("coucou"),
  fileSize(0)
{}

/**
 *	Parallely extract/decrypt file, store clear chunk in clearChunks shared queue.
 **/
void PIRReplyExtraction_internal::extractReply(int aggregated_maxFileSize, shared_queue<char*>* clearChunks)
{
  unsigned long ciphertext_nbr = 1;
  unsigned int data_size, data_size2b;
  char *data, *in_data, *in_data_2b;
  double start = omp_get_wtime();
  uint64_t total_ciphertext_nbr= 0;

  for (unsigned int rec_lvl = pirParams.d ; rec_lvl >= 1 ; rec_lvl--)
  {
    ciphertext_nbr = ceil(static_cast<float>(aggregated_maxFileSize) / static_cast<float>(cryptoMethod.getPublicParameters().getAbsorptionBitsize(0)/GlobalConstant::kBitsPerByte));
#ifdef DEBUG
	std::cout<<"PIRReplyExtraction_internal: First layer ciphertext_nbr="<<ciphertext_nbr<<std::endl;
#endif

    for (unsigned int i = 1 ; i < rec_lvl; i++)
    {
      ciphertext_nbr =  ceil(float(ciphertext_nbr) * float(cryptoMethod.getPublicParameters().getCiphBitsizeFromRecLvl(i)/GlobalConstant::kBitsPerByte) / float(cryptoMethod.getPublicParameters().getAbsorptionBitsize(i)/GlobalConstant::kBitsPerByte));
    }
    total_ciphertext_nbr += ciphertext_nbr;

#ifdef DEBUG
	std::cout<<"PIRReplyextraction: Last layer ciphertext_nbr="<<ciphertext_nbr<<std::endl;
#endif

    char* out_data;

    data_size  = cryptoMethod.getPublicParameters().getCiphBitsizeFromRecLvl(rec_lvl) / GlobalConstant::kBitsPerByte ;
    data_size2b = cryptoMethod.getPublicParameters().getAbsorptionBitsize(rec_lvl-1)/GlobalConstant::kBitsPerByte;
    if (rec_lvl > 1) in_data_2b = (char*) calloc(data_size2b * ciphertext_nbr, sizeof(char));

#ifdef DEBUG
    cout << "PIRReplyExtraction_internal: rec_lvl=" << rec_lvl << " ciphertext_nbr=" << ciphertext_nbr <<  " data_size=" << data_size << " data_size2b=" << data_size2b << endl;
#endif
    cout << "PIRReplyExtraction_internal: Waiting for first replies..." << endl;

    for (unsigned int j = 0 ; j < ciphertext_nbr ; j++)
    {
      data = (rec_lvl == pirParams.d) ? repliesBuffer.pop_front() : in_data+(j*data_size); 
      if (rec_lvl == pirParams.d && j == 0 ) 
      { 
        cout << "PIRReplyExtraction_internal: Starting reply extraction..." << endl;
      }
      out_data = cryptoMethod.decrypt(data, rec_lvl, data_size, data_size2b);
     if (rec_lvl > 1) {
       memcpy(in_data_2b+(data_size2b * j), out_data, data_size2b); 
       free(out_data);
     }
     else 
     {
       clearChunks->push(out_data);
#ifdef CRYPTO_DEBUG
      cout << "PIRReplyExtraction_internal : pushed " << hex;
      for (int k = 0 ; k < data_size2b ; k++)
        cout << (unsigned short) *(out_data+k) << " ";
      cout << dec << endl;
#endif
     }
     if(rec_lvl == pirParams.d) free(data);
    }
    if (rec_lvl < pirParams.d) free(in_data);
    in_data = in_data_2b;
  }
  cout << "PIRReplyExtraction_internal: Reply extraction finished, " << total_ciphertext_nbr <<
    " reply elements decrypted in " << omp_get_wtime() - start << " seconds" << endl;
}


/**
 *	Start extractReply in a thread.
 **/
void PIRReplyExtraction_internal::startExtractReply(int aggregated_maxFileSize, shared_queue<char*>* clearChunks) 
{
  replyThread = boost::thread(&PIRReplyExtraction_internal::extractReply, this, aggregated_maxFileSize, clearChunks);
}

void PIRReplyExtraction_internal::setFileParam(string filename_ ,int fileSize_) 
{
  assert (!filename_.empty());
  assert (fileSize_ >= 0);

  fileSize = fileSize_;
  filename = filename_;
}

PIRReplyExtraction_internal::~PIRReplyExtraction_internal() 
{
  joinThread();

  while (!repliesBuffer.empty())
    free(repliesBuffer.pop_front());

}

/**
 *	Join reply extraction thread if it's possible.
 **/
void PIRReplyExtraction_internal::joinThread() 
{
  if (replyThread.joinable()) replyThread.join();
}

int PIRReplyExtraction_internal::getChosenFileSize()
{
  return fileSize;
}

  //#pragma omp parallel for private(thread_index, buf) shared(index, index_prod)
  //	for (unsigned int i = 0 ; i < paquet_nbr ; i++)
  //	{
  //#pragma omp critical
  //		{
  //			buf = repliesBuffer.pop_front();
  //			thread_index = index++;
  //		}
  //
  //		char *tmp1;
  //		size_t ciph_size = 0, clear_size;
  //		unsigned int exp_factor = pirParams.recLvl;
  //
  //		if (start == 0)
  //#pragma omp critical
  //			start = omp_get_wtime();
  //
  //		for (unsigned int i = 0 ; i < pirParams.d; i++)
  //		{
  //			ciph_size  = cryptoMethod.getPublicParameters().getCiphBitsizeFromRecLvl(pirParams.recLvl - i);
  //			clear_size = cryptoMethod.getPublicParameters().getCiphBitsizeFromRecLvl(pirParams.recLvl - i -1);
  //			tmp1 = buf;
  //
  //			if (exp_factor - i == 1) 
  //				clear_size--;
  //
  //			buf = cryptoMethod.decrypt(tmp1, exp_factor - i, ciph_size, clear_size);
  //
  //			delete[] tmp1;
  //		}
  //
  //		while (thread_index != index_prod)
  //			usleep(5);
  //
  //#pragma omp critical
  //		{
  //			clearChunks.push(buf);
  //			index_prod = thread_index + 1;
  //		}
  //	}
  //printf("\n:: Déchiffrement -> Temp écoulé : %f seconds\n", omp_get_wtime() - start);
