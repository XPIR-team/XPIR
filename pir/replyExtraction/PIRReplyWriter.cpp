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

#include "PIRReplyWriter.hpp"

const std::string PIRReplyWriter::kDefaultFolder("reception");

PIRReplyWriter::PIRReplyWriter(PIRParameters& param, boost::signals2::signal<void (WriteEvent&)> &writeListeners_, boost::signals2::signal<void (MessageEvent&)> &messageListeners_) :
  filePath(kDefaultFolder),
  pirParams(param),
  clearChunks("clear_chunks"), 
  writeListeners(writeListeners_),
  messageListeners(messageListeners_),
  dontWrite(false)
{
  ;//Marco's boost signals fix ! :-D
}

void PIRReplyWriter::setdontWrite(bool newvalue) { dontWrite = newvalue; }

void PIRReplyWriter::writeAggregatedFileSecurely(uint64_t chosenElement, DESC catalog) 
{
  uint64_t firstElement = chosenElement - chosenElement % pirParams.alpha;
  uint64_t lastElement = firstElement + pirParams.alpha - 1;
  uint64_t bytestoskip = 0; 
  char *tmp;
  uint64_t chunkSize = cryptoMethod->getPublicParameters().getAbsorptionBitsize() 
    / GlobalConstant::kBitsPerByte;
  
  WriteEvent event(catalog.getMaxFileSize(), 0);

  if (pirParams.alpha > 1) std::cout << "PIRReplyWriter: Dealing with " << pirParams.alpha << " aggregated files" << std::endl;
  
  if (chosenElement != firstElement) std::cout << "PIRReplyWriter: Skipping aggregated files before chosen element ..." << std::endl;

  // To avoid timing attacks write down all the elements that were aggregated
  for (uint64_t i = firstElement; i <= lastElement; i++)
  {
#ifdef MORE_SIDE_CHANNEL_RESISTANCE
    writeFileSecurely(i, catalog, bytestoskip, event);
#else
    if (i == chosenElement && dontWrite == false) 
    {
      writeFileSecurely(i, catalog, bytestoskip, event);
    }
    else
    {
      bytestoskip += catalog.getMaxFileSize();
    }
#endif
  }

  if (chosenElement != lastElement) std::cout << "PIRReplyWriter: Skipping aggregated files after chosen element ..." << std::endl;

  // Empty clearchunk queue and say we finished
  std::cout << "PIRReplyWriter: Emptying queue ..." << std::endl;
  while (bytestoskip > chunkSize)
  {
    tmp = clearChunks.pop_front();
    free(tmp);
    bytestoskip -= chunkSize;
  }
  
  // ... And remove all of them but the element we are interested in
  if (truncate(std::string(filePath + "/" + catalog.getFileName(chosenElement)).c_str(), 
      catalog.getFileSize(chosenElement)))
  {
    std::cout << "PIRReplyWriter: Unable to truncate retrieved file" << std::endl;
  }
#ifdef MORE_SIDE_CHANNEL_RESISTANCE
  for (uint64_t i = firstElement; i <= lastElement; i++)
  {
    if ( i != chosenElement ) 
    {
      std::remove(std::string(filePath + "/" + catalog.getFileName(i)).c_str());
    }
  }
#endif

  std::cout << "PIRReplyWriter: FINISHED (This should end the client)" << std::endl << std::endl;

}

void PIRReplyWriter::writeFileSecurely(uint64_t element, DESC catalog, uint64_t &bytestoskip, WriteEvent &event)
{
  std::string completeFilename(filePath + "/" + catalog.getFileName(element));
  std::ofstream file(completeFilename.c_str(), ios::out | ios::binary | ios::trunc);

  if (file.good()) 
  {
    char *tmp;
    uint64_t chunkSize = cryptoMethod->getPublicParameters().getAbsorptionBitsize() / GlobalConstant::kBitsPerByte;
    uint64_t leftChars = catalog.getMaxFileSize();
#ifdef DEBUG
    cout << "PIRReplyWriter: Size of the requested file is " << leftChars << endl; 
    cout << "PIRReplyWriter: chunkSize is " << chunkSize << endl;
#endif
    
    while (bytestoskip > chunkSize)
    {
      tmp = clearChunks.pop_front();
      free(tmp);
      bytestoskip -= chunkSize;
      event.addtoWrittenSize(chunkSize);
    }
    writeListeners(event);
    
    while (leftChars != 0) 
    {
      uint64_t writtenchars = 0;
      uint64_t current_bytesskipped = 0; 

      tmp = clearChunks.front();
      if (bytestoskip != 0)
      {
        file.write(tmp+bytestoskip , min(leftChars, chunkSize - bytestoskip));
        writtenchars = min(leftChars, chunkSize - bytestoskip);
        bytestoskip = chunkSize - (writtenchars+bytestoskip);
        if (bytestoskip < 0) std::cout << "PIRReplyWriter: Skipping a negative amount of bytes THIS SHOULD NOT HAPPEN" << std::endl;
      }
      else
      {
        file.write(tmp , min(leftChars, chunkSize));
        writtenchars = min(leftChars, chunkSize);
        if (writtenchars < chunkSize ) bytestoskip = writtenchars;
      }
        
      if (file.fail())
      {
        MessageEvent mEvent(WARNING,"FAIL DURING FILE WRITTING");
        messageListeners(mEvent);
        file.close();
        exit(1);
      }
      file.flush();
      
      // Shall we reuse this clear chunk ?
      if ( bytestoskip == 0 )
      {
        clearChunks.pop_front();
        free(tmp);
      }
      
      event.addtoWrittenSize(writtenchars);
      writeListeners(event);
      leftChars -= writtenchars;

    }
    
    
    file.close();
 
  }
  else
  {
    MessageEvent event(ERROR, "PIRReplyWriter: ERROR WHEN WRITING FILE ! Maybe " + completeFilename + " is write protected?", __FUNCTION__);
    messageListeners(event);
  }
}

void PIRReplyWriter::startFileWritting(uint64_t chosenElement, DESC catalog)
{
  writeThread = boost::thread(&PIRReplyWriter::writeAggregatedFileSecurely,  this, chosenElement, catalog);
}

shared_queue<char*>* PIRReplyWriter::getClearDataQueue()
{
  return &clearChunks;
}

void PIRReplyWriter::setCryptoMethod(HomomorphicCrypto* crypto_method)
{
  cryptoMethod = crypto_method;
}

void PIRReplyWriter::join()
{
  if(writeThread.joinable()) writeThread.join();
}

PIRReplyWriter::~PIRReplyWriter()
{
  join();
}
