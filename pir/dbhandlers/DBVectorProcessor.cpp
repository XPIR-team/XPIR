/* Copyright (C) 2017 Carlos Aguilar Melchor, Joris Barrier, Marc-Olivier Killijian
 *
 * This file is written by Konstantinos Andrikopoulos
 *
 * This file is part of XPIR.
 *
 * XPIR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * XPIR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with XPIR.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include "DBVectorProcessor.hpp"

Element::Element(std::string e_name, uint64_t e_size, char* e_data) :
name(e_name),
data_size(e_size),
data(NULL)
{
    data = (char*) malloc(data_size);
    memcpy(data, e_data, data_size);
}


DBVectorProcessor::DBVectorProcessor(std::vector<element_t>& vector_db) :
elements(vector_db)
{
   maxFileByteSize = 0;

   for (auto e : elements)
   {
       if (e.data_size > maxFileByteSize)
       {
           maxFileByteSize = e.data_size;
       }
   }
}

DBVectorProcessor::~DBVectorProcessor() {}

std::string DBVectorProcessor::getCatalog(const bool typeOfCatalog) {
    std::string buf;
    if(typeOfCatalog) {
        buf = std::to_string((unsigned int)0) + "\n";
        buf += std::to_string(getNbStream()) + "\n";
        for (auto e : elements)
        {
            //auto e = elements[i];
            buf += e.name + "\n" + std::to_string(e.data_size) + "\n";
        }
        return buf;
    }
    else {
        buf = std::to_string((unsigned int)1) + "\n";
        buf += std::to_string(getNbStream());
        buf += std::to_string(getmaxFileBytesize()) + "\n";
        return buf;
    }
}

uint64_t DBVectorProcessor::getNbStream() {
    return elements.size();
}

uint64_t DBVectorProcessor::getmaxFileBytesize() {
    return maxFileByteSize;
}

bool DBVectorProcessor::openStream(uint64_t streamNb, uint64_t requested_offset) {
    if(openStreamOffsets.count(streamNb)) {
        return false;
    }

    char* stream = elements[streamNb].data + requested_offset;
    openStreamOffsets.insert( std::pair<uint64_t, char*>(streamNb, stream));
    return true;
}

uint64_t DBVectorProcessor::readStream(uint64_t streamNb, char * buf, uint64_t size) {
    element_t e = elements[streamNb];
    char* stream = openStreamOffsets[streamNb];
    uint64_t sizeRead = stream - e.data;
    uint64_t sizeRemaining = e.data_size - sizeRead;

    if(sizeRemaining >= size) {
        memcpy(buf, stream, size);
        stream += size;
    }
    else {
        memcpy(buf, stream, sizeRemaining);
        bzero(buf + sizeRemaining, size - sizeRemaining);
        stream += sizeRemaining;
    }

    openStreamOffsets[streamNb] = stream;
    return size;
}

void DBVectorProcessor::closeStream(uint64_t streamNb) {
    openStreamOffsets.erase(streamNb);
}
