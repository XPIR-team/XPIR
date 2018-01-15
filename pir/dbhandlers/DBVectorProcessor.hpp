/* Copyright (C) 2017 Carlos Aguilar Melchor, Joris Barrier, Marc-Olivier Killijian
 *
 * This file is written by Konstantinos Andrikopoulos
 *
 * This file is part of XPIR.
 *
 *  XPIR is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
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

#ifndef DEF_DBVECPROC
#define DEF_DBVECPROC

#include "DBHandler.hpp"

#include <vector>
#include <map>
#include <tuple>
#include <string.h>


struct Element {
    Element(std::string, uint64_t, char*);
    std::string name;
    size_t data_size;
    char *data;
};

typedef Element element_t;

class DBVectorProcessor : public DBHandler
{
private:
    std::vector<element_t>& elements;
    std::map<uint64_t, char*> openStreamOffsets;
    uint64_t maxFileByteSize;

public:
    DBVectorProcessor(std::vector<element_t>& vector_db);
    virtual ~DBVectorProcessor();

    std::string getCatalog(const bool typeOfCatalog);

	uint64_t getNbStream();
	uint64_t getmaxFileBytesize();

	bool openStream(uint64_t streamNb, uint64_t requested_offset);
	uint64_t readStream(uint64_t streamNb, char * buf, uint64_t size);
	void closeStream(uint64_t streamNb);
};

#endif /* DEF_DBVECPROC */
