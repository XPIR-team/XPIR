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

#ifndef DEF_ABSTRACTPUBLICPARAMETERS
#define DEF_ABSTRACTPUBLICPARAMETERS
#include <cstddef>
#include <gmp.h>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

class AbstractPublicParameters {
  protected:
    std::string cryptoName;

  public:
    virtual std::string getSerializedParams(bool shortversion)=0;
    virtual char* getByteModulus()=0;
    virtual unsigned int getAbsorptionBitsize()=0;
    virtual unsigned int getAbsorptionBitsize(unsigned int rec_lvl)=0;
    virtual void setModulus(char* rawPubKey)=0;
    virtual void setMockedPubKey()=0;
    virtual unsigned int getCiphertextBitsize()=0;
    virtual unsigned int getCiphBitsizeFromRecLvl(unsigned int)=0;
    virtual unsigned int getQuerySizeFromRecLvl(unsigned int)=0;
    virtual void computeNewParameters(const std::string& crypto_param_descriptor)=0;
    virtual unsigned int getSerializedModulusBitsize()=0;
    virtual ~AbstractPublicParameters(){};
};

#endif
