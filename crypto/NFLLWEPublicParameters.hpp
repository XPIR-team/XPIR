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


#ifndef DEF_NFLLWEPUBPARAMS
#define DEF_NFLLWEPUBPARAMS

#include <inttypes.h>
#include <string>
#include "AbstractPublicParameters.hpp"
#include <iostream>


class NFLLWE;

class NFLLWEPublicParameters : public AbstractPublicParameters
{
  public :
    NFLLWEPublicParameters();
    NFLLWEPublicParameters(unsigned int modulusBitsize_, unsigned int polyDegree_, int absPCBitsize_);
    void setModulus(char* rawPubKey);
    void setMockedPubKey();

    // Getters
    unsigned int getmodulusBitsize();
    unsigned int getModulusRepresentationBitsize();
    uint64_t* getmoduli();
    unsigned int getpolyDegree();
    unsigned int getSerializedModulusBitsize();
    unsigned int getAbsorptionBitsize();
    unsigned int getAbsorptionBitsize(unsigned int rec_lvl);
    uint64_t getnoiseUB();
    uint64_t getsecurityBits();

    // Setters
    void setmodulus(uint64_t modulus);
    void setpolyDegree(unsigned int polyDegree);
    void setNewParameters(std::string crypto_param_desc);
    void setAbsPCBitsize(int absPCBitsize);
    void setnoiseUB(uint64_t Berr_);
    void setsecurityBits(uint64_t Berr_);
    void setcrypto_container(NFLLWE* c) {crypto_container=c;}

  private :
    //unsigned int modulusBitsize;
    //uint64_t modulus;
    //unsigned int polyDegree;
    int absPerCoordinateBitsize; // Signed as -1 means uninitialized
    uint64_t noise_ub;
    uint64_t securityBits; // Defines the security level
    NFLLWE* crypto_container;


  public:
    char* getByteModulus();
    std::string getSerializedParams(bool shortversion);
    void getParameters();
    unsigned int getCiphertextSize();
    void computeNewParameters(const std::string&);
    void newKeyParameter(unsigned int bitKeySize);
    unsigned int getCiphertextBitsize();
    unsigned int getCiphBitsizeFromRecLvl(unsigned int);
    unsigned int getQuerySizeFromRecLvl(unsigned int);


};
#endif
