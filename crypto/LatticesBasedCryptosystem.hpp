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

#ifndef DEF_LATTICESBASEDCRYPTOSYSTEM
#define DEF_LATTICESBASEDCRYPTOSYSTEM
#include "NFLlib.hpp"
#include "NFLLWEDatatypes.hpp"
#include "HomomorphicCrypto.hpp"

class LatticesBasedCryptosystem : public HomomorphicCrypto {
  public:
    LatticesBasedCryptosystem(const std::string& crypto_name);

    virtual unsigned int getpolyDegree()=0;
    virtual void mulandadd(lwe_cipher rop, lwe_in_data op1, lwe_query op2, lwe_query op2prime, 
        uint64_t current_poly, int rec_lvl)=0;
    virtual void mul(lwe_cipher rop, lwe_in_data op1, lwe_query op2, lwe_query op2prime, 
        uint64_t current_poly, int rec_lvl)=0;

    virtual long setandgetAbsBitPerCiphertext(unsigned int elt_nbr)=0;

    NFLlib& getnflInstance();
    uint64_t* getmoduli();
    uint64_t getsecurityBits();
    unsigned short getnbModuli();
    void setsecurityBits(uint64_t security_bits);
//    virtual void mulandadd(lwe_cipher rop, lwe_in_data op1, lwe_query op2, uint64_t current_poly, 
//        int rec_lvl)=0;
    virtual poly64* deserializeDataNFL(unsigned char **inArrayOfBuffers, uint64_t nbrOfBuffers, 
        uint64_t dataBitsizePerBuffer, uint64_t &polyNumber)=0;

  protected:
    NFLlib nflInstance;
    uint64_t *moduli;
    uint64_t securityBits;
    unsigned short nbModuli;
};

#endif //DEF_LATICEBASEDCRYPTOSYSTEM
