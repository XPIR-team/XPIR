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

#ifndef DEF_NOCRYPTOGRAPHY
#define DEF_NOCRYPTOGRAPHY

#include "../pir/PIRParameters.hpp"
#include "HomomorphicCrypto.hpp"
#include "NoCryptographyPublicParameters.hpp"
#include "AbstractPublicParameters.hpp"


class NoCryptography : public HomomorphicCrypto 
{
  public:
    NoCryptography();
    NoCryptography(const std::string& crypto_name);
    AbstractPublicParameters& getPublicParameters(); 
    unsigned int getAllCryptoParams(std::set<std::string>& crypto_params_set);
    void setNewParameters(const std::string& crypto_param);
    std::string& toString();
    PIRParameters* pirParam;
    NoCryptographyPublicParameters publicParams;
    static unsigned int default_security_bits;
    
    char* encrypt(unsigned int ui, unsigned int );
    char* encrypt(char* data, size_t, unsigned int exponent );
    char* encrypt_perftest();
    char* decrypt(char* cipheredData, unsigned int rec_lvl, size_t, size_t);
    
    unsigned int getCryptoParams(unsigned int k, std::set<std::string>& crypto_params);
    long setandgetAbsBitPerCiphertext(unsigned int elt_nbr);
    std::string getSerializedCryptoParams(bool shortversion=true);
    double estimateAbsTime(std::string crypto_param);
    
};

#endif
