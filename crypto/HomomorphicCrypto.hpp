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

#ifndef DEF_HOMOMORPHICCRYPTO
#define DEF_HOMOMORPHICCRYPTO

#include <string>
#include <set>
#include <boost/algorithm/string.hpp>
#include "../pir/PIRParameters.hpp"
#include "AbstractPublicParameters.hpp"

#ifndef DEF_CRYPTOGRAPHICSYSTEM
#define DEF_CRYPTOGRAPHICSYSTEM

class CryptographicSystem 
{
  public:
    virtual AbstractPublicParameters& getPublicParameters()=0; 
    virtual unsigned int getAllCryptoParams(std::set<std::string>& crypto_params_set)=0;
    virtual void setNewParameters(const std::string& crypto_param)=0;
    virtual ~CryptographicSystem(){};
    virtual std::string& toString()=0;
};

#endif
class HomomorphicCrypto : public CryptographicSystem {
    
protected:
    std::string cryptoName; 
	
    static unsigned int default_security_bits;
public:
    HomomorphicCrypto(const std::string& crypto_name);
    
    virtual char* encrypt(unsigned int ui, unsigned int )=0;
    virtual char* encrypt(char* data, size_t, unsigned int exponent )=0;
    virtual char* encrypt_perftest()=0;
    virtual char* decrypt(char* cipheredData, unsigned int rec_lvl, size_t, size_t)=0;
    
    virtual unsigned int getCryptoParams(unsigned int k, std::set<std::string>& crypto_params)=0;
    virtual long setandgetAbsBitPerCiphertext(unsigned int elt_nbr)=0;
    virtual std::string getSerializedCryptoParams(bool shortversion=true)=0;
    virtual double estimateAbsTime(std::string crypto_param)=0;
    virtual ~HomomorphicCrypto();
    
    double estimatePrecomputeTime(std::string crypto_param) { return 0;}
	uint64_t getCiphertextBytesize() {
		return 	getPublicParameters().getQuerySizeFromRecLvl(0) / 8;

	}
    std::string& toString() { return cryptoName;}
};

#endif
