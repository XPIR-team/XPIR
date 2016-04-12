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

#include "NoCryptography.hpp"

 
NoCryptography::NoCryptography(): 
  HomomorphicCrypto("NoCryptography")
{}

NoCryptography::NoCryptography(const std::string& crypto_name):
  HomomorphicCrypto(crypto_name)
{}

AbstractPublicParameters& NoCryptography::getPublicParameters(){ return publicParams; } 

unsigned int NoCryptography::getAllCryptoParams(std::set<std::string>& crypto_params_set) 
{ 
  crypto_params_set.insert(std::string("NoCryptography"));
  return 1; 
}

void NoCryptography::setNewParameters(const std::string& crypto_param) {}
std::string& NoCryptography::toString() { return cryptoName; }

char* NoCryptography::encrypt(unsigned int ui, unsigned int d )
{ 
  char* charptr = (char *) malloc(sizeof(char)); 
  *charptr =  (char) ui;
  return charptr; 
}

char* NoCryptography::encrypt(char* data, size_t lkjlkj, unsigned int exponent ) 
{
  unsigned int ciphBytesize = publicParams.getCiphertextBitsize()/8;
  char* encrypted_data = (char *) malloc(ciphBytesize);
  memcpy(encrypted_data, data, ciphBytesize);
  return encrypted_data;
}

char* NoCryptography::encrypt_perftest()
{
  return encrypt(1,1);
}

char* NoCryptography::decrypt(char* cipheredData, unsigned int rec_lvl, size_t, size_t)
{
  unsigned int clearBytesize = publicParams.getAbsorptionBitsize()/8;
  char* clear_data = (char *) malloc(clearBytesize);
  memcpy(clear_data, cipheredData, clearBytesize);
  return clear_data;
}

unsigned int NoCryptography::getCryptoParams(unsigned int k, std::set<std::string>& crypto_params)
{
  crypto_params.insert(std::string("NoCryptography"));
  return 1; 
}

long NoCryptography::setandgetAbsBitPerCiphertext(unsigned int elt_nbr)
{
  return 8*1024;
}

std::string NoCryptography::getSerializedCryptoParams(bool shortversion)
{
  return cryptoName;
}

double NoCryptography::estimateAbsTime(std::string crypto_param){
  return 1;
}
//NoCryptography::~NoCryptography() {}

