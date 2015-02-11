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

#include "PaillierPublicParameters.hpp"
#include <iostream>

PaillierPublicParameters::PaillierPublicParameters()
{
  cryptoName = "Paillier";
  securityBits = 0;
}
//PaillierPublicParameters::PaillierPublicParameters(unsigned int bit_key_size) :
//	bitKeySize(bit_key_size),
//	bitAbsSize(bitKeySize - 8),
//	ciphSize(bitKeySize * 2)
//{ 
//  cryptoName = "Paillier";
//  securityBits = 0;
//}

PaillierPublicParameters::~PaillierPublicParameters()
{
}

char* PaillierPublicParameters::getByteModulus() {
  char *key = new char[(getKeyBitsize() / 8 ) + sizeof(int)](); 
  mpz_export(key, NULL, 1, sizeof(char) , 0, 0, *pubkey.getnj(1));
  int init_s = pubkey.getinit_s();
  memcpy(key+(getKeyBitsize()/8), &init_s, sizeof(int));
	return key;
}

unsigned int PaillierPublicParameters::getCiphertextSize() {
	return ciphSize; 
}

unsigned int PaillierPublicParameters::getKeyBitsize() { 
	return bitKeySize;
}

unsigned int PaillierPublicParameters::getAbsorptionBitsize() { 
	return bitAbsSize;
}

unsigned int PaillierPublicParameters::getAbsorptionBitsize(unsigned int rec_lvl)
{
  if (rec_lvl == 0)
  {
    return getAbsorptionBitsize();
  }
  else
  {
    // If we are being used in a recursive scheme take all the bits
    return getKeyBitsize()*(rec_lvl+1);
  }
}

paillier_pubkey* PaillierPublicParameters::getPubKey() {
	return &pubkey;
}

void PaillierPublicParameters::setModulus(char* rawPubKey){
  pubkey.init_key(getKeyBitsize(), rawPubKey);	
}

//lvl starts at 1 !!!
unsigned int PaillierPublicParameters::getCiphBitsizeFromRecLvl(unsigned int lvl) 
{
  return (lvl+pubkey.getinit_s()) * bitKeySize;
}

unsigned int PaillierPublicParameters::getQuerySizeFromRecLvl(unsigned int lvl)
{
  return getCiphBitsizeFromRecLvl(lvl);
}

unsigned int PaillierPublicParameters::getSerializedModulusBitsize()
{
  return  getKeyBitsize() + sizeof(int) * 8;
}

void PaillierPublicParameters::setModulusbits(unsigned int modulusbits_)
{
	bitKeySize  = modulusbits_ ;
	bitAbsSize = (bitKeySize - 8);
	ciphSize    = (bitKeySize * 2);
}

void PaillierPublicParameters::setSecurityBits(unsigned int securitybits_)
{
  securityBits = securitybits_;
} 

void PaillierPublicParameters::computeNewParameters(const std::string& crypto_param_descriptor)
{
  std::vector<std::string> fields;
  boost::algorithm::split(fields, crypto_param_descriptor, boost::algorithm::is_any_of(":"));
  securityBits = (unsigned)atoi (fields[1].c_str());
  unsigned int keySize = (unsigned)atoi(fields[2].c_str());
  setModulusbits(keySize);
}


// Get a serialized version of the parameters
std::string PaillierPublicParameters::getSerializedParams(bool shortversion)
{
  std::string params;
  
  // Name:security:plaintext_modulusbitsize:ciphertext_modulusbitsize:abs_bits
  params = cryptoName + ":" + std::to_string(securityBits) + ":" + std::to_string(bitKeySize) + ":" + std::to_string(ciphSize);
  
  if (!shortversion) params += ":" + std::to_string(bitAbsSize);

  return params;
}


unsigned int PaillierPublicParameters::getCiphertextBitsize()
{
  return getCiphBitsizeFromRecLvl(1);
}

void PaillierPublicParameters::setMockedPubKey()
{
  pubkey.init_key(getKeyBitsize());
}
    

