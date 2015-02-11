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

#include "NFLLWEPublicParameters.hpp"
#include "NFLLWE.hpp"
#include <string.h>
#include "NFLParams.hpp"
using namespace std;

NFLLWEPublicParameters::NFLLWEPublicParameters()//:
  //polyDegree(0)
{
    cryptoName = "LWE";
    crypto_container=nullptr;
    // -1 means uninitialized, 0 means no absorption possible
    absPerCoordinateBitsize=-1;
    noise_ub=0;
}


NFLLWEPublicParameters::NFLLWEPublicParameters(unsigned int modulusBitsize_, unsigned int polyDegree_, int absPCBitsize_):
    absPerCoordinateBitsize(absPCBitsize_)
{
    crypto_container=nullptr;
    cryptoName = "LWE";
}

// Function that sets the modulus and polyDegree from a raw description sent over the network
// Expected structure of rawPubKey (uint modulusREPRESATIONBitsize):(uint modulus)
void NFLLWEPublicParameters::setModulus(char* rawPubKey)
{
    // MOK 05122013 there is no need for setting the public key in NFLLWE
}

void NFLLWEPublicParameters::setMockedPubKey()
{
    // MOK 05122013 there is no need for setting the public key in NFLLWE
}


// Getters
unsigned int NFLLWEPublicParameters::getmodulusBitsize() { return crypto_container->getmodulusBitsize(); }
uint64_t* NFLLWEPublicParameters::getmoduli() { return crypto_container->getmoduli(); }
uint64_t NFLLWEPublicParameters::getnoiseUB() { return noise_ub; }
uint64_t NFLLWEPublicParameters::getsecurityBits() { return crypto_container->getsecurityBits(); }
unsigned int NFLLWEPublicParameters::getpolyDegree() { return crypto_container->getpolyDegree(); }

// Setters
void NFLLWEPublicParameters::setnoiseUB(uint64_t noise_upper_bound) { noise_ub = noise_upper_bound;}
void NFLLWEPublicParameters::setAbsPCBitsize(int bitSize_)
{
	absPerCoordinateBitsize = bitSize_;
	
}
void NFLLWEPublicParameters::setsecurityBits(uint64_t security_bits_) 
{ 
		crypto_container->setsecurityBits(security_bits_);
}
void NFLLWEPublicParameters::setmodulus(uint64_t modulus_)
{
		crypto_container->setmodulus(modulus_);
}
void NFLLWEPublicParameters::setpolyDegree(unsigned int polyDegree_)
{
		crypto_container->setpolyDegree(polyDegree_);
}

unsigned int NFLLWEPublicParameters::getModulusRepresentationBitsize()
{
  // We represent each 60 bit modulus by a 64 bit integer
  return ceil((double)getmodulusBitsize()/kModulusRepresentationBitsize)*kModulusRepresentationBitsize;
}


unsigned int NFLLWEPublicParameters::getSerializedModulusBitsize()
{
  return getModulusRepresentationBitsize();
}


// Expected format of the parameters
// k:polyDegree:modululusBitsize:AbsorptionBitsize
void NFLLWEPublicParameters::setNewParameters(std::string crypto_param_descriptor)
{
  // We want to get rid of public parameter objects so we transfer the most 
  // we can to the crypto object 
  crypto_container->setNewParameters(crypto_param_descriptor);
}


// Get a serialized version of the parameters
std::string NFLLWEPublicParameters::getSerializedParams(bool shortversion)
{
  std::string params;
  
  // Name:security:degree:modulusbitsize
  // WARNING send modulus representation
  params = cryptoName + ":" + std::to_string(getsecurityBits()) + ":" + std::to_string(getpolyDegree()) + ":" + std::to_string(getmodulusBitsize());

  if (!shortversion)
  {
    // Add :abs_per_coordinate if defined or :? otherwise
    if (absPerCoordinateBitsize==-1) params += ":?";
    else params += ":" + std::to_string(absPerCoordinateBitsize);
  }

  return params;
}


char* NFLLWEPublicParameters::getByteModulus()
{
 char* byte_pub_key  = new char[getpolyDegree() * sizeof(uint64_t)](); 
 memcpy(byte_pub_key, &P64, getpolyDegree() * sizeof(uint64_t));
	return byte_pub_key;
}


void NFLLWEPublicParameters::getParameters()
{
}


unsigned int NFLLWEPublicParameters::getAbsorptionBitsize()
{
	return (absPerCoordinateBitsize < 0) ? 0 : getpolyDegree() * absPerCoordinateBitsize;
}

unsigned int NFLLWEPublicParameters::getAbsorptionBitsize(unsigned int i)
{
	return (absPerCoordinateBitsize < 0) ? 0 : getpolyDegree() * absPerCoordinateBitsize;
}

unsigned int NFLLWEPublicParameters::getCiphertextSize()
{
  return getpolyDegree() * 64 * 2 *crypto_container->getnbModuli();
}

unsigned int NFLLWEPublicParameters::getCiphertextBitsize()
{
	return getModulusRepresentationBitsize() * getpolyDegree() * 2 ;
}
unsigned int NFLLWEPublicParameters::getCiphBitsizeFromRecLvl(unsigned int d)
{
	return getCiphertextBitsize() ;
}

unsigned int NFLLWEPublicParameters::getQuerySizeFromRecLvl(unsigned int)
{
	return getCiphertextBitsize();
}

void NFLLWEPublicParameters::computeNewParameters(const std::string& crypto_param_descriptor)
{
  setNewParameters(crypto_param_descriptor);
}

