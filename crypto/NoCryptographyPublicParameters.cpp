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

#include "NoCryptographyPublicParameters.hpp"

using namespace std;

NoCryptographyPublicParameters::NoCryptographyPublicParameters()
{
    cryptoName = "NoCryptography";
}


string NoCryptographyPublicParameters::getSerializedParams(bool shortversion) { return string("NoCryptography"); }
char* NoCryptographyPublicParameters::getByteModulus() { return (char*) malloc(1); }
unsigned int NoCryptographyPublicParameters::getAbsorptionBitsize(){ return 8*1024; }
unsigned int NoCryptographyPublicParameters::getAbsorptionBitsize(unsigned int rec_lvl){ return 8*1024; }
void NoCryptographyPublicParameters::setModulus(char* rawPubKey){}
void NoCryptographyPublicParameters::setMockedPubKey(){}
unsigned int NoCryptographyPublicParameters::getCiphertextBitsize(){ return 8*1024; }
unsigned int NoCryptographyPublicParameters::getCiphBitsizeFromRecLvl(unsigned int){ return 8*1024; }
unsigned int NoCryptographyPublicParameters::getQuerySizeFromRecLvl(unsigned int){ return 64;}//64 bits as an uint64_t
void NoCryptographyPublicParameters::computeNewParameters(const std::string& crypto_param_descriptor){}
unsigned int NoCryptographyPublicParameters::getSerializedModulusBitsize(){ return 0;}


