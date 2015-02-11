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

#ifndef DEF_PAILLIERPUBLICPARAMETERS
#define DEF_PAILLIERPUBLICPARAMETERS

#include <string.h>

#include "AbstractPublicParameters.hpp"
#include "PaillierKeys.hpp"

class PaillierPublicParameters : public AbstractPublicParameters {

	private:
		unsigned int bitKeySize;
		unsigned int byteKeySize;
		unsigned int bitAbsSize;
		unsigned int ciphSize;
    unsigned int securityBits;
		paillier_pubkey pubkey;

	public:
		unsigned int getCiphBitsizeFromRecLvl(unsigned int);
    unsigned int getCiphertextBitsize();
		unsigned int getQuerySizeFromRecLvl(unsigned int);
		PaillierPublicParameters();
		~PaillierPublicParameters();
    void setModulusbits(unsigned int modulusbits);
    void setSecurityBits(unsigned int securitybits);
    void computeNewParameters(const std::string& crypto_param_descriptor);
		char* getByteModulus();
    std::string getSerializedParams(bool shortversion);
		unsigned int getKeyBitsize();
		unsigned int getAbsorptionBitsize();
    unsigned int getAbsorptionBitsize(unsigned int i);

		unsigned int getCiphertextSize();
    unsigned int getSerializedModulusBitsize();
		paillier_pubkey* getPubKey();
		void setModulus(char* pubKey);
    void setMockedPubKey();

		
};

#endif
