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

#ifndef DEF_NFLLWE
#define DEF_NFLLWE

#define SHOUP
//#define TESTSHOUP

#include <omp.h>
#include <inttypes.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "NFLParams.hpp"
#include "NFLlib.hpp"
#include "NFLLWEDatatypes.hpp"
#include "LatticesBasedCryptosystem.hpp"
#include "crypto/HomomorphicCrypto.hpp"
#include "CryptographicSystem.hpp"
#include "NFLLWEPublicParameters.hpp"
#include <string>
#include <cstddef>
#include <gmp.h>

class NFLLWE : public LatticesBasedCryptosystem
{

  public:
    NFLLWEPublicParameters publicParams;
    NFLLWE();
    ~NFLLWE();

    std::string& toString();
    
    unsigned int getpolyDegree();
    poly64* getsecretKey();
  	void recomputeNoiseAmplifiers();
    
    // Setters
    void setmodulus(uint64_t modulus);
    void setpolyDegree(unsigned int polyDegree);
    void setNewParameters(const std::string& crypto_param_descriptor);
    void setNewParameters(unsigned int polyDegree, unsigned int modulusBitsize, int absPCBitsize_);
   
    // Crypto related functions
    long setandgetAbsBitPerCiphertext(unsigned int elt_nbr);
    void enc(lwe_cipher *c, poly64 m);
	  void dec(poly64 m, lwe_cipher *c);	
    char* encrypt(unsigned int ui, unsigned int );
    char* encrypt(char* data, size_t, unsigned int exponent );
    char* encrypt_perftest();
    char* decrypt(char* cipheredData, unsigned int, size_t, size_t);

    // Data importation and exportation
    poly64* deserializeDataNFL(unsigned char **inArrayOfBuffers, uint64_t nbrOfBuffers, 
        uint64_t dataBitsizePerBuffer, uint64_t &polyNumber);
    
    // Functions for PIROptimizer and PIRClient
    std::string getSerializedCryptoParams(bool shortversion);
    unsigned int getCryptoParams(unsigned int k, std::set<std::string>& crypto_params);
    unsigned int getAllCryptoParams(std::set<std::string>& crypto_params);
    AbstractPublicParameters&  getPublicParameters();
    unsigned int findMaxModulusBitsize(unsigned int security_bits, unsigned int poly_degree);
    bool checkParamsSecure(unsigned int security_bits, unsigned int poly_degree, unsigned int p_size);
    double lllOutput(unsigned int n, double& p, double delta);
    double estimateAbsTime(std::string crypto_param);
    double estimatePrecomputeTime(std::string crypto_param);
    unsigned int estimateSecurity(unsigned int n, unsigned int p_size);
	unsigned int getmodulusBitsize();
    
    // **********************************
    // Modular ciphertext manipulation 
    // **********************************

    // Additions
    void add(lwe_cipher rop, lwe_cipher op1, lwe_cipher op2, int d);
    // Fused Multiplications-Additions
    void mulandadd(lwe_cipher rop, lwe_in_data op1, lwe_query op2, int rec_lvl);
    void mulandadd(lwe_cipher rop, lwe_in_data op1, lwe_query op2, uint64_t current_poly, 
        int rec_lvl);
    //Shoup version
	  void mulandadd(lwe_cipher rop, lwe_in_data op1, lwe_query op2, lwe_query op2prime, 
        uint64_t current_poly, int rec_lvl);
	  void mul(lwe_cipher rop, lwe_in_data op1, lwe_query op2, lwe_query op2prime, 
        uint64_t current_poly, int rec_lvl);

    void mulandaddCiphertextNTT(lwe_cipher rop, lwe_in_data op1, lwe_query op2);
    void mulandaddCiphertextNTT(lwe_cipher rop, lwe_in_data op1, lwe_query op2, 
        uint64_t current_poly);


  private:
    // Attributes
    unsigned int oldNbModuli;
    unsigned int polyDegree;
    poly64 *secretKey; // The secret key
    poly64 *secretKeyShoup; // The secret key Shoupified
	uint64_t *Abit_mod,*Abit_mod_shoup;

    void clearSecretKeys();
};

#endif
