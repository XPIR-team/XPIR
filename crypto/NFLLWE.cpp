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

#include "NFLLWE.hpp"
#include <fstream> 
//#define bench
//#define Repetition 10000

void NFLLWE_DEBUG_MESSAGE(const char *s,poly64 p, unsigned int n){
#ifdef CRYPTO_DEBUG
	std::cout<<s;
	NFLlib::print_poly64hex(p,n);
#endif
}

// *********************************************************
// Constructors and initialization
// The constructors are not able to set all the parameters 
// and setNewParameters has to be called afterward, 
// the attribute alreadyInit reflects this uninitialized 
// status
// *********************************************************


NFLLWE::NFLLWE():
    LatticesBasedCryptosystem("LWE"),
    oldNbModuli(0),
    polyDegree(0)
{
  publicParams.setcrypto_container(this);
}


// Expected format of the parameters
// k:polyDegree:modululusBitsize:AbsorptionBitsize
void NFLLWE::setNewParameters(const std::string& crypto_param_descriptor)
{
  unsigned int polyDegree_, aggregatedModulusBitsize_;
  int abspc_bitsize = -1; // We don't know the absorption bit size yet

  std::vector<std::string> fields;
  boost::algorithm::split(fields, crypto_param_descriptor, boost::algorithm::is_any_of(":"));

  setsecurityBits(atoi(fields[1].c_str()));
  polyDegree_ = atoi(fields[2].c_str());
  aggregatedModulusBitsize_ = atoi(fields[3].c_str());
  // Does the fourth parameter exist ? If so set it 
  if (fields.size() >= 5) abspc_bitsize = atoi(fields[4].c_str()); 

  setNewParameters(polyDegree_,aggregatedModulusBitsize_, abspc_bitsize);
}


// The setNewParameters method does the actual parameterization of the crypto object
// it sets the alreadyInit attribute to reflects this
void  NFLLWE::setNewParameters(unsigned int polyDegree_, unsigned int aggregatedModulusBitsize_, int absPCBitsize_)
{
	// Our public parameters need a pointer on us	
  publicParams.setcrypto_container(this);

	// We still need to transfer this two attributes to the crypto_object
	// for the transition towards public parameter elimination
	publicParams.setAbsPCBitsize(absPCBitsize_);

	publicParams.setnoiseUB(5*getsecurityBits()/2);

//#ifdef DEBUG
//  std::cout << "Security bits " << getsecurityBits()<<std::endl;
//  std::cout << "Noise UB " << publicParams.getnoiseUB()<<std::endl;
//#endif

  // We don't use here the polyDegree setter as we would call twice NFLlib init
	polyDegree = polyDegree_;

  nflInstance.setNewParameters(polyDegree_,aggregatedModulusBitsize_);
  clearSecretKeys();
  nbModuli = nflInstance.getnbModuli();
  //used to free memory
  oldNbModuli = nbModuli;
  moduli= nflInstance.getmoduli();


	secretKey = new poly64[nbModuli];
	secretKeyShoup = new poly64[nbModuli];
	Abit_mod = new uint64_t[nbModuli];
	Abit_mod_shoup = new uint64_t[nbModuli];


	// initialize the secret key 
  secretKey[0] = nflInstance.allocBoundedRandomPoly(0,true);
	for (unsigned short currentModulus = 0; currentModulus < nbModuli; currentModulus++) {
	  secretKey[currentModulus] = secretKey[0] + polyDegree*currentModulus;
    secretKeyShoup[currentModulus] = (uint64_t*) calloc(polyDegree,sizeof(uint64_t));
		// compute the Shoup representation of the secret key
		for (unsigned int i=0; i < polyDegree; i++) {
			secretKeyShoup[currentModulus][i]=((uint128_t) secretKey[currentModulus][i] << 64) / moduli[currentModulus];
    }
  }
 recomputeNoiseAmplifiers();
  
}

// *********************************************************
// Getters
// *********************************************************
poly64* NFLLWE::getsecretKey() { return secretKey; }
unsigned int NFLLWE::getpolyDegree() { return polyDegree; }

// *********************************************************
// Setters
// *********************************************************
void NFLLWE::setmodulus(uint64_t modulus_)
{
	// The modulus cannot be set from outside
	std::cout << "Warning(NFLLWE.c): Modulus cannot be set externally." << std::endl;
}
void NFLLWE::setpolyDegree(unsigned int polyDegree_)
{
  polyDegree = polyDegree_;
  nflInstance.setpolyDegree(polyDegree_);
}

// *********************************************************
//         Serialize/Deserialize
// *********************************************************

poly64 *NFLLWE::deserializeDataNFL(unsigned char **inArrayOfBuffers, uint64_t nbrOfBuffers, uint64_t dataBitsizePerBuffer, uint64_t &polyNumber) {
  return nflInstance.deserializeDataNFL(inArrayOfBuffers, nbrOfBuffers, dataBitsizePerBuffer, publicParams.getAbsorptionBitsize()/polyDegree, polyNumber);
}



// *********************************************************
//         Additions and Multiplications of ciphertexts
// *********************************************************


void NFLLWE::add(lwe_cipher rop, lwe_cipher op1, lwe_cipher op2, int d)
{
  nflInstance.addmodPoly(rop.a, op1.a, op2.a);
  nflInstance.addmodPoly(rop.b, op1.b, op2.b);
}

void NFLLWE::mulandadd(lwe_cipher rop, lwe_in_data op1, lwe_query op2, uint64_t current_poly, int rec_lvl)
{
	NFLLWE_DEBUG_MESSAGE("in_data[0].p : ",op1.p[0],4);
  	NFLLWE_DEBUG_MESSAGE("in_data[0].a : ",op2.a,4);
  	NFLLWE_DEBUG_MESSAGE("in_data[0].b : ",op2.b,4);
	
	mulandaddCiphertextNTT(rop, op1, op2, current_poly);
	
	NFLLWE_DEBUG_MESSAGE("out_data[0].a : ",rop.a,4);
  	NFLLWE_DEBUG_MESSAGE("out_data[0].b : ",rop.b,4);
}

// Shoup version
void NFLLWE::mulandadd(lwe_cipher rop, const lwe_in_data op1, const lwe_query op2, const lwe_query op2prime, const uint64_t current_poly, int rec_lvl)
{
  // Don't modify the pointers inside the data or it will be permanent
  poly64 ropa = rop.a, ropb = rop.b, op2a = op2.a, op2b = op2.b, op2primea = op2prime.a, 
         op2primeb = op2prime.b, op1pcurrent = op1.p[current_poly];

	
	 	const unsigned int K = polyDegree;
		const unsigned int md = nbModuli;
	for(unsigned short currentModulus=0;currentModulus<md;currentModulus++) 
  {
	  
 		for (unsigned i = 0; i < K; i++)
		{
			nflInstance.mulandaddShoup(ropa[i],op1pcurrent[i],op2a[i],op2primea[i],moduli[currentModulus]);
		}
 		for (unsigned i = 0; i < K; i++)
		{
			nflInstance.mulandaddShoup(ropb[i],op1pcurrent[i],op2b[i],op2primeb[i],moduli[currentModulus]);
		}
		ropa+=K;
		ropb+=K;
		op1pcurrent+=K;
		op2a+=K;
		op2b+=K;
		op2primea+=K;
		op2primeb+=K;
	}

}

void NFLLWE::mul(lwe_cipher rop, const lwe_in_data op1, const lwe_query op2, const lwe_query op2prime, const uint64_t current_poly, int rec_lvl)
{
  // Don't modify the pointers inside the data or it will be permanent
  poly64 ropa = rop.a, ropb = rop.b, op2a = op2.a, op2b = op2.b, op2primea = op2prime.a, 
         op2primeb = op2prime.b, op1pcurrent = op1.p[current_poly];

	NFLLWE_DEBUG_MESSAGE("in_data[0].p : ",op1.p[current_poly],4);
	NFLLWE_DEBUG_MESSAGE("in_data[0].a : ",op2.a,4);
	NFLLWE_DEBUG_MESSAGE("in_data[0].b : ",op2.b,4);
	NFLLWE_DEBUG_MESSAGE("in_data[0].a' : ",op2prime.a,4);
	NFLLWE_DEBUG_MESSAGE("in_data[0].b' : ",op2.b,4);
	
	for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) 
  {
		for (unsigned i = 0; i < polyDegree; i++)
		{
			ropa[i] = nflInstance.mulmodShoup(op1pcurrent[i],op2a[i],op2primea[i],moduli[currentModulus]);
			ropb[i] = nflInstance.mulmodShoup(op1pcurrent[i],op2b[i],op2primeb[i],moduli[currentModulus]);
		}
		ropa+=polyDegree;
		ropb+=polyDegree;
		op1pcurrent+=polyDegree;
		op2a+=polyDegree;
		op2b+=polyDegree;
		op2primea+=polyDegree;
		op2primeb+=polyDegree;
	}
	NFLLWE_DEBUG_MESSAGE("out_data[0].a : ",rop.a,4);
	NFLLWE_DEBUG_MESSAGE("out_data[0].b : ",rop.b,4);
}

// Same comment as for musAndAddCiphertextNTT we do a simpler version above
void NFLLWE::mulandadd(lwe_cipher rop, lwe_in_data op1, lwe_query op2, int rec_lvl)
{
	NFLLWE_DEBUG_MESSAGE("in_data p: ",op1.p[0],4);
	NFLLWE_DEBUG_MESSAGE("in_data a: ",op2.a,4);
	NFLLWE_DEBUG_MESSAGE("in_data b: ",op2.b,4);
	
  	mulandaddCiphertextNTT(rop, op1, op2);
	
	NFLLWE_DEBUG_MESSAGE("out_data.a : ",rop.a,4);
	NFLLWE_DEBUG_MESSAGE("out_data.b : ",rop.b,4);
}

// Deal just with one polynomial
inline void NFLLWE::mulandaddCiphertextNTT(lwe_cipher rop, lwe_in_data op1, lwe_query op2, uint64_t current_poly)
{
    nflInstance.mulandaddPolyNTT(rop.a, op1.p[current_poly], op2.a);
    nflInstance.mulandaddPolyNTT(rop.b, op1.p[current_poly], op2.b);
}

// Good method but too greedy in memory we start with a simpler one (below)
// Needs to change as we always write in the same rop
void NFLLWE::mulandaddCiphertextNTT(lwe_cipher rop, lwe_in_data op1, lwe_query op2)
{
  for(uint64_t i=0;i<op1.nbPolys;i++)
  {
    nflInstance.mulandaddPolyNTT(rop.a, op1.p[i], op2.a);
    nflInstance.mulandaddPolyNTT(rop.b, op1.p[i], op2.b);
  }
}



//*********************************
// Encryption and decryption 
//********************************* 

// The internal encrypt method
void  NFLLWE::enc(lwe_cipher *c, poly64 m)
{
	bool uniform = true;

	NFLLWE_DEBUG_MESSAGE("Encrypting m: ",m, 4);

	c->a = (poly64) calloc(polyDegree * 2 * nbModuli,  sizeof(uint64_t));
	c->b = c->a + polyDegree * nbModuli;

	// tmpa and tmpb are used to access the nbModuli polynoms of the CRT
	poly64 tmpa = c->a;
	poly64 tmpb = c->b;
	poly64 tmpm = m;

	//   b = (a*s) % f + e * A + m;
    
	// Noise creation
	uint64_t Berr=publicParams.getnoiseUB();
	uint64_t A_bits= publicParams.getAbsorptionBitsize() / publicParams.getpolyDegree();	
	
	// We deal with the nbModuli polynoms at once because the noise is the same size for all of them
	nflInstance.setBoundedRandomPoly(c->b, 2*Berr-1, !uniform);

	NFLLWE_DEBUG_MESSAGE("Noise used: ",c->b, 4);
#ifdef CRYPTO_DEBUG
	std::cout << "NFLLWE: Noise amplifier: " << A_bits << std::endl;
#endif 
		
	// Adjustments and addition to plaintext
	for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) {
		for(unsigned int i=0;i<polyDegree;i++) {
			// e is multiplied by the amplifier A for which we know the size A_bits
			// tmpb[i] = tmpb[i] << (unsigned) A_bits;
			//  std::cout << "noise: " << tmpb[i] << std::endl;
			//std::cout << std::hex << tmpb[i] << " " << std::dec;			
				
			tmpb[i] = nflInstance.mulmodShoup(tmpb[i], Abit_mod[currentModulus],Abit_mod_shoup[currentModulus], moduli[currentModulus]);
      
			// and shifted to be in [-(Berr-1) .. (Berr-1)]
			//tmpb[i] += moduli[currentModulus]-((Berr-1)<<A_bits);   

			// We add the shifted noise to the plaintext
			tmpb[i] = nflInstance.addmod(tmpb[i], tmpm[i], moduli[currentModulus]);
      
			// And reduce the whole if needed 
			if(tmpb[i]>moduli[currentModulus]) tmpb[i]-=moduli[currentModulus];
			
		}
		tmpb+=polyDegree;
		tmpm+=polyDegree;
	}
	tmpb=c->b;

	NFLLWE_DEBUG_MESSAGE("Amplified noise and message: ",c->b, 4);

	// Noise and plaintext are the only things that are not yet in the NTT space
	nflInstance.nttAndPowPhi(c->b);

	// We still have to get a. No NTT needed because uniformly taken
	nflInstance.setBoundedRandomPoly(tmpa, 0, uniform);

#ifdef DEBUG
	poly64 tmp = (poly64) calloc(polyDegree*nbModuli, sizeof(uint64_t));
#endif

	for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) 
	{
    
		// We multiply it by s and add to the previous message and noise
		for (unsigned int i = 0 ; i < polyDegree ; i++)
		{
			nflInstance.mulandaddShoup(tmpb[i],tmpa[i], secretKey[currentModulus][i],
			secretKeyShoup[currentModulus][i], moduli[currentModulus]);

#ifdef DEBUG
			nflInstance.mulandaddShoup(tmp[i+currentModulus*polyDegree], tmpa[i],secretKey[currentModulus][i],secretKeyShoup[currentModulus][i], moduli[currentModulus]);
#endif

#ifdef SHOUP
			tmpa[i]=tmpa[i]%moduli[currentModulus];
			tmpb[i]=tmpb[i]%moduli[currentModulus];
#endif
		}
		tmpa+=polyDegree;
		tmpb+=polyDegree;
	}

	// There is already a ifdef debug inside this function but 
	// tmp is not defined if we are not in debug mode
#ifdef DEBUG 
	NFLLWE_DEBUG_MESSAGE("a*s: ",tmp, 4);
  free(tmp);
#endif

	NFLLWE_DEBUG_MESSAGE("Ciphertext a: ",c->a, 4);
	NFLLWE_DEBUG_MESSAGE("Ciphertext b: ",c->b, 4);
}


void NFLLWE::dec(poly64 m, lwe_cipher *c)
{
	uint64_t A_bits = publicParams.getAbsorptionBitsize() / publicParams.getpolyDegree();
  const uint64_t bitmask = (1ULL<<A_bits) -1; 	
  mpz_t moduliProduct; 
	
  // Get the product of all moduli from the nflInstance object;
  nflInstance.copymoduliProduct(moduliProduct);

	// tmpa and tmpb are used to access the nbModuli polynoms of the CRT
	poly64 tmpa=c->a;
	poly64 tmpb=c->b;
	poly64 tmpm=m;

	for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) {

		// We first get the amplified noise plus message (e*A+m =b-a*S)
		for (unsigned int i=0 ; i < polyDegree; i++) 
    {
			uint64_t temp=0;
			nflInstance.mulandaddShoup(temp, tmpa[i], secretKey[currentModulus][i],
      secretKeyShoup[currentModulus][i], moduli[currentModulus]);
			tmpm[i] = nflInstance.submod(tmpb[i], temp, moduli[currentModulus]);
		}
		tmpa+=polyDegree;
		tmpb+=polyDegree;
    tmpm+=polyDegree;
	}
  tmpm=m;

	// In order to mask the noise bits we need to get out of NTT space through an inverse NTT
  nflInstance.invnttAndPowInvPhi(tmpm);

  NFLLWE_DEBUG_MESSAGE("Amplified noise and message (dec): ",tmpm, 4);
  NFLLWE_DEBUG_MESSAGE("Amplified noise and message (dec): ",tmpm+polyDegree, 4);


  if(nbModuli>1) {
	   	  
	  mpz_t *tmprez=nflInstance.poly2mpz(tmpm);
	  
	  

  	// If e*A+m < p/2 we mask the message bits: bitmask = (1ULL<<A_bits) -1
	  // If e *A+m > p/2 we do a little trick to avoid signed integers and modulus reduction
	  // e[i]= e[i] + 2**61 - p (we replace p by 2**61) and then bitmask the message.
    mpz_t magicConstz;
    mpz_init(magicConstz);
    mpz_ui_pow_ui(magicConstz, 2, (kModulusBitsize + 1) * nbModuli);
    mpz_sub(magicConstz,magicConstz, moduliProduct);
    mpz_t bitmaskz;
    mpz_init(bitmaskz);
    mpz_ui_pow_ui(bitmaskz, 2, A_bits);
    mpz_sub_ui(bitmaskz, bitmaskz, 1);
#ifdef CRYPTO_DEBUG
    gmp_printf("Mask used: %Zx\n",bitmaskz);
#endif	
	  // Shall we prefetch here ?
    mpz_t tmpz;
    mpz_init(tmpz);
    // We need to zero tmpm as export writes nothing on the output for null values
    bzero(tmpm,polyDegree*nbModuli*sizeof(uint64_t));
    for (unsigned int i = 0 ; i < polyDegree ; i++)
	  {
		  //For testing we may do a hardcoded modulus but not always. m[i] = m[i] % modulus;
      mpz_mul_ui(tmpz, tmprez[i], 2UL);
      if (mpz_cmp(tmpz, moduliProduct)==1)// tmprez[i] > moduliProduct / 2
      {
        mpz_add(tmpz, tmprez[i], magicConstz);
        mpz_and(tmprez[i], tmpz, bitmaskz);
      }
      else
      {
        mpz_and(tmprez[i], tmprez[i], bitmaskz);
      }
	
    	// Combien d'uint32 ?
	    int combien = ceil((double)A_bits/32);
	    mpz_export(((uint32_t*)tmpm)+i*combien, NULL, -1, sizeof(uint32_t), 0, 0, tmprez[i]);
	    mpz_clear(tmprez[i]);
	  }
    
    delete[] tmprez;
    mpz_clears(moduliProduct, tmpz, magicConstz, bitmaskz, NULL);

  } else { // nbModuli=1
	
	
		// If e*A+m < p/2 we mask the message bits: bitmask = (1ULL<<A_bits) -1
		// If e*A+m > p/2 we do a little trick to avoid signed integers and modulus reduction
		// e[i]= e[i] + 2**61 - p (we replace p by 2**61) and then bitmask the message.
		const uint64_t magicConst = (1ULL<<61)-moduli[0];// 2**61 - p
		
		// Shall we prefetch here ?
		for (unsigned int i = 0 ; i < polyDegree ; i++)
		{
			//For testing we may do a hardcoded modulus but not always. m[i] = m[i] % modulus;
			tmpm[i] = (tmpm[i] > moduli[0]/2) ? (tmpm[i] + magicConst)& bitmask : tmpm[i] & bitmask;
		}	
	}
}

// MOK is here for the CRT modification

// encrypts a uint (e.g. for producing a equest element with a 0 or a 1)
// does not return a lwe_cipher but the (char*)pointer on two consecutively allocated poly64 (a and b)
char* NFLLWE::encrypt(unsigned int ui, unsigned int d)
{
	if ( ceil(log2(static_cast<double>(ui))) >= publicParams.getAbsorptionBitsize())
	{
		std::cerr << "NFFLWE: The given unsigned int does not fit in " << publicParams.getAbsorptionBitsize() << " bits"<< std::endl;
		ui %= 1<<publicParams.getAbsorptionBitsize();
	}

	lwe_cipher c; 
	poly64 m = (poly64)calloc(nbModuli*polyDegree,sizeof(uint64_t));
	for (unsigned int cm = 0 ; cm < nbModuli ; cm++)
  {
    m[cm*polyDegree]=(uint64_t)ui;
  }
	enc(&c,m);
	free(m);
	return (char*) c.a;
}

char* NFLLWE::encrypt(char* data, size_t s, unsigned int exponent ){
    std::cerr << "char* NFLLWE::encrypt(char* data, size_t, unsigned int exponent) is not implemented"<< std::endl;
    return nullptr;
}

// Do a ciphertext for a plaintext with alternating bits (for performance tests) 
char* NFLLWE::encrypt_perftest()
{
	lwe_cipher c; 
  poly64 m = nflInstance.allocBoundedRandomPoly(0, true);
	enc(&c,m);
	free(m);
	return (char*) c.a;
}

char* NFLLWE::decrypt(char* cipheredData, unsigned int rec_lvl, size_t, size_t)
{
  lwe_cipher ciphertext;
  ciphertext.a = (poly64)cipheredData;
  ciphertext.b = ciphertext.a + nbModuli * polyDegree;
  poly64 clear_data = (poly64) calloc(nbModuli * polyDegree, sizeof(uint64_t));
  unsigned int bits_per_coordinate = publicParams.getAbsorptionBitsize()/polyDegree;
  
#ifdef DEBUG
  std::cout<<"Allocated (bytes): "<<nbModuli * polyDegree * sizeof(uint64_t)<<std::endl;
  std::cout<<"Bits per coordinate: "<<bits_per_coordinate<<std::endl;
#endif

  dec(clear_data, &ciphertext);

  NFLLWE_DEBUG_MESSAGE("Decrypting ciphertext a: ",ciphertext.a, 4);
  NFLLWE_DEBUG_MESSAGE("Decrypting ciphertext b: ",ciphertext.b, 4);
  NFLLWE_DEBUG_MESSAGE("Result: ",clear_data, 4);

  // unsigned char* out_data = (unsigned char*) calloc(nbModuli * polyDegree+1, sizeof(uint64_t));
  // nflInstance.serializeData64 (clear_data, out_data, bits_per_coordinate, polyDegree);

  unsigned char* out_data = (unsigned char*) calloc(bits_per_coordinate*polyDegree/64 + 1, sizeof(uint64_t));
  if (nbModuli == 1)
  {
    nflInstance.serializeData64(clear_data, out_data, bits_per_coordinate, ceil((double)bits_per_coordinate/64)* polyDegree);
  }
  else // nbModuli > 1
  {
    nflInstance.serializeData32 ((uint32_t*)clear_data, out_data, bits_per_coordinate, ceil((double)bits_per_coordinate/32)* polyDegree);
  }
#ifdef DEBUG
  //std::cout<<"Bitgrouped into: "<<out_data<<std::endl;
#endif
  free(clear_data);
  return (char*) out_data;
}


unsigned int NFLLWE::getAllCryptoParams(std::set<std::string>& crypto_params)
{
  unsigned int params_nbr  = 0;
  unsigned int k_array_size = 5;
  unsigned int k[5] = {80, 100, 128, 192, 256};

  for (unsigned int i = 0 ; i < k_array_size ; i++)
  {
    params_nbr += getCryptoParams(k[i], crypto_params);
  }

  return params_nbr;
}


unsigned int NFLLWE::getCryptoParams(unsigned int k, std::set<std::string>& crypto_params)
{
  using namespace std;
  unsigned int p_size, params_nbr = 0;
  string k_str  = to_string(k);

  for (unsigned int degree = kMinPolyDegree ; degree <= kMaxPolyDegree; degree <<= 1)
  {
    string param;
    p_size = findMaxModulusBitsize(k, degree);
    
    // We give a very small margin 59 instead of 60 so that 100:1024:60 passes the test
    //for (unsigned int i = 1; i * 59 <= p_size ; i++)//(p_size > 64) && ((p_size % 64) != 0))
    for (unsigned int i = 1; i * 59 <= p_size && i * 60 <= 240; i++)
    {
      param =  cryptoName + ":" + to_string(estimateSecurity(degree,i*kModulusBitsize)) + ":" + to_string(degree) + ":" + to_string(i*kModulusBitsize) ;
      if (crypto_params.insert(param).second) params_nbr++;
      param = "";
    }
  }

  return params_nbr;
}

void NFLLWE::recomputeNoiseAmplifiers() {
	uint64_t A_bits= publicParams.getAbsorptionBitsize() / publicParams.getpolyDegree();	
	mpz_t tmpz1,tmpz2; 
	mpz_init(tmpz1);
	mpz_init(tmpz2);
	for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) {
		mpz_ui_pow_ui(tmpz2, 2, A_bits);
		mpz_import(tmpz2, 1, 1, sizeof(uint64_t), 0, 0, moduli+currentModulus);
		mpz_mod(tmpz1, tmpz1, tmpz2);
		Abit_mod[currentModulus]=0;
		mpz_export(&Abit_mod[currentModulus], NULL, 1, sizeof(uint64_t), 0, 0, tmpz1);
		Abit_mod_shoup[currentModulus]=((uint128_t) Abit_mod[currentModulus] << 64) / moduli[currentModulus];
	}
  mpz_clears(tmpz1, tmpz2, NULL);
}

unsigned int NFLLWE::estimateSecurity(unsigned int n, unsigned int p_size)
{
  unsigned int estimated_k = 5;//Estimate K can not be too low

  while(!checkParamsSecure(estimated_k,n,p_size)) estimated_k++; 

  return --estimated_k;
}


long NFLLWE::setandgetAbsBitPerCiphertext(unsigned int elt_nbr)
{
    double Berr = static_cast<double>(publicParams.getnoiseUB());
    double nb_sum = elt_nbr;
    double p_size = getmodulusBitsize();
    double nbr_bit = floor(( (p_size - 1) - log2(nb_sum) - log2(Berr) -log2(static_cast<double>(polyDegree))) / 2.0);

    publicParams.setAbsPCBitsize(nbr_bit);
	
	recomputeNoiseAmplifiers();
	
    return long(nbr_bit);
}


unsigned int NFLLWE::findMaxModulusBitsize(unsigned int k, unsigned int n)
{
  unsigned int p_size;
  //p_size can not be too low
  p_size = 10;
  while (!checkParamsSecure(k,n,p_size)) p_size++;

  return --p_size;
}


bool NFLLWE::checkParamsSecure(unsigned int k, unsigned int n, unsigned int p_size)
{
  double p, beta, logBerr = 8, epsi, lll;

  //We take an advantage of 2**(-k/2) and an attack time of 2**(k/2)
  epsi = pow(2, -static_cast<double>(k/2));
  //log(time) = 1.8/ log(delta) − 110 and -80 to compute processor cycles so we take pow(2, k/2) = 1.8/log(delta) - 80
  double delta = pow(2,1.8/(k/2 + 80));

  p    = pow(2, p_size) -  1;
  beta = (p / logBerr) * sqrt(log1p( 1 / epsi) / M_PI);
  lll  = lllOutput(n, p, delta);

  // We love ugly tricks !
  return (lll < beta);// && cout << "beta : " << beta << " p_size : " << p_size << " n :"<< n << " k : "<< k << endl;
}


double NFLLWE::lllOutput(unsigned int n, double& p, double delta)
{
  double m = 2*n + 128;

  //execution log(time) = 1.8/ log(delta) − 110 and -80 to compute processor cycles. We add a margin of 20 so we take k/2 = 1.8/log(delta) - 100
  double lll1 = pow(delta, m) * pow(p, n/m);

  double lll2 = 2 * sqrt(n * log2(p) * log2(delta));
  lll2 = pow(2, lll2);

  return std::min(lll1, lll2);
}

double NFLLWE::estimateAbsTime(std::string crypto_param)
{
  using namespace std;
  vector<string> fields;
  boost::algorithm::split(fields, crypto_param, boost::algorithm::is_any_of(":"));
  unsigned int p_size = (unsigned) atoi(fields[3].c_str());
  double a = (p_size < 64) ? 1 : ceil(static_cast<double>(p_size)/64.0);
  unsigned int degree = (unsigned) atoi(fields[2].c_str());
  double b = degree/1024;

  return 1/(1.75 * pow(10, 5)/(a*b));
}

double NFLLWE::estimatePrecomputeTime(std::string crypto_param)
{
  using namespace std;
  vector<string> fields;
  boost::algorithm::split(fields, crypto_param, boost::algorithm::is_any_of(":"));
  unsigned int p_size = (unsigned) atoi(fields[3].c_str());
  double a = (p_size < 64) ? 1 : ceil(static_cast<double>(p_size)/64.0);
  unsigned int degree = (unsigned) atoi(fields[2].c_str());
  double b = degree/1024;

  return 1/(0.75*pow(10, 5)/(a*b));
}

unsigned int NFLLWE::getmodulusBitsize() {
	return nbModuli*kModulusBitsize;
}

// *********************************************************
// AbstractPublicParameters stuff
// *********************************************************
AbstractPublicParameters& NFLLWE::getPublicParameters()
{
	//This was bug chasing but should not be necessary!
	publicParams.setcrypto_container(this);
  	return publicParams;
}

std::string NFLLWE::getSerializedCryptoParams(bool shortversion)
{
  return publicParams.getSerializedParams(shortversion);
}


NFLLWE::~NFLLWE()
{
  clearSecretKeys();
}


std::string& NFLLWE::toString()
{
  return cryptoName;
}

void NFLLWE::clearSecretKeys()
{
  if(oldNbModuli)
  {
    // secreKey was allocated with a single allocation
    delete[] Abit_mod;
    delete[] Abit_mod_shoup;
    free(secretKey[0]);
    delete[] secretKey;
  }

  if(oldNbModuli)
  {
    for (unsigned int i = 0; i < oldNbModuli; i++) {
      free(secretKeyShoup[i]);
    }
    delete[] secretKeyShoup;
  }

  oldNbModuli = 0;
}


//This main is for benchmarking and tests
// 
// int main(int c,char **v) {
// 	
// // Benchs et correctness enc/dec
// 	 	NFLLWE n;	
// 		n.setNewParameters(1024,64,22);
//  		n.setmodulus(P64);
//  		n.getPublicParameters().computeNewParameters("lwe:80:1024:64:22");
//  	
//  		poly64 p=n.boundedRandomPoly(1024, 1023);
//  		poly64 result=(poly64)calloc(1024,sizeof(uint64_t));
//  		
//  	 	std::cout<<"0-RND polynom: ";n.print_poly64(p,4);std::cout<<std::endl;
//  		
//  		
//  		lwe_cipher cyph;  
//  #ifdef bench
//  		double start     = omp_get_wtime();
//  		for(int i        = 0;i<Repetition;i++) {
//  #endif 
//  			n.enc(&cyph,p);
//  #ifdef bench
//  		
//  			}
//  		double end        = omp_get_wtime();
//  		std::cout<<Repetition/(end-start)<<" chiffre/s"<<std::endl;
//  		
//  		
//  		 start     = omp_get_wtime();
//  		for(int i        = 0;i<Repetition;i++) {
//  #endif 
//  			
//  			n.dec(result,&cyph);
//  #ifdef bench
//  		}
//  		 end        = omp_get_wtime();
//  		std::cout<<Repetition/(end-start)<<" dechiffre/s"<<std::endl;
//  #endif 
// 		NFLLWE_DEBUG_MESSAGE("Encrypted into a",cyph.a,4);
// 		NFLLWE_DEBUG_MESSAGE("Encrypted into b",cyph.b,4);
//  	 	NFLLWE_DEBUG_MESSAGE("1-Encoded-Decoded (but not unshoupified): ",result,4);
//  		
// 	
//  		for(int i = 0;i<1024;i++) {
//  		  if((result[i]%P64)!=(result[i]%P64)) {
//  			std::cout<<"err "<<(p[i])<<" != "<<(result[i]%P64)<<std::endl;
//       	  	exit(1);
//  			break;
//  		  }
//  		}
//  	
// 	std::cout<< "enc/dec test passed"<<std::endl;
// 	
// 	
// 	
// 	// int bytesize=1024*22/8+1;
// 	//  	char* mydata=(char*)calloc(bytesize,1);
// 	//  	for(int i=0;i<bytesize;i++) {
// 	//  		mydata[i]='A'+(i%('Z'-'A'));
// 	//  	}
// 	//  	std::cout<<"Initial data : "<<mydata<<std::endl;
// 	// 
// 	//  	uint64_t bitsize=bytesize*8;
// 	//  	uint64_t nbOfPolys;
// 	//  	Warning, need to tranform this line with the new version of deserializeDataNTT which takes 4 parameters instead of 3 poly64 *mydata_poly=n.deserializeDataNTT((unsigned char*)mydata,bitsize,nbOfPolys);
// 	//  	
// 	//  	std::cout<<"The string has been encoded into ";n.print_poly64hex(*mydata_poly,1024*nbOfPolys);
// 	//   	std::cout<<std::endl<<"nbOfPolys = "<<nbOfPolys<<std::endl;
// 	//  	
// 	//  	// encrypt of poly64 simulation
// 	//  	lwe_cipher *cyphertext=new lwe_cipher[nbOfPolys];
// 	// for(int i=0;i<nbOfPolys;i++) { 
// 	//  		n.enc(&(cyphertext[i]),*(mydata_poly + i*1024));
// 	//  	//std::cout<<"The string has been cyphered into a["<<i<<"]";n.print_poly64(cyphertext[i].a,1024);
// 	//  	//std::cout<<"The string has been cyphered into b["<<i<<"]";n.print_poly64(cyphertext[i].a,1024);
// 	// }
// 	//  
// 	//  	poly64 unciphereddata_poly[nbOfPolys];// = (poly64*)calloc(1024*nbOfPolys,sizeof(uint64_t));
// 	// for(int i=0;i<nbOfPolys;i++) { 	
// 	// 	unciphereddata_poly[i]=(poly64)n.decrypt((char*)(cyphertext[i].a), (unsigned int)0,(size_t) 0,(size_t) 0);
// 	//  	std::cout<<"Decoded into polynom: ";n.print_poly64hex((poly64)unciphereddata_poly[i],128);std::cout<<std::endl;
// 	// }
// 	//  	
// 	//  	
// 	//  	unsigned char* unciphereddata=n.serializeData(unciphereddata_poly[0], nbOfPolys,bitsize, true);
// 	// 
// 	//  	std::cout<<"Decoded into "<<std::hex<<unciphereddata<<std::endl;
// 	//  	std::cout<<"A= "<<std::dec<<(short)'A'<<std::endl<<std::endl;
// 	//  	std::cout<<"G= "<<std::dec<<(short)'G'<<std::endl<<std::endl;
// 	//  
// 	//  
//  	//  
//  	// //lwe_cipher *cypherui;
//  	// //cypherui=(lwe_cipher *)n.encrypt(1,0);
//  	// //unciphereddata_poly = (poly64)n.decrypt((char*)(cypherui->a), (unsigned int)0,(size_t) 0,(size_t) 0);
//  	// char *charptr;
//  	// charptr=n.encrypt(1,0);
//  	// unciphereddata_poly = (poly64)n.decrypt(charptr, (unsigned int)0,(size_t) 0,(size_t) 0);
//  	// std::cout<<"1-Decoded into polynom: ";n.print_poly64(unciphereddata_poly,1024);std::cout<<std::endl;
//  	// 	std::cout<<"Decoded into "<<std::hex<<unciphereddata_poly<<std::endl;
//  }
