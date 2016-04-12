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

#include "PaillierAdapter.hpp"
#include <set>
#include <iostream>
namespace {
const int kMagicConstant1 = 1;
const int kMillion = 1000000;
}

/**
 * Static fonction to feet approximately standards. NIST SP800-57 table 2 page 64
 * Beyond 128 bits of security modules will not be used as getCryptoParams blocks when ciphertext modulus is above 10K bits.
 **/
unsigned int PaillierAdapter::securityToModulus(int security_bits) 
{
  if(security_bits <= 40)
    return 128;
  if (security_bits <= 80)
    return 1024;
  if (security_bits <= 112)
    return 2048;
  if (security_bits <= 128)
    return 3072;
  if (security_bits <= 192)
    return 7680;
  // Unimplemented
  //return 15360;
  return 0;
}

PaillierAdapter::PaillierAdapter() :
  HomomorphicCrypto("Paillier")
{
  initRandomGenerator();
}

/**
 *	PaillierAdapter constructor
 * 	Params:
 * 		- int security_bits : number of security bit wanted (size of the cipher space);
 * 		- int recLvl				: number of recursion level.
 **/
PaillierAdapter::PaillierAdapter(int _security_bits, int init_s ) : 
  HomomorphicCrypto("Paillier"),
  security_bits(_security_bits),
  privateParameters(),
  modulusbits(securityToModulus(_security_bits)),
  publicParameters()
{
  initRandomGenerator();
  publicParameters.setModulusbits(modulusbits);

  keygen(
      modulusbits,
      publicParameters.getPubKey(),
      privateParameters.getPrvKey(),
      init_s );
}

void PaillierAdapter::initRandomGenerator()
{
  mpz_t s;
  mpz_init(s);

  gmp_randinit_default(rand);
  // Get a 1024 bit seed from /dev/urandom
  getRandFromfile( 128 , &s);
  gmp_randseed(rand, s);

  mpz_clear(s);

}

/**
 *	Encrypt an unsigned int and return the encrypted bytes in char*
 *	Params :
 *		- unsigned int ui : the value to encrypt ;
 *		- unsigned int s  : the exponent (recursion lvl).
 *
 *	Return :
 *		- char* : an allocated pointer to the encrypted data.
 **/
  char* 
PaillierAdapter::encrypt(unsigned int ui, unsigned int s) 
{
  unsigned int ciph_size = publicParameters.getKeyBitsize()*(s+1)/8;
  char *tmp, *request = (char*) calloc((ciph_size + 1), sizeof(char));
  size_t n;
  mpz_t pt, ct;

  mpz_init(ct);
  mpz_init_set_ui( pt, ui); // convert the 0 or the 1 to an mpz_t

  enc( this->publicParameters.getPubKey(), pt, s, ct );
#ifdef CRYPTO_DEBUG
  gmp_printf("Created query element: %Zd with args %Zd and %d\n",ct, pt, s);
#endif
  tmp = (char*)mpz_export(NULL, &n, 1, sizeof(char), 0, 0, ct);
  //Padding
  memcpy(request+sizeof(char)*((ciph_size) - n), tmp, n);	

  //Free memory
  mpz_clears(ct, pt, NULL);

  free(tmp);

  return request;
}

// To test decryption performance
char* PaillierAdapter::encrypt_perftest()
{
	return encrypt(1,1);// Decryption performance test work well with this small plaintext
}

/**
 *	Encrypt an char* of size dataSize and return char*
 *	Params :
 *		- char* data 				: data to encrypt ;
 *		- size_t dataSize 	: size of the data ;
 *		- unsigned exponent : the exponent (recursion lvl).
 *	Return :
 *		- char* : an allocated pointer to the encrypted data.
 **/
char* PaillierAdapter::encrypt(char* data, size_t dataSize,  unsigned int exponent) 
{
  char* request;
  mpz_t ct, imported_data;

  mpz_inits(ct, imported_data, NULL);

  mpz_import(imported_data, dataSize, 1, sizeof(char), 0, 0, data);

  enc( this->publicParameters.getPubKey(),
      imported_data,
      exponent,
      ct );

  request = (char*)mpz_export(NULL, NULL, 1, sizeof(char) , 0, 0, ct);

  mpz_clears(ct, imported_data, NULL);

  return request;
}

/**
 *	Decrypt reply sended by the server and return clear bytes in char*.
 *	Params :
 * 		- char* cipheredData 	 : encrypted data ;
 * 		- unsigned int exp_fac : the exponent (recursion lvl) ;
 * 		- size_t ciph_size		 : size of encrypted data ;
 * 		- size_t clear_size		 : size of decrypted data to return.
 * 		Return :
 * 		- char* : an allocated pointer to the decrypted data.
 **/
char* PaillierAdapter::decrypt( char* cipheredData, 
    unsigned int exp_fac, 
    size_t ciph_size , 
    size_t clear_size) 
{
  size_t n = 0; 
  mpz_t ct, res;
  char* tmp;
  mpz_inits(ct, res, NULL);

  mpz_import(ct, ciph_size, 1, sizeof(char), 0, 0, cipheredData);

  dec(publicParameters.getPubKey(),
      privateParameters.getPrvKey(),
      ct, 
      exp_fac,
      res );

#ifdef DEBUG_ENCRYPT
  gmp_printf("PaillierAdapter : Decrypting %Zd into %Zd\n\n", ct, res);
#endif

  tmp = (char*) mpz_export(NULL, &n, 1, sizeof(char) , 0, 0, res);
  
  mpz_clears(res, ct, NULL); //clear content of mpz_t, t

  if( n < clear_size)
  {
    char *pt = (char*)calloc(clear_size+1, sizeof(char)); 
    memcpy(pt+sizeof(char)*((clear_size) - n), tmp, n);
    free(tmp); // malloc'd
    return pt;
  }
  else
  {
    return tmp;
  }
}


PaillierAdapter::~PaillierAdapter() 
{
  gmp_randclear(rand);
}

AbstractPublicParameters& PaillierAdapter::getPublicParameters() 
{
  return publicParameters;
}

/**
 *	Private encrypt methode.
 *	Params :
 *		- paillier_pubkey* pub : Paillier public key pointer ;
 *		- mpz_t m	: data to encrypt ;
 *		- unsigned int s : recursion level ;
 *		- mpz_t c				 : result operande .
 **/
  void 
PaillierAdapter::enc( paillier_pubkey* pub,
    mpz_t m,
    unsigned int s,
    mpz_t c )
{
  mpz_t gm, rns, r;
  mpz_inits( gm, rns, r, NULL);
  unsigned int modulusbits = publicParameters.getKeyBitsize();
  do
  {
    mpz_urandomb(r, rand, (s+1)*modulusbits); 
  } while( mpz_cmp(r, *pub->getnj(s+1)) >= 0 );
  
  mpz_powm(gm, *pub->getg(), m,  *pub->getnj(s+1));
  mpz_powm(rns, r, *pub->getnj(s), *pub->getnj(s+1));
  mpz_mul(c, gm, rns);
  mpz_mod(c,c, *pub->getnj(s+1));

#ifdef DEBUG_ENCRYPT
  mpz_t test;
  mpz_init(test);
  dec(pub, privateParameters.getPrvKey(), c, s, test); 
  gmp_printf("Decrypting %Zd into %Zd\n\n", c, test);
#endif

  mpz_clears(gm, rns, r, NULL);
}

/**
 *	Reimplemented version of decrypt to not use paillier structures. 
 *	Params:
 *		- paillier_prvkey* prv : Paillier private key pointer ;
 *		- mpz_t c : encrypted data ;
 *		- unsigned int s : recursion level ;
 *		- mpz_t i : result operande .
 *  This function is an exact copy of the one in libdj
 **/
/*
	libdj - A library implementing the Damgård-Jurik cryptosystem.

	Copyright 2012 Fred Douglas (fed2@illinois.edu)
	This library is an extension/modification of libpaillier.
	libpaillier is copyright 2006 SRI International (written by John Bethencourt).

	This file is part of libdj.

    libdj is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    libdj is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with libdj.  If not, see <http://www.gnu.org/licenses/>.
*/


  void 
PaillierAdapter::dec( paillier_pubkey* pub,
    paillier_prvkey* prv,
    mpz_t c, 
    unsigned int s,
    mpz_t res)
{
	int i,j,k;
	mpz_t cprime, temp, tempDivisor, i_jCur, i_jPrev;

	mpz_init(cprime);
	
	//cprime = c^key mod n^{s+1}
	mpz_powm(cprime,c,prv->d,*pub->getnj(s+1));
	
	
	
	//the algorithm from extracting i from (1+n)^i described in the paper.
	
	//the algorithm iteratively extracts i mod n, i mod n^2, etc. i mod n starts us off, and is just input - 1 / n.
	mpz_init(temp);
	mpz_init(tempDivisor);
	mpz_init(i_jCur);
	mpz_init_set(i_jPrev,cprime);
	mpz_sub_ui(i_jPrev,i_jPrev,1);
	mpz_divexact(i_jPrev,i_jPrev,*pub->getnj(1));
	mpz_mod(i_jPrev,i_jPrev,*pub->getnj(2));
	
	//just in case s=1; in that case we need this line to actually set i_jcur
	if(s==1)
		mpz_set(i_jCur,i_jPrev);
	
	//extract i_j = i mod n^j given i_{j-1}. this is done by taking (input - 1 mod n^{j+1}) / n , and subtracting
	//from that: ((  (i_{j-1} choose 2)*n  + (i_{j-1} choose 3)*n^2  + ... +  (i_{j-1} choose j)*n^{j-1}  ))   mod n^j
	for(j=2;j<=s;j++)
	{
		//L((1+n)^i) as they call it 
		mpz_mod(i_jCur,cprime,*pub->getnj(j+1));
		mpz_sub_ui(i_jCur,i_jCur,1);
		mpz_divexact(i_jCur,i_jCur,*pub->getnj(1));
		
		mpz_mod(i_jCur,i_jCur,*pub->getnj(j));
		
		//subtract each of the binomial things
		for(k=2;k<=j;k++)
		{
			mpz_bin_ui(temp,i_jPrev,k);
			mpz_mul(temp,temp,*pub->getnj(k-1));
			mpz_mod(temp,temp,*pub->getnj(j));
			mpz_sub(i_jCur,i_jCur,temp);
			mpz_mod(i_jCur,i_jCur,*pub->getnj(j));
		}
		mpz_set(i_jPrev,i_jCur);
	}
	
	//i_jCur is currently the message times the private key.
	mpz_invert(temp, prv->d, *pub->getnj(s));
	mpz_mul(res, i_jCur, temp);
	mpz_mod(res, res, *pub->getnj(s));
	
	//cleanup and return
	mpz_clear(cprime);
	mpz_clear(i_jPrev);
	mpz_clear(i_jCur);
	mpz_clear(temp);
	mpz_clear(tempDivisor);
	

  //  unsigned int kfac, k;
//  mpz_t t1, t2, t3, a, id ;
//  std::cout << "dec called with s="<<s<<std::endl;
//  mpz_inits(t1, t2, t3, a, NULL);
//  mpz_init_set_ui(id,  0);
//
//  mpz_powm(a, c, prv->d, *pub->getnj(s+1));
//
//  for(unsigned int j = 1 ; j <= s; j++)
//  {
//    mpz_set(t1, a);
//    mpz_mod(t1, t1, *pub->getnj(j+1));
//    mpz_sub_ui(t1, t1, 1);
//    mpz_div(t1, t1, *pub->getnj(1));
//    mpz_set(t2, id);
//    kfac = 1;
//
//    for (k = 2; k <= j; k++)
//    {
//      kfac *= k;
//      mpz_sub_ui(id, id, 1);
//      mpz_mul(t2, t2, id);
//      mpz_mod(t2, t2, *pub->getnj(j));
//      mpz_set_ui(t3, kfac);
//
//      mpz_invert(t3, t3, *pub->getnj(j));
//
//      mpz_mul(t3, t3, t2);
//      mpz_mod(t3, t3, *pub->getnj(j));
//      mpz_mul(t3, t3, *pub->getnj(k-1));
//      mpz_mod(t3, t3, *pub->getnj(j));
//      mpz_sub(t1, t1, t3);
//      mpz_mod(t1, t1, *pub->getnj(j));
//    }
//    mpz_set(id, t1);
//  }
//
//  mpz_mul(rop, id, prv->inv_d);
//
//  mpz_mod(rop, rop, *pub->getnj(s));
//
//  mpz_clears(t1, t2, t3, a, id, NULL);
}


  void
PaillierAdapter::keygen( unsigned int modulusbits,
    paillier_pubkey* pub,
    paillier_prvkey* prv,
    unsigned int init_s)
{
  mpz_t p, q;

  mpz_inits(p, q, NULL);

  /* pick random (modulusbits/2)-bit primes p and q */
  do
  {
    get_prime_of_size(p, modulusbits / 2);
    get_prime_of_size(q, modulusbits / 2);

    mpz_mul(*pub->getnj(1), p, q);

  }while(mpz_sizeinbase(*pub->getnj(1), 2) != modulusbits);

  pub->setbits(modulusbits);
  pub->setinit_s(init_s);
  pub->complete_key(init_s+1); /*we don't know if it will be used beyond one level*/

  /* compute the private key lambda = lcm(p-1,q-1) */
  mpz_sub_ui(p, p, 1);
  mpz_sub_ui(q, q, 1);
  mpz_lcm(prv->d, p, q);
  mpz_invert(prv->inv_d, prv->d ,*pub->getnj(init_s));

  /* clear temporary integers */
  mpz_clears(p, q, NULL);
}

/**
 *	Get prime number of given size.
 *	Parms :
 *		- mpz_t rop  : result operand ;
 *		- unsigned int length : length in bit.
 **/
  void
PaillierAdapter::get_prime_of_size(mpz_t rop, unsigned int length)
{

  mpz_t z_p, z_min;
  mpz_init(z_p);
  mpz_init_set_ui(z_min, 1);

  mpz_mul_2exp(z_min, z_min, length - 1);
  do
  {
    mpz_urandomb(z_p, rand, length);
    mpz_add(z_p, z_min, z_p);
    mpz_nextprime(z_p, z_p);
  }
  while(mpz_sizeinbase(z_p, 2) != length);

  mpz_set(rop, z_p);
  mpz_clears(z_p, z_min, NULL);
}

  void 
PaillierAdapter::getRandFromfile(int len, mpz_t* val )
{
  char* p = new char[len + 1]();
  std::ifstream f("/dev/urandom", ios::in | ios::binary);

  for(int i = 0; i < len; i++)
  {
    f.read(p+i, 1);
  }
  mpz_import(*val, len, 1, sizeof(char), 0, 0, p);

  f.close();

  delete[] p;
}

  void
PaillierAdapter::setNewParameters(const std::string& crypto_params)
{
  std::vector<std::string> fields;
  boost::algorithm::split(fields, crypto_params, boost::algorithm::is_any_of(":"));
  unsigned int keySize = (unsigned)atoi(fields[2].c_str());
  security_bits = atoi(fields[1].c_str());
  publicParameters.setModulusbits(keySize);
  publicParameters.setSecurityBits(security_bits);
  // TODO(performance) : Takes too much time when generating the encrypt and decrypt caches
  // An option should be included so that it is a fixed key from a file when called to build caches
  keygen(
      keySize,
      publicParameters.getPubKey(),
      privateParameters.getPrvKey(),
      kMagicConstant1);
}
  double
PaillierAdapter::getDecCost(unsigned int length, unsigned int s)
{
  double  start, result;
  unsigned int cores = omp_get_num_procs(); 
  unsigned int rounds = 2000/pow((length*(s+1)/2048), 3);//Amount of rounds for a one second test on the poor devlopper's computer.
  mpz_t z_r, z_a[cores];
  mpz_init(z_r);

#pragma omp parallel for
  for (unsigned int i = 0 ; i < cores ; i++)
    mpz_init(z_a[i]);

  gmp_randstate_t prng;
  gmp_randinit_default(prng);
  gmp_randseed_ui(prng, time(NULL));

  getRandInteger(z_r, length, prng);

  start = omp_get_wtime();

  for(unsigned int round = 0 ; round < rounds ; round++)
  {
#pragma omp parallel for
    for (unsigned int i = 0 ; i < cores ; i++)
    {
      dec(publicParameters.getPubKey(),
          privateParameters.getPrvKey(),
          z_r, 
          s ,
          z_a[i] );
      usleep(5);
    }
  }
  mpz_clear(z_r);

#pragma omp parallel for
  for (unsigned int i = 0 ; i < cores ; i++)
    mpz_clear(z_a[i]);

  result = omp_get_wtime() - start;

  return result /(cores)/rounds ; 
}

  void 
PaillierAdapter::getRandInteger(mpz_t rop, unsigned int length, gmp_randstate_t prng)
{
  do {
    mpz_urandomb(rop, prng, length);	
  }while(mpz_sizeinbase(rop, 2) != length);
}

unsigned int PaillierAdapter::getAllCryptoParams(std::set<std::string>& crypto_params)
{
  unsigned int params_nbr = 0;
  unsigned int k_array_size = 4;
  // k=256 gives a 15K bit plaintext modulus and a 30K bit ciphertext modulus
  // as this is huge we turn off this proposal by default
  unsigned int k[4] = {80, 100, 128, 192};

  for (unsigned int i = 0 ; i < k_array_size ; i++)
  {
    params_nbr += getCryptoParams(k[i], crypto_params);
  }

  return params_nbr;
}

unsigned int PaillierAdapter::getCryptoParams(unsigned int k, std::set<std::string>& crypto_params)
{
  unsigned int params_nbr = 0;
  string k_str = std::to_string(k);

  unsigned int modulus_size = securityToModulus(k);
  if (modulus_size == 0) return 0;
  // We'll increase params_nbr when the client handles better s
  for (unsigned int s = 1 ; params_nbr < 1 ; s++)
  {
    string param = cryptoName + ":" + k_str + ":" + to_string(modulus_size*s) + ":" + to_string(modulus_size*(s+1)); 
    if(crypto_params.insert(param).second) params_nbr++;
    param = "";
  }
  return params_nbr;
}


long PaillierAdapter::setandgetAbsBitPerCiphertext(unsigned int elt_nbr)
{
  return publicParameters.getAbsorptionBitsize();
}

std::string PaillierAdapter::getSerializedCryptoParams(bool shortversion)
{
  return publicParameters.getSerializedParams(shortversion);
}

double PaillierAdapter::estimateAbsTime(std::string crypto_param)
{
  vector<string> fields;
  boost::algorithm::split(fields, crypto_param, boost::algorithm::is_any_of(":"));
  const double mod_size_s = static_cast<double>(atoi(fields[3].c_str()));
  return 1.0 / (kMillion/(3.0*1016) * pow(2048 / mod_size_s, 3));
}

/**
 *	Compute a modular exponentiation in crypted space who is a operation in clear space.
 *	Params :
 *		- mpz_t res : operation result;
 *		- mpz_t a   : first operand, crypted data ;
 *		- mpz_t b   : second operans, crypted data;
 *		- int		s		: modulus index  .
 **/
void PaillierAdapter::e_add(mpz_t res, mpz_t a, mpz_t b, int s)
{
	mpz_mul(res, a, b);
	mpz_mod(res, res, *publicParameters.getPubKey()->getnj(s));
}

/**
 *	Compute modular multiplication in crypted space who is a multiplicaton in clear space.
 *	Params :
 *		- mpz_t res : operation result;
 *		- mpz_t a   : first operand, crypted data ;
 *		- mpz_t n   : second operand, a constant  ;
 *		- int   s   : modulus index .
 **/
void PaillierAdapter::e_mul_const(mpz_t res, mpz_t a, mpz_t n, int s)
{
	mpz_powm(res, a, n , *publicParameters.getPubKey()->getnj(s));
}
