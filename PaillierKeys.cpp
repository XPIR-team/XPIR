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

#include "PaillierKeys.hpp"
#include <iostream>
#include <cstring>

/*************** PRVKEY**************/
paillier_prvkey::paillier_prvkey(){
  init_key();
}

paillier_prvkey::~paillier_prvkey(){
  clear_key();
}

void paillier_prvkey::init_key(){
  mpz_inits(d, inv_d, NULL);
}

void paillier_prvkey::clear_key()
{
  mpz_clears(d, inv_d, NULL);
}


/*************** PUBKEY**************/
paillier_pubkey::paillier_pubkey() :
  bits(0), 
  init_s(1)
{
  init_key();
}

paillier_pubkey::paillier_pubkey(unsigned int bits, char* rawKey) :
  bits(0), 
  init_s(1)
{
  init_key(bits, rawKey);
}

void paillier_pubkey::init_key() {
  for (int i = 0; i <= MAX_S; i++)
  {
    mpz_init_set_ui(nj[i],1);
  }
  mpz_init_set_ui(g,1);
}

void paillier_pubkey::init_key(unsigned int _bits, char* rawKey) {
  int init_s_;
  bits = _bits;

  init_key();
  
  mpz_import(nj[1], _bits / 8, 1, sizeof(char), 0, 0, rawKey);
  mpz_add_ui(g, nj[1], 1);
  memcpy(&init_s_, rawKey+_bits/8, sizeof(int));
  
  // The client should not be using s above MAX_S
  if (init_s_ >= MAX_S)
  {
    std::cout << "PaillierKeys: WARNING. The client tries to use s>=MAX_S. Setting s=MAX_S-1."<<std::endl;
    init_s = MAX_S-1;
  } 
  else init_s = init_s_;

  for (int i = 2; i <= init_s+1; i++)
  {
    mpz_pow_ui(nj[i], nj[1], i);
  }
}

//mocked key
void paillier_pubkey::init_key(unsigned int key_bit_size)
{
  init_s = 1;
  gmp_randstate_t rand;
  gmp_randinit_default(rand);

  mpz_urandomb(nj[1], rand, key_bit_size); 
  mpz_add_ui(g, nj[1], 1);
  
  for (int i = 2; i <= init_s+1; i++)
	{
		mpz_pow_ui(nj[i], nj[1], i);
	}
  gmp_randclear(rand);
}

inline void paillier_pubkey::init_nj(int i)
{
	mpz_pow_ui(nj[i], nj[1], i);
}


paillier_pubkey::~paillier_pubkey(){
  mpz_clear(g);

  for(int i = 0; i < MAX_S; i++)
    mpz_clear(nj[i]);
}


// Complete nj array up to index s_
void paillier_pubkey::complete_key(unsigned int s_){
  int s = s_;

  // The client should not be using moduli above MAX_S
  if (s > MAX_S)
  {
    std::cerr << "PaillierKeys: WARNING trying to complete keys above MAX_S bounding it to MAX_S" << std::endl;
    s = MAX_S;
  }

  // If g's value has not been initialized do it now
  if (mpz_get_ui(g) == 1 ) mpz_add_ui(g, nj[1], 1);
  
  // Initialize the array's values
  for (unsigned int i = 2; i <= s ; i++){
    // Should we save polar bears ? if (mpz_get_ui(nj[i]) == 1 ) 
    init_nj(i);
  }
}


// Provides the ciphertext modulus key i levels above the s defined in the class
// Ugly to return an mpz_t but the function purpose is to ensure that it has the correct
// value before its reference is returned
mpz_t* paillier_pubkey::getnj(int s_)
{
  int s = s_;

  // The client should not be using moduli above MAX_S
  if (s > MAX_S)
  {
    std::cerr << "PaillierKeys: WARNING trying to get key above MAX_S bounding it to MAX_S" << std::endl;
    s = MAX_S;
  }

  // If the key has been defined do it now
  if (mpz_get_ui(nj[s]) == 1 ) init_nj(s);
  
  return &nj[s];
}

// Simple getters
int paillier_pubkey::getinit_s()
{
  return init_s;
}

mpz_t* paillier_pubkey::getg()
{
  return &g;
}

int paillier_pubkey::getbits()
{
  return bits;
} 

// Simple setters
void paillier_pubkey::setinit_s(int init_s_)
{
  init_s = init_s_;
}

void paillier_pubkey::setbits(int bits_)
{
  bits = bits_;
} 

