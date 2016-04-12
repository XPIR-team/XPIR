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

#ifndef DEF_PAILLIER_STRUCT
#define DEF_PAILLIER_STRUCT
#include <cstddef>
#include <gmp.h>

// Maximum s such that encryption is donne from n^s to n^(s+1)
#define MAX_S 7 

class paillier_prvkey
{
public:
  paillier_prvkey();
  void init_key();
  void clear_key();
  mpz_t d; /* use CRT, d = mod n^s and d = O mod lambda*/
  mpz_t inv_d;
  ~paillier_prvkey();
};

class paillier_pubkey
{
public:
  paillier_pubkey();
  paillier_pubkey(unsigned int bits, char* rawKey);
  ~paillier_pubkey();
  void init_key();
  void init_key(unsigned int bits, char* rawKey);
  void init_key(unsigned int key_bit_size);
  // Complete nj array up to index s
  void complete_key(unsigned int s);
  // Get nj[s] = n^s initializing it if needed
  mpz_t* getnj(int s);
  // Simple getters
  mpz_t* getg();
  int getinit_s();
  int getbits();
  // Simple setters
  void setinit_s(int init_s_);
  void setbits(int bits_);

  void clear_key();
private:
  // Bit-size of the modulus
  int bits; 
  // nj[s] is n^s, nj[0] is therefore 1 (should not be used)
  mpz_t nj[MAX_S+1];  
  // Basic plaintext space is n^init_s and ciphertext space n^(init_s+1)
  unsigned int init_s;
  // Generator, i.e: n+1
  mpz_t g; 
  // Function initializing nj[i] for i>=2 when needed
  void init_nj(int i);
};

#endif
