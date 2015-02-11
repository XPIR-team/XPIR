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
 *
 * This file incorporates work covered by the following copyright and  
 * permission notice:  
 * 
 *   Demonstration code for the paper "Faster arithmetic for number-theoretic
 *   transforms", May 2012
 *
 *   Copyright (C) 2014, David Harvey
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions are met:
 *
 *   * Redistributions of source code must retain the above copyright notice, this
 *     list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef DEF_NFLlib
#define DEF_NFLlib

#define SHOUP

#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <cstddef>
#include <gmp.h>
#include <iostream>
#include "NFLParams.hpp"
#include "prng/fastrandombytes.h"
#include <cstring>
#include <fstream>
#include <cassert>
// SSE/AVX(2) instructions
#include <immintrin.h>
#include <x86intrin.h>

	
using namespace std;

// Polynomials are pointers to arrays of uint64_t coefficients
typedef uint64_t* poly64;

template <size_t Align, class T>
static inline T* malloc_align(const size_t n)
{
       T* ret;
       if(posix_memalign((void**) &ret, Align, n*sizeof(T)))
         std::cout << "NFLlib: WARNING posix_memalign failed" << std::endl;
       return ret;
}

template <size_t Align, class T>
static inline T* malloc_align(const size_t n, T const* const)
{
       return malloc_align<Align, T>(n);
}

template <size_t Align, class T>
static inline T* calloc_align(const size_t n)
{
       T* ret = malloc_align<Align, T>(n);
       memset(ret, 0, n*sizeof(T));
       return ret;
}

template <size_t Align, class T>
static inline T* calloc_align(const size_t n, T const* const)
{
       T* ret = malloc_align<Align, T>(n);
       memset(ret, 0, n*sizeof(T));
       return ret;
}

template <size_t Align, class T, class V>
static inline T* calloc_align_other(const size_t n)
{
       T* ret = malloc_align<Align, T>(n);
       memset(ret, 0, n*sizeof(V));
       return ret;
}

// Number-Theoretic Transform Tools class
class NFLlib  
{

  public:
    // Constructors and initializing functions
    NFLlib();
    void configureNTT(); // Only called from setmodulus and setpolyDegree

    // Fundamental setters (result in a call to initializing functions above) 
    void setmodulus(uint64_t modulus);
    void setpolyDegree(unsigned int polyDegree);
    void setNewParameters(unsigned int polyDegree, unsigned int modulusBitsize);
 
    // Getters
    uint64_t* getmoduli();
    unsigned int getpolyDegree();
    unsigned short getnbModuli();
    void copymoduliProduct(mpz_t dest);
   
    // Random polynomial generation functions
    void setBoundedRandomPoly(poly64 res, uint64_t upperBound_, bool uniform); 
    poly64 allocBoundedRandomPoly(uint64_t upperBound_, bool uniform);

    // Data import and export main functions
    poly64* deserializeDataNFL(unsigned char **inArrayOfBuffers, uint64_t nbrOfBuffers, 
        uint64_t dataBitsizePerBuffer, unsigned bitsPerCoordinate, uint64_t &polyNumber);
    static void serializeData64 (uint64_t* indata, unsigned char* outdata, 
        unsigned int bitsPerChunk, uint64_t nb_of_uint64);
    static void serializeData32 (uint32_t* indata, unsigned char* outdata, unsigned int bitsPerChunk, 
        uint64_t nb_of_uint32);

    // Helper functions
    poly64 allocpoly(bool nullpoly);
    mpz_t* poly2mpz (poly64 p);
    static void print_poly64(poly64 p, unsigned int coeff_nbr);
    static void print_poly64hex(poly64 p, unsigned int coeff_nbr);

    // Destructors memory freeing routines
    ~NFLlib();
    void freeNTTMemory(); // Only called from configureNTT and ~NFLlib
    
    
    // ****************************************************
    // Modular and Polynomial computation inlined functions
    // ****************************************************

    ///////////////////// Operations over uint64_t (polynomial coefficients)
    // Additions
    static uint64_t addmod(uint64_t x, uint64_t y, uint64_t p);
    static uint64_t submod(uint64_t x, uint64_t y, uint64_t p);
    // Multiplications
    static uint64_t mulmod(uint64_t x, uint64_t y, uint64_t p);
    static uint64_t mulmodShoup(uint64_t x, uint64_t y, uint64_t yprime, uint64_t p);
    // Fused Multiplications-Additions
    static void mulandadd(uint64_t &rop, uint64_t x, uint64_t y, uint64_t p);
    static void mulandaddShoup(uint64_t &rop, uint64_t x, uint64_t y, uint64_t yprime, uint64_t p);

    ///////////////////// Operations over polynomials
    // Additions
    void addmodPoly(poly64 rop, poly64 op1, poly64 op2);
    void submodPoly(poly64 rop, poly64 op1, poly64 op2);
    // Multiplications 
    void mulmodPolyNTT(poly64 rop, poly64 op1, poly64 op2);
    void mulmodPolyNTTShoup(poly64 rop, poly64 op1, poly64 op2, poly64 op2prime);
    // Fused Multiplications-Additions
    void mulandaddPolyNTT(poly64 rop, poly64 op1, poly64 op2);
    void mulandaddPolyNTTShoup(poly64 rop, poly64 op1, poly64 op2, poly64 op2prime);
    
    // Number-Theoretic Transorm functions
    static void ntt(uint64_t*  x, const uint64_t*  wtab, const uint64_t*  winvtab, 
        const unsigned K, const uint64_t p);
    static void inv_ntt(uint64_t* const x, const uint64_t* const inv_wtab, const uint64_t* const inv_winvtab,
        const unsigned K, const uint64_t invK, const uint64_t p, const uint64_t* const inv_indexes);
    static void prep_wtab(uint64_t* wtab, uint64_t* winvtab, uint64_t w, unsigned K, uint64_t p);
    void nttAndPowPhi(poly64 op);
    void invnttAndPowInvPhi(poly64);
    
    // Pre-computation functions
    poly64 allocandcomputeShouppoly(poly64 x);

  private: 
    // Says whether the object has been initialized for the int amount of moduli it is set to
    unsigned int alreadyInit;
    // Polynomial attributes
    uint64_t *moduli;
    unsigned short nbModuli;
    unsigned int polyDegree;
    // Chinese Remainder Theorem (CRT) and Inverse CRT related attributes
    mpz_t *liftingIntegers;
    mpz_t moduliProduct;
    
    // NTT and inversse NTT related attributes
    // omega**i values (omega = polydegree-th primitive root of unity), and their inverses
    // phi**i values (phi = square root of omega), and their inverses multiplied by a constant
    // Shoup variables contain redundant data to speed up the process
    // NOTE : omega and derived values are ordered following David Harvey's algorithm
    // w**0 w**1 ... w**(polyDegree/2) w**0 w**2 ... w**(polyDegree/2) w**0 w**4 ...
    // w**(polyDegree/2) etc. (polyDegree values)
    uint64_t **phis, **shoupphis, **invpoly_times_invphis, 
             **shoupinvpoly_times_invphis, **omegas, 
             **shoupomegas, **invomegas, **shoupinvomegas, *invpolyDegree;
    uint64_t* inv_indexes;

    // WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    // ENTERING THE DEN... Internal functions you should never read and we will never comment
    uint64_t* bitsplitter (unsigned char** inDataBuffers, uint64_t nbrOfBuffers, 
        uint64_t bitsPerBuffer, unsigned int bitsPerChunk);
    void internalLongIntegersToCRT(uint64_t* tmpdata, poly64 outdata, uint64_t int_uint64PerChunk,
        uint64_t totalNbChunks) ;
    void bs_finish(poly64 &outdata, uint64_t int_uint64PerChunk, uint64_t polyNumber, 
        uint64_t *splitData, uint64_t nbrOfBuffers, uint64_t bitsPerBuffer, unsigned int bitsPerChunk);
    void bs_loop (unsigned char** inDataBuffers, uint64_t nbrOfBuffers, 
        uint64_t bitsPerBuffer, unsigned int bitsPerChunk, uint64_t *&tmpdata, 
        uint64_t bufferIndex, uint64_t &bitsread, size_t &subchunkIndex);
    uint64_t* bitsplitter_backtoback_internal_test (unsigned char** inDataBuffers, 
        uint64_t nbrOfBuffers, uint64_t bitsPerBuffer, unsigned int bitsPerChunk, 
        uint64_t totalbitsread,uint64_t bitsread, uint64_t* pointer64, unsigned int bitsToRead);
    // END OF THE DEN
};


// We hate C++ inlining
// These function sources were included here to force inlining
//
/* WARNING: For performance, our modular functions often do not do a complete modular reduction.
 * For example when adding two integers x, y modulo an integer p we do not check if they are
 * larger than p or not and only subtract to the result at most p. This is noted x + y lazymod p
 * in our comments. Applications using our functions should ensure that no overflow is possible given
 * our implementation. This is dangerous and a bad programming habit but mandatory for high
 * performance.
 * */


// *********************************************************
//     Operations over uint64_t (polynomial coefficients)
// *********************************************************

// Modular addition plain and simple : add then test if sum is superior 
// to the modulus in which case substract the modulus to the sum
// ASSUMPTION: x + y < 2**64
// OUTPUT: x + y lazymod p
inline uint64_t NFLlib::addmod(uint64_t x, uint64_t y, uint64_t p)
{
    uint64_t z = x + y;
    return z -= ((z >= p) ? p : 0);
}


// Modular subtract trick: if y is smaller than p return x+p-y else return x+p-(y%p)
// OUTPUT: x - y lazymod p
inline uint64_t NFLlib::submod(uint64_t x, uint64_t y, uint64_t p)
{
    return (y < p) ? addmod(x, p-y, p) : addmod(x, p - (y % p), p);
}


// Modular multiplication: trivial (and very costly) approach with complete reduction
// OUTPUT: x * y mod p
inline uint64_t NFLlib::mulmod(uint64_t x, uint64_t y, uint64_t p)
{
  return (uint128_t) x * y % p;
}


// Modular multiplication: much faster alternative
// Works only if yprime = ((uint128_t) y << 64) / p has been pre-computed
// ASSUMPTION: p is at most 62 bits
// OUTPUT: x * y lazymod p (if x and y are below p the result is x * y mod p)
inline uint64_t NFLlib::mulmodShoup(uint64_t x, uint64_t y, uint64_t yprime, uint64_t p)
{
  uint128_t res;
  uint64_t q = ((uint128_t) x * yprime) >> 64;
  res = x * y - q * p;
  res -= ((res>=p) ? p : 0);
  return res;  
}


// Fused Multiply and Addition (FMA) : trivial approach
inline void NFLlib::mulandadd(uint64_t& rop, uint64_t x, uint64_t y, uint64_t p)
{
  rop = (((uint128_t) (x) * y) + rop) % p; 
}


// Fused Multiply and Addition (FMA) : Shoup's precomputed approach again
// Works only if yprime = ((uint128_t) y << 64) / p has been pre-computed
// ASSUMPTION: p is at most 62 bits
// OUTPUT: x * y lazymod p (if x and y are below p the result is below 2p)
inline void NFLlib::mulandaddShoup(uint64_t& rop, const uint64_t x, const uint64_t y, const uint64_t yprime, const uint64_t p)
{
  //mulandadd(rop, x, y, p);
#ifdef DEBUG  
  // This must be before the real computation modifies rop
  uint128_t res = ((uint128_t) x * y + rop) % p;
#endif  
  uint64_t q = ((uint128_t) x * yprime) >> 64;
  rop += x * y - q * p;
  rop -= ((rop>=p) ? p : 0);
#ifdef DEBUG  
  // And this must be after the modification of rop
  if((res%p)!=(rop%p))
  {
    std::cout<<"NFLlib: CRITICAL Shoup multiplication failed (prob. orig. precomputation or bounds)"<<std::endl;
    exit(1);
  }
#endif      
}




// *********************************************************
//     Operations over polynomials
// *********************************************************

// Apply addmod to all the coefficients of a polynomial
inline void NFLlib::addmodPoly(poly64 rop, poly64 op1,poly64 op2) {
  for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) 
  {
    for(unsigned i=0;i<polyDegree;i++)
    {
      rop[i]=addmod(op1[i],op2[i],moduli[currentModulus]);
    }
    rop+=polyDegree;
    op1+=polyDegree;
    op2+=polyDegree;
  }
}


// Apply submod to all the coefficients of a polynomial
inline void NFLlib::submodPoly(poly64 rop, poly64 op1,poly64 op2) {
  for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) 
  {
    for(unsigned i=0;i<polyDegree;i++)
    {
      rop[i]=submod(op1[i],op2[i],moduli[currentModulus]);
    }
    rop+=polyDegree;
    op1+=polyDegree;
    op2+=polyDegree;
  }
}


// Apply mulmod to all the coefficients of a polynomial
// This is a polynomial multiplication mod X**n+1 iff the operands
// have been processed through nttAndPowPhi
inline void NFLlib::mulmodPolyNTT(poly64 rop, poly64 op1, poly64 op2)
{
  for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) {  
    for (unsigned i = 0 ; i < polyDegree ;i++)
    {
      rop[i] = mulmod(op1[i],op2[i],moduli[currentModulus]);
    }
    rop+=polyDegree;
    op1+=polyDegree;
    op2+=polyDegree;
  }
}


// Apply mulmodShoup to all the coefficients of a polynomial (much faster)
// This is a polynomial multiplication mod X**n+1 iff the operands
// have been processed through nttAndPowPhi
// op2prime must be a polynomial with op2 coefficients converted with Shoup's precomputation
inline void NFLlib::mulmodPolyNTTShoup(poly64 rop, poly64 op1,poly64 op2,poly64 op2prime)
{
  for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) {
    //#pragma omp parallel for
    for (unsigned i = 0 ; i < polyDegree ;i++) {
      rop[i] = mulmodShoup(op1[i],op2[i],op2prime[i],moduli[currentModulus]);
    }
    rop+=polyDegree;
    op1+=polyDegree;
    op2+=polyDegree;
    op2prime+=polyDegree;
  }
}


// Same as mulmodPolyNTT but with fused multiplication and addition
inline void NFLlib::mulandaddPolyNTT(poly64 rop, poly64 op1,poly64 op2)
{
  for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) {  
    for (unsigned i = 0 ; i < polyDegree ;i++)
    {
      mulandadd(rop[i],op1[i],op2[i],moduli[currentModulus]);
    }
    rop+=polyDegree;
    op1+=polyDegree;
    op2+=polyDegree;
  }
}


// Same as mulmodPolyNTTShoup but with fused multiplication and addition
inline void NFLlib::mulandaddPolyNTTShoup(poly64 rop, poly64 op1,poly64 op2,poly64 op2prime)
{
  for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) {
  //#pragma omp parallel for
    for (unsigned i = 0 ; i < polyDegree ;i++) {
      mulandaddShoup(rop[i],op1[i],op2[i],op2prime[i],moduli[currentModulus]);
    }
    rop+=polyDegree;
    op1+=polyDegree;
    op2+=polyDegree;
    op2prime+=polyDegree;
  }
}




// *********************************************************
//     Number-Theoretic Functions
// *********************************************************

// *****************************************************************
// NTT functions from David Harvey 
// From : http://web.maths.unsw.edu.au/~davidharvey/papers/fastntt/
// Most of the functions have been modified for our needs
// The associated copyright notice is at the beginning of this file
// *****************************************************************


// Number-Theoretic Transform: replaces a coefficient representation of a polynomial
// by the values of the polynomial on a set of points. Allows to do coefficient wise 
// multiplication of polynomials instead of the usual convolution product
// - x points to the polynomial to which we will apply the NTT
// - wtab is a pre-computed table with the powers of a root of unity
// - winvtab is another pre-computed table with the powers of the inverse of a root of unity
// - K is the log_2 of the degree of the polynomial
// - p is the modulus coefficient operations are done with
// The NTT is computed in-place on x, so x is the output
inline void NFLlib::ntt(uint64_t* x, const uint64_t* wtab, const uint64_t* winvtab,
    const unsigned K, const uint64_t p)
{
  if (K == 1)
    return;

  // special case
  if (K == 2)
  {
    uint64_t u0 = x[0];
    uint64_t u1 = x[1];
    uint64_t t0 = u0 + u1;
    uint64_t t1 = u0 - u1;
    t0 -= (t0 >= 2*p) ? (2*p) : 0;
    t1 += ((int64_t) t1 < 0) ? (2*p) : 0;
    x[0] = t0;
    x[1] = t1;
    return;
  }

  size_t N0 = K;           // size of block
  //size_t M0 = 1;                         // number of blocks

  // log2(N)-2
  const int J = _bit_scan_reverse(K)-2;
  for (int w = 0; w < J; w++) {
	const size_t M = 1 << w; 
	const size_t N = N0 >> w;
//#pragma omp parallel for schedule(static) num_threads(4)
	for (size_t r = 0; r < M; r++) {
		for (size_t i = 0; i < N/2; i++) {
			uint64_t u0 = x[N * r + i];
			uint64_t u1 = x[N * r + i + N/2];

			uint64_t t0 = u0 + u1;
			t0 -= ((t0 >= 2*p) ? (2*p) : 0);

			uint64_t t1 = u0 - u1 + 2*p;

			uint64_t q = ((uint128_t) t1 * winvtab[i]) >> 64;
			uint64_t t2 = t1 * wtab[i] - q * p;

			x[N * r + i] = t0;
			x[N * r + i + N/2] = t2;
		}
	}
    wtab += N/2;
    winvtab += N/2;
  }

  const size_t M = 1 << J;
  // last two layers
  for (size_t r = 0; r < M; r++, x += 4)
  {
    uint64_t u0 = x[0];
    uint64_t u1 = x[1];
    uint64_t u2 = x[2];
    uint64_t u3 = x[3];

    uint64_t v0 = u0 + u2;
    v0 -= (v0 >= 2*p) ? (2*p) : 0;
    uint64_t v2 = u0 - u2;
    v2 += ((int64_t) v2 < 0) ? (2*p) : 0;

    uint64_t v1 = u1 + u3;
    v1 -= (v1 >= 2*p) ? (2*p) : 0;
    uint64_t t = u1 - u3 + 2*p;

    uint64_t q = ((uint128_t) t * winvtab[1]) >> 64;
    uint64_t v3 = t * wtab[1] - q * p;

    uint64_t z0 = v0 + v1;
    z0 -= (z0 >= 2*p) ? (2*p) : 0;
    uint64_t z1 = v0 - v1;
    z1 += ((int64_t) z1 < 0) ? (2*p) : 0;

    uint64_t z2 = v2 + v3;
    z2 -= (z2 >= 2*p) ? (2*p) : 0;
    uint64_t z3 = v2 - v3;
    z3 += ((int64_t) z3 < 0) ? (2*p) : 0;

    x[0] = z0 ;
    x[1] = z1 ;
    x[2] = z2 ;
    x[3] = z3 ;
  }

}

inline static void permut(uint64_t* const y, uint64_t const* const x, uint64_t const* const inv_indexes, size_t K)
{
  for (size_t i = 0; i < K; i += 1)
  {
    y[inv_indexes[i]] = x[i];
  }
}

// Inverse NTT: replaces values representation by the classic coefficient representation
inline void NFLlib::inv_ntt(uint64_t* const x, const uint64_t* const inv_wtab, const uint64_t* const inv_winvtab,
        const unsigned K, const uint64_t invK, const uint64_t p, const uint64_t* const inv_indexes)
{
	// Do we need to align on 32 ?
  uint64_t* y = calloc_align<16, uint64_t>(K+1);

  if (K == 1)
    return;

  // bit-reverse
  permut(y, x, inv_indexes, K);

  ntt(y, inv_wtab, inv_winvtab, K, p);

  // bit-reverse again
  permut(x, y, inv_indexes, K);
  
  free(y);
}


// This function just computes the powers of the root of unity and its inverse
// It is included here just to have all the NTT functions together
inline void NFLlib::prep_wtab(uint64_t* wtab, uint64_t* wtabshoup, uint64_t w, unsigned K, uint64_t p)
{

  while (K >= 2)
  {
    uint64_t wi = 1;     // w^i
    for (size_t i = 0; i < K/2; i++)
    {
      *wtab++ = wi;
      *wtabshoup++ = ((uint128_t) wi << 64) / p;
      wi = mulmod(wi, w, p);
    }
    w = mulmod(w, w, p);
    K /= 2;
  }
}
// *****************************************************************
// END OF NTT functions from David Harvey 
// From : http://web.maths.unsw.edu.au/~davidharvey/papers/fastntt/
// Most of the functions have been modified for our needs
// The associated copyright notice is at the beginning of this file
// *****************************************************************


// In order to have polynomial operations mod X**n + 1 as we want
// we must multiply the polynomial coefficients by phi powers before doing the NTT
inline void NFLlib::nttAndPowPhi(poly64 op) 
{
  for (unsigned int cm = 0 ; cm < nbModuli ; cm++)
  {
    for (unsigned int i = 0; i < polyDegree; i++)
    {
      op[i] = mulmodShoup(op[i], phis[cm][i], shoupphis[cm][i], moduli[cm]);
    }
    ntt(op, omegas[cm], shoupomegas[cm], polyDegree, moduli[cm]);
    op+=polyDegree;
  }
}


// In order to have polynomial operations mod X**n + 1 as we want
// we must multiply the polynomial coefficients by invphi powers before doing the inverse-NTT
inline void NFLlib::invnttAndPowInvPhi(poly64 op) 
{
  for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) 
  {
    inv_ntt(op, invomegas[currentModulus], shoupinvomegas[currentModulus], 
        polyDegree, invpolyDegree[currentModulus] , moduli[currentModulus], inv_indexes);
    for (unsigned int i = 0; i < polyDegree; i++)
    {
      op[i] = mulmodShoup(op[i], invpoly_times_invphis[currentModulus][i],
          shoupinvpoly_times_invphis[currentModulus][i], moduli[currentModulus]);
    }
    op+=polyDegree;
  }
}




// *********************************************************
// Pre-computation functions
// *********************************************************

// Pre-compute quotients for Shoup's multiplication for all the coefficients of a polynomial
inline poly64 NFLlib::allocandcomputeShouppoly(poly64 x)
{
  poly64 res = (poly64) malloc(polyDegree*nbModuli*sizeof(uint64_t));
  poly64 res_orig = res;
  for (unsigned short currentModulus = 0; currentModulus < nbModuli ; currentModulus++)
  {
    for (unsigned int i = 0; i < polyDegree; i++)
    {
      res[i] = ((uint128_t) x[i] << 64) / moduli[currentModulus];
    }
    res+=polyDegree;
    x+=polyDegree;
  }
  return res_orig;
}


// All the other functions are in the cpp file. A more carefull study should be done to see what needs inlining and what not.

#endif
