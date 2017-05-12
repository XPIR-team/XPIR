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

#include "NFLlib.hpp"
#include <x86intrin.h>
	
void DEBUG_MESSAGE(const char *s, poly64 p, unsigned int n){
#ifdef DEBUG
  std::cout<<s;
  NFLlib::print_poly64hex(p,n);
#endif
}

// *********************************************************
//    Constructors and init functions
// *********************************************************
// The constructors are not able to set all the parameters and setNewParameters has to be called
// after them. The attribute alreadyInit reflects this uninitialized status.

NFLlib::NFLlib():
    alreadyInit(0),
    nbModuli(0),
    polyDegree(0)
{
}


// Compute all data needed to do NTT operations
// Depends on the moduli used and polyDegree
void  NFLlib::configureNTT()
{
  uint64_t phi, invphi, omega, invomega, temp;

  // If this is a reconfiguration free memory before reallocating it
  freeNTTMemory();

  // Allocate space for the first dimension of the arrays needed for the NTT and CRT
  phis = (uint64_t **) malloc(nbModuli * sizeof(uint64_t *));
  shoupphis = (uint64_t **) malloc(nbModuli * sizeof(uint64_t *));
  invpoly_times_invphis = (uint64_t **) malloc(nbModuli * sizeof(uint64_t *));
  shoupinvpoly_times_invphis = (uint64_t **) malloc(nbModuli * sizeof(uint64_t *));
  omegas = (uint64_t **) malloc(nbModuli * sizeof(uint64_t *));
  shoupomegas = (uint64_t **) malloc(nbModuli * sizeof(uint64_t *));
  invomegas = (uint64_t **) malloc(nbModuli * sizeof(uint64_t *));
  shoupinvomegas = (uint64_t **) malloc(nbModuli * sizeof(uint64_t *));  
  invpolyDegree = (uint64_t *) malloc(nbModuli * sizeof(uint64_t));
  liftingIntegers = new mpz_t[nbModuli];
  moduli=new uint64_t[nbModuli]();

  // From now on, we have to do everything nbModuli times
  for(unsigned short currentModulus=0;currentModulus<nbModuli;currentModulus++) 
  {

    // Define the moduli we use (useful ?)
    moduli[currentModulus] = P64[currentModulus];
 
    // Allocation of the second dimension of the NTT parameters
    phis[currentModulus] = (uint64_t *) malloc(polyDegree * sizeof(uint64_t));
    shoupphis[currentModulus] = (uint64_t *) malloc(polyDegree * sizeof(uint64_t));
    invpoly_times_invphis[currentModulus] = (uint64_t *) malloc(polyDegree * sizeof(uint64_t));
    shoupinvpoly_times_invphis[currentModulus] = (uint64_t *) malloc(polyDegree * sizeof(uint64_t));
    omegas[currentModulus] = (uint64_t *) malloc((2 * polyDegree) * sizeof(uint64_t));
    shoupomegas[currentModulus] = omegas[currentModulus] + polyDegree;
    invomegas[currentModulus] = (uint64_t *) malloc((2 * polyDegree) * sizeof(uint64_t));
    shoupinvomegas[currentModulus] = invomegas[currentModulus] + polyDegree;
    
    // We start by computing phi 
    // The roots in the array are primitve 2**14-th roots 
    // Squared 14-log2(polyDegree) times they become
    // polyDegree-th roots as required by the NTT
    // But first we get phi = sqrt(omega) squaring them 13-log2(polyDegree) times 
    phi = primitive_roots[currentModulus];
    for (unsigned int i = 0 ; i < 13 - log2(polyDegree) ; i++)
    {
      phi = mulmod(phi, phi, moduli[currentModulus]);
    }
    
    // Now that temp = phi we initialize the array of phi**i values  
    // Initialized to phi**0
    temp = 1;
    for (unsigned int i = 0 ; i < polyDegree ; i++)
    {
      phis[currentModulus][i] = temp;
      shoupphis[currentModulus][i] = ((uint128_t) temp << 64) / moduli[currentModulus];
      // phi**(i+1)
      temp = mulmod(temp, phi, moduli[currentModulus]);
    }
    // At the end of the loop temp = phi**polyDegree

    // Computation of invphi
    // phi**(2*polydegree)=1 -> temp*phi**(polyDegree-1) = phi**(-1)
    invphi = mulmod(temp, phis[currentModulus][polyDegree-1], moduli[currentModulus]); 

    // Computation of the inverse of polyDegree using the inverse of kMaxPolyDegree
    invpolyDegree[currentModulus] = mulmod(invkMaxPolyDegree[currentModulus], 
        kMaxPolyDegree/polyDegree, moduli[currentModulus]);

    // Now we can compute the table invpoly_times_invphis
    temp = invpolyDegree[currentModulus];
    for (unsigned int i = 0 ; i < polyDegree ; i++)
    {
      invpoly_times_invphis[currentModulus][i] = temp;
      shoupinvpoly_times_invphis[currentModulus][i] = ((uint128_t) temp << 64) 
        / moduli[currentModulus];
      // This is invpolyDegree*invphi**(i+1)
      temp = mulmod(temp, invphi, moduli[currentModulus]);
    }

    // For the omegas it is easy, we just use the function of David Harvey modified for our needs
    omega = mulmod(phi, phi, moduli[currentModulus]);
    prep_wtab(omegas[currentModulus], shoupomegas[currentModulus], omega, polyDegree,
      moduli[currentModulus]);

    // And again for the invomegas
    invomega = mulmod(invphi, invphi, moduli[currentModulus]);
    prep_wtab(invomegas[currentModulus], shoupinvomegas[currentModulus], invomega, polyDegree,
        moduli[currentModulus]);

    // Inverse-CRT constants
    mpz_t tmpz, mpz_inverse;
    mpz_init(tmpz);
    mpz_init(mpz_inverse);

    // Compute first the product of all moduli but the current
    mpz_init_set_ui(liftingIntegers[currentModulus],1UL);
    for (unsigned int i = 0 ; i < nbModuli ; i++)
    {
      if (i == currentModulus) continue;
      mpz_import(tmpz, 1, 1, sizeof(uint64_t), 0, 0, &P64[i]); 
      mpz_mul(liftingIntegers[currentModulus], liftingIntegers[currentModulus], tmpz);
    }
    
    // Compute the inverse of the product modulo the current modulus and multiply it with the product
    mpz_import(tmpz, 1, 1, sizeof(uint64_t), 0, 0, &moduli[currentModulus]);
    mpz_invert(mpz_inverse, liftingIntegers[currentModulus], tmpz);
    mpz_mul(liftingIntegers[currentModulus], liftingIntegers[currentModulus], mpz_inverse);
    mpz_clear(tmpz);
    mpz_clear(mpz_inverse);
  }
  
  // Compute  the product of all moduli 
  mpz_t tmpz;
  mpz_init(tmpz);
  mpz_init_set_ui(moduliProduct,1UL);
  for (unsigned int i = 0 ; i < nbModuli ; i++)
  {
    mpz_import(tmpz, 1, 1, sizeof(uint64_t), 0, 0, &moduli[i]);
    mpz_mul(moduliProduct, moduliProduct, tmpz);
  }
  mpz_clear(tmpz);
  
  uint64_t* inv_indexes_tmp = malloc_align<32, uint64_t>(polyDegree);
  inv_indexes = malloc_align<16>(polyDegree, inv_indexes);

  // Compute permutation indexes for inv_ntt
  for (size_t i = 0; i < polyDegree; i++) 
  {
    size_t ii = i, r = 0; 
    for (unsigned h = 1; h < polyDegree; h=h<<1)
    {    
      r = (r << 1) | (ii & 1);
      ii >>= 1;
    }    

    inv_indexes_tmp[i] = r; 
  }

  // Invert the previous permutation. I think we can do better than that :)
  for (size_t i = 0; i < polyDegree; i++) {
      for (size_t j = 0; j < polyDegree; j++) {
          if (inv_indexes_tmp[j] == i) { 
              inv_indexes[i] = j; 
              break;
          }    
      }    
  }
  free(inv_indexes_tmp);
  
  alreadyInit = nbModuli;
}




// *********************************************************************
// Fundamental setters (resulting on a call to the init function above)
// *********************************************************************

// Set or reset the degrees of the polynomials used and the modulus size (which fixes the modulus)
void NFLlib::setNewParameters(unsigned int polyDegree_, unsigned int aggregatedModulusBitsize_)
{
  // We don't use setpolyDegree to avoid configureNTT being called twice
  polyDegree = polyDegree_;

  // We use a special function for setting the modulus as it can be called independently 
  // and requires some processing. This function also ensures that the NTT params are configured.
  setmodulus(aggregatedModulusBitsize_);
}


// Sets the modulus size AND configures NTT parmeters (which fixes a given modulus)
void NFLlib::setmodulus(uint64_t aggregatedModulusBitsize_)
{
  // For the CRT, from the aggregated modulus bitsize, we compute the number of necessary moduli
  if (aggregatedModulusBitsize_ % kModulusBitsize != 0)
  {
    std::cout << "NFLlib: CRITICAL. Modulus of " << aggregatedModulusBitsize_ 
      << " requested but only integer multiples of " << kModulusBitsize << " bits implemented. Exiting ..." << std::endl;
    exit(-1);
  }
  nbModuli=aggregatedModulusBitsize_/kModulusBitsize;
  
  configureNTT();
}


// Sets polyDegree AND configures NTT
void NFLlib::setpolyDegree(unsigned int polyDegree_)
{
  polyDegree = polyDegree_;

  configureNTT();
}




// *********************************************************
// Getters
// *********************************************************

uint64_t* NFLlib::getmoduli() { return moduli; }
unsigned short NFLlib::getnbModuli() { return nbModuli; }
unsigned int NFLlib::getpolyDegree() { return polyDegree; }
void NFLlib::copymoduliProduct(mpz_t dest) { mpz_init_set(dest, moduliProduct); }




// **************************************
// Random polynomial generation functions
// **************************************

// Two modes uniform or bounded (if uniform is false)
// WARNING : The bounded mode only works for a bound 
// below the smaller of the moduli -> we use or for bounded noise

// Allocates and sets a bounded random polynomial in FFT form calling setBoundedRandomPoly
poly64 NFLlib::allocBoundedRandomPoly(uint64_t upperBound_, bool uniform_) 
{
  poly64 res = (poly64)calloc(polyDegree * nbModuli, sizeof(uint64_t));
  setBoundedRandomPoly(res, upperBound_, uniform_);
  return res;
}


// Sets a pre-allocated random polynomial in FFT form
// If uniform = true upperBound is ignored and the coefficients are uniformly random
// ASSUMPTION: if uniform = false upperBound is below all of the moduli used
void NFLlib::setBoundedRandomPoly(poly64 res, uint64_t upperBound_, bool uniform_) 
{
  poly64 rnd, rnd_orig;
  uint64_t mask;

  if (uniform_ == false){
    // In bounded mode upperBound must be below the smaller of the moduli
    for (unsigned int cm = 0 ; cm < nbModuli ; cm++)
    {
      if (upperBound_ >= moduli[cm]) 
      {
        std::cout << "NFLlib: upperBound is larger than the moduli in setBoundedRandomPoly.";
        std::cout << " Unpredictable results ..." << std::endl;
        break;
      }
    }

    // We play with the rnd pointer (in the uniform case), and thus
    // we need to remember the allocated pointer to free it at the end
    rnd_orig = (poly64) malloc(polyDegree * sizeof(uint64_t)); 
    rnd = rnd_orig; 

    // Get some randomness from the PRNG
    fastrandombytes((unsigned char *)rnd, polyDegree * sizeof(uint64_t));

    // upperBound is below the moduli so we create the same mask for all the moduli
    mask=(1ULL <<  (unsigned int)ceil(log2(upperBound_))) -1;    
    
    for(unsigned int i=0;i<polyDegree;i++) {
    
      // First remove the heavy weight bits we dont need
      rnd[i]=(rnd[i]&mask);
    
      // When the random is still too large, reduce it 
      // In order to follow strictly a uniform distribution we should
      // get another rnd but in order to follow the proofs of security
      // strictly we should also take noise from a gaussian ...
      if (rnd[i]>=upperBound_) 
      {
        rnd[i]-=upperBound_;
      }
      for (unsigned int cm = 0 ; cm < nbModuli ; cm++)
      {
        res[polyDegree*cm+i] = rnd[i];
      }
    }
  }
  else // uniform == true
  {  
    // In uniform mode we need randomness for all the polynomials in the CRT
    rnd_orig = (poly64) malloc(polyDegree * nbModuli * sizeof(uint64_t)); 
    // We play with the rnd pointer (in the uniform case), and thus
    // we need to remember the allocated pointer to free it at the end
    rnd = rnd_orig; 
    fastrandombytes((unsigned char *)rnd, polyDegree * nbModuli * sizeof(uint64_t));
  
    for (unsigned int cm = 0 ; cm < nbModuli ; cm++)
    {
      // In the uniform case, instead of getting a big random (within the general moduli),
      // We rather prefer, for performance issues, to get smaller randoms for each module
      // The mask should be the same for all moduli (because they are the same size)
      // But for generality we prefer to compute it for each moduli so that we could have
      // moduli of different bitsize
      
      mask=(1ULL << (int)ceil(log2(moduli[cm]))) -1;    
    
      for(unsigned int i=0;i<polyDegree;i++) 
      {
        // First remove the heavy weight bits we dont need
        rnd[i]=(rnd[i]&mask);
    
        // When the random is still too large, reduce it 
        if (rnd[i]>=moduli[cm]) 
        {
          rnd[i]-=moduli[cm];
        }
        res[i] = rnd[i];
      }
      rnd+=polyDegree;
      res+=polyDegree;
    }
  }
  free(rnd_orig);
}




// *********************************************************
// Data import and export main functions
// *********************************************************

// Takes an array of buffers and:
// 1) Converts them into a set of polynomials with arbitrary large coefficients
// 2) Reduces the polys through CRT to have nbModuli contiguous polys with uint64_t coefficients 
// 3) Does the NTT transform
// - inArrayOfBuffers array of buffers to take the bits from
// - nbrOfBuffers nbr of buffers in the array
// - dataBitsizePerBuffer bits that can be taken from each buffer
// - bitsPerCoordinate bits used to create each coefficient (can be > 64 !)
// - polyNumber set by the function to say how many polynomials are in the returned pointer
poly64 *NFLlib::deserializeDataNFL(unsigned char **inArrayOfBuffers, uint64_t nbrOfBuffers, uint64_t dataBitsizePerBuffer, unsigned bitsPerCoordinate, uint64_t &polyNumber) {

  // We need to handle dataBitsize bits of data per buffer
  // each poly can take publicParams.getAbsorptionBitsize() bits so
  polyNumber = ceil((double)dataBitsizePerBuffer*(double)nbrOfBuffers/(double)(bitsPerCoordinate*polyDegree));

  // The uint64_t arrays are allocated and filled with zeros
  // So that we do not have to pad with zeros beyond the limit
  poly64* deserData = (poly64 *) calloc(polyNumber, sizeof(poly64));

  // bitsplitter does all the hard work WITHOUT using large numbers !
  deserData[0] = bitsplitter(inArrayOfBuffers, nbrOfBuffers, dataBitsizePerBuffer, bitsPerCoordinate);
 
  // We finish the work by applying the NTT transform
#ifdef MULTI_THREAD
  #pragma omp parallel for 
#endif
  for (unsigned int i = 0 ; i < polyNumber ; i++)
  {
    deserData[i] = deserData[0]+i*nbModuli*polyDegree;
#ifndef SIMULATE_PRE_NTT_DATA
    nttAndPowPhi(deserData[i]);
#endif
  }
  
  return deserData;
}

// Serialize an array of poly64 elements into a compact byte buffer
// Takes a set of polynomial coefficients and outputs their concatenation
// - indata points to the polynomial coefficients
// - outdata points to the concatenation obtained
// - bitsPerChunk defines how many bits has each coefficient 
// - nb_of_uint64 defines how many coefficients must be concatenated
// ASSUMPTION: all the polynomials are contiguously allocated
// ASSUMPTION: outdata has allocated one more uint64_t than needed
// ASSUMPTION:  all the coefficients have the same size which is below 56 bits
void NFLlib::serializeData64 (uint64_t* indata, unsigned char* outdata, unsigned int bitsPerChunk, uint64_t nb_of_uint64)
{
  unsigned char *tmppointer;
  uint64_t *pointer64;
  pointer64 = (uint64_t *) outdata;
  uint32_t bitswritten=0;

  // Tricky approach playing with pointers to be able to infinitely add
  // up to 56 bits to any bit string present
  for (uint64_t i = 0 ; i < nb_of_uint64 ;)
  {
    while(bitswritten + bitsPerChunk <= 64)
    {
      *pointer64 |= (*indata++)<<bitswritten; i++;
      bitswritten += bitsPerChunk;
      if(i==nb_of_uint64) break;
    }
    tmppointer = (unsigned char*) pointer64;
    tmppointer+=bitswritten>>3;
    pointer64 = (uint64_t *) (tmppointer);
    bitswritten -=8*(bitswritten>>3);
  }
}


// Serialize an array of poly64 elements into a compact byte buffer
// Takes a set of polynomial coefficients and outputs their concatenation
// - indata points to the polynomial coefficients
// - outdata points to the concatenation obtained
// - bitsPerChunk defines how many bits has each coefficient
// - nb_of_uint64 defines how many coefficients must be concatenated
// ASSUMPTION: all the polynomials are contiguously allocated
// ASSUMPTION: outdata has allocated one more uint64_t than needed
// IMPORTANT NOTE: Unlike in serializeData64 bitsPerChunk can be arbitrarily large. This function considers that large coefficients are retrieved by blocks of varying size up to 32 bit. For example 100-bit coefficients will be retrieved by three blocks of 32 bits and a block of 4 looping that way for each 100-bit coefficient.
void NFLlib::serializeData32 (uint32_t* indata, unsigned char* outdata, unsigned int bitsPerChunk, uint64_t nb_of_uint32){
  unsigned char *tmppointer;
  uint64_t *pointer64;
  pointer64 = (uint64_t *) outdata;
  uint32_t bitswritten=0;

  // See through how many block (=sub-chunk) we will have to loop to get each coefficient (=chunk)
  const double uint32PerChunk = (double)bitsPerChunk/32;
  const uint64_t int_uint32PerChunk = ceil(uint32PerChunk);
  const bool isint_uint32PerChunk = (uint32PerChunk==(double)int_uint32PerChunk);
 
  // Build masks for each sub-chunk
  uint64_t subchunkMasks[int_uint32PerChunk];
  unsigned int subchunkSizes[int_uint32PerChunk];
  
  // Increment with subchunkIndex=((subchunkIndex+1)%int_uint64PerChunk) 
  // and use with subchunkSizes[subchunkIndex]; 
  unsigned int subchunkIndex = 0;  
  for (int i = 0 ; i < int_uint32PerChunk - 1 ; i++)
  {
    subchunkSizes[i]=32;
    subchunkMasks[i] = (1ULL<<32)-1; // const mask for extracting 32 bits
  }
  subchunkSizes[int_uint32PerChunk-1] = bitsPerChunk - 32 * (int_uint32PerChunk - 1);
  subchunkMasks[int_uint32PerChunk-1] = (1ULL<<(subchunkSizes[int_uint32PerChunk-1]))-1;
  
  // Apply the same approach than in serializeData64 but with varying sizes
  for (uint64_t i = 0 ; i < nb_of_uint32 ;)
  {
    while(bitswritten + subchunkSizes[subchunkIndex] <= 64)
    {
      *pointer64 |= ((uint64_t)(*indata++))<<bitswritten; i++;
      bitswritten += subchunkSizes[subchunkIndex];
      subchunkIndex=((subchunkIndex+1)%int_uint32PerChunk);
      if(i==nb_of_uint32) break;
    }
    tmppointer = (unsigned char*) pointer64;
    tmppointer+=bitswritten>>3;
    pointer64 = (uint64_t *) (tmppointer);
    bitswritten -=8*(bitswritten>>3);
  }
}



// *********************************************************
// Helper functions
// *********************************************************

// Allocate a polynomial potentially with all coefficients set to zero if nullpoly = true
poly64 NFLlib::allocpoly(bool nullpoly)
{
  if (nullpoly == true) return (poly64) calloc(polyDegree*nbModuli,sizeof(uint64_t));
  else return (poly64) malloc(polyDegree*nbModuli*sizeof(uint64_t));
}


// Lift a polynomial in CRT representation, into a polynomial with large integer coefficients
mpz_t* NFLlib::poly2mpz(poly64 p)
{
  mpz_t* resultmpz;
  mpz_t* tmpzbuffer;
  resultmpz=new mpz_t[polyDegree];
  tmpzbuffer=new mpz_t[nbModuli];
  for(unsigned i=0;i<polyDegree;i++) {
    mpz_init2(resultmpz[i],192);
  }
  for(int cm = 0; cm < nbModuli;cm++) {
    mpz_init2(tmpzbuffer[cm],192);
  }
  
  for(unsigned i=0;i<polyDegree;i++) {
    mpz_set_ui(resultmpz[i],0UL);
    
    for(int cm = 0; cm < nbModuli;cm++) {
      mpz_import(tmpzbuffer[cm], 1, 1, sizeof(uint64_t), 0, 0, p+i+polyDegree*cm);
      mpz_mul(tmpzbuffer[cm], liftingIntegers[cm], tmpzbuffer[cm]);
      mpz_add(resultmpz[i], resultmpz[i], tmpzbuffer[cm]);
    }
    mpz_mod(resultmpz[i], resultmpz[i], moduliProduct);
  }
  for(int cm = 0; cm < nbModuli;cm++) {
  	mpz_clear(tmpzbuffer[cm]);
  	}
  delete[] tmpzbuffer;
  return resultmpz;
}


// Debug printing function
void NFLlib::print_poly64hex(poly64 p, unsigned int coeff_nbr)
{
 std::cout << "[";
 for (unsigned int i = 0 ; i < coeff_nbr ; i++)
 {
   std::cout << std::hex << (unsigned int)(p[i]>>32)<<(unsigned int) p[i] << " ";
 }
 std::cout << "]" << std::dec << std::endl;
}

// Debug printing function
void NFLlib::print_poly64(poly64 p, unsigned int coeff_nbr)
{
 std::cout << "[";
 for (unsigned int i = 0 ; i < coeff_nbr ; i++)
 {
   std::cout << p[i] << " ";
 }
 std::cout << "]" << std::endl;
}



// *********************************************************
// Destructors and closing or freeing functions
// *********************************************************

// Destructor
NFLlib::~NFLlib()
{
  freeNTTMemory();
}


// Everything is commented out because of hard to predict issues with MacOsX
// Uncomment on a computer with that OS before releasing
void NFLlib::freeNTTMemory(){
  // alreadyInit says how many arrays we have allocated for the NTT
  for (unsigned i = 0 ; i < alreadyInit; i++)
  {
    free(phis[i]);
    free(shoupphis[i]);
    free(invpoly_times_invphis[i]);
    free(shoupinvpoly_times_invphis[i]);
    free(omegas[i]);
    free(invomegas[i]);
    mpz_clear(liftingIntegers[i]);
  
    if (i == alreadyInit - 1)
    {
      free(phis);
      free(shoupphis);
      free(invpoly_times_invphis);
      free(shoupinvpoly_times_invphis);
      free(omegas);
      free(shoupomegas);
      free(invomegas);
      free(shoupinvomegas);
      free(invpolyDegree);
      delete[] liftingIntegers;
      free(inv_indexes);
      mpz_clear(moduliProduct);
      delete[] moduli;
    }
  }

  alreadyInit = 0;
}





















// ****************************************************************************************
// THE DEN: Uncommented howling functions and pointer blood magic. Enter at your own risk.
// ****************************************************************************************

// We define first a back to back funtion to test our bitsplitter function
// If DEBUG_BITSPLIT_B2B the function is used, else it is ignored 

//#define DEBUG_BITSPLIT_B2B
#ifdef DEBUG_BITSPLIT_B2B
#define DEBUG_BITSPLIT
#define B2BTEST(inDataBuffers, nbrOfBuffers, bitsPerBuffer, bitsPerChunk, totalbitsread, bitsread, pointer64, bitsToRead) bitsplitter_backtoback_internal_test (inDataBuffers, nbrOfBuffers, bitsPerBuffer, bitsPerChunk, totalbitsread, bitsread, pointer64, bitsToRead)
#else
#define B2BTEST(inDataBuffers, nbrOfBuffers, bitsPerBuffer, bitsPerChunk, totalbitsread, bitsread, pointer64, bitsToRead) if(0);
#endif

uint64_t* NFLlib::bitsplitter_backtoback_internal_test (unsigned char** inDataBuffers, uint64_t nbrOfBuffers, uint64_t bitsPerBuffer, unsigned int bitsPerChunk, uint64_t totalbitsread,uint64_t bitsread, uint64_t* pointer64, unsigned int bitsToRead)
{
  unsigned int bufferIndex = totalbitsread/bitsPerBuffer;
  uint64_t bitPositionInBuffer = totalbitsread - bufferIndex*bitsPerBuffer;
  uint64_t bytePositionInBuffer = bitPositionInBuffer/8;
  unsigned int bitPositionInByte = bitPositionInBuffer%8;
  if (bitPositionInBuffer+sizeof(uint64_t)*8 > bitsPerBuffer)
  {
    std::cerr << "WARNING: Bitsplit goes beyond buffer size (let's hope it is allocated)" << std::endl;
  }
  if (bitPositionInBuffer + bitsToRead > bitsPerBuffer)
  {
    std::cerr << "CRITICAL: Bitsplit reads AND USES data beyond buffer space" << std::endl;
    std::cerr << "totalbitsread " << totalbitsread << std::endl;
    std::cerr << "bufferIndex " << bufferIndex << std::endl;
    std::cerr << "bitPositionInBuffer " << bitPositionInBuffer << std::endl;
    std::cerr << "bytePositionInBuffer " << bytePositionInBuffer << std::endl;
    std::cerr << "CRITICAL: On B2B test called for " << bitsToRead << " bits" << std::endl;
    exit(-1);
  }
  uint64_t b2bresult = ((*((uint64_t *)(inDataBuffers[bufferIndex]+bytePositionInBuffer)))>>bitPositionInByte) & ((1ULL<<bitsToRead)-1) ;
  uint64_t bitsplitresult = ((*pointer64)>>bitsread) & ((1ULL<<bitsToRead)-1);
  if (b2bresult  != bitsplitresult)
  {
    std::cerr << "CRITICAL: Bitsplit different from back to back function" << std::endl;
    std::cerr << "CRITICAL: Left " << b2bresult << std::endl;
    std::cerr << "CRITICAL: Right " << bitsplitresult << std::endl;
    exit(-1);
  }
  return 0;
}

inline void NFLlib::bs_loop (unsigned char** inDataBuffers, uint64_t nbrOfBuffers, uint64_t bitsPerBuffer, unsigned int bitsPerChunk, uint64_t *&tmpdata, uint64_t bufferIndex, uint64_t &bitsread, size_t &subchunkIndex)
{
  // We redefine the amount of uint64 that can be produced given that we may have already read
  uint64_t bitstoread = bitsPerBuffer - bitsread;
  double nbChunks = (double)(bitsPerBuffer-bitsread)/bitsPerChunk;
  uint64_t int_nbChunks = floor(nbChunks);

  // How many uint64_t are needed to encode a chunk
  const double uint64PerChunk = (double)bitsPerChunk/56;
  const uint64_t int_uint64PerChunk = ceil(uint64PerChunk);
  const bool  isint_uint64PerChunk = (uint64PerChunk==(double)int_uint64PerChunk);
  
  // Compute subchunk sizes and masks
  uint64_t subchunkMasks[int_uint64PerChunk];
  // Increment with subchunkIndex=((subchunkIndex+1)%int_uint64PerChunk) and 
  // use subchunkSizes[subchunkIndex]; 
  unsigned int subchunkSizes[int_uint64PerChunk];
  for (unsigned i = 0 ; i < int_uint64PerChunk - 1 ; i++)
  {
    subchunkSizes[i]=56;
    subchunkMasks[i] = (1ULL<<56)-1; // const mask for extracting 56 bits
  }
  subchunkSizes[int_uint64PerChunk-1] = bitsPerChunk - 56 * (int_uint64PerChunk - 1);
  subchunkMasks[int_uint64PerChunk-1] = (1ULL<<(subchunkSizes[int_uint64PerChunk-1]))-1;
  
#ifdef DEBUG_BITSPLIT

  for (int i = 0 ; i < int_uint64PerChunk ; i++)
  {
    std::cerr<<"bitsplit0 i="<<i<<std::endl;
    std::cerr<<"bitsplit0 subchunkSizes[i]="<<subchunkSizes[i]<<std::endl;
    std::cerr<<"bitsplit0 subchunkMasks[i]="<<subchunkMasks[i]<<std::endl;
  }
#endif

  // We compute how many extra subchunks are available taking into account that we may
  // start this new buffer on the middle of a chunk
  uint64_t supplementalSubchunks = 0;
  uint64_t cumulatedsize = 0;
  for (unsigned i = 0 ; i < int_uint64PerChunk ; i++)
  {
    if (cumulatedsize + subchunkSizes[(subchunkIndex + i) % int_uint64PerChunk] <= (bitsPerBuffer-bitsread)-int_nbChunks*bitsPerChunk)        
    {
      supplementalSubchunks++;
      cumulatedsize += subchunkSizes[(subchunkIndex + i) % int_uint64PerChunk];
    }
    else
    {
      break;
    }
  }
  uint64_t totalSubChunks=int_nbChunks*int_uint64PerChunk + supplementalSubchunks;
#ifdef DEBUG_BITSPLIT
  std::cerr<<"bitsplit1 nbrOfBuffers="<<nbrOfBuffers<<std::endl;
  std::cerr<<"bitsplit1 nbChunks="<<nbChunks<<std::endl;
  std::cerr<<"bitsplit1 int_nbChunks="<<int_nbChunks<<std::endl;
  std::cerr<<"bitsplit1 uint64PerChunk="<<uint64PerChunk<<std::endl;
  std::cerr<<"bitsplit1 int_uint64PerChunk="<<int_uint64PerChunk<<std::endl;
  std::cerr<<"bitsplit1 isint_uint64PerChunk="<<isint_uint64PerChunk<<std::endl;
  std::cerr<<"bitsplit1 totalSubChunks="<<totalSubChunks<<std::endl;
  std::cerr<<"bitsplit1 cumulatedsize " << cumulatedsize << std::endl;
  std::cerr<<"bitsplit1 bitsPerBuffer " << bitsPerBuffer << std::endl;
  std::cerr<<"bitsplit1 bitsread " << bitsread << std::endl;
  std::cerr<<"bitsplit1 nextsubchunksize " <<  subchunkSizes[(subchunkIndex)]<< std::endl;
  std::cerr<<"bitsplit1 nextsubchunk " <<  subchunkIndex<< std::endl;
#endif
  
  unsigned  char *tmppointer;
  uint64_t *pointer64;
  pointer64 = (uint64_t *) inDataBuffers[bufferIndex];
  uint64_t bitsremaining=0; 

  
  
  // Loop over the subchunks in the current buffer
  for (uint64_t i = 0 ; i < totalSubChunks ; )
  {
    // Get up to 64bits
    while (bitsread + subchunkSizes[subchunkIndex] <= 64)
    {
      *tmpdata = ((*pointer64)>>bitsread) & subchunkMasks[subchunkIndex];
      tmpdata++;i++;
      bitsread += subchunkSizes[subchunkIndex];
      subchunkIndex= (subchunkIndex+1 == int_uint64PerChunk ? 0 : subchunkIndex + 1);
      if(i==totalSubChunks) break;
    }
    if (bitstoread > 128)
    {
      unsigned shift = bitsread >>3; 
      tmppointer = (unsigned char*) pointer64;
      tmppointer += shift;
      pointer64 = (uint64_t *) (tmppointer);
      bitstoread-= shift<<3;
      bitsread -= shift<<3;
    } 
    else
    {
      tmppointer = (unsigned char*) pointer64;
      while ((64 - bitsread < subchunkSizes[subchunkIndex]) && bitstoread > 0)
      {
        tmppointer++;
        bitsread -= 8;
        bitstoread -= 8;
      }
      pointer64 = (uint64_t *) (tmppointer);
    }
  }

  // If there is a last partial subchunk in this buffer, read it part from this buffer and part 
  // from next buffer if available
  bitsremaining = (uint64_t) round((nbChunks-int_nbChunks)*bitsPerChunk - cumulatedsize );
#ifdef DEBUG_BITSPLIT
    std::cout<<"bitsplit2 bitsremaining (should be <56)="<<bitsremaining<<std::endl;
#endif
  if (bitsremaining !=0) 
  {
    size_t shift=(64-(bitsread+bitsremaining))/8;
    bitsread+=(shift<<3);
    tmppointer = (unsigned char*) pointer64;
    tmppointer-=shift; 
    pointer64 = (uint64_t *) (tmppointer);
    *tmpdata = ((*pointer64)>>bitsread) & ((1ULL<<bitsremaining)-1);
    // If there is another buffer to deal with, finish the current tmpdata uint64_t 
    if (bufferIndex < nbrOfBuffers - 1)  
    {
      pointer64 = (uint64_t *) inDataBuffers[bufferIndex+1];
      *tmpdata |= ((*pointer64)<<bitsremaining) & subchunkMasks[subchunkIndex];
      // We restart bitsread to the bits read in the new buffer
      bitsread = subchunkSizes[subchunkIndex] - bitsremaining;  
      subchunkIndex= (subchunkIndex+1 == int_uint64PerChunk ? 0 : subchunkIndex + 1);
      tmpdata++;
    }
  }
  else
  {
    bitsread = 0;
  }
}


inline void NFLlib::bs_finish(poly64 &outdata, uint64_t int_uint64PerChunk, uint64_t polyNumber, uint64_t* splitData, uint64_t nbrOfBuffers, uint64_t bitsPerBuffer, unsigned int bitsPerChunk)
{
  if(int_uint64PerChunk>1) {
    outdata=(poly64) calloc(polyNumber*nbModuli*polyDegree + 1,sizeof(uint64_t));   
    internalLongIntegersToCRT( splitData, outdata,   int_uint64PerChunk, ceil(((double)bitsPerBuffer*nbrOfBuffers)/bitsPerChunk));
    free(splitData);
  }
  else 
  {
    if (nbModuli > 1)
    {
      outdata=(poly64) calloc(polyNumber*nbModuli*polyDegree + 1,sizeof(uint64_t));   
      for (unsigned i = 0 ; i < polyNumber ; i++)
      {
        for (int cm = 0 ; cm < nbModuli ; cm++)
        {
          memcpy(outdata + i*polyDegree*nbModuli + cm*polyDegree, 
              splitData + i*polyDegree, polyDegree*sizeof(uint64_t));
        }
      }
	  free(splitData);
    }
    else
    {
      outdata=splitData;
    }
  }
}


// This function does all the hard work of deserializeDataNFL
// 1) Converts input into a set of polynomials with arbitrary large coefficients
// 2) Reduces the polys through CRT to have nbModuli contiguous polys with uint64_t coefficients 
// Do not try to understand it, it is a nightmare, we won't try to explain it :)
uint64_t* NFLlib::bitsplitter (unsigned char** inDataBuffers, uint64_t nbrOfBuffers, uint64_t bitsPerBuffer, unsigned int bitsPerChunk)
{
  // If you don't need to change me don't try to understand me
  // If you need to change me, build me again from scratch :)

  // How many uint64_t are needed to encode a chunk
  const double uint64PerChunk = (double)bitsPerChunk/56;
  const uint64_t int_uint64PerChunk = ceil(uint64PerChunk);
  // How many polynomials are needed to encode the data
  uint64_t polyNumber = 
    ceil((double)bitsPerBuffer*(double)nbrOfBuffers/(double)(bitsPerChunk*polyDegree));
  uint64_t* splitData =
    (uint64_t*)(calloc(polyNumber*polyDegree*int_uint64PerChunk+1,sizeof(uint64_t)));
  uint64_t* tmpdata=splitData;

  uint64_t bitsread=0;
  size_t subchunkIndex=0;
  
  // Loop over the buffers
  for (uint64_t h = 0 ; h < nbrOfBuffers ; h++)
  { 
    bs_loop (inDataBuffers, nbrOfBuffers, bitsPerBuffer, bitsPerChunk, 
        tmpdata, h, bitsread, subchunkIndex);
  }
 
  poly64 outdata;

  bs_finish(outdata, int_uint64PerChunk, polyNumber, splitData, nbrOfBuffers, bitsPerBuffer, bitsPerChunk);
  
  return outdata;
}



// Subroutine for bitsplitter, the daemonic function. This is the function that allows
// us to circumvect GMP
void NFLlib::internalLongIntegersToCRT(uint64_t* tmpdata, poly64 outdata, uint64_t int_uint64PerChunk, uint64_t totalNbChunks) 
{
  uint64_t* outdataPtr=outdata;
  uint64_t* indataPtr=tmpdata;
  uint64_t multiplier[nbModuli][int_uint64PerChunk];

  uint64_t* chunkParts[int_uint64PerChunk];

  for(int cm=0;cm<nbModuli;cm++) 
  {
    multiplier[cm][0]=1;
    for(unsigned j=1;j<int_uint64PerChunk;j++) 
    {
      multiplier[cm][j] = mulmod(multiplier[cm][j-1],1ULL<<56,moduli[cm]);  
    }
  }

  for(unsigned i=0;i<totalNbChunks;i++) 
  {
    for(unsigned j=0;j<int_uint64PerChunk;j++) 
    {
      for(int cm=0;cm<nbModuli;cm++) 
      {
      // set to zero before computation if not calloc'd
      *(outdataPtr+cm*polyDegree) += mulmod(*(indataPtr+j), multiplier[cm][j],moduli[cm]);  
      }
    }
    indataPtr+=int_uint64PerChunk;
    outdataPtr++;
    if((i+1)%polyDegree==0) 
    {
      outdataPtr += polyDegree*(nbModuli-1);
    }
  }    
}


// ****************************************************************************************
// YOU EXITED THE DEN ALIVE! Either you jumped a lot of lines,
// or you became completely insane ... in the latter case ...
// Welcome to the group!
// ****************************************************************************************
