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

#ifndef PARAMS_H
#define PARAMS_H
#include <inttypes.h> 

// signed and unsigned 128-bit types
typedef int int128_t __attribute__((mode(TI)));
typedef unsigned int uint128_t __attribute__((mode(TI)));

// The moduli used in each 64 bit block (60 bits long each)
extern const unsigned int kMaxNbModuli;
extern const uint64_t P64[];
extern const unsigned int kModulusBitsize;
extern const unsigned int kModulusRepresentationBitsize;

// A primitive 2**14 root of unity for each one of the moduli
extern const uint64_t primitive_roots[];

// Inverses of kMaxPolyDegree (for the other degrees it can be derived easily)
// for the different moduli
extern const uint64_t invkMaxPolyDegree[];

// Polynomial related data
extern const unsigned int kMinPolyDegree;
extern const unsigned int kMaxPolyDegree;

#endif

