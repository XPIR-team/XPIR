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

#ifndef PARAMS_C
#define PARAMS_C

#include "NFLParams.hpp"

// signed and unsigned 128-bit types
typedef int int128_t __attribute__((mode(TI)));
typedef unsigned int uint128_t __attribute__((mode(TI)));

// The moduli used in each 64 bit block (60 bits long each)
// They are of the form p = 2**61 - i*2**14 + 1 for increasing i
const unsigned int kMaxNbModuli = 17;
const uint64_t P64[kMaxNbModuli] = { 2305843009213317121ULL,  2305843009213120513ULL,  2305843009212694529ULL,2305843009212399617ULL,  2305843009211662337ULL,  2305843009211596801ULL, 2305843009211400193ULL,  2305843009210580993ULL,  2305843009210515457ULL, 2305843009210023937ULL,  2305843009209057281ULL,  2305843009208795137ULL, 2305843009208713217ULL,  2305843009208123393ULL,  2305843009207468033ULL, 2305843009206976513ULL,  2305843009206845441ULL};
const unsigned int kModulusBitsize = 60;
const unsigned int kModulusRepresentationBitsize = 64; 

// A primitive 2**14 root of unity for each one of the moduli
const uint64_t primitive_roots[kMaxNbModuli] = {1187132827279672845ULL, 753478288701480417ULL, 717492273700781398ULL, 1965831349000139179ULL, 2188074825493578002ULL, 232964966911499637ULL, 395215708149128273ULL, 241993652537061162ULL, 370764889455808762ULL, 1598724739621481826ULL, 610870591815427803ULL, 1880981984715338041ULL, 1389374064030612248ULL, 1508488901479281364ULL, 2111887543055470169ULL, 2020509703903079203ULL, 1120240637471598814ULL};

// Inverses of kMaxPolyDegree (for the other degrees it can be derived easily)
// for the different moduli
const uint64_t invkMaxPolyDegree[kMaxNbModuli] = {2305561534236606511ULL, 2305561534236409927ULL, 2305561534235983995ULL, 2305561534235689119ULL, 2305561534234951929ULL, 2305561534234886401ULL, 2305561534234689817ULL, 2305561534233870717ULL, 2305561534233805189ULL, 2305561534233313729ULL, 2305561534232347191ULL, 2305561534232085079ULL, 2305561534232003169ULL, 2305561534231413417ULL, 2305561534230758137ULL, 2305561534230266677ULL, 2305561534230135621ULL};

// Polynomial related data
const unsigned int kMinPolyDegree = 512;
const unsigned int kMaxPolyDegree = 8192;

#endif

