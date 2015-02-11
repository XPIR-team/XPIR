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

#include <iostream>
#include <string>
#define BOOST_TEST_DYN_LINK        // this is optional
#define BOOST_TEST_MODULE TestCryptoLWE   // specify the name of your test module
#include <boost/test/included/unit_test.hpp>
#include "../crypto/NTTLWE.hpp"
#include "../crypto/BasicLWE.hpp"

/* Functions checked :
 *  - computeNewParameters
 *  - getSerializedCryptoParams
 *  - getAbsBitPerCipherText
 *  - getSerializedModulusBitsize
 *  - encryit ui
 *  - decrypt encrypted ui
 *  - getCryptoParams
 *  - getAllCryptoParams
 *  - estimateAbsTime
 *  - getCiphSize
 *  - getAbsorptionBitsize
 */

NTTLWE lwe_to_test;
BasicLWE lwe_ref;
std::string serialized_params = "LWE:80:1024:60:8"; 
bool shortversion = true; 

BOOST_AUTO_TEST_CASE(test_nttlwe_computeNewParameters)   
{
  lwe_to_test.getPublicParameters().computeNewParameters(serialized_params);
  lwe_ref.getPublicParameters().computeNewParameters(serialized_params);
  auto params_test = lwe_to_test.getSerializedCryptoParams(shortversion);
  static auto params_ref = lwe_ref.getSerializedCryptoParams(shortversion);

  BOOST_CHECK_EQUAL(params_test, params_ref);
}

BOOST_AUTO_TEST_CASE(test_nttlwe_getAbsBitPerCipherText)   
{
  lwe_to_test.getPublicParameters().computeNewParameters(serialized_params);
  lwe_ref.getPublicParameters().computeNewParameters(serialized_params);
  auto tested = lwe_to_test.setandgetAbsBitPerCiphertext(1, 8);
  static auto expected = lwe_ref.setandgetAbsBitPerCiphertext(1, 8);
  BOOST_CHECK_EQUAL(tested, expected);
}

BOOST_AUTO_TEST_CASE(test_nttlwe_getSerializedModulusBitsize)   
{
  lwe_to_test.getPublicParameters().computeNewParameters(serialized_params);
  lwe_ref.getPublicParameters().computeNewParameters(serialized_params);
  auto tested = lwe_to_test.getPublicParameters().getSerializedModulusBitsize();
  static auto expected = lwe_ref.getPublicParameters().getSerializedModulusBitsize();
  BOOST_CHECK_EQUAL(tested, expected);
}

BOOST_AUTO_TEST_CASE(test_nttlwe_getQuerySizeFromRecLvl)   
{
  lwe_to_test.getPublicParameters().computeNewParameters(serialized_params);
  lwe_ref.getPublicParameters().computeNewParameters(serialized_params);

  for (unsigned int i = 0 ; i < 4 ; i++)
  {
    auto tested = lwe_to_test.getPublicParameters().getQuerySizeFromRecLvl(i);
    static auto expected = lwe_ref.getPublicParameters().getQuerySizeFromRecLvl(i);
    BOOST_CHECK_EQUAL(tested, expected);
  }
}

BOOST_AUTO_TEST_CASE(test_nttlwe_encryptui_decrypt)   
{
  lwe_to_test.getPublicParameters().computeNewParameters(serialized_params);
  unsigned int to_cipher = 100;
  char* ciphered_data = lwe_to_test.encrypt(to_cipher, 1);
  char* clear_data = lwe_to_test.decrypt(ciphered_data, 0, 0,0);

  BOOST_CHECK_EQUAL(to_cipher, (unsigned int) clear_data[0]);
  free(ciphered_data);
  free(clear_data);
}

BOOST_AUTO_TEST_CASE(test_nttlwe_getCryptoParams)   
{
  lwe_to_test.getPublicParameters().computeNewParameters(serialized_params);
  lwe_ref.getPublicParameters().computeNewParameters(serialized_params);

  unsigned int security_bits = 80;
  std::vector<std::string> params_expected, params_tested;

  auto params_nbr = lwe_ref.getCryptoParams(security_bits, params_expected);
  lwe_to_test.getCryptoParams(security_bits, params_tested);

  for (unsigned int i = 0; i < params_nbr ; i++)
  {
    BOOST_CHECK_EQUAL(params_tested[i], params_expected[i]);
  }
}


BOOST_AUTO_TEST_CASE(test_nttlwe_getAllCryptoParams)   
{
  lwe_to_test.getPublicParameters().computeNewParameters(serialized_params);
  lwe_ref.getPublicParameters().computeNewParameters(serialized_params);

  std::vector<std::string> params_expected, params_tested;

  auto params_nbr = lwe_ref.getAllCryptoParams(params_expected);
  lwe_to_test.getAllCryptoParams(params_tested);

  for (unsigned int i = 0; i < params_nbr ; i++)
  {
    BOOST_CHECK_EQUAL(params_tested[i], params_expected[i]);
  }
}

BOOST_AUTO_TEST_CASE(test_nttlwe_estimateAbsTime)   
{
  lwe_to_test.getPublicParameters().computeNewParameters(serialized_params);
  lwe_ref.getPublicParameters().computeNewParameters(serialized_params);

  std::vector<std::string> params_expected, params_tested;

  auto params_nbr = lwe_ref.getAllCryptoParams(params_expected);
  lwe_to_test.getAllCryptoParams(params_tested);

  for (unsigned int i = 0; i < params_nbr ; i++)
  {
    BOOST_CHECK_EQUAL(lwe_to_test.estimateAbsTime(params_tested[i]), lwe_ref.estimateAbsTime(params_expected[i]));
  }
}

BOOST_AUTO_TEST_CASE(test_nttlwe_getAbsorptionBitsize)   
{
  lwe_to_test.getPublicParameters().computeNewParameters(serialized_params);
  lwe_ref.getPublicParameters().computeNewParameters(serialized_params);

  for (auto i : {1, 2, 3, 7, 10})
  {
    BOOST_CHECK_EQUAL(lwe_to_test.getPublicParameters().getAbsorptionBitsize(i), lwe_ref.getPublicParameters().getAbsorptionBitsize(i));
  }
}

BOOST_AUTO_TEST_CASE(test_nttlwe_getCiphBitsizeFromRecLvl)   
{
  lwe_to_test.getPublicParameters().computeNewParameters(serialized_params);
  lwe_ref.getPublicParameters().computeNewParameters(serialized_params);

  for (auto i : {1, 2, 3, 7, 10})
  {
    BOOST_CHECK_EQUAL(lwe_to_test.getPublicParameters().getCiphBitsizeFromRecLvl(i), lwe_ref.getPublicParameters().getCiphBitsizeFromRecLvl(i));
  }
}

//BOOST_AUTO_TEST_CASE(test_nttlwe_pseudo_PIR)
//{
//  unsigned int poly_nbr = 10, poly_degree = 1024, indexOmega, minPolyDegree = 512;
//  NTTLWE lwe;
//  lwe.getPublicParameters().computeNewParameters("lwe:80:1024:64:10");
//  lwe_in_data data_to_absorb[poly_nbr];
//  poly64 poly;
//
//  for (indexOmega = 0; (1<<indexOmega) * minPolyDegree < poly_degree; indexOmega++);
//
//  for (unsigned int i = 0 ; i < poly_nbr ; i++) {
//    //poly  = lwe.boundedRandomPoly(poly_degree, 1<<9);
//    poly = (poly64) calloc(poly_degree, sizeof(uint64_t));
//    poly[0] = 1;
//    data_to_absorb[i].p  = &poly;
//    data_to_absorb[i].nbPolys = 1;
//  }
//   /* for (unsigned int j = 0 ; j < poly_degree ; j++) lwe.mulmod(data_to_absorb[i].p[0][j], phi_modP64[indexOmega][j], P64);
//
//    lwe.ntt_new(data_to_absorb[i].p[0], lwe.wtab, lwe.winvtab, poly_degree, P64);
//  }*/
//
//  lwe_query query_elements[poly_nbr];
//  poly64 zero = (poly64) calloc(poly_degree, sizeof(uint64_t));
//
//  for (unsigned int i = 0 ; i <poly_nbr - 1; i++) lwe.enc((lwe_cipher*)(&query_elements[i]), zero);
//
//  poly64 one = zero;
//  one[0] = static_cast<uint64_t>(1);
//
//  lwe.enc((lwe_cipher*)(&query_elements[poly_nbr - 1]), one);
//
//  free(one);
//
//  lwe_cipher result;
//  result.a = (poly64) calloc(poly_degree * 2, sizeof(uint64_t)); 
//  result.b = result.a + poly_degree;
//
//  for (unsigned int i = 0 ; i < poly_nbr ; i++) lwe.mulandadd(result, data_to_absorb[i], query_elements[i], 0);
//  
//  poly64 clear_data_poly = (poly64) calloc(poly_degree,sizeof(uint64_t));
//
//  lwe.dec(clear_data_poly, &result);
//
//  for (unsigned int i = 0 ; i < 10/*poly_degree*/ ; i++)
//  {
//    std::cout << clear_data_poly[i] << " ";
//  }
//  std::cout << std::endl;
//  BOOST_CHECK_EQUAL(clear_data_poly[0], data_to_absorb[poly_nbr - 1].p[0][0]);
//
//  for (unsigned int i = 0 ; i < poly_nbr; i++)
//  {
////    free(&data_to_absorb[i].p);
//    free(query_elements[i].a);
//  }
//  free(result.a);
//  free(clear_data_poly);
//  
//}
/*
 * Benchmark here.
 */
//BOOST_AUTO_TEST_CASE(test_nttlwe_intern_functions)
//{
//  int rounds = 100;
//  NTTLWE n;
//  n.setNewParameters(1024,64,10);
//  n.setmodulus(P64);
//  poly64 p=n.boundedRandomPoly(1024, P64);
//  poly64 p2=(poly64)calloc(1024,sizeof(uint64_t));
//  lwe_cipher cyph;//=(lwe_cipher *)calloc(1024*2,sizeof(uint64_t));
//
//  double start = omp_get_wtime();
//  for(int i = 0 ; i < rounds ; i++) {
//    n.enc(&cyph,p);
//    free(cyph.a);
//  }
//  double end = omp_get_wtime();
//
//  std::cout<<rounds/(end - start)<<" chiffre/s"<<std::endl;
//
//  start = omp_get_wtime();
//  for (int i = 0 ; i < rounds ; i++) {
//    n.dec(p2,&cyph);
//  }
//  end = omp_get_wtime();
//
//  std::cout<<rounds/(end-start)<<" dechiffre/s"<<std::endl;
//
//  free(p2);
//  free(p);
//
// for(int i        = 0;i<1024;i++) {
// 	if(p2[i]!=p[i]) {
// 		std::cout<<"err"<<std::endl;
// 	}
// }

//  start = omp_get_wtime();
//  for(int i = 0;i<100000;i++) {
//    poly64 p=n.boundedRandomPoly(1024, P64);
//    // for(int j=0;j<1024;j++) {
//    // 		if(p[j]>P64) {
//    // 			std::cerr<<"ERROR"<<std::endl;
//    // 			exit(-1);
//    // 		}
//    // 		//std::cout<<(unsigned long long)p[j]<<" ";
//    // 		}
//    // 		//std::cout<<std::endl;
//  }
//
//  end        = omp_get_wtime();
//  long nbULL        = 100000*1024;
//  double nbULLpeSec = nbULL/(end-start);
//  std::cout<<nbULLpeSec<<"ULL/s"<<std::endl;
//  std::cout<<nbULLpeSec*8*8/1000000000<<" Gbps with ULL mod P"<<std::endl;
//  unsigned char rndbuffer[1024*sizeof(uint64_t)];
//  start     = omp_get_wtime();
//  for(int i        = 0;i<100000;i++) {
//    fastrandombytes((unsigned char *)rndbuffer, 1024);
//  }
//  end        = omp_get_wtime();
//  std::cout<<100000*1024*8/(end-start)/1000000000<<" Gbps with raw chars"<<std::endl;
//}
