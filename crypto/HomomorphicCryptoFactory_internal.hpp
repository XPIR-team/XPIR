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

#ifndef DEF_HOMOMORPHIC_FACTORY
#define DEF_HOMOMORPHIC_FACTORY

#include <string>
#include "PaillierAdapter.hpp"
#include "NFLLWE.hpp"
#include "NoCryptography.hpp"

class HomomorphicCryptoFactory_internal
{
public:
  static HomomorphicCrypto* getCrypto(std::string cryptoType);
  static HomomorphicCrypto* getCryptoMethod(std::string cryptoType);
  static void getAllCryptoSystems(std::vector<HomomorphicCrypto*>& crypto_sys_vec);
  static void getOneCryptoSystem(std::vector<HomomorphicCrypto*>& crypto_sys_vec, std::string crypto_system_desc);
  static const unsigned int crypto_method_nbr;
  static const std::vector<std::string> crypto_method_name_vec;
};

#endif
