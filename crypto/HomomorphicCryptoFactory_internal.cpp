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

#include "HomomorphicCryptoFactory_internal.hpp"


const unsigned int HomomorphicCryptoFactory_internal::crypto_method_nbr = 2;
const std::vector<std::string>  HomomorphicCryptoFactory_internal::crypto_method_name_vec = [] //use lambda function to initialise the vector
{
  std::vector<std::string> v;
  v.push_back("Paillier");
  v.push_back("LWE");
  return v;
}();

HomomorphicCrypto* HomomorphicCryptoFactory_internal::getCrypto(std::string cryptoType)
{
  HomomorphicCrypto* h;
  
  if (cryptoType == "Paillier")
  {
    h = new PaillierAdapter();
  }
  else if (cryptoType == "LWE")
  {
    h = new NFLLWE();
  }
  else if(cryptoType == "NoCryptography")
  {
    h = new NoCryptography();
  }
  else
  {
    std::cerr << "HomomorphicCryptoFactory_internal: Warning, unrecognized cryptosystem. Returning NULL" << std::endl;
    h = NULL;
  }
  
  return h;
}
  
HomomorphicCrypto* HomomorphicCryptoFactory_internal::getCryptoMethod(std::string crypto_system_desc)
{
  std::vector<std::string> fields;
  boost::algorithm::split(fields, crypto_system_desc, boost::algorithm::is_any_of(":"));
  HomomorphicCrypto* h;

  h = getCrypto(fields[0]);

  h->setNewParameters(crypto_system_desc);
  return h;
}

  
// Returns a vector with all the cryptosystems THAT WE WANT TO BE USED
void HomomorphicCryptoFactory_internal::getAllCryptoSystems(std::vector<HomomorphicCrypto*>& crypto_sys_vec)
{
  HomomorphicCrypto* crypto_ptr = new PaillierAdapter();
  crypto_sys_vec.push_back(crypto_ptr);
  crypto_ptr = new NFLLWE();
  crypto_sys_vec.push_back(crypto_ptr);
}

void HomomorphicCryptoFactory_internal::getOneCryptoSystem(std::vector<HomomorphicCrypto*>& crypto_sys_vec, std::string crypto_system_desc)
{
  std::vector<std::string> fields;
  boost::algorithm::split(fields, crypto_system_desc, boost::algorithm::is_any_of(":"));
  
  HomomorphicCrypto* h;
  h = getCrypto(fields[0]);
  crypto_sys_vec.push_back(h);
}


