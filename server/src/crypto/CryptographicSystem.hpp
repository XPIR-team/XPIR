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

#ifndef DEF_CRYPTOGRAPHICSYSTEM
#define DEF_CRYPTOGRAPHICSYSTEM
#include <set>
#include <string>
#include "AbstractPublicParameters.hpp"
#include <boost/algorithm/string.hpp>

class CryptographicSystem 
{
  public:
    virtual AbstractPublicParameters& getPublicParameters()=0; 
    virtual unsigned int getAllCryptoParams(std::set<std::string>& crypto_params_set)=0;
    virtual void setNewParameters(const std::string& crypto_param)=0;
    virtual ~CryptographicSystem(){};
    virtual std::string& toString()=0;
};

#endif
