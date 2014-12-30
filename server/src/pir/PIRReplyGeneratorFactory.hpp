/* Copyright (C) 2014 Carlos Aguilar Melchor, Joris Barrier, Marc-Olivier Killijian
 * This file is part of XPIRe.
 *
 *  XPIRe is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  XPIRe is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XPIRe.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DEF_PIRREPLYGENERATORFACTORY
#define DEF_PIRREPLYGENERATORFACTORY

#include <string>
#include <vector>

#include "PIRParameters.hpp"
#include "../server/DBHandler.hpp"
#include "PIRReplyGeneratorGMP.hpp"
#include "PIRReplyGeneratorNFL.hpp"
#include "PIRReplyGeneratorTrivial.hpp"

class PIRReplyGeneratorFactory
{
  public:
    static GenericPIRReplyGenerator* getPIRReplyGenerator(const std::string& genName, vector <std::string>& database_, PIRParameters& param, DBHandler* db);
    static GenericPIRReplyGenerator* getPIRReplyGenerator(const std::string& genName);
};

#endif
