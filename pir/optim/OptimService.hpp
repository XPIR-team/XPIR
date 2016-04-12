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

#ifndef DEF_OPTIM_SERVICES
#define DEF_OPTIM_SERVICES

#include <map>
#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>

#include "OptimVars.hpp"
#include "crypto/HomomorphicCrypto.hpp"
#include "crypto/HomomorphicCryptoFactory_internal.hpp"

using namespace std;

class OptimService
{
  protected:
    static const std::string testValuesFileName;
    static const std::string getCurrentTime();

  public:
    static const std::string folderName;
    static const std::string fileName;
    static const std::string absFileExtension;
    static const std::string encFileExtension;
    static const std::string decFileExtension;

    static bool fileOutdated(std::string crypto_name, std::string extension);

    static int readEntireFile(std::string& file_content, const std::string& file_path);
    static void getAllOptimData(std::vector<FixedVars>& fixed_vars_vec, std::string testValuesFileName);
    static int readOptimData(map<std::string, double>& values, const std::string& file_path);
    static int writeOptimData(double encrypt_time, double decrypt_time, std::string crypto_params_desc, std::string crypto_name);
    static int writeOptimDataBuffer(const std::string& buffer, const std::string& file_path);
    static int readTestValues(unsigned int i, FixedVars& vars, std::string testValuesFileName);
    
    static void gotoLine(std::ifstream& file, unsigned int num);
    static int getNumberOfExperiences(std::string testValuesFileName);
    static unsigned int getNumberOfLines(std::ifstream& f);

    static int verifyOptimData(set<string> crypto_params_set, const std::string& fenc_path, const std::string& fdec_path);
 
    static void writeHeadFile(unsigned int i, FixedVars& fixedVars);
    static void writeTestCurrentResult(unsigned int alpha, unsigned int alphaMul, unsigned int d, unsigned int i, OptimVars& vars);
    static void writeTestCurrentResult(unsigned int alpha_min, unsigned int alpha_max, unsigned int alpha_curr, unsigned int a_inf_bound, unsigned int a_sup_bound, unsigned int d, unsigned int i, OptimVars& vars);
    static void writeFootFile(unsigned int i);

    static void writeMessage(unsigned int i, std::string const& message);
    static void writeConfigFile(unsigned int alpha, unsigned int alphaMul, unsigned int d, unsigned int exp_nbr);

  static void writeLWEFile(unsigned int order, unsigned int p_size, unsigned int exp_nbr);
};

#endif
