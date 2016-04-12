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

#include "ServerService.hpp"

const unsigned int ServerService::kPIRParamsOptions = 2;

int ServerService::writePIRParameters(const PIRParameters& pir_params, const std::string& file_path)
{
  std::ofstream file(file_path);

  if(!file.is_open()) return 1;

  file << "d " << pir_params.d << std::endl;
  file << "dim ";

  for (unsigned int i = 0; i < pir_params.d ; i++) file << std::to_string(pir_params.n[i]) << " ";

  file << std::endl;
  file << pir_params.crypto_params << std::endl;

  file.close();

  return 0;
}


int ServerService::readPIRParameters(PIRParameters& pir_params, const std::string& file_path)
{
  std::ifstream file(file_path);

  if (!file.is_open()) return 1;

  std::string line;
  std::vector<std::string> fields;

  std::getline(file, line);
  boost::algorithm::split(fields, line, boost::algorithm::is_any_of(" "));

  pir_params.d = atoi(fields.at(1).c_str());

  if (file.eof()) return 1;

  std::getline(file, line);
  boost::algorithm::split(fields, line, boost::algorithm::is_any_of(" "));

  unsigned int i;
  for (i = 0 ; i < pir_params.d && i < fields.size() - 1 ; i++)
  {
    pir_params.n[i] = atoi(fields.at(i+1).c_str());
  }

  if (i != pir_params.d) return 1;
  
  std::getline(file, line);

  pir_params.crypto_params = line;

  file.close();
  return 0;
}
