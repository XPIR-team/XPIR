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

#include "OptimService.hpp"

const std::string OptimService::folderName = "exp/";
const std::string OptimService::fileName = "preCompute";
const std::string OptimService::absFileExtension = ".abs";
const std::string OptimService::decFileExtension = ".dec";
const std::string OptimService::encFileExtension = ".enc";

const std::string OptimService::getCurrentTime()
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];

    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

void OptimService::getAllOptimData(std::vector<FixedVars>& fixed_vars_vec, std::string testValuesFileName)
{
  unsigned int experience_nbr = getNumberOfExperiences(testValuesFileName);

  for(unsigned int i = 0 ; i < experience_nbr ; i++)
  {
    FixedVars vars;
    readTestValues(i, vars, testValuesFileName);
    fixed_vars_vec.push_back(vars);
  }
}

int OptimService::getNumberOfExperiences(std::string testValuesFileName)
{
  std::string line;
  std::ifstream f(testValuesFileName);
  int i = 0;

  if(!f.is_open())
    return -1;

  while (std::getline(f, line)) 
  {  
    if(!(line.c_str()[0] == '#')) i++;//jump over commented line
  }
  f.close();
  return --i;
}

int OptimService::readTestValues(unsigned int i, FixedVars& vars, std::string testValuesFileName)
{
  std::ifstream	f(testValuesFileName, std::ios::in);
  std::string line;
  std::vector<std::string> fields;

  if (f.is_open())
  {
    for (unsigned k = 0 ; k < i+1 ; k++)
    {
      std::getline(f, line);
      if(line.c_str()[0] == '#') i++;//jump over commented line
    }
    boost::algorithm::split(fields, line, boost::algorithm::is_any_of(" "));

    vars.n 			  = atoi(fields[0].c_str());
    vars.l 			  = atoi(fields[1].c_str()); 
    vars.Tupc 		= vars.Tdos = static_cast<double>(atol(fields[2].c_str())); 
    vars.Tdoc 	  = vars.Tups = static_cast<double>(atol(fields[3].c_str()));
    vars.k 			  = atoi(fields[4].c_str()); 
    vars.alphaMax = atoi(fields[5].c_str());
    vars.dMax     = atoi(fields[6].c_str()); 
  }
  else
  {
    return 1;
  }

  f.close();

  return 0;
}

void OptimService::writeHeadFile(unsigned int i, FixedVars& fixedVars)
{
  std::ofstream file(std::string("exp/exp"+ std::to_string(i)).c_str(), std::ios::out);

  file << "# " << std::string("exp" + std::to_string(i)) << " " << getCurrentTime() << std::endl;
  file << "#Fixed Param : n " 	 << fixedVars.n 	 << ", l " 	  << fixedVars.l << ", Tupc "  << fixedVars.Tupc << ", Tdoc " << fixedVars.Tdoc  <<", k " << fixedVars.k << std::endl;
  file << "#Bound Param : " <<  ", alphaMax " << fixedVars.alphaMax << ", dMax " << fixedVars.dMax << std::endl;
  file << "#1:d 2a:alpha_min\t 2b:alpha_max\t 2c:current_best_alpha\t2d:alpha_lowbound\t2e:alpha_upbound\t3:GenQ    \t 4:SendQ  \t 5:GenR \t 6:SendR  \t 7:DecR \t 8:Total Time (pot. pipelined)" << std::endl;

  file.close();
}

int OptimService::readEntireFile(std::string& file_content, const std::string& file_path)
{
  std::ifstream f(file_path);

  if (!f.is_open()) 
  {
    f.close();
    return 1;
  }

  getline(f, file_content, (char)EOF);

  f.close();

  return 0;
}

int OptimService::readOptimData(map<std::string,double>& values, const std::string& file_path)
{
  std::ifstream	f(file_path);

  if(!f.good()) 
  {
    f.close();
    return 1;
  }

  std::string line;
  std::vector<std::string> fields;

  while (std::getline(f, line))
  {
    boost::algorithm::split(fields, line, boost::algorithm::is_any_of(" "));
    values[fields.at(0)] = atof(fields.at(1).c_str());
  }
  f.close();
  return 0;
}


/**
 * Write test result into file exp/exp{$exp_nbr}
 * Params are self explanatory
 **/
void OptimService::writeTestCurrentResult(unsigned int alpha_min, unsigned int alpha_max, unsigned int alpha_curr, unsigned int a_inf_bound, unsigned int a_sup_bound, unsigned int d, unsigned int exp_nbr, OptimVars& vars)
{
  // Open output file exp/exp{$exp_nbr}
  std::ofstream file(std::string("exp/exp"+ std::to_string(exp_nbr)).c_str(), std::ios::out | std::ios::app );

  // Try to output double values always with the same amount of decimals
  file.setf( std::ios::fixed, std:: ios::floatfield );

  // Output test result in a line
  file << d << "\t" << alpha_min << "\t" << alpha_max << "\t" << alpha_curr << "\t" << a_inf_bound << " \t" << a_sup_bound << "\t" << vars.getGenQ() << "   \t " << vars.getSendQ() << "   \t " << vars.getGenR() << "\t" << vars.getSendR() << "   \t " << vars.getDecR() << "\t" << vars.getValue()  << std::endl;

  file.close();
}

int OptimService::writeOptimData(double encrypt_time, double decrypt_time, std::string crypto_params_desc, std::string crypto_name)
{
  std::ofstream	fdec(folderName + fileName + crypto_name + decFileExtension, ios::app);
  std::ofstream fenc(folderName + fileName + crypto_name + encFileExtension, ios::app);
  fdec.setf( std::ios::fixed, std:: ios::floatfield );
  fenc.setf( std::ios::fixed, std:: ios::floatfield );

  if (!fdec.good() || !fenc.good()) 
  {
    fenc.close();
    fdec.close();
    return 1;
  }

  fenc << crypto_params_desc << " " << encrypt_time << std::endl;
  fdec << crypto_params_desc << " " << decrypt_time << std::endl;

  return 0;
}

int OptimService::writeOptimDataBuffer(const std::string& buffer, const std::string& file_path)
{
  std::ofstream	f(file_path);

  if (!f.good())
  {
    f.close();
    return 1;
  }

  f << buffer;
  
  f.close();
  return 0;
}


void OptimService::gotoLine(std::ifstream& file, unsigned int num)
{
  if (num > 0)
  {
    file.seekg(std::ios::beg);

    for (unsigned int i = 0 ; i < num ; ++i)
    {
      file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
  }
}

void OptimService::writeFootFile(unsigned int i)
{  
  std::ofstream file(std::string("exp/exp"+ std::to_string(i)).c_str(), std::ios::out | std::ios::app );

  file << std::endl <<"#End " << getCurrentTime() << std::endl;

  file.close();
}

void OptimService::writeTestCurrentResult(unsigned int alpha, unsigned int alphaMul, unsigned int d, unsigned int i, OptimVars& vars)
{
  std::ofstream file(std::string("exp/exp"+ std::to_string(i)).c_str(), std::ios::out | std::ios::app );

  file.setf( std::ios::fixed, std:: ios::floatfield );

  file << d << "\t" << alpha << " \t" << alphaMul << "\t" << vars.getGenQ() << "   \t " << vars.getSendQ() << "   \t " << vars.getGenR() << "\t" << vars.getSendR() << "   \t " << vars.getDecR() << "\t" << vars.getValue() << std::endl;

  file.close();
}

void OptimService::writeMessage(unsigned int i, std::string const& message)
{
  std::ofstream file(std::string("exp/exp"+ std::to_string(i)).c_str(), std::ios::out | std::ios::app );

  file << message << std::endl;

  file.close();
}

void OptimService::writeConfigFile(unsigned int alpha, unsigned int alphaMul, unsigned int d, unsigned int exp_nbr)
{
  std::ofstream file(std::string("configFile" + std::to_string(exp_nbr)).c_str(), std::ios::out);

  file << "alpha\t"  << alpha    << endl;
  file << "alphaM\t" << alphaMul << endl;

  file.close();
}

void OptimService::writeLWEFile(unsigned int order, unsigned int p_size, unsigned int exp_nbr)
{
  std::ofstream file(std::string("LWEFile" + std::to_string(exp_nbr)).c_str(), std::ios::out);

  file << "degree\t" << order  << endl;
  file << "p_size\t" << p_size  << endl;

  file.close();
}

// Ugly test can cause segfault if file is somehow too big
// Should try to build a cache dictionary and test whether the
// desired crypto_params have been filled
int OptimService::verifyOptimData(set<string> crypto_params_set, const std::string& fenc_path, const std::string& fdec_path)
{
  string line_dec, line_enc;
  unsigned int params_nbr = crypto_params_set.size();
  int return_value = 1;
  std::ifstream	fdec(fenc_path);
  std::ifstream fenc(fdec_path);

  if (!fenc.is_open() || !fdec.is_open()) 
  {
    return 0;
  }
  if(getNumberOfLines(fdec) < params_nbr || getNumberOfLines(fenc) < params_nbr)
  {
    fdec.close();
    fenc.close();
    std::ofstream fdec_in(fdec_path);
    std::ofstream fenc_in(fenc_path);
    fdec_in.close();
    fdec_in.close();
    return_value = 0;
  }
  fdec.close();
  fenc.close();
  return return_value;
}

unsigned int OptimService::getNumberOfLines(std::ifstream& f)
{
  unsigned int i = 0;
  string line;

  while (std::getline(f, line)) 
  {  
    if(!(line.c_str()[0] == '#')) i++;//jump over commented line
  }
  return i;
}

// Returns true if optimization file does not exist or is outdated
bool OptimService::fileOutdated(std::string crypto_name, std::string extension)
{
  map<string, double> cache;
  std::string file_path(OptimService::folderName + OptimService::fileName + crypto_name 
      + extension); 

  // Try to open and read the file
  // If it fails suppose that it is because the file does not exist
  if(readOptimData(cache, file_path)) 
  {
    std::cout << "OptimService: Could not access cache file" << std::endl;
    return true;
  }

  // Get a set with all the crypto parameters of the requested cryptosystem
  CryptographicSystem* crypto_ptr = HomomorphicCryptoFactory_internal::getCrypto(crypto_name);
  std::set<std::string> crypto_params_set;
  crypto_ptr->getAllCryptoParams(crypto_params_set);

  // Try to find each crypto_param in the cache and remove it 
  for (auto crypto_param : crypto_params_set)
  {
    // If there is an element missing in the cache file is outdated
    if (cache.erase(crypto_param) == 0) 
    {
      std::cout << "OptimService: "<< crypto_param << " not found in the cache" << std::endl;
      delete crypto_ptr;
      return true;
    }
  }

  // If some values in the cache do not correspond to a crypto_param file is outdated
  if(!cache.empty()) 
  {
    std::cout << "OptimService: " << extension << " cache has too many entries" << std::endl;
    delete crypto_ptr;
    return true;
  }
  delete crypto_ptr;
  return false;
}

