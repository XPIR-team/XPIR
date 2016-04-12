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

#include "PIRReplyGeneratorFactory.hpp"

GenericPIRReplyGenerator* PIRReplyGeneratorFactory::getPIRReplyGenerator(const std::string& gen_name, PIRParameters& param, DBHandler *db)
{
  GenericPIRReplyGenerator* generator;
  
  if ((gen_name == "LWE"))
	{
		generator = new PIRReplyGeneratorNFL_internal( param, db); 
	}
	else if (gen_name == "Paillier")
	{
    generator = new PIRReplyGeneratorGMP( param, db);
	}
  else if (gen_name == "NoCryptography")
  {
    generator = new PIRReplyGeneratorTrivial( param, db);
  }
  else
  {
    std::cout << "PIRReplyGeneratorFactory: Could not find an appropriate generator for "<< gen_name << ", returning NULL" << std::endl;
    generator = NULL;
  }

  return generator;
}

// Function dedicated to the optimizer. Should not be called with NoCryptography
GenericPIRReplyGenerator* PIRReplyGeneratorFactory::getPIRReplyGenerator(const std::string& gen_name)
{
  GenericPIRReplyGenerator* generator;
 
  if ((gen_name == "LWE"))
	{
	  generator = new PIRReplyGeneratorNFL_internal(); 
	}
	else if (gen_name == "Paillier")
	{
    generator = new PIRReplyGeneratorGMP();
	}
  else
  {
    std::cout << "PIRReplyGeneratorFactory (optim function): Could not find an appropriate generator for "<< gen_name << ", returning NULL" << std::endl;
    generator = NULL;
  }

  return generator;
}
