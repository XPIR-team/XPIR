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

#include "ToolsBox.hpp"

namespace tools
{
  const std::string getCurrentTime()
  {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];

    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
  }

 // char* ZZToChar(const ZZ& z)
 // {
 //   unsigned long l;
 //   conv(l,z);
 //   char* c = new char[sizeof(long)]();
 //   memcpy(c, &l, sizeof(long));
 //   return c;
 // }

 // void charToZZ_pX(char* c, size_t s, ZZ_pX& z)
 // {
 //   z.SetLength(s);

 //   for (unsigned int i = 0 ; i < s ; i++)
 //   {
 //     conv(z[i], (long)c[i]);
 //   }
 // }
}

