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

#ifndef DEF_PIRCLIENTCLI
#define DEF_PIRCLIENTCLI

#include <boost/signals2.hpp>
#include <string>
#include <iostream>

#include "PIRView.hpp"

#define RESET_COLOR "\033[m"
#define RED "\033[31m"
#define BOLD "\033[1;30m"
#define ORANGE "\033[33m"
#define DEFAULT_BOX_SIZE 24

class PIRController;

class PIRViewCLI : public PIRView
{
	public:
		PIRViewCLI(PIRController* controller_);
		
		void messageUpdate(MessageEvent& event);
		void catalogUpdate(CatalogEvent& event);
		void writeUpdate(WriteEvent& 		 event);

		void getUserInputRetry();
		void getUserInputFile(int maxValue);
};

#endif
