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

#ifndef DEF_PIRCLIENTLOG
#define DEF_PRICLIENTLOG

#include <fstream>

#include "PIRView.hpp"
#include "../common/ToolsBox.hpp"

class PIRClientLog : public PIRView
{
	private:
		std::ofstream file;

		void writeMessage(std::string, MessageEvent& event);
	
	public:
		PIRClientLog(PIRController* controller_);
		~PIRClientLog();

		void messageUpdate(MessageEvent& event);
		void catalogUpdate(CatalogEvent& event);
		void writeUpdate(WriteEvent& 		 event);
};

#endif
