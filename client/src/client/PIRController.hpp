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

#ifndef DEF_PIRCONTROLLER
#define DEF_PIRCONTROLLER

#include <vector>

#include "PIRViewCLI.hpp"
#include "PIRClientLog.hpp"

class PIRClientSimple;

class PIRController
{
	private:
		PIRClientSimple& 	model_ptr;
		PIRViewCLI 				view_cli;
		PIRClientLog			logger;

		std::vector<boost::signals2::connection> connectionMessage, connectionMenu, connectionWrite;
	
	public:
			PIRController(PIRClientSimple& model_ptr_);
			~PIRController();

			void addListenersToModel();
			void notifyClientChoice(char c);
			void notifyClientChoice(int  c);

};

#endif
