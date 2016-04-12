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

#ifndef DEF_MESSAGEEVENT
#define DEF_MESSAGEEVENT

#include <string>
#include "pir/ClientDefines.hpp"

/**
 * Used to send a message (like error message, info message etc )
 * of type message_type to the view.
 **/

class MessageEvent {
		private :
			message_type mtype;
			std::string message, info;

		public :
			MessageEvent();
			MessageEvent(message_type m_type);
			MessageEvent(std::string message_);	
			MessageEvent(message_type mtype_, std::string message_);
			MessageEvent(message_type mtype_, std::string message_, std::string info);	

			message_type getMessageType();
			const std::string& getMessage();
			const std::string& getInfo();
			void  setMessage(std::string message_);
};

#endif
