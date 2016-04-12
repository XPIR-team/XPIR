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

#include "MessageEvent.hpp"

MessageEvent::MessageEvent():
	mtype(DEFAULT)
{}

MessageEvent::MessageEvent(message_type mtype_):
	mtype(mtype_)
{}

MessageEvent::MessageEvent(std::string message_):
	mtype(DEFAULT),
	message(message_)
{}

MessageEvent::MessageEvent(message_type mtype_, std::string message_):
	mtype(mtype_),
	message(message_),
	info(message_)
{}

MessageEvent::MessageEvent(message_type mtype_, std::string message_, std::string info_):
	mtype(mtype_),
	message(message_),
	info(info_)
{}

message_type MessageEvent::getMessageType()
{
	return mtype;
}

const std::string& MessageEvent::getMessage()
{
	return message;
}

const std::string& MessageEvent::getInfo()
{
	return info;
}
void MessageEvent::setMessage(std::string message_)
{
	message = message_;
}
