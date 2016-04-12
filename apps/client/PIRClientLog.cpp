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

#include "PIRClientLog.hpp"

PIRClientLog::PIRClientLog(PIRController* controller_) :
	PIRView(controller_),
	file("clientError.log", std::ios::app)
{}

PIRClientLog::~PIRClientLog()
{
	file.close();
}

const std::string PIRClientLog::getCurrentTime()
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];

    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

/**
 *	Called by the model when a information as to be displayed.
 *	Param :
 *		- MessageEvent& event : the message to write. 
 **/
void PIRClientLog::messageUpdate(MessageEvent& event)
{
	std::string m;
	switch(event.getMessageType())
	{
		case WARNING :
			{
				m = "Warning";
				break;
			}
		case ERROR :
			{
				m = "Message";
				break;
			}
		case CHOICE:
			{
				m = "Choice ";
				break;
			}
		case DEFAULT:
			{
				m = "Default";
				break;
			}
		case RETRY:
			{
				m = "Retry";
				break;
			}
	}
	writeMessage(m, event);
}
/**
 *	Called by the model when catalog is entierely received.
 *	Param :
 *		- CatalogEvent& event : new catalog event to write.
 **/
void PIRClientLog::catalogUpdate(CatalogEvent& event)
{
	file << "Proposed files : " << getCurrentTime() << std::endl;
	file << "---------------  " << std::endl;
	const	std::vector<std::string>& fileList = event.getCatalog();

	for (unsigned int i = 0 ; i < fileList.size() ; i++)
	{
		file << fileList.at(i) << " ";
		if (i+1 % 10 == 0)
			file << std::endl;
	}
	file << std::endl;
	file << "---------------  " << std::endl;
	std::flush(file);
}
void PIRClientLog::writeUpdate(WriteEvent& event)
{
	event.getSizeToWrite();
}

void PIRClientLog::writeMessage(std::string type, MessageEvent& event)
{
	file << getCurrentTime() << " " << type << " : " << event.getMessage() << " : " << event.getInfo() << std::endl;	
	std::flush(file);
}

