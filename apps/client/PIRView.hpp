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

#ifndef DEF_PIRVIEW
#define DEF_PIRVIEW

#include <vector>
#include <string>

#include <boost/signals2.hpp>

#include "pir/events/MessageEvent.hpp"
#include "pir/events/CatalogEvent.hpp"
#include "pir/events/WriteEvent.hpp"


class PIRController;

class PIRView
{
	private:
		 boost::signals2::connection m_connection, m_connectionMenu;
	
	protected:
		 PIRController* controller;

	public:
		PIRView(PIRController* controller_);
		virtual ~PIRView();

		virtual void messageUpdate(MessageEvent& event)=0;
		virtual void catalogUpdate(CatalogEvent& event)=0;
		virtual void writeUpdate(WriteEvent& 		 event)=0;

};

#endif
