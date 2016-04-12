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

#include "PIRController.hpp"
#include "PIRClient.hpp"

/**
 *	Class constructor.
 *	Param :
 *		- pirClient_ptr model_ptr_ : PIRClient shared_ptr.
 **/
PIRController::PIRController(PIRClientSimple& model_ptr_):
  model_ptr(model_ptr_),
  view_cli(this),
  logger(this)
{
  addListenersToModel();
}

/**
 *	Add listners (attribute of PIRController) to the model (model_ptr).
 **/
void PIRController::addListenersToModel()
{
  connectionMessage.push_back(model_ptr.addMessageListener(boost::bind(&PIRView::messageUpdate, &view_cli, _1)));
  connectionMessage.push_back(model_ptr.addMessageListener(boost::bind(&PIRView::messageUpdate, &logger, _1)));

  connectionMenu.push_back(model_ptr.addMenuListener(boost::bind(&PIRView::catalogUpdate, &view_cli, _1))); 
  connectionMenu.push_back(model_ptr.addMenuListener(boost::bind(&PIRView::catalogUpdate, &logger, _1)));

  connectionWrite.push_back(model_ptr.addWriteListener(boost::bind(&PIRView::writeUpdate, &view_cli, _1))); 
}


/**
 *	Notify model of a client choice.
 *	Param :
 *		- int c : int input.
 **/
void PIRController::notifyClientChoice(int c)
{
  model_ptr.setChosenElement(c);
}

PIRController::~PIRController()
{
  for (auto &messageView : connectionMessage)
    messageView.disconnect();

  for (auto &menuView : connectionMenu)
    menuView.disconnect();

  for (auto &writeView : connectionWrite)
    writeView.disconnect();
}
