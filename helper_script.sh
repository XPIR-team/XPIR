#!/bin/bash
#/* Copyright (C) 2014 Carlos Aguilar Melchor, Joris Barrier, Marc-Olivier Killijian
# * This file is part of XPIR.
# *
# *  XPIR is free software: you can redistribute it and/or modify
# *	it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation, either version 3 of the License, or
# *  (at your option) any later version.
# *
# *  XPIR is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with XPIR.  If not, see <http://www.gnu.org/licenses/>.
#*/

wget --no-check-certificate https://github.com/XPIR-team/XPIR-dependencies/raw/master/dependencies.tgz
tar zxf dependencies.tgz
rm dependencies.tgz
mkdir local

CONFIGURE="./configure CFLAGS=-I$PWD/local/include LDFLAGS=-L$PWD/local/lib --prefix=$PWD/local/"

  cd dependencies/mpfr-3.1.2 && $CONFIGURE && make && make install
  cd ../..
  cd dependencies/gmp-6.0.0 && $CONFIGURE --enable-cxx && make && make check && make install
  cd ../..
  LOCAL_PATH="$PWD/local/ "
  # Boostrap the build module
  cd dependencies/boost; 
  ./bootstrap.sh
   ./bjam --prefix=$LOCAL_PATH install 
  if [ `uname` = "Darwin" ]
    then 
	cd ../../local/lib/
	install_name_tool -change libboost_system.dylib `pwd`/libboost_system.dylib libboost_thread.dylib
  fi
  cd ../..
