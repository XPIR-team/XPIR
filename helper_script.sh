#/* Copyright (C) 2014 Carlos Aguilar Melchor, Joris Barrier, Marc-Olivier Killijian
# * This file is part of XPIRe.
# *
# *  XPIRe is free software: you can redistribute it and/or modify
# *	it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation, either version 3 of the License, or
# *  (at your option) any later version.
# *
# *  XPIRe is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with XPIRe.  If not, see <http://www.gnu.org/licenses/>.
#*/

CONFIGURE="./configure CFLAGS=-I$PWD/local/include LDFLAGS=-L$PWD/local/lib --prefix=$PWD/local/"

if [ "$1" = "mpfr" ] && [ -d dependencies/mpfr-3.1.2 ] 
then
  cd dependencies/mpfr-3.1.2 && $CONFIGURE && make && make install
  exit $?
fi

if [ $1 = "gmp" ] && [ -d dependencies/gmp-6.0.0 ]
then
  cd dependencies/gmp-6.0.0 && $CONFIGURE && make && make check && make install
  exit $?
fi
                            
if [ "$1" = "mpfrclean" ] && [ -d dependencies/mpfr-3.1.2 ] 
then
  cd dependencies/mpfr-3.1.2 && $CONFIGURE && make distclean
  exit $?
fi

if [ $1 = "gmpclean" ] && [ -d dependencies/gmp-6.0.0 ]
then
  cd dependencies/gmp-6.0.0 && $CONFIGURE && make distclean
  exit $?
fi
                            
if [ $1 = "boost" ] && [ -d dependencies/boost ]
then
	LOCAL_PATH="$PWD/local/ "
  # Boostrap the build module
  cd dependencies/boost; 
	./bootstrap.sh
  # Use the resulting bjam to compile everything 
	./bjam --prefix=$LOCAL_PATH install 
	#fix linking issue on OSX
	if [ `uname` == "Darwin" ]
	then 
		cd ../../local/lib/
		install_name_tool -change libboost_system.dylib `pwd`/libboost_system.dylib libboost_thread.dylib
	fi
  exit $?
fi

