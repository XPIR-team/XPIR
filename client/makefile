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

all : client

client :
	@make -f genMakefile  MODULES='client pir crypto crypto/prng events common optim .' EXEC='PIRClient'

client_debug :
	@make -f genMakefile  MODULES='client pir crypto crypto/prng events common optim .' EXEC='PIRClient' DEBUG='1'

tests :
	@make -f genMakefile  MODULES='crypto crypto/prng events common tests' EXEC='Tests' LDLIBS='boost_unit_test_framework ntl mpfr  gmp boost_system boost_program_options boost_thread boost_signals'

clean :
	rm -rf build/* 

