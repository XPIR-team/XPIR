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

all : server

optim :
	@make -f genMakefile  MODULES='pir server crypto crypto/prng optim optim_main common' EXEC='PIROptimizer'

optim_debug :
	@make -f genMakefile  MODULES='pir server crypto crypto/prng optim optim_main common' EXEC='PIROptimizer' DEBUG='1'

server :
	@make -f genMakefile  MODULES='crypto crypto/prng optim main_server server pir common' EXEC='PIRServer'

server_debug :
	 @make -f genMakefile  MODULES='main_server server pir crypto crypto/prng common optim' EXEC='PIRServer' DEBUG='1'
database :
	@/bin/bash makedb.sh 1 4

tests :
	@make -f genMakefile  MODULES='server pir optim crypto crypto/prng common tests' EXEC='Tests' LDLIBS='boost_unit_test_framework ntl mpfr gmp boost_system boost_program_options boost_thread boost_signals' DEBUG='1'

clean :
	rm -rf build/*


