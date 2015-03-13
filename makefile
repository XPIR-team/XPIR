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

all: boost gmp mpfr database per_user.mk client server

boost:
	@echo "#############################################"
	@echo "Building boost..."
	@echo "#############################################"
	@bash helper_script.sh boost
	@echo "\n\n"
gmp:
	@echo "#############################################"
	@echo "Building gmp..."
	@echo "#############################################"
	@bash helper_script.sh gmp
	@echo "\n\n"
mpfr:
	@echo "#############################################"
	@echo "Building mpfr..."
	@echo "#############################################"
	@bash helper_script.sh mpfr
	@echo "\n\n"
database:
	@echo "#############################################"
	@echo "Building server databases..."
	@echo "#############################################"
	@if [ -d server/check.repo ]; then echo "Directory server/check.repo exists, skipping this step..."; else cd server; ./mkdb-correctness.sh; ln -s check.repo/db-1280000-10 db; fi
	@echo "\n\n"
per_user.mk:
	@echo "#############################################"
	@echo "Ensuring per_user.mk exists..."
	@echo "#############################################"
	@if [ ! -f per_user.mk ]; then touch per_user.mk; fi 
	@echo "\n\n"
client:
	@echo "#############################################"
	@echo "Building client..."
	@echo "#############################################"
	@cd client; make
	@echo "\n\n"
server:
	@echo "#############################################"
	@echo "Building server..."
	@echo "#############################################"
	@cd server; make
	@echo "\n\n"

.PHONY: clean boost client server   

clean:
	@if [ -d dependencies/boost ]; then cd dependencies/boost; ./bootstrap.sh; ./bjam --clean; cd ../..; fi
	@bash helper_script.sh gmpclean
	@bash helper_script.sh mpfrclean
	@cd client; make clean; cd ..
	@cd server; make clean; cd ..
