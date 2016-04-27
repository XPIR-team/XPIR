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


ONE_KBIT=1024
HUNDRED_KBIT=102400
ONE_MBIT=1024000
TEN_MBIT=10240000
HUNDRED_MBIT=102400000
ONE_GBIT=1024000000

#files: 1kbits, 100kbits, 10mbits 1gbit
#bases: 1Mbits, 10M, 100M, 1G, 10G

if [[ -e check.repo ]]
then
  echo "Check repo exists, not rebuilding it"
  echo "Remove it manually if you want it rebuilt"
  exit 
fi

mkdir check.repo
cp -r ../client/exp check.repo/
cd check.repo
mkdir reception

if [[ -f source.random ]]
then
	echo "Reusing the source.random file"
else
	echo "Creating the source.random file"
  dd if=/dev/urandom of=source.random count=`python -c"print( int($ONE_GBIT / (1024000) ) );"` bs=128000
fi

echo "Creating the databases .."
for DB in $ONE_MBIT $TEN_MBIT $HUNDRED_MBIT $ONE_GBIT
do
  for L in $HUNDRED_KBIT $TEN_MBIT $ONE_GBIT $ONE_KBIT 
	do
    N=`python -c"print( int($DB / $L) );"`
		
    # Only build directories with 1000 or less files
		if [[ ( $DB -gt $L ) && ( $N -le 1000 )]]; 
			then
			
			#First we need to obtain the appropriate parameters for N and L fixed
	
      L_KBIT=`python -c"print( int($L / 1024) );"`
      L_BYTE=`python -c"print( int($L / 8) );"`
				
			rm -fr db
			mkdir db
			
			dd if=source.random of=db/test1 count=$L_KBIT bs=128
			
      cd db
			for ((  i = 2 ;  i <= $N ; i++  ))
				do
					ln -s test1 test$i
				done
      cd ..

			mv db db-$L_BYTE-$N
    fi      
	done
  
  # For tests with more than 1000 files use a single file with split_file option
  rm -fr db
	mkdir db
  dd if=source.random of=db/test1 count=`python -c"print( int($DB / 1024000) );"` bs=128000
  mv db db-`python -c"print( int($DB / 8) );"`
done
rm -f source.random
cd ..
