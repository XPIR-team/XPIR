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


########### CONSTANTS ##############
MAX_REC=3
MIN_ALPHA=0		    # Set MIN and MAX to 1 to forbid aggregation
MAX_ALPHA=0		   
NO_REREAD=1
NO_PIPELINE=0
VERBOSE=0

TEST_PAILLIER=0
TEST_NOCRYPTOGRAPHY=1
TEST_LWE=1

REMOTE=0
IP=169.254.4.46
USER=marco

ONE_KBIT=1024
HUNDRED_KBIT=102400
ONE_MBIT=1024000
TEN_MBIT=10240000
HUNDRED_MBIT=102400000
ONE_GBIT=1024000000

#files: 1kbits, 100kbits, 10mbits 1gbit
#bases: 1Mbits, 10M, 100M, 1G 10G

########### SUBROUTINES ##############
deal_with_options() 
{
    if [[ ( $NO_REREAD == 1 ) ]];
    then
	S_OPTION="-z"
    fi
    if [[ ( $NO_PIPELINE == 1 ) ]];
    then
	S_OPTION=$S_OPTION" --no-pipeline"
	C_OPTION="--no-pipeline"
	echo "Mode --no-pipeline selected"
    fi
}

do_a_test() 
{
    rm -f reception/* 2> /dev/null
	if [[ $VERBOSE == 2 ]]; then
		echo $BASE_DIR/../server/pir_server $S_OPTION
	    echo $BASE_DIR/../client/pir_client -r $PARAM $C_OPTION $@ -c 
	fi
  if [[ $N -le 1000 ]];
  then 
    $BASE_DIR/../server/pir_server $S_OPTION > /tmp/checkpirserver.stdout 2>/tmp/checkpirserver.stderr &
  else
    $BASE_DIR/../server/pir_server $S_OPTION -s $N > /tmp/checkpirserver.stdout 2>/tmp/checkpirserver.stderr &
  fi
  PID=$!
  sleep 1
  $BASE_DIR/../client/pir_client -r $PARAM $C_OPTION $@ -c  > /tmp/checkpirclient.stdout 2> /tmp/checkpirclient.stderr 
}

exploit_results()
{
    FILE_RETRIEVED=`ls reception`
    MD5_R1=`sha1sum reception/$FILE_RETRIEVED 2>/dev/null |cut -d\  -f1`
    MD5_DB=`dd if=db/test1 bs=1 count=$L_BYTE 2>/dev/null |sha1sum 2>/dev/null |cut -d\  -f1`
    
    if [[ ( $NO_REREAD -eq 1 ) ]];
    then
						#	to check no-reread-database, do it a second time
	rm -f reception/* 2> /dev/null
	$BASE_DIR/../client/pir_client -r $PARAM $C_OPTION $@ -c  >> /tmp/checkpirclient.stdout 2>> /tmp/checkpirclient.stderr 
	MD5_R2=`sha1sum reception/* 2>/dev/null |cut -d\  -f1`
	if [[ $FILE_RETRIEVED != "" && ($MD5_DB == $MD5_R1) && ($MD5_DB == $MD5_R2) ]]; then
	    CORRECT=1;
	else
	    CORRECT=0;
	fi
    else
	if [[ $FILE_RETRIEVED != "" && ($MD5_DB == $MD5_R1) ]]; then
	    CORRECT=1;
	else
	    CORRECT=0;
	fi
    fi

    if [[ $CORRECT == 1 ]]; then
	    echo -e "$DB:$L:$PARAM \033[32mCORRECT\033[m"
    else
	    echo -e "$DB:$L:$PARAM \033[31m*************** NOT CORRECT **********\033[m"
	    if [[ $VERBOSE -ge 1 ]]; then
	      echo "Database : check.repo/db-$L_BYTE-$N"
	      echo "Server : $BASE_DIR/../server/pir_server $S_OPTION"
	      echo "Client : $BASE_DIR/../client/pir_client -r $PARAM $C_OPTION $@ -c  "
	      echo "*************** Server stdout **********"
        cat /tmp/checkpirserver.stdout
	      echo "*************** Server stderr **********"
        cat /tmp/checkpirserver.stderr
	      echo "*************** Client stdout **********"
        cat /tmp/checkpirclient.stdout
	      echo "*************** Client stderr **********"
        cat /tmp/checkpirclient.stderr
	      echo "hit <enter> to continue";read
	    fi    
    fi
    
    (kill $PID >/dev/null 2>/dev/null)
    # Notify when waiting for kill
    # Use ANSI escape sequences to stay on the same line
    while [[ `ps -ef|grep pir_server|wc -l` -ne 1 ]] ; do 
      echo -e "Could not kill pir_server, waiting ..."
      echo -e "\033[2A"
      sleep 1
    done
    # Use ANSI escape sequences again to erase and reuse the line
    echo  "                                     "
    echo -e "\033[2A"
    rm -f /tmp/checkpir* >/dev/null 2>/dev/null
}


########### MAIN ##############

echo -e "##########################################################################"
echo -e "This tool tests that pir_server and pir_client run correctly and that an"
echo -e "element can be retrieved without errors. You should obtain CORRECT or "
echo -e "\"Skipping test...\" for all tests. THE FIRST TEST CAN BE QUITE LONG if"
echo -e "performance caches need to be built (first run for the server or client)"
echo -e "##########################################################################"

killall -9 pir_server >/dev/null 2>/dev/null; sleep 1

# Notify when waiting for kill
# Use ANSI escape sequences to stay on the same line
while [[ `ps -ef|grep pir_server|wc -l` -ne 1 ]]; do 
  echo -e "Could not kill pir_server, waiting ..."
  echo -e "\033[2A"
  sleep 1
done
# Use ANSI escape sequences again to erase and reuse the line
echo  "                                     "
echo -e "\033[2A"

deal_with_options

BASE_DIR=$PWD

cd check.repo


# Paillier tests only for small databases
if [[ TEST_PAILLIER -eq 1 ]]; then
  echo -e "\nPaillier tests\n#################\n"
  for DB in $ONE_MBIT #$TEN_MBIT 
    do
	for L in $ONE_KBIT $HUNDRED_KBIT $TEN_MBIT $ONE_GBIT 
	do
		
    N=`python -c"print(int($DB / $L));"`
    L_BYTE=`python -c"print(int($L / 8 ));"`
	    
	    if [[ ( $DB -gt $L ) && ( -f  db-$L_BYTE-$N/test1 ) ]]; 
	    then
		
		rm -fr db
		mkdir reception 2> /dev/null
		mkdir exp 2> /dev/null
		
		
		ln -s db-$L_BYTE-$N db
	  echo Checking db-$L_BYTE-$N 

		for QP in "80:1024:2048"
		#for QP in "80:1024:2048:1016"
		do
		    for REC in `eval echo {1..$MAX_REC}`
		    do	
		    # TODO use alpha (aggregation does not work yet)
			for ALPHA in `eval echo {$MIN_ALPHA..$MAX_ALPHA}`
			do	
			    PARAM="Paillier:$QP --reclvl $REC --alpha $ALPHA "
			    
			    do_a_test
			    exploit_results
			done
		    done
		done
	    fi
	done
    done
fi

echo -e "\n\nTests\n#################\n"
for DB in $ONE_MBIT $TEN_MBIT $HUNDRED_MBIT $ONE_GBIT 
do
    for L in $ONE_KBIT $HUNDRED_KBIT $TEN_MBIT $ONE_GBIT  
    do
      N=`python -c"print(int($DB / $L));"`
      L_BYTE=`python -c"print(int($L / 8) );"`
      DB_BYTE=`python -c"print(int($DB / 8));"`

  # If N <= 10000 a database with different files must exist 
  # If not a database with a single file to split must exist
  if [[ ( ( $N -le 10000 ) && ( -f  db-$L_BYTE-$N/test1 ) ) || ( $N -gt 1000 ) && ( -f db-$DB_BYTE/test1 ) ]]; 
	then
	    
	    rm -fr db
	    mkdir reception 2> /dev/null
	    mkdir exp 2> /dev/null
	    
	    if [[ $N -le 1000 ]]; 
      then
	      ln -s db-$L_BYTE-$N db
	      echo Checking db-$L_BYTE-$N 
      else
	      ln -s db-$DB_BYTE db
	      echo Checking db-$DB_BYTE with split_value=$N 
      fi

	    if [[ TEST_NOCRYPTOGRAPHY -eq 1 ]]; then
	        # Test No Cryptography
		PARAM="NoCryptography"
		do_a_test
		exploit_results
	    fi

	    # Test LWE
	    for QP in "4096:180" "2048:120" "1024:60"
	    #for QP in "180:73" "120:43" "60:13"
	    do
		#TODO use python to compute absorption #math.floor(($Q-math.ceil(math.log($SEC/2,2))-math.ceil(math.log($N,2))-math.ceil(math.log($DEG,2)))/2)
		for REC in `eval echo {1..$MAX_REC}`
		do	
		    for ALPHA in `eval echo {$MIN_ALPHA..$MAX_ALPHA}`
		    do	
			if [[ TEST_LWE -eq 1 ]]; then
			  PARAM="LWE:.*:$QP --reclvl $REC --alpha $ALPHA "
        CIPH_SIZE=`echo $QP|tr : \*`"*2"
        QUERY_SIZE="$CIPH_SIZE*$REC*$N.0**(1/$REC.0)"
        DB_FFT_SIZE="6*$DB"
        if [[ ( `python -c "print(($QUERY_SIZE+$DB_FFT_SIZE)/10**9 > 40);"` == "True" ) || ( `python -c "print(($CIPH_SIZE*$N+$QUERY_SIZE)/10**9 > 40);"` == "True" ) ]]; then
          echo "Skipping tests requiring more than 5Gbytes RAM"
        else
			    do_a_test
			    exploit_results
        fi
			fi
		    done
		done
	    done
	fi
    done
done
cd ..


