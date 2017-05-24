# -*- coding: utf-8 -*-

from estimator import *  # Import lwe estimator of Martin Albrecht

from sage.rings.real_mpfr import RRtoRR  # Import RealField

import numpy as np
import os
import re
import datetime


print "Estimate the complexity of solving LWE with XPIR parameters\n"


# Precision of 100 (high)
RR = RealField(100)


# Get the current path
pathScript = os.getcwd()
# Initialize the path of the data file
path = pathScript + "/../_build/apps/server/exp/preComputeLWE.abs"


# Chech that the data file exists
if (os.path.isfile(path)):
    
    print "Please wait...\n"
    
    
    # Open the data file which contains the XPIR parameters
    with open(path, "r") as cryptoParamsFile: 
        # Read the file
        cryptoParams = cryptoParamsFile.read()
        # Put the paramters in a list
        cryptoParamsList = filter(None, re.split("[\n :]+", cryptoParams))
        cryptoParamsFile.close()
    
    # Initialize n and q lists
    n = []
    q = []

    # Add the parameters of the data file in the lists
    for i in range (2, len(cryptoParamsList), 6):
    
        n.append(int(cryptoParamsList[i]))
        q.append(int(cryptoParamsList[i + 1]))

    print "XPIR parameters loaded\n"


    # Initialize the list of numbers of bits of security
    nbrBits = []

    # Compute the number of bits for each parameters with the Martin Albrecht algortihm
    for i in range (len(n)):
    
        security = estimate_lwe(n[i], RR(80 / RR((2 ** q[i])) ) , 2 ** q[i], skip=("mitm", "bkw", "arora-gb"))
    
        # Add the result to the list
        nbrBits.append(int(np.log2(min(security['sis']['bkz2'], security['dec']['bkz2'], security['kannan']['bkz2']))))
    
        print "estimate parameters ", i + 1, " : done"
    

    # Open 'security_estimations.txt' file, if it does not exitst, it creates it
    paramsSecure = open("security_estimations.txt", 'w')

    # Write in the file the current date & hour
    paramsSecure.write(str(datetime.datetime.now()))

    paramsSecure.write("\n\n")
    paramsSecure.write("n:q:nbrBits\n")

    # Write the parameters and the numbers of bits of security in the file
    for i in range (len(n)):
    
        paramsSecure.write(str(n[i]))
        paramsSecure.write(":")
        paramsSecure.write(str(q[i]))
        paramsSecure.write(":")
        paramsSecure.write(str(nbrBits[i]))
        paramsSecure.write("\n")
 
    # Close the file 
    paramsSecure.close()

    print "\nResults of the estimation written\n\nScript finished !"
    
    
else:
    print "ERROR : data file that contains XPIR parameters does not find... Please run XPIR for the first time before !"


