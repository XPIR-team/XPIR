# -*- coding: utf-8 -*-


from sage.rings.real_mpfr import RRtoRR  # Import RealField from SageMath


import numpy as np
import os
import re
import datetime


import sys

# Get the current path
pathScript = os.getcwd()
# Initialize the path of Martin Albrecht script
pathModule = pathScript + "/lwe-estimator"
sys.path.append(pathModule)
# Import lwe estimator of Martin Albrecht
from estimator import *  


print "Estimate the complexity of solving LWE with XPIR parameters\n"


# Precision of 100 (high)
RR = RealField(100)


# Initialize the path of the NFLParams.cpp file
pathNFLParameters = pathScript + "/NFLParams.cpp"

# Chech that the file exits
if (os.path.isfile(pathNFLParameters)):
    
    
    print "Please wait...\n"
    
    
    # Open the data file which contains NFL parameters
    with open(pathNFLParameters) as paramsFile:
        
        # Check all lines
        for line in paramsFile:

            # Check the line that contains kMinPolyDegree
            if 'const unsigned int kMinPolyDegree' in line:
                
                # Find the index of caracters before and after kMinPolyDegree
                index1 = line.find('=')
                index2 = line.find('\n', index1+1)
                
                # Set kMinPolyDegree
                kMinPolyDegree = int(line[index1 + 2 : index2 - 1])
            
            # Check the line that contains kMaxPolyDegree
            if 'const unsigned int kMaxPolyDegree' in line:

                # Find the index of caracters before and after kMaxPolyDegree
                index1 = line.find('=')
                index2 = line.find('\n', index1+1)
                
                # Set kMaxPolyDegree
                kMaxPolyDegree = int(line[index1 + 2 : index2 - 1])
                
            # Check the line that contains kMaxAggregatedModulusBitsize    
            if 'const unsigned int kModulusBitsize' in line:

                # Find the index of caracters before and after kModulusBitsize
                index1 = line.find('=')
                index2 = line.find('\n', index1+1)
                
                # Set kModulusBitsize
                kModulusBitsize = int(line[index1 + 2 : index2 - 1])
                
            # Check the line that contains kMaxAggregatedModulusBitsize
            if 'const unsigned int kMaxAggregatedModulusBitsize' in line:

                # Find the index of caracters before and after kMaxAggregatedModulusBitsize
                index1 = line.find('=')
                index2 = line.find('\n', index1+1)
                
                # Set kMaxAggregatedModulusBitsize
                kMaxAggregatedModulusBitsize = int(line[index1 + 2 : index2 - 1])
 
    
    # Initialize the path of the NFLLWESecurityEstimated.hpp file
    pathNFLLWESecurityEstimatedHPP = pathScript + "/../NFLLWESecurityEstimated.hpp"
    
    # Open NFLLWESecurityEstimated.hpp, if it does not exist, it will create it
    paramsSecure = open(pathNFLLWESecurityEstimatedHPP, 'w')
    paramsSecure.write('#pragma once\n')
    paramsSecure.write("#include <string>\n")
    paramsSecure.write('\n')
    paramsSecure.write("using namespace std;\n")
    paramsSecure.write('\n')
    paramsSecure.write('string securityParameters = "') 
    
    
    # Initialize the number of estimations
    i =0
    
    # Scan n from kMinPolyDegree to kMaxPolyDegree
    for log2n in range(int(np.log2(kMinPolyDegree)), int(np.log2(kMaxPolyDegree)) + 1, 1):
        n = 2 ** log2n
        
        # Scan log2q from kModulusBitsize to kMaxAggregatedModulusBitsize
        for log2q in range(kModulusBitsize, kMaxAggregatedModulusBitsize + 1, 60):
            
            # Increment the number of estimations
            i += 1
            
            # Compute the number of bits for each parameters with the Martin Albrecht algortihm
            security = estimate_lwe(n, RR(80 / RR((2 ** log2q)) ) , 2 ** log2q, skip=("mitm", "bkw", "arora-gb"))
            # Select the security and return the number of bits
            nbrBits = int(np.log2(min(security['sis']['bkz2'], security['dec']['bkz2'], security['kannan']['bkz2'])))
            
            # Write security parameters
            paramsSecure.write(str(n))
            paramsSecure.write(":")
            paramsSecure.write(str(log2q))
            paramsSecure.write(":")
            paramsSecure.write(str(nbrBits))
            paramsSecure.write('\\n')
            paramsSecure.write('\\')
            paramsSecure.write('\n')
            
            # Print to the user the number of the last estimation done
            print "estimate parameters ", i, " : done"
            
            
    paramsSecure.write('";\n')
    # Close the NFLLWESecurityEstimated.hpp file
    paramsSecure.close()
    
    
    print "\nResults of the estimation written\n\nScript finished !"
            
            
else:
    # Error if XPIR file not found
    print "ERROR : XPIR files not found !"
