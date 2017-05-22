# -*- coding: utf-8 -*-

from estimator import * # import lwe estimator of Martin Albrecht

from sage.rings.real_mpfr import RRtoRR

import numpy as np
import os
import re
import datetime


print "Estimate the complexity of solving LWE with XPIR parameters\n\nPlease wait...\n"


RR = RealField(100)


pathScript = os.getcwd()
path = pathScript + "/../_build/apps/server/exp/preComputeLWE.abs"

with open(path, "r") as cryptoParamsFile: 
     
    cryptoParams = cryptoParamsFile.read()
    cryptoParamsList = filter(None, re.split("[\n :]+", cryptoParams))
    cryptoParamsFile.close()
    

n = []
q = []

for i in range (2, len(cryptoParamsList), 6):
    
    n.append(int(cryptoParamsList[i]))
    q.append(int(cryptoParamsList[i + 1]))

print "XPIR parameters loaded\n"


nbrBits = []

for i in range (len(n)):
    
    security = estimate_lwe(n[i], RR(80 / RR((2 ** q[i])) ) , 2 ** q[i], skip=("mitm", "bkw", "arora-gb"))
    
    nbrBits.append(int(np.log2(min(security['sis']['bkz2'], security['dec']['bkz2'], security['kannan']['bkz2']))) + 1)
    
    print "estimate parameters ", i + 1, " : done"
    

paramsSecure = open("security_estimations.txt", 'w')

paramsSecure.write(str(datetime.datetime.now()))

paramsSecure.write("\n\n")
paramsSecure.write("n:q:nbrBits\n")

for i in range (len(n)):
    
    paramsSecure.write(str(n[i]))
    paramsSecure.write(":")
    paramsSecure.write(str(q[i]))
    paramsSecure.write(":")
    paramsSecure.write(str(nbrBits[i]))
    paramsSecure.write("\n")
        
paramsSecure.close()

print "\nResults of the estimation written\n\nScript finished !"
   
    


