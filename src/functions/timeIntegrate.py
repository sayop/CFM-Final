import numpy as np
from solutionMethod import *

import time

def timeIntegrate(inputDict):
   Cr      = float(inputDict['Courant'])
   imax    = int(inputDict['iDim'])
   jmax    = int(inputDict['jDim'])
   gamma   = float(inputDict['gamma'])
   Rgas    = float(inputDict['gasConst'])

   # start to count time for calculting computation performance
   start = time.clock()

   a = computeSpeedOfSound(flowVars.T,gamma,Rgas)

   #
   # time elapsed:
   elapsedTime = (time.clock() - start)
   print "## Elapsed time: ", elapsedTime



