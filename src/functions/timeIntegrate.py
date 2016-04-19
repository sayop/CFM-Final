import numpy as np
from solutionMethod import *

import time

def timeIntegrate(inputDict):
   tStart  = 0.0
   imax    = int(inputDict['iDim'])
   jmax    = int(inputDict['jDim'])
   maxIter = int(inputDict['maxIter'])
   nIterWrite  = int(inputDict['nIterWrite'])
   iVisc   = int(inputDict['iVisc'])

   # start to count time for calculting computation performance
   start = time.clock()

   #
   # Time Marching:
   #
   print '=============================================='
   print '# Time integration starts at t = %s' % tStart
   print '=============================================='
   t = tStart
   nIter = 0
   while True:
      nIter += 1
      # Find time step that may stabilize the solution with given Courant number
      dt = computeTimeStep(inputDict,imax,jmax)

      # update flux vector for inviscid and viscid terms
      updateFluxVectors(iVisc)
      
      # update Q vector: Q vector contains finite differenced state vectors to update state vector in time
      updateQvector(imax,jmax)

      t += dt
      print "|- nIter = %s" % nIter, ", t = %.5e" % t, ", dt = %.5e" % dt

      #if (nIter >= maxIter or resNorm <= residualMin): break
      if (nIter >= maxIter): break

   #
   # time elapsed:
   elapsedTime = (time.clock() - start)
   print "## Elapsed time: ", elapsedTime



