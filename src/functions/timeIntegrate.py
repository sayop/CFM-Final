import numpy as np
from solutionMethod import *
from boundaryCondition import *
from post import *

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

      # update state vector from Q vector
      updateStateVector(dt)

      # compute primitive variables from state vector elements: update only interior points
      updatePrimitiveVariables(inputDict,imax,jmax)

      # update boundary condition
      updateBC(imax,jmax)

      
      t += dt
      print "|- nIter = %s" % nIter, ", t = %.5e" % t, ", dt = %.5e" % dt, ", Tmax = %.5e" % flowVars.T.max(), ", Tmin = %.5e" % flowVars.T.min(), ", Pmax = %.5e" % flowVars.P.max()

      if (nIter % nIterWrite == 0):
         plotContour(domainVars.x, domainVars.y, flowVars.U, flowVars.V, nIter)

      #if (nIter >= maxIter or resNorm <= residualMin): break
      if (nIter >= maxIter): break

   #
   # time elapsed:
   elapsedTime = (time.clock() - start)
   print "## Elapsed time: ", elapsedTime



