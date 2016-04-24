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
   nonDim  = int(inputDict['nonDim'])

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
      # plot initial condition
      if nIter == 0: 
         if nonDim == 1: dimensionalize(1, 1)
         plotUmagContour(domainVars.x, domainVars.y, flowVars.U, flowVars.V, nIter)
         if nonDim == 1: nondimensionalize(1, 1)

      nIter += 1
      # Find time step that may stabilize the solution with given Courant number
      dt = computeTimeStep(inputDict,imax,jmax)

      # update flux vector for inviscid and viscid terms
      updateFluxVectors(inputDict,imax,jmax,iVisc)

      # update Q vector: Q vector contains finite differenced state vectors to update state vector in time
      updateQvector(imax,jmax,iVisc)

      # update state vector from Q vector
      integrateStateVector(dt)

      # compute primitive variables from state vector elements: update only interior points
      updatePrimitiveVariables(inputDict,imax,jmax)

      # update boundary condition
      if nonDim == 1: dimensionalize(0, 1)
      updateBC(inputDict,imax,jmax)
      if nonDim == 1: nondimensionalize(0, 1)
      
      t += dt
      print "|- nIter = %s" % nIter, ", t = %.5e" % t, ", dt = %.5e" % dt, ", Tmax = %.5e" % flowVars.T.max(), ", Tmin = %.5e" % flowVars.T.min(), ", Pmax = %.5e" % flowVars.P.max(), ", Umax = %.5e" % flowVars.U.max(), ", Vmax = %.5e" % flowVars.V.max()

      if (nIter % nIterWrite == 0):
         if nonDim == 1: dimensionalize(1, 1)
         plotStreamLine(domainVars.x, domainVars.y, flowVars.U, flowVars.V, nIter)
         plotUmagContour(domainVars.x, domainVars.y, flowVars.U, flowVars.V, nIter)
         #field = 'P'
         #plotContour(domainVars.x, domainVars.y, flowVars.P, nIter, field)
         field = 'T'
         plotContour(domainVars.x, domainVars.y, flowVars.T, nIter, field)
         if nonDim == 1: nondimensionalize(1, 1)

      #if (nIter >= maxIter or resNorm <= residualMin): break
      if (nIter >= maxIter): break

   #
   # time elapsed:
   elapsedTime = (time.clock() - start)
   print "## Elapsed time: ", elapsedTime



