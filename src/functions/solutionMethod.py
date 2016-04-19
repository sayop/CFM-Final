from variables import *

import numpy as np

def computeTimeStep(inputDict,imax,jmax):
   Cr      = float(inputDict['Courant'])
   gamma   = float(inputDict['gamma'])
   Rgas    = float(inputDict['gasConst'])
   dx = domainVars.dx
   dy = domainVars.dy

   # calculate speed of sound
   a = computeSpeedOfSound(flowVars.T,gamma,Rgas)

   Uconv = abs(flowVars.U) + a
   Vconv = abs(flowVars.V) + a
   dtAll = Cr / (Uconv / dx + Vconv / dy)

   # minimum dt should be found from the interior points
   dtMin = dtAll[1:imax-1,1:jmax-1].min()

   return dtMin

def computeSpeedOfSound(temp,gamma,Rgas):

   a = (gamma * Rgas * temp) ** (0.5)

   return a


def updateFluxVectors(iVisc):

   #State vector
   FDM.phi[0] = flowVars.rho
   FDM.phi[1] = flowVars.rho * flowVars.U
   FDM.phi[2] = flowVars.rho * flowVars.V
   FDM.phi[3] = flowVars.rho * flowVars.et

   # flux vector in x-direction
   FDM.Fi[0] = flowVars.rho * flowVars.U
   FDM.Fi[1] = flowVars.rho * flowVars.U * flowVars.U + flowVars.P
   FDM.Fi[2] = flowVars.rho * flowVars.U * flowVars.V
   FDM.Fi[3] = flowVars.rho * flowVars.U * flowVars.et

   # flux vector in y-direction
   FDM.Gi[0] = flowVars.rho * flowVars.V
   FDM.Gi[1] = flowVars.rho * flowVars.V * flowVars.U
   FDM.Gi[2] = flowVars.rho * flowVars.V * flowVars.V + flowVars.P
   FDM.Gi[3] = flowVars.rho * flowVars.V * flowVars.et

def updateQvector(imax,jmax):

   f = centralFiniteDifference(-FDM.Fi[0], 'x')

#   for n in range(4):
#      FDM.Q[n] = 


def centralFiniteDifference(phi, direction):
   imax = len(phi[:,0])
   jmax = len(phi[0,:])
   dx = domainVars.dx
   dy = domainVars.dy
   f = np.zeros((imax,jmax))

   # x-derivative
   if direction == 'x':
      f = phi[

   return f

