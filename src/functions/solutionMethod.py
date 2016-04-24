from variables import *

import numpy as np

def computeTimeStep(inputDict,imax,jmax):
   Cr      = float(inputDict['Courant'])
   gamma   = float(inputDict['gamma'])
   Rgas    = float(inputDict['gasConst'])
   beta    = float(inputDict['beta'])
   dx = domainVars.dx
   dy = domainVars.dy

   # calculate speed of sound
   if beta > 0.0:
      a = 1.0 / np.sqrt(beta)
   else:
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

def populateStateVector(inputDict,imax,jmax):
   Rgas  = float(inputDict['gasConst'])
   gamma = float(inputDict['gamma'])
   beta  = float(inputDict['beta'])
   # transform primative variables to elements of state vector
   #State vector
   # each elements are only effective in interior node points
   if beta == 0: 
      FDM.phi[0] = flowVars.rho
      FDM.phi[1] = flowVars.rho * flowVars.U
      FDM.phi[2] = flowVars.rho * flowVars.V
      FDM.phi[3] = flowVars.rho * flowVars.et
   elif beta > 0.0:
      FDM.phi[0] = flowVars.P
      FDM.phi[1] = flowVars.U
      FDM.phi[2] = flowVars.V
      FDM.phi[3] = flowVars.et


def updateFluxVectors(inputDict,imax,jmax,iVisc):
   beta  = float(inputDict['beta'])
   
   if beta > 0.0:
      FDM.Fi[0] = flowVars.U / beta
      FDM.Fi[1] = flowVars.U * flowVars.U + flowVars.P
      FDM.Fi[2] = flowVars.U * flowVars.V
      FDM.Fi[3] = flowVars.U * flowVars.et

      FDM.Gi[0] = flowVars.V / beta
      FDM.Gi[1] = flowVars.V * flowVars.U
      FDM.Gi[2] = flowVars.V * flowVars.V + flowVars.P
      FDM.Gi[3] = flowVars.V * flowVars.et

   else:
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

   # update diffusive flux
   if iVisc == 1:
      Tau, Qj = computeDiffusiveTransport(imax,jmax)
      if beta > 0.0:
         Tau = Tau / flowVars.rhoRef
         Qj  = Qj  / flowVars.rhoRef

      FDM.Fv[0] = np.zeros((imax,jmax))
      FDM.Fv[1] = Tau[:,:,0,0]
      FDM.Fv[2] = Tau[:,:,0,1]
      FDM.Fv[3] = -Qj[:,:,0] + flowVars.U * Tau[:,:,0,0] + flowVars.V * Tau[:,:,0,1]

      FDM.Gv[0] = np.zeros((imax,jmax))
      FDM.Gv[1] = Tau[:,:,1,0]
      FDM.Gv[2] = Tau[:,:,1,1]
      FDM.Gv[3] = -Qj[:,:,1] + flowVars.U * Tau[:,:,1,0] + flowVars.V * Tau[:,:,1,1]



def computeDiffusiveTransport(imax,jmax):
   from transport import sutherland
   # gas transport properties
   mu, k = sutherland(flowVars.T)
   # UiXj: Xj derivative of Ui
   UiXj = np.zeros((imax,jmax,2,2))

   for i in range(2):
      if i == 0: Ui = flowVars.U
      if i == 1: Ui = flowVars.V
      for j in range(2):
         if j == 0: direction = 'x'
         if j == 1: direction = 'y'
         UiXj[:,:,i,j] = centralFiniteDifference(Ui, direction)

   # UiXi: divergence of velocity vector
   # UiXi = dUdX + dVdY
   UiXi = np.zeros((imax,jmax))
   for j in range(jmax):
      for i in range(imax):
         UiXi[i,j] = UiXj[i,j,0,0] + UiXj[i,j,1,1]

   # evaluate viscous stress tensor
   Tau = np.zeros((imax,jmax,2,2))
   for i in range(2):
      for j in range(2):
         # Kronecker delta function
         if i == j:
            kDelta = 1.0
         else:
            kDelta = 0.0
         # Bulk viscosity coefficient: Lambda = -2/3 * mu
         Lambda = -2.0 / 3.0 * mu[i,j]
         Tau[:,:,i,j] = mu[i,j] * (UiXj[:,:,i,j] + UiXj[:,:,j,i]) + kDelta * Lambda * UiXi


   # evaluate heat flux
   Qj = np.zeros((imax,jmax,2))
   for j in range(2):
      if j == 0: direction = 'x'
      if j == 1: direction = 'y'
      Qj[:,:,j] = k * centralFiniteDifference(-flowVars.T, direction)

   return Tau, Qj

def updateQvector(imax,jmax,iVisc):

   for n in range(4):
      FDM.Q[n] = np.zeros((imax,jmax))
      FDM.Q[n] += centralFiniteDifference(-FDM.Fi[n], 'x')
      FDM.Q[n] += centralFiniteDifference(-FDM.Gi[n], 'y')
      if iVisc == 1:
         FDM.Q[n] += centralFiniteDifference(FDM.Fv[n], 'x')
         FDM.Q[n] += centralFiniteDifference(FDM.Gv[n], 'y')


def integrateStateVector(dt):
   for n in range(4):
      FDM.phi[n] += dt * FDM.Q[n]

def updatePrimitiveVariables(inputDict,imax,jmax):
   Rgas  = float(inputDict['gasConst'])
   gamma = float(inputDict['gamma'])
   beta  = float(inputDict['beta'])
   if beta == 0:
      # Cv: constant volume specific heat
      Cv    = Rgas / (gamma - 1.0)
      # update only interior points
      flowVars.rho[1:imax-1,1:jmax-1] = FDM.phi[0][1:imax-1,1:jmax-1]
      flowVars.U[1:imax-1,1:jmax-1]   = FDM.phi[1][1:imax-1,1:jmax-1] / flowVars.rho[1:imax-1,1:jmax-1]
      flowVars.V[1:imax-1,1:jmax-1]   = FDM.phi[2][1:imax-1,1:jmax-1] / flowVars.rho[1:imax-1,1:jmax-1]
      flowVars.et[1:imax-1,1:jmax-1]  = FDM.phi[3][1:imax-1,1:jmax-1] / flowVars.rho[1:imax-1,1:jmax-1]
      # internal energy
      flowVars.ei[1:imax-1,1:jmax-1] = flowVars.et[1:imax-1,1:jmax-1] - 0.5 * (flowVars.U[1:imax-1,1:jmax-1] ** 2 - flowVars.V[1:imax-1,1:jmax-1] ** 2)
      flowVars.T[1:imax-1,1:jmax-1]   = flowVars.ei[1:imax-1,1:jmax-1] / Cv
      flowVars.P[1:imax-1,1:jmax-1]   = flowVars.rho[1:imax-1,1:jmax-1] * Rgas * flowVars.T[1:imax-1,1:jmax-1]
   else:
      flowVars.P[1:imax-1,1:jmax-1] = FDM.phi[0][1:imax-1,1:jmax-1]
      flowVars.U[1:imax-1,1:jmax-1] = FDM.phi[1][1:imax-1,1:jmax-1]
      flowVars.V[1:imax-1,1:jmax-1] = FDM.phi[2][1:imax-1,1:jmax-1]


def centralFiniteDifference(phi, direction):
   imax = len(phi[:,0])
   jmax = len(phi[0,:])
   dx = domainVars.dx
   dy = domainVars.dy
   f = np.zeros((imax,jmax))

   # x-derivative
   if direction == 'x':
      #f[1:imax-1,1:jmax-1] = 0.5 * (phi[2:imax,1:jmax-1] - phi[0:imax-2,1:jmax-1]) / dx 
      f[1:imax-1,:] = 0.5 * (phi[2:imax,:] - phi[0:imax-2,:]) / dx 
      # Forward
      #f[0,:] = (phi[1,:] - phi[0,:]) / dx
      f[0,:] = 0.5 * (-3.0 * phi[0,:] + 4.0 * phi[1,:] - phi[2,:]) / dx
      # Backward
      #f[imax-1,:] = (phi[imax-1,:] - phi[imax-2,:]) / dx
      f[imax-1,:] = 0.5 * (3.0 * phi[imax-1,:] - 4.0 * phi[imax-2,:] + phi[imax-3,:]) / dx

   elif direction == 'y':
      #f[1:imax-1,1:jmax-1] = 0.5 * (phi[1:imax-1,2:jmax] - phi[1:imax-1,0:jmax-2]) / dy
      f[:,1:jmax-1] = 0.5 * (phi[:,2:jmax] - phi[:,0:jmax-2]) / dy
      # Forward
      #f[:,0] = (phi[:,1] - phi[:,0]) / dy
      f[:,0] = 0.5 * (-3.0 * phi[:,0] + 4.0 * phi[:,1] - phi[:,2]) / dy
      # Backward
      #f[:,jmax-1] = (phi[:,jmax-1] - phi[:,jmax-2]) / dy
      f[:,jmax-1] = 0.5 * (3.0 * phi[:,jmax-1] - 4.0 * phi[:,jmax-2] + phi[:,jmax-3]) / dy
      

   return f

