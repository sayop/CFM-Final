import numpy as np
from variables import *

def updateBC(inputDict,imax,jmax):
   Rgas  = float(inputDict['gasConst'])
   beta  = float(inputDict['beta'])
   gamma = float(inputDict['gamma'])
   Poutflow = float(inputDict['Poutflow'])
   Cv    = Rgas / (gamma - 1.0)
   
   # Left boundary: adiabatic, no-slip wall
   # Do NOT update on velocity fields
   # pressure gradient zero normal to the wall
   flowVars.U[0,:] = 0.0
   flowVars.V[0,:] = 0.0
   # dP/dN = 0
   flowVars.P[0,:]   = flowVars.P[1,:]
   #Adiabatic wall
   flowVars.T[0,:]   = flowVars.T[1,:]
   # For compressible solution with beta = 0,
   if beta == 0: 
      flowVars.rho[0,:] = flowVars.P[0,:] / (Rgas * flowVars.T[0,:])
 
   flowVars.ei[0,:]  = Cv * flowVars.T[0,:]
   flowVars.et[0,:]  = flowVars.ei[0,:] + 0.5 * (flowVars.U[0,:] ** 2 + flowVars.V[0,:] ** 2)

   # Right boundary: outflow boundary
   flowVars.U[imax-1,0:jmax] = 0.0
   flowVars.V[imax-1,0:jmax] = 0.0
   flowVars.P[imax-1,0:jmax] = flowVars.P[imax-2,0:jmax]
   flowVars.T[imax-1,0:jmax] = flowVars.T[imax-2,0:jmax]
   if beta == 0:
      flowVars.rho[imax-1,0:jmax] = flowVars.P[imax-1,0:jmax] / (Rgas * flowVars.T[imax-1,0:jmax])

   flowVars.ei[imax-1,0:jmax]  = Cv * flowVars.T[imax-1,0:jmax]
   flowVars.et[imax-1,0:jmax]  = flowVars.ei[imax-1,0:jmax] + 0.5 * (flowVars.U[imax-1,0:jmax] ** 2 + flowVars.V[imax-1,0:jmax] ** 2)

   # Bottom boundary: adiabatic, no-slip wall
   # Do NOT update on velocity fields
   # pressure gradient zero normal to the wall
   flowVars.U[:,0] = 0.0
   flowVars.V[:,0] = 0.0
   flowVars.P[:,0]   = flowVars.P[:,1]
   flowVars.T[:,0]   = flowVars.T[:,1]
   if beta == 0:
      flowVars.rho[:,0] = flowVars.P[:,0] / (Rgas * flowVars.T[:,0])

   flowVars.ei[:,0]  = Cv * flowVars.T[:,0]
   flowVars.et[:,0]  = flowVars.ei[:,0] + 0.5 * (flowVars.U[:,0] ** 2 + flowVars.V[:,0] ** 2)

   # Upper boundary: adiabatic, no-slip wall
   # Do NOT update on velocity fields
   # pressure gradient zero normal to the wall
   flowVars.U[:,jmax-1] = flowVars.Ujet
   flowVars.V[:,jmax-1] = 0.0
   flowVars.P[:,jmax-1] = flowVars.P[:,jmax-2]
   flowVars.T[:,jmax-1] = flowVars.T[:,jmax-2]
   if beta == 0:
      flowVars.rho[:,jmax-1] = flowVars.P[:,jmax-1] / (Rgas * flowVars.T[:,jmax-1])

   flowVars.ei[:,jmax-1]  = Cv * flowVars.T[:,jmax-1]
   flowVars.et[:,jmax-1]  = flowVars.ei[:,jmax-1] + 0.5 * (flowVars.U[:,jmax-1] ** 2 + flowVars.V[:,jmax-1] ** 2)


def updateBC_backup(inputDict,imax,jmax):
   Rgas  = float(inputDict['gasConst'])
   gamma = float(inputDict['gamma'])
   Poutflow = float(inputDict['Poutflow'])
   Cv    = Rgas / (gamma - 1.0)
   
   # Left boundary: adiabatic, no-slip wall
   # Do NOT update on velocity fields
   # pressure gradient zero normal to the wall
   flowVars.U[0,:] = 0.0
   flowVars.V[0,:] = 0.0
   # dP/dN = 0
   flowVars.P[0,:]   = flowVars.P[1,:]
   #Adiabatic wall
   flowVars.T[0,:]   = flowVars.T[1,:]
   flowVars.rho[0,:] = flowVars.P[0,:] / (Rgas * flowVars.T[0,:])
 
   # Set nozzle inlet condition
   updateNozzleInletBC(inputDict,imax,jmax,flowVars.Ujet,flowVars.jetTemp)
  
   
   flowVars.ei[0,:]  = Cv * flowVars.T[0,:]
   flowVars.et[0,:]  = flowVars.ei[0,:] + 0.5 * (flowVars.U[0,:] ** 2 + flowVars.V[0,:] ** 2)

   # Right boundary: outflow boundary
   flowVars.U[imax-1,0:jmax]   = 2.0 * flowVars.U[imax-2,0:jmax] - flowVars.U[imax-3,0:jmax]
   flowVars.V[imax-1,0:jmax]   = 2.0 * flowVars.V[imax-2,0:jmax] - flowVars.V[imax-3,0:jmax]
   flowVars.P[imax-1,0:jmax]   = Poutflow
   flowVars.T[imax-1,0:jmax]   = 2.0 * flowVars.T[imax-2,0:jmax] - flowVars.T[imax-3,0:jmax]
   flowVars.rho[imax-1,0:jmax] = flowVars.P[imax-1,0:jmax] / (Rgas * flowVars.T[imax-1,0:jmax])

   flowVars.ei[imax-1,0:jmax]  = Cv * flowVars.T[imax-1,0:jmax]
   flowVars.et[imax-1,0:jmax]  = flowVars.ei[imax-1,0:jmax] + 0.5 * (flowVars.U[imax-1,0:jmax] ** 2 + flowVars.V[imax-1,0:jmax] ** 2)

   # Bottom boundary: adiabatic, no-slip wall
   # Do NOT update on velocity fields
   # pressure gradient zero normal to the wall
   flowVars.U[:,0] = 0.0
   flowVars.V[:,0] = 0.0
   flowVars.P[:,0]   = flowVars.P[:,1]
   flowVars.T[:,0]   = flowVars.T[:,1]
   flowVars.rho[:,0] = flowVars.P[:,0] / (Rgas * flowVars.T[:,0])

   flowVars.ei[:,0]  = Cv * flowVars.T[:,0]
   flowVars.et[:,0]  = flowVars.ei[:,0] + 0.5 * (flowVars.U[:,0] ** 2 + flowVars.V[:,0] ** 2)

   # Upper boundary: adiabatic, no-slip wall
   # Do NOT update on velocity fields
   # pressure gradient zero normal to the wall
   flowVars.U[:,jmax-1] = 0.0
   flowVars.V[:,jmax-1] = 0.0
   flowVars.P[:,jmax-1]   = flowVars.P[:,jmax-2]
   flowVars.T[:,jmax-1]   = flowVars.T[:,jmax-2]
   flowVars.rho[:,jmax-1] = flowVars.P[:,jmax-1] / (Rgas * flowVars.T[:,jmax-1])

   flowVars.ei[:,jmax-1]  = Cv * flowVars.T[:,jmax-1]
   flowVars.et[:,jmax-1]  = flowVars.ei[:,jmax-1] + 0.5 * (flowVars.U[:,jmax-1] ** 2 + flowVars.V[:,jmax-1] ** 2)


def updateNozzleInletBC(inputDict,imax,jmax,Ujet,jetTemp):
   jetRadius = float(inputDict['jetRadius'])
   Pinflow = float(inputDict['Pinflow'])
   ymin = float(inputDict['ymin'])
   ymax = float(inputDict['ymax'])
   Rgas = float(inputDict['gasConst'])
   yNozOrig = 0.5 * (ymin + ymax)
   yminNoz = yNozOrig - jetRadius
   ymaxNoz = yNozOrig + jetRadius


   for j in range(jmax):
      if (domainVars.y[j] >= yminNoz) and (domainVars.y[j] <= ymaxNoz):
         flowVars.U[0,j] = Ujet
         flowVars.V[0,j] = 0.0
         flowVars.T[0,j] = jetTemp
         flowVars.rho[0,j] = 2.0*flowVars.rho[1,j] - flowVars.rho[2,j]
         flowVars.P[0,j] = flowVars.rho[0,j] * Rgas * flowVars.T[0,j]

         #flowVars.U[0,j] = 2.0 * flowVars.U[1,j] - flowVars.U[2,j]
         #flowVars.V[0,j] = 2.0 * flowVars.V[1,j] - flowVars.V[2,j]
         #flowVars.T[0,j] = jetTemp
         #flowVars.P[0,j] = Pinflow
         #flowVars.rho[0,j] = Pinflow  / (Rgas * jetTemp)
