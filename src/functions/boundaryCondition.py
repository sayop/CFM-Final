import numpy as np
from variables import *

def updateBC(imax,jmax):
   # Left boundary: adiabatic, no-slip wall
   # Do NOT update on velocity fields
   # pressure gradient zero normal to the wall
   #flowVars.U[0,:] = 0.0
   #flowVars.V[0,:] = 0.0
   flowVars.rho[0,:] = flowVars.rho[1,:]
   flowVars.P[0,:]   = flowVars.P[1,:]
   flowVars.T[0,:]   = flowVars.T[1,:]
   flowVars.et[0,:]  = flowVars.et[1,:]
   
  
   # Right boundary: outflow boundary
   flowVars.rho[imax-1,0:jmax] = flowVars.rho[imax-2,0:jmax]
   flowVars.U[imax-1,0:jmax]   = flowVars.U[imax-2,0:jmax]
   flowVars.V[imax-1,0:jmax]   = flowVars.V[imax-2,0:jmax]
   flowVars.P[imax-1,0:jmax]   = flowVars.P[imax-2,0:jmax]
   flowVars.T[imax-1,0:jmax]   = flowVars.T[imax-2,0:jmax]
   flowVars.et[imax-1,0:jmax]  = flowVars.et[imax-2,0:jmax]


   # Bottom boundary: adiabatic, no-slip wall
   # Do NOT update on velocity fields
   # pressure gradient zero normal to the wall
   #flowVars.U[0,:] = 0.0
   #flowVars.V[0,:] = 0.0
   flowVars.rho[:,0] = flowVars.rho[:,1]
   flowVars.P[:,0]   = flowVars.P[:,1]
   flowVars.T[:,0]   = flowVars.T[:,1]
   flowVars.et[:,0]  = flowVars.et[:,1]


   # Upper boundary: adiabatic, no-slip wall
   # Do NOT update on velocity fields
   # pressure gradient zero normal to the wall
   #flowVars.U[0,:] = 0.0
   #flowVars.V[0,:] = 0.0
   flowVars.rho[:,jmax-1] = flowVars.rho[:,jmax-2]
   flowVars.P[:,jmax-1]   = flowVars.P[:,jmax-2]
   flowVars.T[:,jmax-1]   = flowVars.T[:,jmax-2]
   flowVars.et[:,jmax-1]  = flowVars.et[:,jmax-2]




def updateNozzleInletBC(inputDict,imax,jmax,Ujet,jetTemp):
   jetRadius = float(inputDict['jetRadius'])
   ymin = float(inputDict['ymin'])
   ymax = float(inputDict['ymax'])
   yNozOrig = 0.5 * (ymin + ymax)
   yminNoz = yNozOrig - jetRadius
   ymaxNoz = yNozOrig + jetRadius


   for j in range(jmax):
      if (domainVars.y[j] >= yminNoz) and (domainVars.y[j] <= ymaxNoz):
         flowVars.U[0,j] = Ujet
         flowVars.T[0,j] = jetTemp
  
