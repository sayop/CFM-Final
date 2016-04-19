import numpy as np
from variables import *

def updateVelocityBC(imax,jmax):
   # update boundary values of u and v in dimensionalized form
   # Left boundary
   flowVars.U[0,:] = 0.0
   flowVars.V[0,:] = 0.0
   # Bottom boundary (No-slip)
   flowVars.U[:,0] = 0.0
   flowVars.V[:,0] = 0.0
   # Upper boundary (No-slip)
   flowVars.U[:,jmax-1] = 0.0
   flowVars.V[:,jmax-1] = 0.0



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
  
