import numpy as np
from variables import *

def updateDirichletBC(inputDict,imax,jmax,Ujet):
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
