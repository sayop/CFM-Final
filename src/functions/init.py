import numpy as np
from variables import *
from transport import sutherland
from boundaryCondition import *

def initSimulationVars(inputDict):
   print '# Initializing flow variables...'
   imax = int(inputDict['iDim'])
   jmax = int(inputDict['jDim'])

   Pinit = float(inputDict['Pinit'])
   Tinit = float(inputDict['Tinit'])
   Uinit = float(inputDict['Uinit'])
   Vinit = float(inputDict['Vinit'])
   Rgas  = float(inputDict['gasConst'])
   gamma = float(inputDict['gamma'])

   # Reference parameters defined at jet exit
   jetTemp = float(inputDict['jetTemp'])
   flowVars.jetTemp = jetTemp
   muRef, kRef = sutherland(jetTemp)
   print '# Reference gas viscosity = ', muRef
   print '# Reference gas thermal conductivity = ', kRef

   jetRe = float(inputDict['jetRe'])
   flowVars.jetRe = jetRe
   Dref  = 2.0 * float(inputDict['jetRadius'])
   RHOjet = Pinit / (Rgas * jetTemp)
   Ujet  = muRef * jetRe / (RHOjet * Dref)
   flowVars.Ujet = Ujet

   print '# Gas jet velocity = ', Ujet
   print '# Jet diameter = ', Dref
   print '# Jet Reynolds no = ', jetRe


   flowVars.P    = Pinit * np.ones((imax,jmax))
   flowVars.T    = Tinit * np.ones((imax,jmax))
   flowVars.U    = Uinit * np.ones((imax,jmax))
   flowVars.V    = Vinit * np.ones((imax,jmax))

   flowVars.rho  = flowVars.P / (Rgas * flowVars.T)

   # energy per unit mass
   # ei: internal energy per unit mass
   ei = flowVars.P / flowVars.rho / (gamma - 1.0)
   flowVars.et   = ei + 0.5 * (flowVars.U ** 2 + flowVars.V ** 2)

   updateBC(inputDict,imax,jmax)
