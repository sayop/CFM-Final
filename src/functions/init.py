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

   RHOinit = Pinit / (Rgas * Tinit)

   flowVars.rho  = RHOinit * np.ones((imax,jmax))
   flowVars.P    = Pinit * np.ones((imax,jmax))
   flowVars.T    = Tinit * np.ones((imax,jmax))
   flowVars.U    = Uinit * np.ones((imax,jmax))
   flowVars.V    = Vinit * np.ones((imax,jmax))
   # energy per unit mass
   flowVars.et   = flowVars.P / (gamma - 1.0) / flowVars.rho + 0.5 * (flowVars.U ** 2 + flowVars.V ** 2)

   # Reference parameters defined at jet exit
   muRef, kRef = sutherland(Tinit)
   print '# Reference gas viscosity = ', muRef
   print '# Reference gas thermal conductivity = ', kRef

   jetRe = float(inputDict['jetRe'])
   flowVars.jetRe = jetRe
   Dref  = 2.0 * float(inputDict['jetRadius'])
   Ujet  = muRef * jetRe / (RHOinit * Dref)
   flowVars.Ujet = Ujet

   print '# Gas jet velocity = ', Ujet
   print '# Jet diameter = ', Dref
   print '# Jet Reynolds no = ', jetRe

   # Populate boundary values that is specified by user in inputs.in
   updateVelocityBC(imax,jmax)

   # update jet nozzle inlet BC
   updateNozzleInletBC(inputDict,imax,jmax,Ujet)
