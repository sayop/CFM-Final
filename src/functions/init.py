import numpy as np
from variables import *
from transport import sutherland
from boundaryCondition import *
from solutionMethod import *

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
   # For Artificial Compressibility Method
   beta  = float(inputDict['beta'])
   
   # Cv: constant volume specific heat
   Cv    = Rgas / (gamma - 1.0)

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
   # For ACM, convert dynamic pressure into kinematic pressure
   if beta > 0.0:
      # Since ACM is for incompressible solution, rhoRef is used for constant density.
      # ACM does NOT update density. INCOMPRESSIBLE solution.
      flowVars.rhoRef = Pinit / (Rgas * Tinit)
      flowVars.P = flowVars.P / flowVars.rhoRef
   else:
      flowVars.rho  = flowVars.P / (Rgas * flowVars.T)

   # compute energy per unit mass from EOS relation
   flowVars.ei   = Cv * flowVars.T
   flowVars.et   = flowVars.ei + 0.5 * (flowVars.U ** 2 + flowVars.V ** 2)

   updateBC(inputDict,imax,jmax)

   # state vector is populated with given initial condition at very first beginning.
   # Then, state vector will updated during the time integration process only.
   populateStateVector(inputDict,imax,jmax)
