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
   nonDim = int(inputDict['nonDim'])   

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
   Djet  = 2.0 * float(inputDict['jetRadius'])
   RHOjet = Pinit / (Rgas * jetTemp)
   Ujet  = muRef * jetRe / (RHOjet * Djet)
   flowVars.Ujet = Ujet

   # Set reference scales
   domainVars.Lref = Djet
   flowVars.Uref = Ujet
   flowVars.RHOref = RHOjet
   flowVars.Tref = jetTemp
   flowVars.MUref = RHOjet * Ujet * Djet
   flowVars.Kref  = flowVars.MUref * Ujet ** 2  / jetTemp
   # dimensionless gas constant
   flowVars.Rref  = Ujet ** 2 / jetTemp

   print '# Gas jet velocity = ', Ujet
   print '# Jet diameter = ', Djet
   print '# Jet Reynolds no = ', jetRe


   flowVars.P    = Pinit * np.ones((imax,jmax))
   if beta == 0: flowVars.T    = Tinit * np.ones((imax,jmax))
   if beta > 0:  flowVars.T    = jetTemp * np.ones((imax,jmax))
   flowVars.U    = Uinit * np.ones((imax,jmax))
   flowVars.V    = Vinit * np.ones((imax,jmax))
   flowVars.rho  = flowVars.P / (Rgas * flowVars.T)
   # compute energy per unit mass from EOS relation
   if beta == 0:
      flowVars.ei   = flowVars.P / (flowVars.rho * (gamma - 1.0))
      flowVars.et   = flowVars.ei + 0.5 * (flowVars.U ** 2 + flowVars.V ** 2)
   # if ACM, change P (dynamic pressure) to kinematic pressure and rho (density) will never be updated.
   if beta > 0:
      flowVars.P = flowVars.P / flowVars.RHOref

   updateBC(inputDict,imax,jmax)

   if nonDim == 1:
      nondimensionalize(inputDict,1,1)

   # state vector is populated with given initial condition at very first beginning.
   # Then, state vector will updated during the time integration process only.
   populateStateVector(inputDict,imax,jmax)
