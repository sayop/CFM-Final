import numpy as np
from variables import *

def updateBC(inputDict,imax,jmax):
   Rgas  = float(inputDict['gasConst'])
   beta  = float(inputDict['beta'])
   gamma = float(inputDict['gamma'])
   Poutflow = float(inputDict['Poutflow'])
   
   # Left boundary: adiabatic, no-slip wall
   # Do NOT update on velocity fields
   # pressure gradient zero normal to the wall
   flowVars.U[0,:] = 0.0
   flowVars.V[0,:] = 0.0
   flowVars.P[0,:] = flowVars.P[1,:]
   # For compressible solution with beta = 0,
   if beta == 0: 
      flowVars.T[0,:]   = flowVars.T[1,:]
      flowVars.rho[0,:] = flowVars.P[0,:] / (Rgas * flowVars.T[0,:])
      flowVars.ei[0,:]  = flowVars.P[0,:] / (flowVars.rho[0,:] * (gamma - 1.0))
      flowVars.et[0,:]  = flowVars.ei[0,:] + 0.5 * (flowVars.U[0,:] ** 2 + flowVars.V[0,:] ** 2)

   # Set nozzle inlet condition
   updateNozzleInletBC(inputDict,imax,jmax, 'Left')


   # Right boundary: adiabatic, no-slip wall
   flowVars.U[imax-1,:] = 0.0
   flowVars.V[imax-1,:] = 0.0
   flowVars.P[imax-1,:] = flowVars.P[imax-2,:]
   if beta == 0:
      flowVars.T[imax-1,:] = flowVars.T[imax-2,:]
      flowVars.rho[imax-1,:] = flowVars.P[imax-1,:] / (Rgas * flowVars.T[imax-1,:])
      flowVars.ei[imax-1,:]  = flowVars.P[imax-1,:] / (flowVars.rho[imax-1,:] * (gamma - 1.0))
      flowVars.et[imax-1,:]  = flowVars.ei[imax-1,:] + 0.5 * (flowVars.U[imax-1,:] ** 2 + flowVars.V[imax-1,:] ** 2)
   # Right boundary: outflow boundary
   #flowVars.U[imax-1,:] = 2.0 * flowVars.U[imax-2,:] - flowVars.U[imax-3,:]
   #flowVars.V[imax-1,:] = 2.0 * flowVars.V[imax-2,:] - flowVars.V[imax-3,:]
   ##if beta > 0: flowVars.P[imax-1,:] = Poutflow / flowVars.RHOref
   #if beta > 0: flowVars.P[imax-1,:] = Poutflow / flowVars.RHOref
   #if beta == 0:
   #   flowVars.P[imax-1,:] = Poutflow
   #   flowVars.T[imax-1,:] = 2.0 * flowVars.T[imax-2,:] - flowVars.T[imax-3,:]
   #   flowVars.rho[imax-1,:] = flowVars.P[imax-1,:] / (Rgas * flowVars.T[imax-1,:])
   #   flowVars.ei[imax-1,:]  = flowVars.P[imax-1,:] / (flowVars.rho[imax-1,:] * (gamma - 1.0))
   #   flowVars.et[imax-1,:]  = flowVars.ei[imax-1,:] + 0.5 * (flowVars.U[imax-1,:] ** 2 + flowVars.V[imax-1,:] ** 2)

   # Set nozzle inlet condition
   updateNozzleInletBC(inputDict,imax,jmax, 'Right')

   # Bottom boundary: adiabatic, no-slip wall
   # Do NOT update on velocity fields
   # pressure gradient zero normal to the wall
   #flowVars.U[:,0] = 0.0
   #flowVars.V[:,0] = 0.0
   #flowVars.P[:,0] = flowVars.P[:,1]
   #if beta == 0:
   #   flowVars.T[:,0]   = flowVars.T[:,1]
   #   flowVars.rho[:,0] = flowVars.P[:,0] / (Rgas * flowVars.T[:,0])
   #   flowVars.ei[:,0]  = flowVars.P[:,0] / (flowVars.rho[:,0] * (gamma - 1.0))
   #   flowVars.et[:,0]  = flowVars.ei[:,0] + 0.5 * (flowVars.U[:,0] ** 2 + flowVars.V[:,0] ** 2)
   #
   # Symmetry boundary condition
   #
   flowVars.U[:,0] = flowVars.U[:,1]
   flowVars.V[:,0] = 0.0
   flowVars.P[:,0] = flowVars.P[:,1]
   if beta == 0:
      flowVars.T[:,0]   = flowVars.T[:,1]
      flowVars.rho[:,0] = flowVars.P[:,0] / (Rgas * flowVars.T[:,0])
      flowVars.ei[:,0]  = flowVars.P[:,0] / (flowVars.rho[:,0] * (gamma - 1.0))
      flowVars.et[:,0]  = flowVars.ei[:,0] + 0.5 * (flowVars.U[:,0] ** 2 + flowVars.V[:,0] ** 2)


   # Upper boundary: adiabatic, no-slip wall
   # Do NOT update on velocity fields
   # pressure gradient zero normal to the wall
   #flowVars.U[:,jmax-1] = 0.0
   #flowVars.V[:,jmax-1] = 0.0
   #flowVars.P[:,jmax-1] = flowVars.P[:,jmax-2]
   #if beta == 0:
   #   flowVars.T[:,jmax-1] = flowVars.T[:,jmax-2]
   #   flowVars.rho[:,jmax-1] = flowVars.P[:,jmax-1] / (Rgas * flowVars.T[:,jmax-1])
   #   flowVars.ei[:,jmax-1]  = flowVars.P[:,jmax-1] / (flowVars.rho[:,jmax-1] * (gamma - 1.0))
   #   flowVars.et[:,jmax-1]  = flowVars.ei[:,jmax-1] + 0.5 * (flowVars.U[:,jmax-1] ** 2 + flowVars.V[:,jmax-1] ** 2)
   # 
   # Upper boundary: open outflow
   #
   flowVars.U[:,jmax-1] = 2.0 * flowVars.U[:,jmax-2] - flowVars.U[:,jmax-3]
   flowVars.V[:,jmax-1] = 2.0 * flowVars.V[:,jmax-2] - flowVars.V[:,jmax-3]
   if beta > 0: flowVars.P[:,jmax-1] = 2.0 * flowVars.P[:,jmax-2] - flowVars.P[:,jmax-3]
   #if beta > 0: flowVars.P[:,jmax-1] = Poutflow / flowVars.RHOref
   if beta == 0:
      flowVars.P[:,jmax-1] = Poutflow
      flowVars.T[:,jmax-1] = 2.0 * flowVars.T[:,jmax-2] - flowVars.T[:,jmax-3]
      flowVars.rho[:,jmax-1] = flowVars.P[:,jmax-1] / (Rgas * flowVars.T[:,jmax-1])
      flowVars.ei[:,jmax-1]  = flowVars.P[:,jmax-1] / (flowVars.rho[:,jmax-1] * (gamma - 1.0))
      flowVars.et[:,jmax-1]  = flowVars.ei[:,jmax-1] + 0.5 * (flowVars.U[:,jmax-1] ** 2 + flowVars.V[:,jmax-1] ** 2)

 

def updateNozzleInletBC(inputDict,imax,jmax,nozzleDir):
   jetRadius = float(inputDict['jetRadius'])
   Pinflow = float(inputDict['Pinflow'])
   ymin = float(inputDict['ymin'])
   ymax = float(inputDict['ymax'])
   Rgas = float(inputDict['gasConst'])
   beta = float(inputDict['beta'])
   #yNozOrig = 0.5 * (ymin + ymax)
   yNozOrig = ymin
   yminNoz = yNozOrig - jetRadius
   ymaxNoz = yNozOrig + jetRadius


   for j in range(jmax):
      if (domainVars.y[j] >= yminNoz) and (domainVars.y[j] <= ymaxNoz):
         if nozzleDir == 'Left':
            flowVars.U[0,j] = flowVars.Ujet
            flowVars.V[0,j] = 0.0
            #if beta > 0: flowVars.P[0,j] = Pinflow / flowVars.RHOref
            if beta > 0: flowVars.P[0,j] = 2.0 * flowVars.P[1,j] - flowVars.P[2,j]
            if beta == 0:
                #flowVars.P[0,j] = Pinflow
                flowVars.P[0,j] = 2.0 * flowVars.P[1,j] - flowVars.P[2,j]
                flowVars.T[0,j] = jetTemp
                flowVars.rho[0,j] = flowVars.P[0,j] / (Rgas * flowVars.T[0,j])
                flowVars.ei[0,j]  = flowVars.P[0,j] / (flowVars.rho[0,j] * (gamma - 1.0))
                flowVars.et[0,j]  = flowVars.ei[0,j] + 0.5 * (flowVars.U[0,j] ** 2 + flowVars.V[0,j] ** 2)

         elif nozzleDir == 'Right':
            flowVars.U[imax-1,j] = -flowVars.Ujet
            flowVars.V[imax-1,j] = 0.0
            #if beta > 0: flowVars.P[0,j] = Pinflow / flowVars.RHOref
            if beta > 0: flowVars.P[imax-1,j] = 2.0 * flowVars.P[imax-2,j] - flowVars.P[imax-3,j]
            if beta == 0:
                #flowVars.P[0,j] = Pinflow
                flowVars.P[imax-1,j] = 2.0 * flowVars.P[imax-2,j] - flowVars.P[imax-3,j]
                flowVars.T[imax-1,j] = jetTemp
                flowVars.rho[imax-1,j] = flowVars.P[imax-1,j] / (Rgas * flowVars.T[imax-1,j])
                flowVars.ei[imax-1,j]  = flowVars.P[imax-1,j] / (flowVars.rho[imax-1,j] * (gamma - 1.0))
                flowVars.et[imax-1,j]  = flowVars.ei[imax-1,j] + 0.5 * (flowVars.U[imax-1,j] ** 2 + flowVars.V[imax-1,j] ** 2)

