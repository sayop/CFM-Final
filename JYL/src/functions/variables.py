import numpy as np

class domainVars:
   x = []
   y = []
   dx = 0.0
   dy = 0.0
   Lref = 0.0

class flowVars:
   Ujet = 0.0
   jetTemp = 0.0
   jetRe = 0.0     # Reynolds number
   RHOref = 0.0
   Uref   = 0.0
   Tref   = 0.0
   MUref  = 0.0
   Kref   = 0.0
   Rref   = 0.0

   rho = []
   P   = []
   U   = []
   V   = []
   T   = []
   et  = []
   ei  = []


class FDM:
   Q   = [None] * 4
   phi = [None] * 4
   Fi  = [None] * 4
   Gi  = [None] * 4
   Fv  = [None] * 4
   Gv  = [None] * 4
   
