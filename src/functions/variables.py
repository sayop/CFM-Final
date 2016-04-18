import numpy as np

class domainVars:
   x = []
   y = []
   dx = 0.0
   dy = 0.0
   Lref = 0.0

class flowVars:
   Ujet = 0.0
   jetRe = 0.0     # Reynolds number
   rho = []
   P   = []
   U   = []
   V   = []
   T   = []
   Et  = []


class FDM:
   Fi  = [] * 4
   Gi  = [] * 4
   Fv  = [] * 4
   Gv  = [] * 4
   
