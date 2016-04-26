def sutherland(temp):

   C1 = 1.458e-6
   C2 = 110.4
   C3 = 2.495e-3
   C4 = 194.0

   mu = C1 * temp ** (1.5) / (temp + C2)
   k  = C3 * temp ** (1.5) / (temp + C4)

   return mu, k
 
