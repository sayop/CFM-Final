from variables import *

import numpy as np

def computeSpeedOfSound(temp,gamma,Rgas):

   a = (gamma * Rgas * temp) ** (0.5)

   return a
