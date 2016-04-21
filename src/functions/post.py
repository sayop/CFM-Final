import numpy as np
import matplotlib.pyplot as plt
import csv

def plotContour(x,y,U,V,nIter):
   from matplotlib import mlab, cm
   cmap = cm.PRGn

   pltFile = 'contour_vel_%5.5d.png' % int(nIter)
   x = np.asarray(x)
   y = np.asarray(y)
   #U = np.swapaxes(U,1,0)
   #V = np.swapaxes(V,1,0)
   #xi, yi = np.meshgrid(x, y)
   #zi = np.swapaxes(phi,1,0)

   Umag = np.sqrt( U ** 2 + V ** 2 )
   Umag = np.swapaxes(Umag,0,-1)
   phiMin = Umag.min()
   phiMax = Umag.max()
   #nLevels = 10
   #dPhi = (phiMax - phiMin) / nLevels
   #levels = np.arange(phiMin, phiMax*1.001, dPhi)

   #plt.contour(xi,yi,zi,levels)
   plt.imshow(Umag, vmin=phiMin, vmax=phiMax, extent=[x.min(), x.max(), y.min(), y.max()])
   plt.colorbar()

   plt.xscale('linear')
   plt.yscale('linear')
   plt.xlabel('x', fontsize=18)
   plt.ylabel('y', fontsize=18)
   #plt.grid(True)
   ax = plt.gca()
   xlabels = plt.getp(ax, 'xticklabels')
   ylabels = plt.getp(ax, 'yticklabels')
   plt.setp(xlabels, fontsize=15)
   plt.setp(ylabels, fontsize=15)

   fig = plt.gcf()
   fig.set_size_inches(6,5)
   plt.tight_layout()
   plt.savefig(pltFile, format='png')
   plt.close()

   print "%s DONE!!" % (pltFile)
   plt.show()


def plotStreamLine(x,y,U,V,nIter):
   pltFile = 'streamLine_%5.5d.png' % int(nIter)
   x = np.asarray(x)
   y = np.asarray(y)
   U = np.swapaxes(U,1,0)
   V = np.swapaxes(V,1,0)

   strm = plt.streamplot(x,y,U,V, color='k', density=1, linewidth=1)

   plt.axis([x.min(), x.max(), y.min(), y.max()])
   plt.xscale('linear')
   plt.yscale('linear')
   plt.xlabel('x [m]', fontsize=18)
   plt.ylabel('y [m]', fontsize=18)
   plt.grid(True)
   ax = plt.gca()
   xlabels = plt.getp(ax, 'xticklabels')
   ylabels = plt.getp(ax, 'yticklabels')
   plt.setp(xlabels, fontsize=10)
   plt.setp(ylabels, fontsize=10)

   fig = plt.gcf()
   fig.set_size_inches(6,5)
   plt.tight_layout()
   plt.savefig(pltFile, format='png')
   plt.close()

   print "%s DONE!!" % (pltFile)
   plt.show()


def plotTempContour(x,y,T,nIter):
   from matplotlib import mlab, cm
   cmap = cm.PRGn

   pltFile = 'contour_Temp_%5.5d.png' % int(nIter)
   x = np.asarray(x)
   y = np.asarray(y)

   phiMin = T.min()
   phiMax = T.max()

   plt.imshow(T, vmin=phiMin, vmax=phiMax, extent=[x.min(), x.max(), y.min(), y.max()])
   plt.colorbar()

   plt.xscale('linear')
   plt.yscale('linear')
   plt.xlabel('x', fontsize=18)
   plt.ylabel('y', fontsize=18)
   #plt.grid(True)
   ax = plt.gca()
   xlabels = plt.getp(ax, 'xticklabels')
   ylabels = plt.getp(ax, 'yticklabels')
   plt.setp(xlabels, fontsize=15)
   plt.setp(ylabels, fontsize=15)

   fig = plt.gcf()
   fig.set_size_inches(6,5)
   plt.tight_layout()
   plt.savefig(pltFile, format='png')
   plt.close()

   print "%s DONE!!" % (pltFile)
   plt.show()
