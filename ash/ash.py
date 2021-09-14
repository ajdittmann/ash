import numpy as np

def ash1d(vals, nbins=20, nshifts=10, weights=None, padVals=False):
  vmax = np.max(vals)
  vmin = np.min(vals)
  dx = (vmax-vmin)/nbins
  ds = dx/nshifts

  ashgrid = np.linspace(vmin, vmax, (nbins+2)*nshifts)
  ashden = np.zeros(ashgrid.shape)
  for i in range(nshifts):
    hist_range = (vmin + i*ds, vmax + i*ds - dx)
    hist, edges = np.histogram(vals, nbins+1, range=hist_range, density=False, weights = weights)
    histdat = np.repeat(hist, nshifts, axis=0)
    if (padVals):
      histdat = np.append( np.ones(i)*hist[0], histdat, axis=0)
      histdat = np.append( histdat, np.ones(nshifts-i)*hist[-1], axis=0)
    else:
      histdat = np.append( np.zeros(i), histdat, axis=0)
      histdat = np.append( histdat, np.zeros(nshifts-i), axis=0)
    ashden += histdat/np.sum(histdat*(vmax - vmin - dx)/((nbins+1)*nshifts))
  ashden = ashden/nshifts
  ashind = np.where(ashden>0)
  ashgrid = ashgrid[ashind]
  ashden = ashden[ashind]
  return ashgrid, ashden


def ash2d(xvals, yvals, nbins=20, nshifts=10, weights=None, padVals=False):
  xmax = np.max(xvals)
  xmin = np.min(xvals)
  dx = (xmax-xmin)/nbins
  dsx = dx/nshifts

  ymax = np.max(yvals)
  ymin = np.min(yvals)
  dy = (ymax-ymin)/nbins
  dsy = dy/nshifts

  nsg = (nbins+2)*nshifts
  xgrid = np.linspace(xmin, xmax, (nbins+2)*nshifts)
  ygrid = np.linspace(xmin, xmax, (nbins+2)*nshifts)
  ashmesh = np.meshgrid(xgrid, ygrid)
  ashden = np.zeros((nsg, nsg))

  for i in range(nshifts):  #x loop
    for j in range(nshifts):  #y loop
      histrange = ( (xmin + i*dsx, xmax + i*dsx -dx), (ymin + j*dsy, ymax + j*dsy - dy) )
      hist, xedge, yedge = np.histogram2d(xvals, yvals, nbins+1, range=histrange, density=False, weights=weights)
      histdat = np.repeat(hist, nshifts, axis=0)
      histdat = np.repeat(histdat, nshifts, axis=1)
      if (padVals):
        appx1 = np.ones((i,(nbins+1)*nshifts))
        appx2 = np.ones((nshifts - i,(nbins+1)*nshifts))
        appy1 = np.ones(((nbins+2)*nshifts, j))
        appy2 = np.ones(((nbins+2)*nshifts, nshifts-j))
        for k in range(i):
          appx1[k,:] = histdat[0,:]
        for k in range(nshifts-i):
          appx2[k,:] = histdat[-1,:]
        histdat = np.append( appx1, histdat, axis=0)
        histdat = np.append( histdat, appx2, axis=0)
        for k in range(j):
          appy1[:,k] = histdat[:,0]
        for k in range(nshifts-j):
          appy2[:,k] = histdat[:,-1]
        histdat = np.append( appy1, histdat, axis=1)
        histdat = np.append( histdat, appy2, axis=1)
      else:
        histdat = np.append( np.zeros((i,(nbins+1)*nshifts)), histdat, axis=0)
        histdat = np.append( histdat, np.zeros((nshifts-i, (nbins+1)*nshifts)), axis=0)
        histdat = np.append( np.zeros(((nbins+2)*nshifts,j)), histdat, axis=1)
        histdat = np.append( histdat, np.zeros(( (nbins+2)*nshifts,nshifts-j)), axis=1)
      ashden += histdat/np.sum(histdat*(xmax - xmin - dx)*(ymax - ymin - dy)/((nbins+1)*nshifts)**2)

  ashden /= nshifts

  return ashmesh, ashden

