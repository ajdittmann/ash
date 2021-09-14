import numpy as np

def ash1d(vals, nbins=20, nshifts=10, weights=None):
  vmax = np.max(vals)
  vmin = np.min(vals)
  dx = (vmax-vmin)/nbins
  ds = dx/nshifts

  ashgrid = np.linspace(vmin, vmax, (nbins+2)*nshifts) + 0.5*(vmax-vmin)/(nbins+1)
  ashden = np.zeros(ashgrid.shape)
  for i in range(nshifts):
    hist_range = (vmin + i*ds, vmax + i*ds)
    hist, edges = np.histogram(vals, nbins+1, range=hist_range, density=True, weights = weights)
    histdat = np.repeat(hist, nshifts, axis=0)
    histdat = np.append( np.zeros(i), histdat, axis=0)
    ashden += np.append( histdat, np.zeros(nshifts-i), axis=0)
  ashden = ashden/nshifts
  ashind = np.where(ashden>0)
  ashgrid = ashgrid[ashind]
  ashden = ashden[ashind]
  return ashgrid, ashden


def ash2d(xvals, yvals, nbins=20, nshifts=10, weights=None):
  xmax = np.max(xvals)
  xmin = np.min(xvals)
  dx = (xmax-xmin)/nbins
  dsx = dx/nshifts

  ymax = np.max(yvals)
  ymin = np.min(yvals)
  dy = (ymax-ymin)/nbins
  dsy = dy/nshifts

  nsg = (nbins+2)*nshifts
  xgrid = np.linspace(xmin, xmax, (nbins+2)*nshifts) + 0.5*(xmax-xmin)/(nbins+1)
  ygrid = np.linspace(xmin, xmax, (nbins+2)*nshifts) + 0.5*(ymax-ymin)/(nbins+1)
  ashmesh = np.meshgrid(xgrid, ygrid)
  ashden = np.zeros((nsg, nsg))

  for i in range(nshifts):  #x loop
    for j in range(nshifts):  #y loop
      histrange = ( (xmin + i*dsx, xmax + i*dsx), (ymin + j*dsy, ymax + j*dsy) )
      hist, xedge, yedge = np.histogram2d(xvals, yvals, nbins+1, range=histrange, density=True, weights=weights)
      histt = np.repeat(hist, nshifts, axis=0)
      histt = np.repeat(histt, nshifts, axis=1)
      histt = np.append( np.zeros((i,(nbins+1)*nshifts)), histt, axis=0)
      histt = np.append( histt, np.zeros((nshifts-i, (nbins+1)*nshifts)), axis=0)
      histt = np.append( np.zeros(((nbins+2)*nshifts,j)), histt, axis=1)
      histt = np.append( histt, np.zeros(( (nbins+2)*nshifts,nshifts-j)), axis=1)
      ashden += histt
  ashden /= nshifts

  return ashmesh, ashden

