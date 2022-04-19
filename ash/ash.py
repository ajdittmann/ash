import numpy as np


def ash1d(vals, nbins=20, nshifts=10, weights=None):
  vmax = np.max(vals)
  vmin = np.min(vals)
  L = vmax-vmin
  h = L/(nbins)
  d = h/nshifts

  kgrid = np.linspace(vmin, vmax, nbins*nshifts + 1)
  khist, edges = np.histogram(vals, bins=kgrid, weights=weights)
  values = np.zeros(nbins*nshifts)
  for k in range(nbins*nshifts):
    for i in range(1-nshifts, nshifts):
      if k+i<0: value = 0
      elif k+i>=nbins*nshifts: value = 0
      else: value = khist[k+i]
      values[k] += value*(1 - np.abs(i)/nshifts)

  values = values*nshifts/(h*np.sum(values))
  return kgrid[:-1]+0.5*d, values

def ash2d(xvals, yvals, nbins=20, nshifts=10, weights=None):
  xmax = np.max(xvals)
  xmin = np.min(xvals)
  hx = (xmax-xmin)/nbins
  dx = hx/nshifts

  ymax = np.max(yvals)
  ymin = np.min(yvals)
  hy = (ymax-ymin)/nbins
  dy = hy/(nshifts)

  xgrid = np.linspace(xmin, xmax, nbins*nshifts + 1)
  ygrid = np.linspace(ymin, ymax, nbins*nshifts + 1)

  ashmesh = np.meshgrid(xgrid[:-1]+dx*0.5, ygrid[:-1] + dy*0.5)

  xyhist, xedges, yedges = np.histogram2d(xvals,yvals, bins=(xgrid, ygrid), weights=weights)
  values = np.zeros((nbins*nshifts, nbins*nshifts))
  for kx in range(nbins*nshifts):
    for ky in range(nbins*nshifts):
      for ix in range(1-nshifts, nshifts):
        for iy in range(1-nshifts, nshifts):
          value = 0.0
          if (ix+kx >= 0) and (iy+ky >= 0) and (ix+kx < nshifts*nbins) and (iy+ky < nshifts*nbins): value = xyhist[kx+ix, ky+iy]
          values[kx,ky] += value*(1 - np.abs(ix)/nshifts)*(1 - np.abs(iy)/nshifts)

  values = values/(dx*dy*np.sum(values))
  return ashmesh, values.T

