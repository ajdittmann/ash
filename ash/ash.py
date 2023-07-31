import numpy as np


def ash1d(vals, nbins=20, nshifts=10, weights=None, periodic=False, range=None):
  if range is None:
    vmax = np.max(vals)
    vmin = np.min(vals)
  else:
    vmax = range[1]
    vmin = range[0]
  L = vmax-vmin
  h = L/(nbins)
  d = h/nshifts

  N = nbins*nshifts
  kgrid = np.linspace(vmin, vmax, N + 1)
  khist, edges = np.histogram(vals, bins=kgrid, weights=weights)
  values = np.zeros(N)
  for k in np.arange(N):
    for i in np.arange(1-nshifts, nshifts):
      if periodic:
        value = khist[(k+i)%N]
      else:
        if k+i<0: value = 0
        elif k+i>=N: value = 0
        else: value = khist[k+i]
      values[k] += value*(1 - np.abs(i)/nshifts)

  values = values*nshifts/(h*np.sum(values))
  return kgrid[:-1]+0.5*d, values

def ash2d(xvals, yvals, nbins=20, nshifts=10, weights=None, periodic=False, xrange=None, yrange=None):
  if xrange is not None:
    xmax = xrange[1]
    xmin = xrange[0]
  else:
    xmax = np.max(xvals)
    xmin = np.min(xvals)
  hx = (xmax-xmin)/nbins
  dx = hx/nshifts

  if yrange is not None:
    ymax = yrange[1]
    ymin = yrange[0]
  else:
    ymax = np.max(yvals)
    ymin = np.min(yvals)
  hy = (ymax-ymin)/nbins
  dy = hy/(nshifts)

  N = nbins*nshifts
  xgrid = np.linspace(xmin, xmax, N + 1)
  ygrid = np.linspace(ymin, ymax, N + 1)

  ashmesh = np.meshgrid(xgrid[:-1]+dx*0.5, ygrid[:-1] + dy*0.5)

  xyhist, xedges, yedges = np.histogram2d(xvals,yvals, bins=(xgrid, ygrid), weights=weights)
  values = np.zeros((N, N))
  for kx in range(N):
    for ky in range(N):
      for ix in range(1-nshifts, nshifts):
        for iy in range(1-nshifts, nshifts):
          value = 0.0
          if periodic: value = xyhist[(kx+ix)%N,(ky+iy)%N]
          else:
            if (ix+kx >= 0) and (iy+ky >= 0) and (ix+kx < N) and (iy+ky < N): value = xyhist[kx+ix, ky+iy]
          values[kx,ky] += value*(1 - np.abs(ix)/nshifts)*(1 - np.abs(iy)/nshifts)

  values = values/(dx*dy*np.sum(values))
  return ashmesh, values.T

