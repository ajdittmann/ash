import numpy as np


def ash1d(vals, nbins=20, nshifts=10, weights=None):
  vmax = np.max(vals)
  vmin = np.min(vals)
  dx = (vmax-vmin)/(nbins)
  ds = dx/nshifts

  #need to make sure these are correct
  ashgrid = np.linspace(vmin+ds/2, vmax-ds/2, nbins*(nshifts))

  ashden = np.zeros(ashgrid.shape)
  for i in range(nshifts):
    bins = vmin + dx*(1+np.arange(nbins-1)) + ds*(i - nshifts + 1) + ds*(nshifts//2)
    bins = np.append(vmin, bins)
    bins = np.append(bins, vmax)
    hist, edges = np.histogram(vals, bins, density=True, weights = weights)
    repeats = np.ones(nbins, dtype=int)*nshifts
    repeats[0] = i + 1 + nshifts//2
    repeats[-1] = nshifts*2 - i - 1 - nshifts//2
    histdat = np.repeat(hist, repeats, axis=0)
    ashden += histdat
  ashden /= (nshifts)
  return ashgrid, ashden


def ash2d(xvals, yvals, nbins=20, nshifts=10, weights=None):
  xmax = np.max(xvals)
  xmin = np.min(xvals)
  dx = (xmax-xmin)/(nbins)
  dsx = (xmax-xmin)/(nshifts*nbins)


  ymax = np.max(yvals)
  ymin = np.min(yvals)
  dy = (ymax-ymin)/(nbins)
  dsy = dy/(nshifts)

  nsg = (nbins)*(nshifts)

  xgrid = np.linspace(xmin+dsx/2, xmax-dsx/2, nbins*nshifts)
  ygrid = np.linspace(ymin+dsy/2, ymax-dsy/2, nbins*nshifts)
  ashmesh = np.meshgrid(xgrid, ygrid)
  ashmesh[0] = ashmesh[0].T
  ashmesh[1] = ashmesh[1].T
  ashden = np.zeros((nsg, nsg))

  for i in range(nshifts): #x loop
    xbins = xmin + dx*(1+np.arange(nbins-1)) + dsx*(i - nshifts + 1) + dsx*(nshifts//2)
    xbins = np.append(xmin, xbins)
    xbins = np.append(xbins, xmax)
    for j in range(nshifts): #y loop
      ybins = ymin + dy*(1+np.arange(nbins-1)) + dsy*(j - nshifts + 1) + dsy*(nshifts//2)
      ybins = np.append(ymin, ybins)
      ybins = np.append(ybins, ymax)

      hist, xedge, yedge = np.histogram2d(xvals, yvals, bins=(xbins, ybins), density=True, weights=weights)
      hist[(np.isfinite(hist)==0)]=0.0

      xreps = np.ones(nbins, dtype=int)*nshifts
      xreps[0] = i + 1 + nshifts//2
      xreps[-1] = nshifts*2 - i - 1 - nshifts//2

      yreps = np.ones(nbins, dtype=int)*nshifts
      yreps[0] = j + 1 + nshifts//2
      yreps[-1] = nshifts*2 - j - 1 - nshifts//2

      histdat = np.repeat(hist, xreps, axis=0)
      histdat = np.repeat(histdat, yreps, axis=1)

      ashden += histdat
  ashden /= nshifts**2

  return ashmesh, ashden
