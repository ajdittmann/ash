# ash
Simple average shifted histograms (ASHs) in python (1D and 2D).

---------------------------
### Basic usage
After installation, the module can be imported using
```
import ash
```
Given some values (`vals`), a number of bins (`nbins`), and number of shifts (`nshifts`), you can calculate a 1D ASH using
```
bins, heights = ash.ash1d(vals, nbins, nshifts)
```
Similarly, given values along two dimensions, (`xvals`) and (`yvals`), you can calculate a 2D ASH using
```
binsGrid, heights2d = ash.ash2d(xvals, yvals, nbins, nshifts)
```
By default, the samples are assumed to be weighted equally. If not, weights for each sample can be passed using the optional `weights` argument. If the data are peridic (such as values of an angle from 0 to 2&pi), periodic boundary conditions can be specified by setting `periodic=True`. 
The upper and lower limits for the histograms can be specified using the `range` argument to `ash1d` (e.g. `range=(xmin, xmax)` or the `xrange` and `yrange` arguments to ash2d. Otherwise, the minimum and maximum data values will be used. 

Once computed, ASHs are straightforward to plot, e.g.
```
import matplotlib.pyplot as plt
plt.plot(bins, heights)
...
plt.pcolormesh(binsGrid[0], binsGrid[1], heights2d)
```
### ash in action 
Sampling from a mixture of three Gaussian distributions:
![Gaussian Mixture](https://github.com/ajdittmann/ash/blob/master/exampleGaussianMixture.png)

Sampling from a mixture of two 2D Gaussian distributions:
![2D Gaussian Mixture](https://github.com/ajdittmann/ash/blob/master/example2DGaussianMixture.png)

### References
* [Scott 1985](http://www.stat.cmu.edu/~rnugent/PCMI2016/papers/ScottASH.pdf)
* [Anderson et al. 2016](https://pubs.acs.org/doi/10.1021/acs.chemmater.6b03430)
### Installation
First, clone the repository using,
`
git clone https://github.com/ajdittmann/ash.git
`.
Then, install with
`
pip install -e .
`
or 
`
python setup.py install
`
.
