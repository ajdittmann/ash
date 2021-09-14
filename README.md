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
Optional arguments include `weights`, which defaults to equally-weighted samples, and `padVals` which can be set to `True` to pad arrays with the first and last histogram values rather than zeroes. 
These can then be plotted, e.g.
```
import matplotlib.pyplot as plt
plt.plot(bins, heights)
...
plt.pcolormesh(binsGrid[0], binsGrid[1], heights2d)
```
### ASH in action 
Sampling from a mixture of three Gaussian distributions:
![Gaussian Mixture](https://github.com/ajdittmann/ash/blob/master/exampleGaussianMixture.png)

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
python setup.py install
`
or 
`
pip install -e .
`.
