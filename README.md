FastCorr
========
This is a fast implementation of a j_0 spherical hankel transform 
in order to create matter-matter correlation functions from a power spectrum.
The implementation is based off of the Ogata 2005 and its implementation
in the hankel.py package written by Steven Murray. This calculation
takes advantage of the fact that the matter-matter correlation
function is a hankel transformation of the power spectrum
and the j_0 spherical bessel function. This allows the roots and
weights to be computed more quickly than in the generalized
hankel transformation algorithms.

Dependencies
------------
* numpy
* GSL

You must have a path setup to to gsl/include called GSLI and
a path to gsl/lib called GSLL.

Installation
------------
From the FastCorr directory, run
```
python setup.py install
```

And if you care about keeping the root directory clean
```
python setup.py clean
```

Usage
-------
Please look at examples/example.py for an example of how to run.

In order to calculate, for instance, a correlation function all
you need to do is
```python
xi = fastcorr.calc_corr(R,k,P)
```
where R is an array of radii, k is an array of wavenumbers,
and P is an array of the power spectrum.

The code can be used to calculate either the correlation function,
Xi_2 or Xi_4. On an intel ASUS X550L with an intel i5 these take
approximately 0.05 seconds each.

Running the extended_example.py code produces the following

![alt text](https://github.com/tmcclintock/FastCorr/blob/master/figures/figure_1.png)