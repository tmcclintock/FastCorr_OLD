FastCorr
========
This is a fast implementation of a j_0 spherical hankel transform 
in order to create matter-matter correlation functions from a power spectrum.
The implementation is based off of the Ogata 2005 and its implementation
in the hankel.py package written by Steven Murray. This calculation
takes advantage of the fact that the matter-matter correlation
function is a hankel transformation of the power spectrum
times the j_0 spherical bessel function. This code achieves
a speed up by calculating the roots and
weights more quickly than in the generalized
hankel transformation algorithms, since they can
be written purely as sin(x) and cos(x) calls.

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

Accuracy
---------
The algorith employed by this module contains two variables that 
control the precision: the number of Bessel function roots (N),
and the step size (h). N will control the number of integrand 
evaluations, while h controls how far out in k-space the algorithm
extends to.

The **default** behavior of the module is tuned such that
it reproduces the matter-matter correlation function
from CAMB reproduced by Eduardo Rozo to better than 1%.
While still running less than 0.05 seconds for one thousand
radii.

If you care about probing **large scales** you need to make
h smaller. If you care about smoothness on **all scales** then
you need to need to increase N. Note that the performance
of the algorithm scales linearly with N and is largely 
unaffected by h. Furthermore, it (anecdotally) appears
that when N and h are powers of 2 the algorithm performs
slightly better.