# FastCorr
This is a fast implementation of a j_0 spherical hankel transform 
in order to create matter-matter correlation functions from a power spectrum.
The implementation is based off of the Ogata 2005 and its implementation
in the hankel.py package written by Steven Murray. This calculation
takes advantage of the fact that the matter-matter correlation
function is a hankel transformation of the power spectrum
and the j_0 spherical bessel function. This allows the roots and
weights to be computed more quickly than in the generalized
hankel transformation algorithms.

# Dependencies
numpy
GSL

# Installation
From the FastCorr directory, run
python setup.py install

And if you care about keeping the root directory clean
python setup.py clean

# Running
Please look at examples/example.py for an example of how to run.