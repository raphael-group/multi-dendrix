#!/usr/bin/python

# Import the Multi-Dendrix MCMC C module from the build directory
import sys, os
for filename in os.listdir('multi_dendrix_mcmc/build/'):
	sys.path.append( 'multi_dendrix_mcmc/build/%s' % filename )
from multidendrixmodule import mcmc