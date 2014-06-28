#!/usr/bin/python

from distutils.core import setup, Extension

compile_args = ['-g', '-O0']

module1 = Extension('multidendrixmodule',
        extra_compile_args = compile_args,
        sources = ['multidendrixmodule.c'],
        )
setup (name = 'Multi-Dendrix MCMC',
        version = '1.0',
        description = 'Runs Multi-Dendrix MCMC.',
        ext_modules = [module1]
        )
