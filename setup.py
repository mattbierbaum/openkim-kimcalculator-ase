#!/usr/bin/env python

from distutils.core import setup
import kimcalculator

setup(name='kimcalculator',
      version=kimcalculator.__version__,
      description='ASE Calculator using KIM potentials',
      author=kimcalculator.__author__,
      py_modules=['kimcalculator'], 
      license='CDDL'
    )
