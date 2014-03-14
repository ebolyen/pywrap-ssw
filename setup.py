#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2014--, biocore team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__version__ = '0.1.0-dev'

from setuptools import find_packages, setup
from distutils.command.build_py import build_py
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

classes = """
    Development Status :: 1 - Planning
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

long_description = """A simple cython wrapper for `mengyao/Complete-Striped-Smith-Waterman-Library` on github."""

setup(name='pywrap-ssw',
      cmdclass={'build_py': build_py,
                'build_ext': build_ext},
      version=__version__,
      license='BSD',
      description='Cython wrapper for SSW',
      long_description=long_description,
      author="Evan Bolyen",
      author_email="ebolyen@gmail.com",
      maintainer="Evan Bolyen",
      maintainer_email="ebolyen@gmail.com",
      url='https://github.com/biocore/bipy',
      packages=find_packages(),
      install_requires=['numpy >= 1.5.1', 'cython >= 0.20.1'],
      extras_require={'test': ["nose >= 0.10.1", "pep8"],
                      'doc': ["Sphinx >= 1.1", "sphinx-bootstrap-theme"]},
      classifiers=classifiers,
      ext_modules = [Extension("ssw_wrapper", ["ssw_wrapper.pyx", "ssw.c"],
                               include_dirs=[np.get_include()])])

