#!/usr/bin/env python
import re
import ast
from glob import glob

from setuptools import find_packages, setup


classes = """
    Development Status :: 5 - Production/Stable
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

description = ('Downloading 16S data from SRA using Deblur')

long_description = description

version = '0.0.1'

setup(name='GetData',
      version=version,
      description=description,
      long_description=long_description,
      long_description_content_type='text/markdown',
      author="Amnon",
      author_email="jamietmorton@gmail.com",
      maintainer="Amnon",
      maintainer_email="jamietmorton@gmail.com",
      packages=find_packages(),
      scripts=glob('*.py'),
      install_requires=[],
      classifiers=classifiers,
      url="https://github.com/amnona/GetData",
      zip_safe=False)
