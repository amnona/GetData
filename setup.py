#!/usr/bin/env python
from setuptools import setup


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

description = 'Downloading 16S data from SRA using Deblur'

setup(
    name='GetData',
    version='0.0.1',
    description=description,
    long_description=description,
    long_description_content_type='text/markdown',
    author='Amnon',
    author_email='jamietmorton@gmail.com',
    maintainer='Amnon',
    maintainer_email='jamietmorton@gmail.com',
    packages=['GetData'],
    py_modules=['process_experiment'],
    scripts=['process_experiment.py'],
    install_requires=['numpy'],
    entry_points={
        'console_scripts': [
            'process-experiment=process_experiment:cli',
        ],
    },
    classifiers=classifiers,
    url='https://github.com/amnona/GetData',
    zip_safe=False,
    python_requires='>=3.8',
)
