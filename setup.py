# -*- coding: utf-8 -*-

from setuptools import setup

for line in open('pyelectro/__init__.py'):
    if line.startswith("__version__"):
        version = line.split("=")[1].strip()[1:-1]

setup(
    name = "pyelectro",
    version = version,
    packages = ['pyelectro'],
    author = "Michael Vella, Padraig Gleeson",
    author_email = "mv333@cam.ac.uk, p.gleeson@gmail.com",
    description = "A Python library for analysis of electrophysiological data",
    license = "BSD",
    url='https://github.com/NeuralEnsemble/pyelectro',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 2.6',
        'Topic :: Scientific/Engineering']
)
