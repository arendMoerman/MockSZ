import os
import pathlib

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext as build_ext_orig

setup(
    name='MockSZ',
    license="MIT",
    version='0.0.1',
    author="Arend Moerman",
    install_requires = ["numpy", "matplotlib", "scipy"],
    package_dir = {'': 'src'},
    packages=['MockSZ'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    python_requires='>=3.8',
)

