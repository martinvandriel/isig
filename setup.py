#!/bin/env python
# -*- coding: utf-8 -*-
u"""
Instaseis Sismosphere Input Generator.

Produces meshes and relevant input to create instaseis databases using normal
modes.

:copyright:
    Martin van Driel, Martin@vanDriel.de
:license:
    None
"""
from setuptools import find_packages, setup
from setuptools.extension import Extension

import inspect
import os
import sys


# Import the version string.
path = os.path.join(os.path.abspath(os.path.dirname(inspect.getfile(
    inspect.currentframe()))), "isig")
sys.path.insert(0, path)
from version import get_git_version


DOCSTRING = __doc__.strip().split("\n")


def get_package_data():
    """
    Returns a list of all files needed for the installation relativ to the
    "isig" subfolder.
    """
    filenames = []
    # The lasif root dir.
    root_dir = os.path.join(os.path.dirname(os.path.abspath(
        inspect.getfile(inspect.currentframe()))), "isig")
    # Recursively include all files in these folders:
    folders = [os.path.join(root_dir, "tests", "data")]
    for folder in folders:
        for directory, _, files in os.walk(folder):
            for filename in files:
                # Exclude hidden files.
                if filename.startswith("."):
                    continue
                filenames.append(os.path.relpath(
                    os.path.join(directory, filename),
                    root_dir))
    filenames.append("RELEASE-VERSION")
    return filenames


# Hack to prevent build_ext from trying to append "init" to the export symbols.
class finallist(list):
    def append(self, object):
        return


class MyExtension(Extension):
    def __init__(self, *args, **kwargs):
        Extension.__init__(self, *args, **kwargs)
        self.export_symbols = finallist(self.export_symbols)



INSTALL_REQUIRES = ["numpy",
                    "pymesher",
                    "future",
                    "flake8>=2",
                    "pytest"]

# Add argparse and ordereddict for Python 2.6. Both are standard library
# packages for Python >= 2.7.
if sys.version_info[:2] == (2, 6):
    INSTALL_REQUIRES.extend(["argparse", "ordereddict"])
# Add mock for Python 2.x. Starting with Python 3 it is part of the standard
# library.
if sys.version_info[0] == 2:
    INSTALL_REQUIRES.append("mock")

setup_config = dict(
    name="isig",
    version=get_git_version(),
    description=DOCSTRING[0],
    long_description="\n".join(DOCSTRING[2:]),
    author=u"Martin van Driel",
    author_email="Martin@vanDriel.de",
    url="http://github.com/martinvandriel/isig",
    packages=find_packages(),
    license="",
    platforms="OS Independent",
    install_requires=INSTALL_REQUIRES,
    # this is needed for "pip install instaseis==dev"
    # download_url=("https://github.com/krischer/instaseis/zipball/master"
    #               "#egg=instaseis=dev"),
    classifiers=[
        # complete classifier list:
        # http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: MacOS",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.6",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics"
        ],
)

if __name__ == "__main__":
    setup(**setup_config)
