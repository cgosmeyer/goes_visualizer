#!/usr/bin/env python

import os
from setuptools import find_packages
from setuptools import setup

setup(name = 'goes_visualizer',
      description = "Visualization tool to display the performance of each chip in a GOES image assessment.",
      author = 'C. Gosmeyer, B. Tan',
      url = '',
      packages = find_packages(),
      include_package_data=True
)