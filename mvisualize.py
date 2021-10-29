#!/usr/bin/env python

"""
Visualization tool to display vector arrows from mutliple measurement files on a 
single GOES image or land/water mask.

Authors:

    C.M. Gosmeyer, B. Tan (2021)

"""

import argparse
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime
import glob
import gzip
import matplotlib as mpl
import matplotlib.pyplot as plt
import metpy
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import os
import xarray

from goes_visualizer.utils import round_dt2minute, round_dt2nearest
from goes_visualizer.visualize import Visualizer


def parse_args():
    """ Parses command line arguments.
    """
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--i', dest='image_file', default='',
        action='store', required=True,
        help="L1B image file and path.")
    parser.add_argument('--b', dest='band2extract', default=1,
        action='store', required=True,
        help="The band to plot. E.g., '3'.")
    parser.add_argument('--s', dest='scene2extract', default=None,
        action='store', required=False,
        help="[OPTIONAL] The scene(s) to plot as a single scene 'HHMM' or a range of scenes 'HHMM HHMM'. By default all scenes plotted.")
    parser.add_argument('--vmax', dest='vmax', default=0.4,
        action='store', required=False,
        help="[OPTIONAL] The max to stretch in range [0,1]. Default of 0.4.")
    parser.add_argument('--overlay', dest='overlay_l1b', default=False,
        action='store_true', required=False,
        help="[OPTIONAL] Overlay the L1B image.")
    parser.add_argument('--l', dest='chip_file', default='',
        action='store', required=False,
        help="[OPTIONAL] File containing list of chip names, one chip name per line.")
    parser.add_argument('--m', dest='measurement_files', default=None,
        action='store', required=False,
        help="[OPTIONAL] Measurement files and paths. If this is selected, it overrides the data specification option (--d). \n \
            NOTE: Do not mix satellites and metrics. The satellite and metric labels will be taken from the first file after they are sorted alphabetically.")
    parser.add_argument('--d', dest='dataspec', default=None,
        action='store', required=False,
        help="[OPTIONAL] The satellite, metric, coverage, and the date range of measurement files to search from MMDDYYYY to MMDDYYYY. E.g., 'G17 NAV FULL 07182020 12312020'. \n \
            This will be overriden if a text file listing specific measurement files (--m) is provided.")

    args = parser.parse_args()

    return args

if __name__=="__main__":

    args = parse_args()

    MV = MVisualizer(image_file=args.image_file, 
        band2extract=int(args.band2extract), scene2extract=args.scene2extract,
        vmax=args.vmax, overlay_l1b=args.overlay_l1b, chip_file=args.chip_file,
        measurement_files=args.measurement_files, dataspec=args.dataspec)
    MV.visualize()
