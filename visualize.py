#!/usr/bin/env python

"""
Visualization tool to display vector arrows from a single measurement files on a GOES
image or land/water mask.

Authors:

    C.M. Gosmeyer, B. Tan (2021)
    
"""

from goes_visualizer.utils import round_dt2minute, round_dt2nearest
from goes_visualizer.visualize import Visualizer

def parse_args():
    """ Parses command line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--m', dest='measurement_file', default='',
        action='store', required=False,
        help="[OPTIONAL] Measurement file and path. If left out, only the image will be plotted.")
    parser.add_argument('--i', dest='image_file', default='',
        action='store', required=True,
        help="L1B image file and path.")
    parser.add_argument('--b', dest='band2extract', default=2,
        action='store', required=True,
        help="The band to plot. E.g., '3'.")
    parser.add_argument('--s', dest='scene2extract', default=None,
        action='store', required=False,
        help="[OPTIONAL] The scene to plot as HHMM-MMDDYYYY. E.g., '1810-07182020'. By default all scenes plotted.")
    parser.add_argument('--vmax', dest='vmax', default=0.4,
        action='store', required=False,
        help="[OPTIONAL] The max to stretch in range [0,1]. Default of 0.4.")
    parser.add_argument('--overlay', dest='overlay_l1b', default=False,
        action='store_true', required=False,
        help="[OPTIONAL] Overlay the L1B image.")
    parser.add_argument('--l', dest='chip_file', default='',
        action='store', required=False,
        help="[OPTIONAL] File containing list of chip names, one chip name per line.")

    args = parser.parse_args()

    return args

if __name__=="__main__":

    args = parse_args()
   
    V = Visualizer(image_file=args.image_file, measurement_file=args.measurement_file, 
         scene2extract=args.scene2extract, band2extract=int(args.band2extract), 
         vmax=args.vmax, overlay_l1b=args.overlay_l1b, chip_file=args.chip_file)
    V.visualize()

