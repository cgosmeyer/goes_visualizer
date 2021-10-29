# goes_visualization

## About

GOES-R is a series of NASA/NOAA geostationary satellites that image Earth's weather, oceans, and environment.

### visualize.py: Visualize

Visualization tool to display the performance of each chip in a single GOES-R image assessment.

1. Extracts each chip from selected scene from the "measurement" file and geolocates them by matching to a chip from either of the chip database files "nav_chipdb.csv" or "other_chipdb.csv".
2. Displays the selected GOES image.
3. Optionally can overlay a land/water mask on the GOES image.
4. Optionally plots the selected (closest matching) scene from a single "measurement" file as vector arrows on top of the GOES image.

### mvisualize.py: Multi-Visualize

Plots vector arrows from a scene or an averaged range of (closest-matching) scenes from multiple "measurement" files.

## Installation (for the scientist)

Create a new [anaconda](https://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/20/conda/) environment with Python 3.6.7. For example,

```
conda create -n yourenvname python=3.6.7 anaconda
```

Then activate the environment and do the following to install the required packages.

```
   conda install cartopy
   conda install netCDF4
   conda install xarray=0.15.0
   conda install -c conda-forge metpy
   conda install -c conda-forge basemap
   conda install -c conda-forge pint=0.9
```

Finally install the visualizer package itself in your anaconda environment by the following command.

```
python setup.py
```

## How to Use

See the documentation in `visualizer.py` and `mvisualizer.py`. The command line options can be found by

```
python visualize.py --h

python mvisualize.py --h
```

## Examples