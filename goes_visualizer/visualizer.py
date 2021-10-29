
"""
Visualizer classes for GOES-R series.

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
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import os
import xarray


class Visualizer(object):
    def __init__(self, image_file, measurement_file, band2extract, scene2extract=None, 
        vmax=0.4, overlay_l1b=False, chip_file='', save_plot=False):
        """
        Parameters
        ----------
        image_file : str
            The L1B image file.
        measurement_file : str
            The measurement file.
        band2extract : int
            The band to extract.
        scene2extract : str
            The scene to extract. E.g., 1810-07182020, meaning scene falling during 
            18:10 on 07/18/2021.
        vmax : int
            The max to stetch. Larger->less contrast.
        overlay_l1b : {True, False}
            Whether to overlay the L1B image. By default shows the generaric 
            land/ocean map.
        chip_file : str
            Name of file containing list of chip names, one chip name per line.
        save_plot : {True, False}
            Whether to save the plot or just show it.
        """
        self.image_file = image_file
        self.measurement_file = measurement_file
        self.band2extract = band2extract
        self.scene2extract = scene2extract
        self.vmax = float(vmax)
        self.overlay_l1b = overlay_l1b
        self.chip_file = chip_file
        self.save_plot = save_plot
        self.scene = ''
        self.nir_flg = False

        if self.measurement_file != '':
            # Extract satellite name
            self.sat = self.measurement_file.split('/')[-1].split('_')[0]
            # Extract the metric type
            self.metric = self.measurement_file.split('/')[-1].split('_')[1]
            # Find coverage
            if 'CONUS' in self.measurement_file:
                self.coverage = 'CONUS'
            else:
                self.coverage = 'FULL'
        else:
            self.sat = ''
            self.metric = ''
            self.coverage = ''

        # Build band name
        if self.band2extract/10 < 1:
            self.band = '0' + str(self.band2extract)
        else:
            self.band = str(self.band2extract)

    def extract_geoloc(self):
        """ Extract the geolocation information for the band of interest from the 
        appropriate Chip DB file.
        """
        # Extract the input date and time 
        if self.scene2extract != None:
            date = datetime.datetime.strptime(self.scene2extract.split('-')[1], '%m%d%Y')
            time = datetime.datetime.strptime(self.scene2extract.split('-')[0], '%H%M')
            date_time = datetime.datetime.strptime(self.scene2extract, '%H%M-%m%d%Y')
        else:
            date = 0
            time = 1

        # If metric is BBR, need unzip the measurements file
        if self.metric == 'BBR':
            with gzip.open(self.measurement_file) as f:
                measure_df = pd.read_csv(self.measurement_file)
        else:
            measure_df = pd.read_csv(self.measurement_file)

        # Create a datetime column.
        activity_date = np.array(measure_df['ACTIVITY_DATE1'])
        activity_time = np.array(measure_df['ACTIVITY_TIME_1'])
        measure_df['DATETIME'] = [datetime.datetime.strptime(activity_date[j]+'_'+activity_time[j], 
            '%m-%d-%Y_%H:%M:%S') for j in range(len(activity_time))]

        # Round the user-inputted time to nearest scene (date/time) in measurement file
        if self.scene2extract != None:
            t = pd.DataFrame(measure_df, columns = ['DATETIME'])
            t_df = pd.DataFrame.drop_duplicates(t)
            t_df = t_df.reset_index()
            df_sort = t_df.iloc[(t_df['DATETIME']-date_time).abs().argsort()[:1]]

            self.scene = df_sort['DATETIME'].iloc[0].strftime('%H:%M')

            # Issue warning message if the requested scene is not in range of file.
            # (in that case, extract either first or last scene)
            if not(date_time >= measure_df['DATETIME'].iloc[0] and date_time <= measure_df['DATETIME'].iloc[-1]):
                print("--WARNING: Requested scene ({}) falls outside measurement file. Using closest scene ({}) instead.--"\
                    .format(self.scene2extract, df_sort['DATETIME'].iloc[0].strftime('%H%M-%m%d%Y')))
                # Set "not in range" flag
                self.nir_flg = True
            else:
                print("--Plotting closest scene in file ({})--"\
                    .format(df_sort['DATETIME'].iloc[0].strftime('%m/%d/%Y %H:%M')))

            # Extract the band of interest and scene (date/time) of interest.
            measure_df = measure_df[measure_df['BAND_NUM'] == self.band2extract]\
                [measure_df['DATETIME'] == df_sort['DATETIME'].iloc[0]]
        else:
            self.scene = 'All'
            # Extract the band of interest.
            measure_df = measure_df[measure_df['BAND_NUM'] == self.band2extract]

        print("Scene: ", self.scene)

        # Read the Chip DB file, depending on the metric
        exe_path = os.path.dirname(os.path.realpath(__file__))
        if self.metric == 'NAV':
            chipdb_df = pd.read_csv(os.path.join(exe_path, 'data', 'other_chipdb.csv'))
            # Remove all columns from chip db except for LANDMARK_S24, ORIGLAT_R, ORIGLON_R.
            chipdb_new = chipdb_df[['LANDMARK_S24', 'NEWLAT_R', 'NEWLON_R']].copy()
            # Rename columns
            chipdb_new = chipdb_new.rename(columns={"LANDMARK_S24":"chip", "NEWLAT_R":"lat", "NEWLON_R":"lon"})
        else:
            chipdb_df = pd.read_csv(os.path.join(exe_path, 'data', 'nav_chipdb.csv'))
            # Remove all columns from chip db except for LANDMARK_S24, ORIGLAT_R, ORIGLON_R.
            chipdb_new = chipdb_df[['name_S24', 'lat_R', 'lon_R']].copy()
            # Rename columns
            chipdb_new = chipdb_new.rename(columns={"name_S24":"chip", "lat_R":"lat", "lon_R":"lon"})

        # Remove all duplicate rows from Chip DB.
        chipdb_new = chipdb_new.drop_duplicates()
        chipdb_new = chipdb_new.reset_index()

        # Pull out columns to speed up search in for loop
        origlat_r = chipdb_new["lat"]
        origlon_r = chipdb_new["lon"]
        landmark_s24 = np.array(chipdb_new["chip"])
        chip_name = np.array(measure_df['CHIP_NAME'])

        # Match chip names from the Chip DB file to those in measurements file in order to match rows in the
        # measurements file to latitudes and longitudes.
        lat_arr = []
        lon_arr = []

        # Extract chip names, if specified
        if self.chip_file != '':
            chip_list = self.extract_chips()
            print("--Only user-specified chips will be plotted: {}--".format(chip_list))
        else:
            chip_list = chip_name

        # Match chip name from measurements file to chip in Chip DB file in order to
        # extract the corresponding lat/lon.
        # If user specifies a chip list, retain only those chips.
        for i in range(len(measure_df)):
            if (chip_name[i] in landmark_s24) and (chip_name[i] in chip_list):
                lat = np.array(origlat_r[chipdb_new["chip"] == chip_name[i]])
                lon = np.array(origlon_r[chipdb_new["chip"] == chip_name[i]])
                if len(lat) > 0:
                    lat_arr.append(lat[0])
                    lon_arr.append(lon[0])
                else:
                    lat_arr.append(0)
                    lon_arr.append(0)
            else:
                lat_arr.append(0)
                lon_arr.append(0)

        # Append lat and lon arrays to measurement dataframe
        measure_df['Lat'] = lat_arr
        measure_df['Lon'] = lon_arr

        measure_df = measure_df[(measure_df["Lat"] != 0)]

        print("Number of vectors: ", len(measure_df["Lat"]))

        return measure_df

    def extract_chips(self):
        """
        """
        chip_list = []
        with open(self.chip_file) as f:
            for line in f:
                chip_list.append(line.strip('\n'))
        return chip_list

    def visualize(self):
        """ Visualize the offsets as vector field on either L1B map or generic
        world map.
        """
        # Remove path to get just filename for parsing purposes
        image_file = self.image_file.split('/')[-1]

        # Extract mode
        mode = image_file.split('_')[1].split('-')[3][:2]

        # Extract geographic coverage
        # Based on coverage, set the orientation for the plot colorbar
        coverage = image_file.split('-')[2].strip('Rad')
        if coverage == 'C':
            coverage = 'CONUS'
            orientation = 'horizontal'
        elif coverage == 'F':
            coverage = 'FULL'
            orientation = 'vertical'
        else:
            ## Say all others should be treated as "FULL" would, for now
            coverage = 'FULL'
            orientation = 'vertical'           

        # Extract satellite from image
        sat = image_file.split('_')[2]

        # Search for the Scan start in the file name
        start = (image_file[image_file.find("s")+1:image_file.find("_e")])
        start_formatted = start[0:4] + " Day " + start[4:7] + " - " + start[7:9] + ":" + \
            start[9:11] + ":" + start[11:13] + "." + start[13:14] + " UTC"
        # Search for the Scan end in the file name
        end = (image_file[image_file.find("e")+1:image_file.find("_c")])
        end_formatted = end[0:4] + " Day " + end[4:7] + " - " + end[7:9] + ":" + end[9:11] + \
            ":" + end[11:13] + "." + end[13:14] + " UTC"

        # Open the file using the NetCDF4 library
        nc = Dataset(self.image_file)

        # Determine the lon_0
        geo_extent = nc.variables['geospatial_lat_lon_extent']
        lon_0 = geo_extent.geospatial_lon_center
        lat_0 = 0

        print("Measurement file satellite: ", self.sat)
        print("Measurement file metric: ", self.metric)
        print("Measurement file band: ", self.band)
        print("Measurement file coverage: ", self.coverage)

        print("Image satellite: ", sat)
        print("Image coverage: ", coverage)
        print("Image start: ", start)
        print("Image end: ", end)

        # Import the measurements dataframe
        if self.measurement_file != '':
            measure_df = self.extract_geoloc()
        else:
            print("No measurement file supplied.")
        # Extract the Brightness Temperature values from the NetCDF
        if 'Rad' in image_file:
            image_kwd = 'Rad'
        elif 'ACMF' in image_file:
            image_kwd = 'BCM'

        data = nc.variables[image_kwd][:]

        geos = ccrs.Geostationary(central_longitude=lon_0, satellite_height=35786023.0, sweep_axis='x')

        # Start figure
        fig=plt.figure(figsize=(12, 8))
        ax=fig.add_axes([0.1,0.1,0.8,0.8], projection=geos)

        open_image = xarray.open_dataset(self.image_file)
        image_data = open_image.metpy.parse_cf(image_kwd)
        image_x = image_data.x
        image_y = image_data.y

        # Set the axis bounds.
        if coverage == 'CONUS':
            ax.set_extent([image_x.min(), image_x.max(), image_y.min(), image_y.max()], crs=geos)
            info_text='cyan'
        elif coverage == 'FULL':
            ax.set_global()
            info_text='k'

        # Overlay the L1B data
        if self.overlay_l1b:
            # De-normalize the vmax from range [0,1] to natural range
            min_range = float(nc.variables[image_kwd].valid_range[0])
            max_range = float(nc.variables[image_kwd].valid_range[1])
            vmax = self.vmax*(max_range - min_range)

            if coverage == 'CONUS':
                vmax = vmax/3.5

            # Plot L1B data
            # Note: Increasing vmax lowers contrast. Vmax=small->black; Vmax=large->white
            ax.imshow(open_image[image_kwd][:], origin='upper', cmap='gray', transform=geos, vmax=vmax,
                extent=(image_x.min(), image_x.max(), image_y.min(), image_y.max()))

            # Draw coatlines, country borders, lakes, and grid
            # See https://scitools.org.uk/cartopy/docs/v0.14/matplotlib/feature_interface.html
            ax.coastlines(linewidth=0.9, linestyle='solid', color='green')
            ax.add_feature(cfeature.BORDERS, linewidth=0.9, linestyle='solid', 
                facecolor='none', edgecolor='green')
            ax.add_feature(cfeature.LAKES, linewidth=0.9, linestyle='solid', 
                facecolor='none', edgecolor='green')
            ax.gridlines(linewidth=0.3, color='white')

        # If no image file selected to overlay, draw ocean and land
        else:
            ax.stock_img()

            # Draw the coastlines, countries, parallels and meridians
            ax.coastlines(linewidth=0.9, linestyle='solid', color='black')
            ax.add_feature(cfeature.BORDERS, linewidth=0.9, linestyle='solid', 
                facecolor='none', edgecolor='black')
            ax.add_feature(cfeature.LAKES, linewidth=0.9, linestyle='solid', 
                facecolor='skyblue', edgecolor='black')
            ax.add_feature(cfeature.RIVERS, linewidth=0.9, linestyle='solid', 
                facecolor='none', edgecolor='skyblue')
            ax.gridlines(linewidth=0.3, color='white')

        # Add a title to the plot
        plt.title(self.sat + " ABI L1B Band " + self.band + " Scene " + \
            self.scene + " Metric " + self.metric + "\n" + coverage + \
            " Scan from " + start_formatted + " to " + end_formatted)

        # Read some variables from the NetCDF header in order to use it in the plot
        center = str(geo_extent.geospatial_lon_center)
        west = str(geo_extent.geospatial_westbound_longitude)
        east = str(geo_extent.geospatial_eastbound_longitude)
        north = str(geo_extent.geospatial_northbound_latitude)
        south = str(geo_extent.geospatial_southbound_latitude)

        # Close netCDF file when finished
        nc.close()
        nc = None  

        # Put the information retrieved from the header in the final image
        plt.text(0.01, 0.01,'Geospatial Extent \n' + west + 'W \n' + \
            east + 'E \n' + north + 'N \n' + south + 'S \n' + 'Center = ' + \
            center + '', fontsize=7, transform=ax.transAxes, color=info_text)

        # Start time to be printed large on image
        start_time = start[7:9] + ":" + start[9:11] + ":" + start[11:13]
        plt.text(0.78, 0.88, start_time, fontsize=24, transform=ax.transAxes, color='red')

        if self.nir_flg:
            plt.text(0.01, 0.94,"WARNING: Selected scene \n{} \nnot in measurement file"\
                .format(self.scene2extract), color='red', fontsize=8, transform=ax.transAxes)

        if self.measurement_file != '':
            # Project the coordinates from measurements dataframe
            x = np.array(measure_df['Lon'])
            y = np.array(measure_df['Lat'])

            # Generate the vectors
            delta_ew = np.array(measure_df['DELTA_EW'])
            delta_ns = np.array(measure_df['DELTA_NS'])
            # Calculate magnitudes so can colorize
            mag = (delta_ew**2 + delta_ns**2)**(0.5)
            # Normalize the arrows
            delta_ew_norm = delta_ew/np.sqrt(delta_ew**2 + delta_ns**2)
            delta_ns_norm = delta_ns/np.sqrt(delta_ew**2 + delta_ns**2)

            # Draw the vectors
            ax.quiver(x, y, delta_ew_norm, delta_ns_norm, mag, width=0.003, 
                cmap='jet', transform=ccrs.PlateCarree())

            # Insert the colorbar
            # Source: https://www.geeksforgeeks.org/matplotlib-pyplot-colorbar-function-in-python/
            norm = mpl.colors.Normalize(vmin=min(mag), vmax=max(mag))
            cmap = plt.get_cmap('jet')
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            plt.colorbar(sm, orientation=orientation, label='Shift Magnitude, urad')

        if 'ACMF' in image_file:
            # Plot the chips as red dots.
            exe_path = os.path.dirname(os.path.realpath(__file__))
            chipdb_df = pd.read_csv(os.path.join(exe_path, 'data', 'nav_chipdb.csv'))
            # Remove all columns from MutliSpecDB except for LANDMARK_S24, ORIGLAT_R, ORIGLON_R.
            chipdb_new = chipdb_df[['LANDMARK_S24', 'NEWLAT_R', 'NEWLON_R']].copy()
            # Rename columns
            chipdb_new = chipdb_new.rename(columns={"LANDMARK_S24":"chip", "NEWLAT_R":"lat", "NEWLON_R":"lon"})
            chipdb_new = chipdb_new.drop_duplicates()
            chipdb_new = chipdb_new.reset_index()

            plt.plot(chipdb_new["lon"], chipdb_new["lat"], color='red', marker='o', 
                linestyle='None', markersize=1.5, transform=ccrs.PlateCarree())

        # Show or save the plot
        if save_plot:
            plt.savefig('vplot.png', bbox_inches='tight')
        else:
            plt.show()
        plt.close()


class MVisualizer(Visualizer):
    def __init__(self, image_file, band2extract, scene2extract,
        vmax, overlay_l1b, chip_file, save_plot, measurement_files, dataspec):
        """
        Parameters
        ----------
        image_file : str
            The L1B image file.
        band2extract : int
            The band to extract.
        vmax : int
            The max to stetch. Larger->less contrast.
        overlay_l1b : {True, False}
            Whether to overlay the L1B image. By default shows the generaric 
            land/ocean map.
        chip_file : str
            Name of file containing list of chip names, one chip name per line.
        save_plot : {True, False}
            Whether to save the plot or just show it.
        measurement_files : str
            File containing list (one per line) of measurement file names.
        dataspec : str
            The range of dates in which to search for measurement files. 
        """
        measurement_file = None
        super().__init__(image_file, measurement_file, band2extract, scene2extract, 
            vmax, overlay_l1b, chip_file, save_plot)

        # Build band name
        if self.band2extract/10 < 1:
            self.band = '0' + str(self.band2extract)
        else:
            self.band = str(self.band2extract)

        if measurement_files != None:
            self.measurement_files = self.extract_from_file(measurement_files)
            # Sort so that files are in order of datetime (unless files are in different locations...)
            self.measurement_files = sorted(self.measurement_files)

            print("Measurement files: ", self.measurement_files)

            # Use the first file to determine the satellite and metric and start date
            # Use the last file to determien end date
            self.sat = self.measurement_files[0].split('/')[-1].split('_')[0]
            self.metric = self.measurement_files[0].split('/')[-1].split('_')[1]
            self.start_range = datetime.datetime.strptime(self.measurement_files[0]\
                .split('/')[-1].split('_')[4].split('.')[0] \
                + '-' + self.measurement_files[0].split('/')[-1].split('_')[3], '%j-%Y')
            self.end_range = datetime.datetime.strptime(self.measurement_files[-1]\
                .split('/')[-1].split('_')[4].split('.')[0] \
                + '-' + self.measurement_files[-1].split('/')[-1].split('_')[3], '%j-%Y')

            if 'CONUS' in self.measurement_files[0]:
                self.coverage = 'CONUS'
            else:
                self.coverage = 'FULL'

            print("Measurement file satellite: ", self.sat)
            print("Measurement file metric: ", self.metric)
            print("Measurement file band:", self.band)
            print("Measurement file coverage: ", self.coverage)       
            print("Measurement file start date: ", self.start_range)
            print("Measurement file end date: ", self.end_range)

        elif dataspec != None:
            print("dataspec: ", dataspec)
            try:
                self.sat = dataspec.split(' ')[0].upper()
                self.metric = dataspec.split(' ')[1].upper()
                self.coverage = dataspec.split(' ')[2].upper()
                self.start_range = datetime.datetime.strptime(dataspec.split(' ')[3], '%m%d%Y') 
                self.end_range = datetime.datetime.strptime(dataspec.split(' ')[4], '%m%d%Y') 

                self.measurement_files = self.searchforfiles()
                print("Measurement files: ", self.measurement_files)

                if self.measurement_files == []:
                    print("Error! No measurement files found.")
                else:
                    print("Measurement file satellite: ", self.sat)
                    print("Measurement file metric: ", self.metric)
                    print("Measurement file band:", self.band)
                    print("Measurement file coverage: ", self.coverage)
                    print("Measurement file start date: ", self.start_range)
                    print("Measurement file end date: ", self.end_range)
            except:
                print("Error! Data specification needs to be in format 'AAA BBB CCC MMDDYYYY MMDDYYYY', where AAA can be G16 or G17; BBB can be FFR, NAV, BBR or WIFR; and CCC can be FUL or CON")
        else:
            print("Error! Please provide either file listing measurement files (--m) or a data specification (satellite, metric, coverage, and date range) to search for measurement files (--d).")

    def extract_geoloc(self, measurement_file):
        """ Extract the geolocation information for the band of interest from the 
        appropriate Chip DB file.
        """
        # Extract the input date and time 
        if self.scene2extract != None:
            print("User-requested starting scene: ", self.scene2extract.split(' ')[0])
            print("User-requested ending scene: ", self.scene2extract.split(' ')[-1])
            start_time = datetime.datetime.strptime(self.scene2extract.split(' ')[0], '%H%M')
            end_time = datetime.datetime.strptime(self.scene2extract.split(' ')[-1], '%H%M')

        # Check if file nseeds to be unzipped
        if 'gz' in measurement_file:
            with gzip.open(measurement_file) as f:
                measure_df = pd.read_csv(measurement_file)
        else:
            measure_df = pd.read_csv(measurement_file)

        # Create a datetime column.
        activity_date = np.array(measure_df['ACTIVITY_DATE1'])
        activity_time = np.array(measure_df['ACTIVITY_TIME_1'])

        measure_df['DATETIME'] = [datetime.datetime.strptime(activity_time[j], '%H:%M:%S') for j in range(len(activity_time))] 

        # Round the user-inputted time to nearest scene (date/time) in measurement file
        if self.scene2extract != None and start_time != end_time:
            t_df = pd.DataFrame(measure_df, columns = ['ACTIVITY_TIME_1'])
            t_df['DATETIME'] = [datetime.datetime.strptime(i, '%H:%M:%S') for i in t_df['ACTIVITY_TIME_1']]
            time_sorted = t_df.sort_values(by='DATETIME')

            # Find the start and ending date and then form a datetime in order to get the range the user wants
            df_sort_start = t_df.iloc[(t_df['DATETIME']-start_time).abs().argsort()[:1]]
            df_sort_end = t_df.iloc[(t_df['DATETIME']-end_time).abs().argsort()[:1]]

            self.scene = df_sort_start['ACTIVITY_TIME_1'].iloc[0] + ' to ' + df_sort_end['ACTIVITY_TIME_1'].iloc[0]

            # Extract the band of interest and scene (date/time) of interest.
            print("--WARNING using closest found scenes as the bounds {}.".format(self.scene))
            measure_df = measure_df[measure_df['BAND_NUM'] == self.band2extract]\
                [(measure_df['DATETIME'] >= df_sort_start['DATETIME'].iloc[0]) & (measure_df['DATETIME'] <= df_sort_end['DATETIME'].iloc[0])]

        elif self.scene2extract != None and start_time == end_time:
            t = pd.DataFrame(measure_df, columns = ['DATETIME'])
            t_df = pd.DataFrame.drop_duplicates(t)
            t_df = t_df.reset_index()
            df_sort = t_df.iloc[(t_df['DATETIME']-start_time).abs().argsort()[:1]]

            self.scene = df_sort['DATETIME'].iloc[0].strftime('%H:%M')

            # Issue warning message if the requested scene is not in range of file.
            # (in that case, extract either first or last scene)
            if not(start_time >= measure_df['DATETIME'].iloc[0] and start_time <= measure_df['DATETIME'].iloc[-1]):
                print("--WARNING: Requested scene ({}) falls outside measurement file. Using closest scene ({}) instead.--"\
                    .format(self.scene2extract, df_sort['DATETIME'].iloc[0].strftime('%H%M-%m%d%Y')))
                # Set "not in range" flag
                self.nir_flg = True
            else:
                print("--Plotting closest scene in file ({})--".format(df_sort['DATETIME'].iloc[0].strftime('%m/%d/%Y %H:%M')))

            # Extract the band of interest and scene (date/time) of interest.
            measure_df = measure_df[measure_df['BAND_NUM'] == self.band2extract]\
                [measure_df['DATETIME'] == df_sort['DATETIME'].iloc[0]]
        else:
            self.scene = 'All'
            # Extract the band of interest.
            measure_df = measure_df[measure_df['BAND_NUM'] == self.band2extract]

        print("Scene: ", self.scene)
        #print("measure_df: ", measure_df)

        # Read the Chip DB file, depending on the metric
        exe_path = os.path.dirname(os.path.realpath(__file__))
        if self.metric == 'NAV':
            chipdb_df = pd.read_csv(os.path.join(exe_path, 'MultiSpecDB_Jan20_2018_v1_2018_Feb13_150634.csv'))
            # Remove all columns from MutliSpecDB except for LANDMARK_S24, ORIGLAT_R, ORIGLON_R.
            chipdb_new = chipdb_df[['LANDMARK_S24', 'NEWLAT_R', 'NEWLON_R']].copy()
            # Rename columns
            chipdb_new = chipdb_new.rename(columns={"LANDMARK_S24":"chip", "NEWLAT_R":"lat", "NEWLON_R":"lon"})
        else:
            chipdb_df = pd.read_csv(os.path.join(exe_path, 'EvalLoc_2018_Feb26.csv'))
            # Remove all columns from MutliSpecDB except for LANDMARK_S24, ORIGLAT_R, ORIGLON_R.
            chipdb_new = chipdb_df[['name_S24', 'lat_R', 'lon_R']].copy()
            # Rename columns
            chipdb_new = chipdb_new.rename(columns={"name_S24":"chip", "lat_R":"lat", "lon_R":"lon"})

        # Remove all duplicate rows from Chip DB.
        chipdb_new = chipdb_new.drop_duplicates()
        chipdb_new = chipdb_new.reset_index()

        # Pull out columns to speed up search in for loop
        origlat_r = chipdb_new["lat"]
        origlon_r = chipdb_new["lon"]
        landmark_s24 = np.array(chipdb_new["chip"])
        chip_name = np.array(measure_df['CHIP_NAME'])

        # Match chip names from the Chip DB file to those in measurements file in order to match rows in the
        # measurements file to latitudes and longitudes.
        lat_arr = []
        lon_arr = []

        # Extract chip names, if specified
        if self.chip_file != '':
            chip_list = self.extract_from_file(self.chip_file)
            print("--Only user-specified chips will be plotted: {}--".format(chip_list))
        else:
            chip_list = chip_name

        # Match chip name from measurements file to chip in Chip DB file in order to
        # extract the corresponding lat/lon.
        # If user specifies a chip list, retain only those chips.
        for i in range(len(measure_df)):
            if (chip_name[i] in landmark_s24) and (chip_name[i] in chip_list):
                lat = np.array(origlat_r[chipdb_new["chip"] == chip_name[i]])
                lon = np.array(origlon_r[chipdb_new["chip"] == chip_name[i]])
                if len(lat) > 0:
                    lat_arr.append(lat[0])
                    lon_arr.append(lon[0])
                else:
                    lat_arr.append(0)
                    lon_arr.append(0)
            else:
                lat_arr.append(0)
                lon_arr.append(0)

        # Append lat and lon arrays to measurement dataframe
        measure_df['Lat'] = lat_arr
        measure_df['Lon'] = lon_arr

        measure_df = measure_df[(measure_df["Lat"] != 0)]

        return measure_df

    def extract_from_file(self, filename):
        """
        """
        item_list = []
        with open(filename) as f:
            for line in f:
                item_list.append(line.strip('\n'))
        return item_list

    def searchforfiles(self):
        """ Creates list of measurement files in the given range for satellite and metric.
        """
        measurement_files = []

        # Calculate how many days are between first and last.
        ndates = (self.end_range - self.start_range).days

        # Loop over number of days and add day each iteration to start date and then form the filename to glob for it
        for d in range(ndates):
            # Add 'd' days to the start_date
            newdate = self.start_range + datetime.timedelta(d) 
            # Convert the new date to year and day of year
            year = newdate.year
            startofyear = datetime.datetime(year=year, month=1, day=1)
            days_since_startofyear = (newdate - startofyear).days
            # Use year and day of year to construct a string to search for measurement file
            if 'F' in self.coverage:
                search_path = ''
            elif 'C' in self.coverage:
                search_path = ''           
            search_string = '{}_{}_measurements_{}_{}.*.csv*'.format(self.sat, 
                self.metric, year, days_since_startofyear)
            # Append the found file to list of files
            measurement_files += glob.glob(os.path.join(search_path, search_string))

        return measurement_files

    def visualize(self):
        """ Visualize the offsets as vector field on either L1B map or generic
        world map.
        """
        measure_df = self.build_measurement_df()

        # Remove path to get just filename for parsing purposes
        image_file = self.image_file.split('/')[-1]

        # Extract mode
        mode = image_file.split('_')[1].split('-')[3][:2]

        # Extract geographic coverage from image
        # Based on coverage, set the orientation for the plot colorbar
        coverage = image_file.split('-')[2].strip('Rad')
        if coverage == 'C':
            coverage = 'CONUS'
            orientation = 'horizontal'
        elif coverage == 'F':
            coverage = 'FULL'
            orientation = 'vertical'

        # Extract satellite from image
        sat = image_file.split('_')[2]  

        # Search for the Scan start in the file name
        start = (image_file[image_file.find("s")+1:image_file.find("_e")])
        start_formatted = start[0:4] + " Day " + start[4:7] + " - " + start[7:9] + ":" + \
            start[9:11] + ":" + start[11:13] + "." + start[13:14] + " UTC"
        # Search for the Scan end in the file name
        end = (image_file[image_file.find("e")+1:image_file.find("_c")])
        end_formatted = end[0:4] + " Day " + end[4:7] + " - " + end[7:9] + ":" + end[9:11] + \
            ":" + end[11:13] + "." + end[13:14] + " UTC"

        # Open the file using the NetCDF4 library
        nc = Dataset(self.image_file)

        # Determine the lon_0
        geo_extent = nc.variables['geospatial_lat_lon_extent']
        lon_0 = geo_extent.geospatial_lon_center
        lat_0 = 0

        print("Image satellite: ", sat)
        print("Image coverage: ", coverage)
        print("Image start: ", start)
        print("Image end: ", end)

        # Extract the Brightness Temperature values from the NetCDF
        data = nc.variables['Rad'][:]

        geos = ccrs.Geostationary(central_longitude=lon_0, satellite_height=35786023.0, sweep_axis='x')

        # Start figure
        fig=plt.figure(figsize=(12, 8))
        ax=fig.add_axes([0.1,0.1,0.8,0.8], projection=geos)

        open_image = xarray.open_dataset(self.image_file)
        image_data = open_image.metpy.parse_cf('Rad')
        image_x = image_data.x
        image_y = image_data.y

        # Set the axis bounds.
        if coverage == 'FULL':
            ax.set_global()
            info_text='k'

        elif coverage == 'CONUS':
            ax.set_extent([image_x.min(), image_x.max(), image_y.min(), image_y.max()], crs=geos)
            info_text='cyan'

        # Overlay the L1B data
        if self.overlay_l1b:
            # De-normalize the vmax from range [0,1] to natural range
            min_range = float(nc.variables['Rad'].valid_range[0])
            max_range = float(nc.variables['Rad'].valid_range[1])
            vmax = self.vmax*(max_range - min_range)

            if coverage == 'CONUS':
                vmax = vmax/3.5

            # Plot L1B data
            # Note: Increasing vmax lowers contrast. Vmax=small->black; Vmax=large->white
            ax.imshow(open_image['Rad'][:], origin='upper', cmap='gray', transform=geos, vmax=vmax,
                extent=(image_x.min(), image_x.max(), image_y.min(), image_y.max()))

            # Draw coatlines, country borders, lakes, and grid
            # See https://scitools.org.uk/cartopy/docs/v0.14/matplotlib/feature_interface.html
            ax.coastlines(linewidth=0.9, linestyle='solid', color='green')
            ax.add_feature(cfeature.BORDERS, linewidth=0.9, linestyle='solid', 
                facecolor='none', edgecolor='green')
            ax.add_feature(cfeature.LAKES, linewidth=0.9, linestyle='solid', 
                facecolor='none', edgecolor='green')
            ax.gridlines(linewidth=0.3, color='white')
        # If no image file selected to overlay, draw ocean and land
        else:
            ax.stock_img()

            # Draw the coastlines, countries, parallels and meridians
            #bmap.drawparallels(np.arange(-90.0, 90.0, 10.0), linewidth=0.3, color='white')
            #bmap.drawmeridians(np.arange(0.0, 360.0, 10.0), linewidth=0.3, color='white')
            ax.coastlines(linewidth=0.9, linestyle='solid', color='black')
            ax.add_feature(cfeature.BORDERS, linewidth=0.9, linestyle='solid', 
                facecolor='none', edgecolor='black')
            ax.add_feature(cfeature.LAKES, linewidth=0.9, linestyle='solid', 
                facecolor='skyblue', edgecolor='black')
            ax.add_feature(cfeature.RIVERS, linewidth=0.9, linestyle='solid', 
                facecolor='none', edgecolor='skyblue')
            ax.gridlines(linewidth=0.3, color='white')

        # Add a title to the plot
        plt.title(self.sat + " ABI L1B Band " + self.band + " Metric " + self.metric + " from " + \
            self.start_range.strftime('%Y Day %j') + " to " + self.end_range.strftime('%Y Day %j') + \
             "\n" + coverage + " Scan from " + start_formatted + " to " + end_formatted)

        # Read some variables from the NetCDF header in order to use it in the plot
        center = str(geo_extent.geospatial_lon_center)
        west = str(geo_extent.geospatial_westbound_longitude)
        east = str(geo_extent.geospatial_eastbound_longitude)
        north = str(geo_extent.geospatial_northbound_latitude)
        south = str(geo_extent.geospatial_southbound_latitude)

        # Close netCDF file when finished
        nc.close()
        nc = None  

        # Put the information retrieved from the header in the final image
        plt.text(0.01, 0.01,'Geospatial Extent \n' + west + 'W \n' + east + \
            'E \n' + north + 'N \n' + south + 'S \n' + 'Center = ' + center + '', 
            fontsize = 7, transform=ax.transAxes, color=info_text)

        if self.nir_flg:
            plt.text(0.01, 0.94,"WARNING: Selected scene \n{} \nnot in measurement file"\
                .format(self.scene2extract), color='red', fontsize=8, transform=ax.transAxes)

        # Project the coordinates from measurements dataframe
        x = np.array(measure_df['Lon'])
        y = np.array(measure_df['Lat'])

        # Generate the vectors
        delta_ew = np.array(measure_df['AVG_EW'])
        delta_ns = np.array(measure_df['AVG_NS'])
        # Calculate magnitudes so can colorize
        mag = (delta_ew**2 + delta_ns**2)**(0.5)
        # Normalize the arrows
        delta_ew_norm = delta_ew/np.sqrt(delta_ew**2 + delta_ns**2)
        delta_ns_norm = delta_ns/np.sqrt(delta_ew**2 + delta_ns**2)

        # Draw the vectors
        ax.quiver(x, y, delta_ew_norm, delta_ns_norm, mag, width=0.003, 
            cmap='jet', transform=ccrs.PlateCarree())

        # Insert the colorbar
        # Source: https://www.geeksforgeeks.org/matplotlib-pyplot-colorbar-function-in-python/
        norm = mpl.colors.Normalize(vmin=min(mag), vmax=max(mag))
        cmap = plt.get_cmap('jet')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, orientation=orientation, label='Mean Shift Magnitude, urad')

        # Show or save the plot
        if save_plot:
            plt.savefig('vplot.png', bbox_inches='tight')
        else:
            plt.show()

    def build_measurement_df(self):
        """
        """
        measure_df = pd.DataFrame([])

        for measurement_file in self.measurement_files:

            # Import the measurements dataframe
            df = self.extract_geoloc(measurement_file)

            # Append dataframes.
            measure_df = measure_df.append(df, ignore_index = True)

        # Calculate the mean EW and NS for each chip
        measure_df['AVG_EW'] = measure_df.groupby(['CHIP_NAME','BAND_NUM'])['DELTA_EW'].transform('mean')
        measure_df['AVG_NS'] = measure_df.groupby(['CHIP_NAME','BAND_NUM'])['DELTA_NS'].transform('mean')

        # Remove the duplicates within each "band_num">"chip_name" subgroups
        measure_df = measure_df.drop_duplicates(['CHIP_NAME', 'BAND_NUM']) 

        return measure_df
