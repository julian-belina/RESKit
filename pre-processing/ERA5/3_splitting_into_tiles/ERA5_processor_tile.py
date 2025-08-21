# coding: utf-8
# @authors: s-ryberg, s-ishmam

#%%
from scipy.interpolate import RectBivariateSpline, interp1d
from datetime import datetime, timedelta
from functools import reduce
import numpy as np
import geokit as gk
import shutil
import sys
import os
import netCDF4 as nc
import argparse
VERSION = "0.0.2"

#%%
####### Command line inputs

parser = argparse.ArgumentParser(
    description='Process ERA5 data from raw Global data into yearly files over a given tile and year')
parser.add_argument('year', type=str, help='Target weather year')
parser.add_argument(  'xi', type=int, help='X-Index of tile')
parser.add_argument(  'yi', type=int, help='Y-Index of tile')
parser.add_argument('zoom', type=int, help='Zoom level')

args = parser.parse_args()

#%%
#ERA5_TOP_DIR = os.path.dirname(os.path.dirname(__file__))
ERA5_TOP_DIR = "/data/gears/weather/ERA5/"
# if os.path.isdir("/data2"):
#     print("USING InfiniBand")
#     ERA5_TOP_DIR = ERA5_TOP_DIR.replace("/data/","/data2/")

MONTHS = np.arange(1, 13)


#%%
####### Clip raw data files to temp location

# Set Targets
target_id = f"z{args.zoom}.x{args.xi}.y{args.yi}.y{args.year}"

target_topdir = os.path.join(ERA5_TOP_DIR, "processed")
source_year = args.year

# Set other constants
source_topdir = os.path.join(ERA5_TOP_DIR, "raw")
source_group = "reanalysis-era5-single-levels"

#%%
# make target directory structure

target_dir = target_topdir
for tmp in [args.zoom, args.xi, args.yi, args.year]:
    target_dir = os.path.join(target_dir, str(tmp))
    # if not os.path.isdir(target_dir):
    #     try:
    #         os.mkdir(target_dir)
    #     except:
    #         pass

target_temp_dir = os.path.join(target_dir, "temp")

# if os.path.isdir(target_dir):
#     shutil.rmtree(target_dir)

if not os.path.isdir(target_dir):
    print("Creating target directory!")
    os.mkdir(target_dir)

if not os.path.isdir(target_temp_dir):
    print("Creating target temporary directory!")
    os.mkdir(target_temp_dir)

#%%
# check if all files exist

target_file_variables = [
    "100m_wind_speed.processed",
    # "100m_wind_direction.processed",
    # "10m_wind_speed.processed",
    # "10m_wind_direction.processed",
    # "boundary_layer_height",
    # "forecast_surface_roughness",
    # "total_sky_direct_solar_radiation_at_surface.processed",
    # "surface_solar_radiation_downwards.processed",
    # "2m_temperature",
    # "2m_dewpoint_temperature",
    # "surface_pressure",
]

target_pathlist = []

for variable in target_file_variables:    
    target_file = f"{source_group}.{target_id}.{variable}.nc"
    target_path = os.path.join(target_dir, target_file)
    target_pathlist.append(target_path)

for target_path in target_pathlist:
    assert not os.path.isfile(target_path), f"{target_path} \nFile exists! All processing skipped!"


# In[5]:
# Write note into processed directory
LOG_FILE = "README.md"
with open(os.path.join(target_dir, LOG_FILE), 'w') as fo:
    fo.write("Begun processing at: " + str(datetime.now()) + "\n")
    fo.write("Processor version: " + VERSION + "\n")
    fo.write("Desired year: " + source_year + "\n")
    fo.write("Zoom level: " + str(args.zoom) + "\n")
    fo.write("X-Index: " + str(args.xi) + "\n")
    fo.write("Y-Index: " + str(args.yi) + "\n")

#%%

extent_single_levels = gk.Extent.fromTile(args.xi, args.yi, args.zoom).castTo(gk.srs.EPSG4326).pad(2).xXyY
print("EXTENT:", extent_single_levels)

with open(os.path.join(target_dir, LOG_FILE), 'a') as fo:
    fo.write("Extent:" + reduce(lambda a,
                                b: f"{a},{b}", extent_single_levels) + "\n")


#%%

def clip_dataset(variable, month, target_dir=target_temp_dir, lon_lat_box=extent_single_levels):
    """Clips a netCDF dataset to the spatial domain given by 'lon_lat_box' and deposits it into 'target_dir'"""

    source_name = f"{source_group}.{source_year}{month:02d}.{variable}.nc"
    print("  ", source_name)
    source = os.path.join(source_topdir, source_year,
                          f"{month:02d}", source_name)

    target = os.path.join(target_temp_dir, source_name)

    r = os.system(
        f"cdo sellonlatbox,{lon_lat_box[0]},{lon_lat_box[1]},{lon_lat_box[2]},{lon_lat_box[3]} {source} {target}")
    if not r == 0:
        raise RuntimeError("File Clipping Failed:", variable, month)


#%%

source_variables = [
    # VARIABLE NAME                                ,
    ### "friction_velocity",
    "100m_v_component_of_wind",
    "100m_u_component_of_wind",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "boundary_layer_height",
    "forecast_surface_roughness",
    "total_sky_direct_solar_radiation_at_surface",
    "surface_solar_radiation_downwards",
    "2m_temperature",
    "2m_dewpoint_temperature",
    "surface_pressure",
]
for variable in source_variables:
    print("CLIPPING:", variable)
    for month in MONTHS:
        clip_dataset(variable, month)
        # continue



#%%
####### Do processing

def copy_basic(raw_ds, target_ds):
    # Copy NC Attributes
    for attr in raw_ds.ncattrs():
        target_ds.setncattr(attr, getattr(raw_ds, attr))

    # Copy Dimensions
    for dim in raw_ds.dimensions:
        target_ds.createDimension(
            raw_ds.dimensions[dim].name, raw_ds.dimensions[dim].size)

    # Copy lat, lon, time
    for variable in ['latitude', 'longitude', 'time']:
        target_ds.createVariable(
            variable, raw_ds[variable].dtype, raw_ds[variable].dimensions)
        for attr in raw_ds[variable].ncattrs():
            target_ds[variable].setncattr(
                attr, getattr(raw_ds[variable], attr))
        target_ds[variable][:] = raw_ds[variable][:]


#%%
####### Process Irradiance files

def process_irradiance(variable, shortname, month,):
    file_name = f"{source_group}.{source_year}{month:02d}.{variable}.nc"
    print("PROCESSING:", file_name)
    raw_ds = nc.Dataset(os.path.join(target_temp_dir, file_name))

    target_file_name_pr = f"{source_group}.{source_year}{month:02d}.{variable}.processed.nc"
    target_ds = nc.Dataset(os.path.join(
        target_temp_dir, target_file_name_pr), mode='w')

    try:
        copy_basic(raw_ds, target_ds)

        # Compute radiation over individual hours
        ssrd_raw = raw_ds[shortname][:]

        no_data = 2**16-1
        max_data = 2**16-2
        ssrd_final = ssrd_raw / 3600  # Convert to W/m2

        # Save to output file
        v = target_ds.createVariable(
            shortname, "uint16", ("time", "latitude", "longitude"), fill_value=no_data)
        if hasattr(v, "standard_name"):
            v.setncattr('standard_name',
                        raw_ds[shortname].standard_name+"_processed")
        v.setncattr(
            'long_name', raw_ds[shortname].long_name+'. Processed to hourly irradiance')
        v.setncattr('units', 'W m**-2')
        v.setncattr('scale_factor',  1400 / max_data)
        v.setncattr('missing_value',  no_data)
        v.set_auto_scale(True)
        v[:] = ssrd_final

        # Update file history
        target_ds.history += str(datetime.today()) + \
            "Processed by Shitab Ishmam"

        # Close files
        raw_ds.close()
        target_ds.close()
    except Exception as e:
        raw_ds.close()
        target_ds.close()
        raise e

for month in MONTHS:
    process_irradiance("surface_solar_radiation_downwards", "ssrd", month)
    process_irradiance(
        "total_sky_direct_solar_radiation_at_surface", "fdir", month)

#%%
####### Process Wind Speed

# Save to netCDF
def make_windvar_nc(target_name, data, height, is_winddir=False):
    reference_ds = f"{source_group}.{source_year}{month:02d}.100m_u_component_of_wind.nc"
    ref_ds = nc.Dataset(os.path.join(target_temp_dir, reference_ds))

    target_ds = nc.Dataset(os.path.join(
        target_temp_dir, target_name), mode='w')

    try:
        copy_basic(ref_ds, target_ds)

        # Make main variable
        if is_winddir:
            var_name = f"wd{height}"
            standard_name = f"wind_direction_at_{height}m"
            long_name = f"Total wind direction at {height} m. Processed from ERA5:u100,v100"
            max_value = 360
            units = "degrees"
        else:
            var_name = f"ws{height}"
            standard_name = f"wind_speed_at_{height}m"
            long_name = f"Total wind speed at {height} m. Processed from ERA5:u10,v10,u100,v100,blh,fsr and ERA5-Land:u10v10"
            max_value = 50
            units = "m s-1"

        no_data = 2**16-1
        max_data = 2**16-2
        v = target_ds.createVariable(
            var_name, "uint16", ("time", "latitude", "longitude"), fill_value=no_data)
        v.setncattr('standard_name', standard_name)
        v.setncattr('long_name', long_name)
        v.setncattr('units', units)
        v.setncattr('scale_factor',  max_value / max_data)
        v.setncattr('missing_value',  no_data)
        v.set_auto_scale(True)
        v[:] = data

        target_ds.close()
        ref_ds.close()
    except Exception as e:
        target_ds.close()
        ref_ds.close()
        raise e


def process_windspeed(height, month):
    # Open files and copy basic parameters
    filename_u = f"{source_group}.{source_year}{month:02d}.{height}m_u_component_of_wind.nc"
    filename_v = f"{source_group}.{source_year}{month:02d}.{height}m_v_component_of_wind.nc"

    print("PROCESSING:", filename_u, filename_v)

    ds_u = nc.Dataset(os.path.join(target_temp_dir, filename_u))
    ds_v = nc.Dataset(os.path.join(target_temp_dir, filename_v))

    # Extract data
    data_u = ds_u[f'u{height}'][:]
    data_v = ds_v[f'v{height}'][:]

    data_uv = np.sqrt(np.power(data_u, 2) + np.power(data_v, 2))

    # V is "towards the North", and U is "towards the east", so order should be good...
    data_uvdir = np.arctan2(data_v, data_u)
    data_uvdir = np.degrees(data_uvdir)
    sel = data_uvdir < 0
    data_uvdir[sel] = data_uvdir[sel]+360

    # Close datasets
    ds_u.close()
    ds_v.close()

    # Write output datasets
    filename = f"{source_group}.{source_year}{month:02d}.{height}m_wind_speed.processed.nc"
    make_windvar_nc(filename, data_uv, height, is_winddir=False)

    filename = f"{source_group}.{source_year}{month:02d}.{height}m_wind_direction.processed.nc"
    make_windvar_nc(filename, data_uvdir, height, is_winddir=True)


for month in MONTHS:
    process_windspeed(100, month)
    process_windspeed(10, month)


# %%
####### Combine files to yearly

def concatenate_year(variable, shortname, max_value, min_value=0, max_data=2**16-2, min_data=0, do_scaling=True,
                     fill_value=2**16-1, dtype='uint16', chunksizes=None, zlib=None, tolerance=1e-1):
    source_files = f"{source_group}.{source_year}{{:02d}}.{variable}.nc"
    print("Concatenating:", source_files)

    time = []
    time_unit = None
    lats = None
    lons = None
    data = []
    for i in MONTHS:
        source_path = os.path.join(target_temp_dir, source_files.format(i))
        ds = nc.Dataset(source_path)
        try:
            if time_unit is None:
                time_unit = ds['time'].units
                lats = ds['latitude'][:]
                lons = ds['longitude'][:]
            else:
                assert time_unit == ds['time'].units
                assert (lats == ds['latitude'][:]).all()
                assert (lons == ds['longitude'][:]).all()

            time.append(ds['time'][:])
            data.append(ds[shortname][:])
        except Exception as e:
            ds.close()
            raise e
        ds.close()

    time = np.concatenate(time)
    data = np.concatenate(data)

#     print(data.min(), data.max(), data.mean())
#     return

    ######################
    # Make output dataset
    target_file = f"{source_group}.{target_id}.{variable}.nc"
    target_path = os.path.join(target_dir, target_file)

    reference_ds = nc.Dataset(source_path)
    target_ds = nc.Dataset(target_path, mode='w')

    try:
        # Copy NC Attributes
        for attr in reference_ds.ncattrs():
            target_ds.setncattr(attr, getattr(reference_ds, attr))

        # Copy Dimensions
        target_ds.createDimension('time', time.shape[0])
        target_ds.createDimension('latitude', lats.shape[0])
        target_ds.createDimension('longitude', lons.shape[0])

        # Copy lat, lon, time
        for variable in ['latitude', 'longitude', 'time']:
            target_ds.createVariable(
                variable, reference_ds[variable].dtype, reference_ds[variable].dimensions)
            for attr in reference_ds[variable].ncattrs():
                target_ds[variable].setncattr(
                    attr, getattr(reference_ds[variable], attr))
            if variable == "time":
                target_ds[variable][:] = time
            else:
                target_ds[variable][:] = reference_ds[variable][:]

        # Make data
        datavar = target_ds.createVariable(shortname, dtype, ('time', 'latitude', 'longitude'),
                                           fill_value=fill_value, chunksizes=chunksizes, zlib=zlib)
        for attr in reference_ds[shortname].ncattrs():
            if attr not in ['_FillValue', 'scale_factor', 'add_offset']:
                target_ds[shortname].setncattr(
                    attr, getattr(reference_ds[shortname], attr))

        datavar.setncattr('missing_value',  fill_value)
        if do_scaling == True:
            scale_factor = (max_value - min_value) / (max_data - min_data)
            add_offset = max_value - scale_factor*max_data

            datavar.setncattr('scale_factor',  scale_factor)
            datavar.setncattr('add_offset',  add_offset)
            datavar.set_auto_scale(True)

        datavar[:] = data

        target_ds.close()
        reference_ds.close()
    except Exception as e:
        target_ds.close()
        reference_ds.close()
        raise e

    # Test written values
    if tolerance is not None:
        target_ds = nc.Dataset(target_path)
        data_written = target_ds[shortname][:]

        try:
            assert np.isclose(data_written, data, atol=tolerance).all()
        except AssertionError as e:
            target_ds.close()
            raise RuntimeError(
                "Written values are not equal (approximately) to the desired target values")
        target_ds.close()


PARAMS = dict(zlib=True, chunksizes=(256, 16, 16))

concatenate_year("100m_wind_speed.processed",
                 'ws100',  max_value=50,    **PARAMS)
concatenate_year("100m_wind_direction.processed",
                 'wd100',  max_value=360,   **PARAMS)
concatenate_year("10m_wind_speed.processed",
                 'ws10',  max_value=50,    **PARAMS)
concatenate_year("10m_wind_direction.processed",
                 'wd10',  max_value=360,   **PARAMS)
concatenate_year("boundary_layer_height",           'blh',
                 max_value=12000, tolerance=1, **PARAMS)
concatenate_year("forecast_surface_roughness",      'fsr',  dtype='f4',
                 max_value=None, do_scaling=False, fill_value=-1, **PARAMS)
concatenate_year("total_sky_direct_solar_radiation_at_surface.processed",
                 'fdir',  max_value=1400, **PARAMS)
concatenate_year("surface_solar_radiation_downwards.processed",
                 'ssrd',  max_value=1400, **PARAMS)
concatenate_year("2m_temperature",                  't2m',
                 min_value=100,      max_value=500,    **PARAMS)
concatenate_year("2m_dewpoint_temperature",         'd2m',
                 min_value=100,      max_value=500,    **PARAMS)
concatenate_year("surface_pressure",                 'sp',
                 min_value=50000,  max_value=130000, tolerance=5, **PARAMS)


# In[30]:

shutil.rmtree(target_temp_dir)

with open(os.path.join(target_dir, LOG_FILE), 'a') as fo:
    fo.write("Finished processing at:" + str(datetime.now()) + "\n")
