information, documented by s.chen@fz-juelich.de on 2025-07-24

**To download necessary variables of ERA5 data for running RESKit model**

# ERA5 data set information page
Please check the general information of ERA5 data so to have a basic idea about this data:
https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview

# Download ERA5 data using API
info: (1) this is the officially recommended way for downloading ERA5 data (2) depending on how busy the system is for sure, usually, the larger your requested data is, the longer your waiting time. So try to tailor your request as really needed

## Get the API download script
1. At the download page of the ERA5 data (https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download), specify the product type (should be reanalysis, change it only when you know why), variable, year, month, day, time, geographical area, data format (choose NetCDF4 as RESKit works with this format), and download format (suggest to choose Zip)

2. You have to log in to be able to download

3. Click the "show API request code", you should now see the API download script that was tailored by yourself just now.

4. Copy this code to be your API download script

5. please add the following code line to the API download script, since RESKit at the moment is not compatiable with the latest ERA5 version from ECMWF

    ```json
    "data_format": "netcdf_legacy",
    ```

## Example API download script
For your convenience, we provide an example API download scirpt: ERA5_API_download.py, containing variables that are necessary for running RESKit model solar and wind workflows

## Download ERA5 data using the API download script
follow this link (https://cds.climate.copernicus.eu/how-to-api) to see how you need to set up your python environment to successfully download the ERA5 data specified by your API download script