information, documented by s.chen@fz-juelich.de on 2025-07-24

**To explain why we need to post-process the downloaded ERA5 data in order to run the RESKit model**

# extract every individual variables from the downloaded bulk
It is optional in the sense of running RESKit model. We perform it for the sake of easy management.

# change for solar variables
Downloaded direct and downward solar varaibles from ERA5 are in the hourly total unit (J m-2). RESKit needs instantaneous value (W m-2), i.e., solar irradiance. Therefore we divide the downloaded ones with the total seconds in one hour. 

Due to the current approaches used in solar PV and CSP calculation, the time shift of one hour towards future is necessary.

# change for wind variables
RESKit expects ready-to-use wind speed, so we need to calculate wind speed based on its components u and v that are downloaded from the host website.

# example processing script: variable_preparation.sh
see the corresponding script variable_preparation.sh for detials of post-processing ERA5 data. FYI, knowledge about Linux tool/command CDO is necessary for executing and understanding this script.

In the near future, we may consider wrapping these command via a python interface using python package - python-cdo.

