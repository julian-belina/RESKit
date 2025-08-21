information, documented by s.chen@fz-juelich.de on 2025-08-21

**To split ERA5 data set into tiles**

# why this is important?
In order to run the RESKit model using ERA5 data, the splitting process is not necessary. But when your study area is a large area like global and when you have thousands of placements to be calculated, the weather data reading part will feel like forever. By splitting the ERA5 data into tiles, RESKit could locate percisely the placements to a certain tile and read the associated tile only for these placements. This would save lots of weather data reading time compared to ERA5 data without splitting.

# how to split ERA5 data set into tiles
0. This procedure was initially created by s-ryberg and s-ishmam, as one can see in the head comments of the python script. Original location is at [jugit](https://jugit.fz-juelich.de/iek-3/groups/global-systems/ishmam/weather_data_processing/-/tree/master?ref_type=heads).

1. job calculation should start by specifing years in submit -> slurm_process_set.sh and execute it, where the ERA5_processor_tile.py script will be then called.

