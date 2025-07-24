# bash, documented by s.chen@fz-juelich.de on 2025-07-24
# customized ERA5 variables to run RESKit model

root_dir="/where/is/your/downloaded/ERA5/data/"
output_dir="/your/output/folder/"

consider_time="2018" # depending on your need, also feel free to change the output file name

# sp
v_era5="sp"
input=${root_dir}/data.nc
output=${output_dir}/ERA5_NER5.1_${v_era5}_${consider_time}_inst_hourly.nc

cdo selname,${v_era5} ${input} ${output}

# t2m
v_era5="t2m"
input=${root_dir}/data.nc
output=${output_dir}/ERA5_NER5.1_${v_era5}_${consider_time}_inst_hourly.nc

cdo selname,${v_era5} ${input} ${output}


# d2m
v_era5="d2m"
input=${root_dir}/data.nc
output=${output_dir}/ERA5_NER5.1_${v_era5}_${consider_time}_inst_hourly.nc

cdo selname,${v_era5} ${input} ${output}


# ssrd
v_era5="ssrd"
n_v_era5="ssrd_t_adj"
input=${root_dir}/data.nc
output1=${output_dir}/ERA5_NER5.1_${v_era5}_${consider_time}_inst_hourly.nc
output2=${output_dir}/ERA5_NER5.1_${v_era5}_${consider_time}_inst_hourly.t_adjusted.nc

unit="W m**-2"

cdo -L selname,${v_era5} -divc,3600 -setattribute,${v_era5}@units="${unit}" ${input} ${output1}
cdo -L chname,${v_era5},${n_v_era5} -selname,${v_era5} -shifttime,+1hour -divc,3600 -setattribute,${v_era5}@units="${unit}" ${input} ${output2}


# fdir
v_era5="fdir"
n_v_era5="fdir_t_adj"
input=${root_dir}/data.nc
output1=${output_dir}/ERA5_NER5.1_${v_era5}_${consider_time}_inst_hourly.nc
output2=${output_dir}/ERA5_NER5.1_${v_era5}_${consider_time}_inst_hourly.t_adjusted.nc

unit="W m**-2"

cdo -L selname,${v_era5} -divc,3600 -setattribute,${v_era5}@units="${unit}" ${input} ${output1}
cdo -L chname,${v_era5},${n_v_era5} -selname,${v_era5} -shifttime,+1hour -divc,3600 -setattribute,${v_era5}@units="${unit}" ${input} ${output2}


# fsr
v_era5="fsr"
input=${root_dir}/data.nc
output=${output_dir}/ERA5_NER5.1_${v_era5}_${consider_time}_inst_hourly.nc

cdo selname,${v_era5} ${input} ${output}


# blh
v_era5="blh"
input=${root_dir}/data.nc
output=${output_dir}/ERA5_NER5.1_${v_era5}_${consider_time}_inst_hourly.nc

cdo selname,${v_era5} ${input} ${output}


# ws100
v_era5="ws100"
input=${root_dir}/data.nc
output=${output_dir}/ERA5_NER5.1_${v_era5}_${consider_time}_inst_hourly.nc

cdo expr,${v_era5}="sqrt(u100*u100+v100*v100)" ${input} ${output}


# ws10
v_era5="ws10"
input=${root_dir}/data.nc
output=${output_dir}/ERA5_NER5.1_${v_era5}_${consider_time}_inst_hourly.nc

cdo expr,${v_era5}="sqrt(u10*u10+v10*v10)" ${input} ${output}

