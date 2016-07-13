#!/bin/bash

#echo "01_return_raw_r_0.dat"
#cd /scratch/projects/reconstruct2contourtiler/output01
#time /bin/bash ./mesh_and_convert.sh
#
#echo "02_return_raw_r_4.dat"
#cd /scratch/projects/reconstruct2contourtiler/output02
#time /bin/bash ./mesh_and_convert.sh

#echo "03_return_interpolated_curvature_gain_1E2_r_0.dat"
#cd /scratch/projects/reconstruct2contourtiler/output03
#time /bin/bash ./mesh_and_convert.sh

echo "04_return_interpolated_curvature_gain_1E2_r_4.dat"
cd /scratch/projects/reconstruct2contourtiler/output04
time /bin/bash ./mesh_and_convert.sh

echo "05_return_interpolated_curvature_gain_0E2_r_0.dat"
cd /scratch/projects/reconstruct2contourtiler/output05
time /bin/bash ./mesh_and_convert.sh

echo "06_return_interpolated_curvature_gain_0E2_r_4.dat"
cd /scratch/projects/reconstruct2contourtiler/output06
time /bin/bash ./mesh_and_convert.sh

echo "07_splined_curvature_gain_1E2_r_0.dat"
cd /scratch/projects/reconstruct2contourtiler/output07
time /bin/bash ./mesh_and_convert.sh

echo "08_splined_curvature_gain_1E2_r_4.dat"
cd /scratch/projects/reconstruct2contourtiler/output08
time /bin/bash ./mesh_and_convert.sh

echo "09_splined_curvature_gain_0E2_r_0.dat"
cd /scratch/projects/reconstruct2contourtiler/output09
time /bin/bash ./mesh_and_convert.sh

echo "10_splined_curvature_gain_0E2_r_4.dat"
cd /scratch/projects/reconstruct2contourtiler/output10
time /bin/bash ./mesh_and_convert.sh
