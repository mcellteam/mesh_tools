#!/bin/bash

prog='/home/jkinney/src/mesh_tools/reconstruct2contourtiler/reconstruct2contourtiler'

echo "01_return_raw_r_0.dat"
mkdir /scratch/projects/reconstruct2contourtiler/output01
time $prog -I cube -I cyl -I grid -I domain1 \
-i /scratch/paper1_spillover/synapses_CA1/contours_from_kristen \
-o /scratch/projects/reconstruct2contourtiler/output01 \
--min_section=60 \
--max_section=160 \
--curvature_gain=1E2 \
--proximity_gain=3E0 \
--min_point_per_contour=0 \
--return_raw_contour_points \
--print_detailed_info \
&>01_return_raw_r_0.dat

echo "02_return_raw_r_4.dat"
mkdir /scratch/projects/reconstruct2contourtiler/output02
time $prog -I cube -I cyl -I grid -I domain1 \
-i /scratch/paper1_spillover/synapses_CA1/contours_from_kristen \
-o /scratch/projects/reconstruct2contourtiler/output02 \
--min_section=60 \
--max_section=160 \
--curvature_gain=1E2 \
--proximity_gain=3E0 \
--min_point_per_contour=4 \
--return_raw_contour_points \
--print_detailed_info \
&>02_return_raw_r_4.dat

echo "03_return_interpolated_curvature_gain_1E2_r_0.dat"
mkdir /scratch/projects/reconstruct2contourtiler/output03
time $prog -I cube -I cyl -I grid -I domain1 \
-i /scratch/paper1_spillover/synapses_CA1/contours_from_kristen \
-o /scratch/projects/reconstruct2contourtiler/output03 \
--min_section=60 \
--max_section=160 \
--curvature_gain=1E2 \
--proximity_gain=3E0 \
--min_point_per_contour=0 \
--return_interpolated_raw_points \
--print_detailed_info \
&>03_return_interpolated_curvature_gain_1E2_r_0.dat

echo "04_return_interpolated_curvature_gain_1E2_r_4.dat"
mkdir /scratch/projects/reconstruct2contourtiler/output04
time $prog -I cube -I cyl -I grid -I domain1 \
-i /scratch/paper1_spillover/synapses_CA1/contours_from_kristen \
-o /scratch/projects/reconstruct2contourtiler/output04 \
--min_section=60 \
--max_section=160 \
--curvature_gain=1E2 \
--proximity_gain=3E0 \
--min_point_per_contour=4 \
--return_interpolated_raw_points \
--print_detailed_info \
&>04_return_interpolated_curvature_gain_1E2_r_4.dat

echo "05_return_interpolated_curvature_gain_0E2_r_0.dat"
mkdir /scratch/projects/reconstruct2contourtiler/output05
time $prog -I cube -I cyl -I grid -I domain1 \
-i /scratch/paper1_spillover/synapses_CA1/contours_from_kristen \
-o /scratch/projects/reconstruct2contourtiler/output05 \
--min_section=60 \
--max_section=160 \
--curvature_gain=0E2 \
--proximity_gain=3E0 \
--min_point_per_contour=0 \
--return_interpolated_raw_points \
--print_detailed_info \
&>05_return_interpolated_curvature_gain_0E2_r_0.dat

echo "06_return_interpolated_curvature_gain_0E2_r_4.dat"
mkdir /scratch/projects/reconstruct2contourtiler/output06
time $prog -I cube -I cyl -I grid -I domain1 \
-i /scratch/paper1_spillover/synapses_CA1/contours_from_kristen \
-o /scratch/projects/reconstruct2contourtiler/output06 \
--min_section=60 \
--max_section=160 \
--curvature_gain=0E2 \
--proximity_gain=3E0 \
--min_point_per_contour=4 \
--return_interpolated_raw_points \
--print_detailed_info \
&>06_return_interpolated_curvature_gain_0E2_r_4.dat

echo "07_splined_curvature_gain_1E2_r_0.dat"
mkdir /scratch/projects/reconstruct2contourtiler/output07
time $prog -I cube -I cyl -I grid -I domain1 \
-i /scratch/paper1_spillover/synapses_CA1/contours_from_kristen \
-o /scratch/projects/reconstruct2contourtiler/output07 \
--min_section=60 \
--max_section=160 \
--curvature_gain=1E2 \
--proximity_gain=3E0 \
--min_point_per_contour=0 \
--print_detailed_info \
&>07_splined_curvature_gain_1E2_r_0.dat

echo "08_splined_curvature_gain_1E2_r_4.dat"
mkdir /scratch/projects/reconstruct2contourtiler/output08
time $prog -I cube -I cyl -I grid -I domain1 \
-i /scratch/paper1_spillover/synapses_CA1/contours_from_kristen \
-o /scratch/projects/reconstruct2contourtiler/output08 \
--min_section=60 \
--max_section=160 \
--curvature_gain=1E2 \
--proximity_gain=3E0 \
--min_point_per_contour=4 \
--print_detailed_info \
&>08_splined_curvature_gain_1E2_r_4.dat

echo "09_splined_curvature_gain_0E2_r_0.dat"
mkdir /scratch/projects/reconstruct2contourtiler/output09
time $prog -I cube -I cyl -I grid -I domain1 \
-i /scratch/paper1_spillover/synapses_CA1/contours_from_kristen \
-o /scratch/projects/reconstruct2contourtiler/output09 \
--min_section=60 \
--max_section=160 \
--curvature_gain=0E2 \
--proximity_gain=3E0 \
--min_point_per_contour=0 \
--print_detailed_info \
&>09_splined_curvature_gain_0E2_r_0.dat

echo "10_splined_curvature_gain_0E2_r_4.dat"
mkdir /scratch/projects/reconstruct2contourtiler/output10
time $prog -I cube -I cyl -I grid -I domain1 \
-i /scratch/paper1_spillover/synapses_CA1/contours_from_kristen \
-o /scratch/projects/reconstruct2contourtiler/output10 \
--min_section=60 \
--max_section=160 \
--curvature_gain=0E2 \
--proximity_gain=3E0 \
--min_point_per_contour=4 \
--print_detailed_info \
&>10_splined_curvature_gain_0E2_r_4.dat

#--deviation_threshold=.005 -S 0.5 \
#$prog -i /scratch/paper1_spillover/synapses_CA1/contours_from_kristen \
