#ifndef SCALE_H
#define SCALE_H

int scale_contour( VertType *lines_ary,short ** new_group_ary, 
				  short *new_group_ind_ary, int group_num);

int scale_back_triangle(VertType *p);

int find_bounding_box(double box[2][2],VertType *lines_ary,
					  short ** new_group_ary, short *new_group_ind_ary, int group_num);

int calculate_factor(double box1[2][2], double box2[2][2], double x_ratio,
					 double y_ratio);

int scale_shift( double xshift, double yshift, VertType *lines_ary0, 
				VertType *lines_ary1, short ** new_group_ary0, short *new_group_ind_ary0, int group_num0,
				short ** new_group_ary1, short *new_group_ind_ary1, int group_num1);

#endif
