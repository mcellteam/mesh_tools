#ifndef IMAGE2CONT2D_H
#define IMAGE2CONT2D_H

void save_contours_pts_format(char *s, int slice_num);

void save_marching_contours_pts_format(char *s, int slice_num);

void set_boundary_zero(float *fary, int dimx, int dimy);

void get_a_slice(float *fary, int slice_num, float (*valfun)(int, int, int),
				 void (*slicefun)(int, int));

void shift_group_sequence(short *group_tab, int tab_size, int origin);

void reverse_group_direction(short *group_tab, int ind);

void free_pos_ary();

int initialization(int begin_slice_num,int dimx, int dimy, int dimz, double *orig, 
				   double *voxel_unit, double thres );

void allocate_image_memory(int area);

void free_image_memory();

int process_one_slice(float *ary);

int read_and_approximate(int k, float (*valfun)(int, int, int),
						 void (*slicefun)(int, int));

int count_triangles_one_slice(int k, float (*valfun)(int, int, int),
							  void (*slicefun)(int, int));

#endif
