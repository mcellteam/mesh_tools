#ifndef TILE_UTIL_H
#define TILE_UTIL_H

void free_group_ary(char **p);

void pass_tiling_parameters(VertType *lines_ary0, VertType *lines_ary1, 
							short ** new_group_ary0, short *new_group_ind_ary0, 
							int group_num0,	short ** new_group_ary1, 
							short *new_group_ind_ary1, int group_num1);

int push_boundary_array(int mode,short group0, short ind0, char tile_dir,short group1, 
						short ind1);

int pop_boundary_array(int mode, int *group0, int *ind0, int *tile_dir,int *group1,
					   int *ind1);

int find_tiling_direction_from_array(int mode, int group0, int group1, int ind0, 
									 int ind1);

void security_check_untiled_stack(int mode);

int is_pt_in_LS(VertType *p1, TileType* tile);

int is_triangle_legal(VertType *p0, VertType *p1, VertType *p2,TileType *tile0, 
					  TileType *tile1);

int is_line_legal(int mode, VertType *p0, VertType *p1,TileType *tile0, TileType *tile1);

void find_best_tiling_group_to_group(int topdown, VertType *lines_ary0, 
									 VertType *lines_ary1, short *new_group_ary0, 
									 short ary_ind0, int group0,short *new_group_ary1, 
									 short ary_ind1, int group1, TileType *best_tiling_tab0,
									 TileType *best_tiling_tab1);

void save_best_tile(char *s, VertType *lines_ary0, VertType *lines_ary1,
					short * group_ary, short num, short ** other_group_ary,
					TileType *tile_table);

int calc_best_tiling_table(VertType *lines_ary0, VertType *lines_ary1, 
						   short ** new_group_ary0, short *new_group_ind_ary0, 
						   int group_num0, short ** new_group_ary1, 
						   short *new_group_ind_ary1, int group_num1,
						   short *group_connection_table);

int find_start_tile_point(int topdown, int group, int index, int *group0, int *group1,
						  int *ind0, int *ind1);

int find_best_start_tile_point(int group0, int ind0, int *ret_group1, int *ret_ind1);

Linklist *record_break_contour(int topdown, int group, int ind0);

int break_contour_if_necessary( VertType *lines_ary0, VertType *lines_ary1,
							   short ** new_group_ary0, short *new_group_ind_ary0, 
							   int group_num0,short ** new_group_ary1,
							   short *new_group_ind_ary1, int group_num1,int *lines_ind0,
							   int *lines_ind1, int lines_ary_size0, int lines_ary_size1);

#endif
