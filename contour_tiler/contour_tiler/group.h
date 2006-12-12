#ifndef GROUP_H
#define GROUP_H

void find_bound_box();

void ReorderData(int begin,int end,int xyz, double value,int *end1, int *begin1);

void put_into_subcubes(int begin, int end, double a1, double a2, double b1, double b2);

void verify_sorted_data();

int  Neighbors(double *point,int *neiborindex);

void remove_redundant_vertice();

void label_group(int offset, short *vert_group_ary, VertType *vert_ary,
				 short **group_ary, short *group_ind_ary, int group_num);

int tri_region_growing(char *slice_group[2], int start_index,
					   short *vert_group_ary);

int trace_contours(VertType *vert_ary, int vert_ind, short **group_ary, 
				   short *group_ind_ary);

int trace_untiled_contours(VertType *vert_ary, int complete_ind, int unused_ind, 
						   int vert_ind, short **group_ary, short *group_ind_ary);

#endif
