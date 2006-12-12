#ifndef CORRESPOND_H
#define CORRESPOND_H

void print_group_relation_table(int group_num0, int group_num1);

void build_polarity_table();

void build_hierachy();

int check_if_two_contours_intersect(VertType *lines_ary0, VertType *lines_ary1,
									short *group_tab0, short num0, short *group_tab1,
									int num1,int index_ary_size, short *index_ary0, 
									short *index_ary1, int *ret_num);

void check_intersection_all_contours();

void check_projection_table();

void build_projection_table();

int have_a_negative_projection_vertex(int topdown, int group);

int have_a_positive_projection_vertex(int topdown, int group);

int have_a_vertex_nec_equal_to_another_contour(int topdown, int g0, int g1);

int check_condition2(int g0, int g1);

int is_the_nec_of_different_slice(int topdown,int g0,int g1);

int check_condition3(int topdown, int g0, int g1);

int check_condition3(int topdown, int g0, int g1);

void    make_all_group_counter_clockwise();

void build_contours_relation(VertType *lines_ary0, VertType *lines_ary1,
							 short ** group_ary0, short *group_ind_ary0, int group_num0,
							 short ** group_ary1, short *group_ind_ary1, int group_num1);

void build_contours_levels(VertType *lines_ary[2],short ** group_ary[2], 
						   short *group_ind_ary[2], int group_num[2],
						   short *level_ary[2]);

#endif
