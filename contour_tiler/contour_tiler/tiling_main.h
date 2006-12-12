#ifndef TILING_H
#define TILING_H

void debug_pt();

void draw_tiled_triangle(VertType *p1, VertType *p2, VertType *p3);

void draw_tiled_triangle_no_mesh(VertType *p1, VertType *p2, VertType *p3);

int check_and_process_optimum_triangle_pair(VertType *lines_ary0, VertType *lines_ary1,
											int group0, int group1,int ind0, int ind1, 
											short group_num0, short group_num1,
											short *group_ary0, short *group_ary1,
											TileType *best_tiling_table0, 
											TileType *best_tiling_table1,int *ret_pt0, 
											int *ret_pt1,int *ret_edge_ind0,
											int *ret_edge_ind1);

int check_and_process_optimum_triangle(int topdown,VertType *lines_ary0,
									   VertType *lines_ary1, int dir, int group0, 
									   int group1, int ind0, int ind1, short group_num0,
									   short *group_ary0, short *group_ary1,
									   TileType *best_tiling_table0, 
									   TileType *best_tiling_table1, int *new_pt);

int check_and_process_regular_triangle(int topdown,VertType *lines_ary0,
									   VertType *lines_ary1, int dir, int group0, 
									   int group1, int ind0, int ind1, short group_num0,
									   short *group_ary0, short *group_ary1,
									   TileType *best_tiling_table0, 
									   TileType *best_tiling_table1, double *dist,
									   int *edge_ret);

int check_and_process_shortest_regular_triangle(int *tile_dir, VertType *lines_ary0, 
												VertType *lines_ary1, int group0, 
												int group1, int ind0, int ind1, 
												short group_num0, short group_num1,
												short *group_tab0, short *group_tab1,
												TileType *top_best_tiling_table, 
												TileType *bot_best_tiling_table,
												int *new_pt0, int *new_pt1);

int is_both_edge_used(TileType *best_tiling_tab, int ind0, short group_num);

void do_best_tiling(int group0, int group1, int ind0, int ind1, VertType *lines_ary0, 
					VertType *lines_ary1, short * group_tab0, short group_num0,
					short * group_tab1, short group_num1);

void do_overlapping_pt_tiling(int group0, int group1, int ind0, int ind1,
							  VertType *lines_ary0, VertType *lines_ary1,
							  short * group_tab0, short group_num0, short * group_tab1,
							  short group_num1);

int debug_check(VertType *p1, double x, double y);

void do_partial_tiling(int group0, int group1, int ind0, int ind1, int tile_dir,
					   VertType *lines_ary0, VertType *lines_ary1, short * group_tab0, 
					   short group_num0, short * group_tab1, short group_num1);

void triangulate_quad( VertType *lines_ary0, VertType *lines_ary1, int group0, 
					  int ind0, int new_ind0, int kk0, int kk1, int group1, int ind1, 
					  int new_ind1, int kk2, int kk3,TileType *best_tiling_table0,
					  TileType *best_tiling_table1);

void do_cross_tiling( VertType *lines_ary0, VertType *lines_ary1,short ** new_group_ary0,
					 short *new_group_ind_ary0, int group_num0, short ** new_group_ary1, 
					 short *new_group_ind_ary1, int group_num1);

int do_tiling(VertType *lines_ary0, VertType *lines_ary1, short ** new_group_ary0, 
			  short *new_group_ind_ary0, int group_num0, short ** new_group_ary1, 
			  short *new_group_ind_ary1, int group_num1);

int christiansen(VertType *lines_ary0, VertType *lines_ary1,short ** group_ary0, 
				 short *group_ind_ary0, int group_num0, short ** group_ary1, 
				 short *group_ind_ary1, int group_num1);

#endif

