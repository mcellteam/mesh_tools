#ifndef LEGAL_REGION_H
#define LEGAL_REGION_H

void  calc_group_legal_equations(VertType *lines_ary, short * new_group_tab, 
								 short group_ind, TileType *best_tiling_table);

void make_group_legal_equation_tables( VertType *lines_ary0, VertType *lines_ary1,
									  short ** new_group_ary0, short *new_group_ind_ary0,
									  int group_num0, short ** new_group_ary1, 
									  short *new_group_ind_ary1, int group_num1,
									  TileType **best_tiling_table[2]);

#endif
