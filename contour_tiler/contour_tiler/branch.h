#ifndef BRANCH_H
#define BRANCH_H

int get_untiled_triangle_number();

int store_unused_contour(VertType **untiled_vert_ary, int *complete_index);

int store_boundary_tiling(VertType **untiled_vert_ary, int index);

int get_convex_triangle_number();

void  draw_convex_triangles(PolygonStruct **ret_poly_ary);

void draw_untiled_contours(VertType *untiled_vert_ary, short *untiled_group_ary, 
						   short num);

int triangulate_with_center(VertType *vert_ary, short *vert_index_ary, 
							int poly_num, VertType *vt);

PolygonStruct ** check_is_a_triangle(VertType *vert_ary, short * group_ary, 
									 short poly_num);

PolygonStruct ** check_is_2_triangles(VertType *vert_ary, short * group_ary, 
									  short poly_num);

int collect_untiled_contours(VertType *lines_ary0, VertType *lines_ary1,
							 short ** new_group_ary0, short *new_group_ind_ary0,
							 int group_num0,short ** new_group_ary1,
							 short *new_group_ind_ary1, int group_num1,
							 VertType *untiled_vert_ary,short ** untiled_group_ary,
							 short *untiled_group_ind_ary);

#endif
