#ifndef DECOM_UTIL_H
#define DECOM_UTIL_H


int is_convex(VertType *,short *,int);

void calc_vector_from_cir(VertType *, VertType *,short *, int, int);
//void free_polygons_ary(PolygonStruct **);//KLC

void put_a_poly_into_poly_structure(PolygonStruct**, VertType *, int *, int ,
									short *, int);

double dot_product_2d(VertType *, VertType *);

int  calc_norm_vector_from_2vp(VertType *, VertType *, VertType *);

int is_center_line_not_cross_poly(VertType *, VertType *, short *, int );

int check_3_sharp_triangle_angle(VertType *, VertType *, VertType *);

int is_angle_not_sharp(VertType *, short *, int , VertType *);

void find_polygon_center(VertType *, short *, int , VertType *);

int find_quad_center(VertType *, short *, int , VertType *);

int is_line_cross_polygon(int , int , VertType *, short *, int );

void swap_int(int *a, int *b);

double calc_edge_distance_2d_of_a_poly(VertType *vert_ary, short *vert_ind_ary,
									   int poly_num, double *dist_ary);

#endif


