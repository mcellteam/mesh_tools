#ifndef CROSS_H
#define CROSS_H

void calc_vector_from_2vp(VertType *ret_v, VertType *p0, VertType *p1);

int is_polygon_self_crossed(VertType *vert_ary, short *vert_ind_ary, int poly_num);

//int is_convex(short *vert_ind_ary,int poly_num, VertType *ret_vt);
//void calc_vector_from_cir(VertType *ret_v, short *vert_ind_ary, int index);
//void free_polygons_ary(PolygonStruct **polygon_ary);
//KLC as this conflicts with the existing functions in decom_util.h
#endif
