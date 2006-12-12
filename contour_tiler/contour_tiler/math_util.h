#ifndef MATH_UTIL_H
#define MATH_UTIL_H

#include "common.h"


double calculate_area(VertType *lines_ary, short *group_ary, short group_size);

int normalize_vector(VertType *v1);

void bubble_sort(double *, int);

double line_distance(VertType *, VertType *);

double line_distance_2d(VertType *, VertType *);

double outer_product_2d(double *, double *);

int calc_line_equation(VertType *, VertType *, double *, double *, double *);

int is_inside_triangle(double *, double [3][2]);

int is_pt_on_line(VertType *, VertType *, VertType *);

int check_bounding_boxes_intersection(VertType *, VertType *, VertType *, VertType *, 
									  int);

int find_intersection(VertType *, VertType *, VertType *, VertType *, VertType *t,
					  double *, double *);

int solve_2_equations(double, double, double, double, double, double,
					  double *, double *);

int find_intersection_3D(VertType *, VertType *, VertType *, VertType *, VertType *, 
						 double *);

int is_inside_contour(VertType *, VertType *, short *, short , double *, int *);

int is_contour_CCW(VertType *lines_ary, short *group_tab, short ary_ind);

double calc_angle(VertType *u, VertType *v);

void draw_small_square(FILE *, double, double, double);

void draw_square(FILE *, double, double, double);

void draw_grid(FILE *, double, double, double);

void draw_line_seg(FILE *, double, double, double, double, double, double);

void draw_vert2vert(FILE *, VertType *,VertType *);

void draw_triangle(FILE *, VertType *,  VertType *, VertType *);

void draw_triangle_FP(FILE *, VertType *,  VertType *, VertType *);

void draw_triangle_mode_FP(FILE *, VertType *,  VertType *, VertType *);

#endif
