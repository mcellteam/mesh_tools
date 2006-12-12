#define MESHALYZER_VERSION "meshalyzer Version 1.00  9/05/2006\n"

#include "vector.h"

/* object types */
#define META_OBJ 1
#define POLY 2
#define NURBS 3

#define EPS_D 1.0e-12

typedef unsigned char byte;

struct object {
        struct object *next;            /**< ptr to next sibling object */
        struct object *parent;          /**< ptr to parent meta object */
        struct object *first_child;     /**< ptr to first child object */
        struct object *last_child;      /**< ptr to last child object */
        char *name;
        byte object_type;
        void *contents;                 /**< ptr to actual physical object */
};

struct object_list {
	struct object *object;
	struct object_list *next;
};

struct vertex {
        int vertex_index;
        int vertex_count;
        int merged_index;
        double projection;
        struct vector3 *vertex;
        struct vector3 *normal;
};

struct vertex_list {
        struct vertex *vertex;
        struct vertex_list *next;
};

struct surface_mesh {
        int n_polys;
        int n_verts;
        int max_verts;
        int merged_verts;
        double volume;
        struct polygon_list *polygon_list;
        struct vertex_list *vertex_list;
        struct vertex **vertex_array;
        struct vector3 llf;
        struct vector3 urb;
};

struct polygon_list {
	struct polygon *polygon;
	struct polygon_list *next;
};

struct element_data {
	int n_verts;
	int *vertex_index;
};

struct polygon {
        int polygon_index;
        double area;
        double aspect_ratio;
        int n_verts;
        struct vertex **vertex_array; 
};

struct nurbs {
	double upts;
	double vpts;
	double uorder;
	double vorder;
	struct double_list *uknots;
	struct double_list *vknots;
	struct ctlpt_list *ctlpts;
};

struct ctlpt {
	byte rational;
	double w;
	double x;
	double y;
	double z;
};

struct ctlpt_list {
	struct ctlpt *ctlpt;
	struct ctlpt_list *next;
};

struct histogram
{
  u_int size;
  double *bin_lower_limit;
  u_int *bin_count;
};

struct double_list {
	double val;
	struct double_list *next;
};

int compare_verts(const void *p1, const void *p2);
int distinguishable(double a,double b,double eps);
void merge_duplicate_verts(struct surface_mesh *smp, double epsilon);
void output_object();
void output_nurbs();
double polygon_area(struct polygon *pop);
double mesh_vol(struct surface_mesh *smp);
void bin_value(double value, struct histogram *histp);
int analyze_surface_mesh(struct surface_mesh *smp);
double triangle_aspect_ratio(struct polygon *pop);
void output_surface_mesh(struct surface_mesh *smp);

#ifndef DEBUG
void no_printf();
#endif
