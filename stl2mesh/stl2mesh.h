#define STL2MESH_VERSION "stl2mesh Version 1.00  11/11/2002\n"

#include "vector.h"

/* object types */
#define META_OBJ 1
#define POLY 2
#define NURBS 3

#define EPS_D 1.0e-12

typedef unsigned char byte;

/*
struct object {
        char *name;
        int object_type;
        void *contents;
	struct object *parent;
};
*/

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
        int merged_index;
        double projection;
        struct vector3 *vertex;
        struct vector3 *normal;
};

struct vertex_list {
        struct vertex *vertex; 
        struct vertex_list *next; 
};

struct polyhedron {
	int n_polys;
	struct polygon_list *polygon_list;
	int n_verts;
	int merged_verts;
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
        int n_verts;
        struct vertex_list *vertex_list; 
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

struct double_list {
	double val;
	struct double_list *next;
};

struct subvolume {
        struct vertex_list_list *vertex_list_list;
        double x1,x2;
        double y1,y2;
        double z1,z2;
};

struct volume {
        int n_x_partitions; /* number of elements in x_partitions array */
        int n_y_partitions; /* number of elements in y_partitions array */
        int n_z_partitions; /* number of elements in z_partitions array */
        double *x_partitions; /* x partition boundaries: -infin,...,+infin */
        double *y_partitions; /* y partition boundaries: -infin,...,+infin */
        double *z_partitions; /* z partition boundaries: -infin,...,+infin */
        int n_subvol; /* total number of subvolumes */
        int n_x_subvol; /* number of subvolumes along x dimension */
        int n_y_subvol; /* number of subvolumes along y dimension */
        int n_z_subvol; /* number of subvolumes along z dimension */
        struct subvolume *subvol; /* array of subvolume structs */
};

int compare_verts(const void *p1, const void *p2);
int distinguishable(double a,double b,double eps);
void merge_duplicate_verts(struct polyhedron *php, double epsilon);
void output_object();
void output_nurbs();
void output_polyhedron(struct polyhedron *php);

#ifndef DEBUG
void no_printf();
#endif
