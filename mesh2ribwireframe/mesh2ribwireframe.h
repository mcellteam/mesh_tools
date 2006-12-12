#define MESH2RIBWIREFRAME_VERSION "mesh2ribwireframe Version 1.00  5/31/2000\n"

#include "vector.h"

/* object types */
#define META_OBJ 1
#define POLY 2
#define NURBS 3

#define EDGE 1

#define HASHSIZE 1001

typedef unsigned char byte;

struct object {
        char *name;
        int object_type;
        void *contents;
	struct object *parent;
};

struct meta_object {
	struct object_list *first_child;
	struct object_list *last_child;
};

struct object_list {
	struct object *object;
	struct object_list *next;
};

struct polyhedron {
	int n_polys;
	struct polygon_list *polygon_list;
	int n_verts;
	struct vector3 **vertex;
        struct vertex_list *unique_vertex;
	struct element_data *element_data;
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
	int vertex_index[4];
};

struct vertex_list {
        int vertex_index;
        int vertex_count;
        struct vector3 *vertex;
        struct vector3 *normal;
        struct vertex_list *next; 
};

struct edge {
        struct sym_table *sym;
	int vert_1;
	int vert_2;
};

struct edge_list {
	struct edge *edge;
        struct edge_list *next;
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

struct sym_table {
        unsigned short sym_type; /*type of symbol stored --
                                   EDGE ...*/
        char *name;              /*name of symbol*/
        void *value;             /*ptr to stored value*/
        struct sym_table *next;  /*next symbol in symbol table*/
};

void output_object();
void output_nurbs();
void output_polyhedron();
struct sym_table **init_symtab();
struct sym_table *retrieve_sym();
struct sym_table *store_sym();

#ifndef DEBUG
void no_printf();
#endif
