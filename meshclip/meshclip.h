#define MESHCLIP_VERSION "meshclip Version 1.00  6/17/98\n"

#include <limits.h>
#include "vector.h"

/* object types */
#define META_OBJ 1
#define POLY 2
#define NURBS 3

/* wall characteristics */
#define RFLCT 1 
#define SINK 2
#define TRANSP 3
#define ND 4
#define SUBVOL 5
#define VOL 6

/* ray tracing polygon collision status */
#define ENTERING -1
#define NO_COLLISION 0
#define LEAVING 1

/* polygon/vertex intersection inclusion status */
#define FULLY_OUTSIDE 0
#define FULLY_INSIDE 1
#define ON_EDGE 2

/* box sides */
#define TP 0
#define BOT 1
#define FRNT 2
#define BCK 3
#define LFT 4
#define RT 5
#define ALL_SIDES INT_MAX

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
	struct vertex_list **vertex_list_array;
        struct vertex_list *unique_vertex;
        struct vector3 llf;
        struct vector3 urb;
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
        byte polygon_status;
        struct vector3 *vertex_array[3]; 
        struct vertex_list *vertex_list_array[3]; 
};

struct vertex_list {
        int vertex_index;
        int vertex_count;
        int fully_outside_index;
        int fully_inside_index;
        int on_edge_index;
        byte vertex_status;
        byte fully_outside_member;
        byte fully_inside_member;
        byte on_edge_member;
        struct vector3 *vertex;
        struct vector3 *normal;
        struct vertex_list *next; 
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

/*Declarations for walls in model*/
struct wall {
	byte *wall_type;                 /*array of wall type -- 
	                                   node, sink, transparent, reflective,
					   subvol: indexed by
					   ligand_info.type*/
	byte side;                       /*side of compartment for this wall
	                                   TP BOT FRNT BCK LFT RT*/
	byte n_vert;                     /*number of vertices in wall*/ 
	struct vector3 **vert;           /*array of ptrs to xyz loc of
	                                   vertices of wall*/
	struct vector3 **vert_normal;    /*array of ptrs to vertex normals
	                                   of wall*/
	double *length;                  /*length of each side of wall*/
	double length_first;             /*length of first side*/
	double length_last;              /*length of last side*/
	double r_length_first;           /*reciprocal of length of first side*/
	double r_length_last;            /*reciprocal of length of last side*/
	struct vector3 normal;           /*normal vector of wall*/
	double d;	                 /*d value for point normal form*/
	byte projection;                 /*direction of 2D projection 0,1,2*/
	struct wall *next_wall;          /*ptr to next wall*/
};

struct wall_list {
	struct wall *wall;
	struct wall_list *next;
};

struct subvolume {
	struct wall_list *wall_list;
	struct wall *walls[6];
        double x1,x2;
        double y1,y2;
        double z1,z2;
};

struct volume {
	struct vector3 corner[8]; /* corners of the world */
	int n_x_partitions; /* number of elements in x_partitions array */
	int n_y_partitions; /* number of elements in y_partitions array */
	int n_z_partitions; /* number of elements in z_partitions array */
	double *x_partitions; /* x partition boundaries: -infin,...,+infin */
	double *y_partitions; /* y partition boundaries: -infin,...,+infin */
	double *z_partitions; /* z partition boundaries: -infin,...,+infin */
        struct wall **x_walls; /* array of wall struct ptrs for x partitions */
        struct wall **y_walls; /* array of wall struct ptrs for y partitions */
        struct wall **z_walls; /* array of wall struct ptrs for z partitions */
	int n_subvol; /* total number of subvolumes */
	int n_x_subvol; /* number of subvolumes along x dimension */
	int n_y_subvol; /* number of subvolumes along y dimension */
	int n_z_subvol; /* number of subvolumes along z dimension */
	struct subvolume *subvol; /* array of subvolume structs */
};


void output_object();
void output_nurbs();
void output_polyhedron();

#ifndef DEBUG
void no_printf();
#endif
