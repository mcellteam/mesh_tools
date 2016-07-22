#ifndef TETGEN2FV_H
#define TETGEN2FV_H

#define TETGEN2FV_VERSION "tetgen2FV Version 1.00  4/14/2005\n"

#include <sys/types.h>
#include "vector.h"
#include "hash.h"

/* object types */
#define META_OBJ 1
#define POLY 2
#define NURBS 3

/* node/edge/tet types */
#define INTERIOR 0
#define INTERFACE 1
#define BOUNDARY 2
#define VIRTUAL 3
#define PARENT 4
#define CHILD 5

/* input modes */
#define NETGEN_MODE 1
#define TETGEN_MODE 2

typedef unsigned char byte;

struct object
{
        char *name;
        int object_type;
        void *contents;
	struct object *parent;
};


struct meta_object
{
	struct object_list *first_child;
	struct object_list *last_child;
};


struct object_list
{
	struct object *object;
	struct object_list *next;
};


struct polyhedron
{
	int n_polys;
	struct polygon_list *polygon_list;
	int n_nodes;
	struct node **node_array;
        struct node_list *unique_node;
	struct element_data *element_data;
};


struct element_data
{
	int n_nodes;
	int *node_index;
};


struct polygon
{
  u_int polygon_index;
  u_int voronoi_index;
  u_int srfnum;
  u_int bcnum;
  u_int domin;
  u_int domout;
  u_int nnodes;
  u_int ntets;
  u_int node_index[3];
  struct tet *shared_tet[2];
  struct vector3 cent;
};


struct polygon_list
{
  struct polygon *polygon;
  struct polygon_list *next;
};


struct voronoi_facet
{
  double area;
  u_int nnodes;
  u_int *node_index;
  struct voronoi_facet *next;  /* facets through virtual edges are */
                               /* composed of multiple pieces */
};


/*
struct voronoi_facet_list
{
  struct voronoi_facet *facet;
  struct voronoi_facet_list *next;
};
*/


struct voronoi_volume_element
{
  double volume;
  u_int nedges;
  u_int nedges_tot;
/*
  u_int ntets;
  u_int nfacets;
  struct voronoi_facet_list *facet_head;
*/
};


struct node 
{
  double x; /* start of struct node must be compatible with struct vector3 */
  double y;
  double z;
  u_int node_index;
  u_int voronoi_index;
  byte node_type;
  u_int nedges;
  u_int ntets;
  u_int nfacets;
  u_int nboundaryfaces;
  u_int ninterfaceneighbors;
  u_int srfnum;
  u_int material_count;
  struct material_list *material_head;
  struct node_material_info *node_material;
  struct voronoi_volume_element *volume_element;
  struct edge_list *edge_head;
  struct polygon_list *boundary_head;
};


struct node_list
{
  struct node *node;
  struct node_list *next;
};


struct node_material_info {
  int matnum;
  u_int child_index;
  double volume;
  u_int nedges;
  u_int ntets;
  u_int nfacets;
  struct edge_list *edge_head;
};


struct edge
{
  double length;
  u_int ntets;
  u_int node_index[2];
  u_int voronoi_index;
  u_int material_count;
  int matnum;  /* material of an INTERIOR edge */
  byte edge_type;
  struct voronoi_facet *facet;
  struct edge_material_info *edge_material;
  struct tet_list *tet_head;
  struct polygon *boundary_face[2];
  struct vector3 *cent;
};


struct edge_list        
{
  struct edge *edge;
  struct edge_list *next;
};


struct edge_material_info
{
  int matnum;
  u_int ntets;
  struct tet_list *tet_head;
  struct voronoi_facet *facet;
};


struct tet
{
  u_int tet_index;
  u_int node_index[4];
  int matnum;
  u_int material_count;
  int srfnum;
  struct vector3 cent;
  double r;
  byte delaunay;
  byte tet_type;
};


struct tet_list
{
  struct tet *tet;
  struct tet_list *next;
};


struct tet_mesh
{
  struct node_list *node_head;
  struct node **nodes;
  struct polygon_list *boundary_head;
  struct tet_list *tet_head;
  struct hash_table **edge_hashtab;
  struct vector3 **voronoi_nodes;
  u_int ntets;
  u_int nnodes;
  u_int nvoronoinodes;
  u_int world_boundary_nodes;
  u_int interface_nodes;
  u_int nedges;
  u_int nboundarynodes;
  u_int nboundaryedges;
  u_int nboundaryfaces;
  u_int voronoi_boundary_nedges_tot;
  u_int voronoi_facet_nedges_tot;
  u_int voronoi_cap_nedges_tot;
  u_int voronoi_nfacets_tot;
  u_int voronoi_facet_ntriangles_tot;
  u_int voronoi_cap_ntriangles_tot;
  u_int nmaterials;
  ub4 hashsize;
  ub4 hashmask;
};


struct material_list {
	int matnum;
        struct material_list *next;
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


void output_object();
void output_nurbs();
void output_polyhedron();

#ifndef DEBUG
void no_printf();
#endif

#endif
