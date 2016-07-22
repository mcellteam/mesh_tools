#ifndef VORONOI_H
#define VORONOI_H

#include "tetgen2FV.h"

int construct_voronoi_mesh(struct tet_mesh *tet_mesh);

int construct_node_boundary_face_list(struct tet_mesh *tet_mesh);

int construct_edge_list(struct tet_mesh *tet_mesh);

int construct_material_info(struct tet_mesh *tet_mesh);

int output_voronoi_element(struct tet_mesh *tet_mesh, struct node *node);

int output_voronoi_element_by_material(struct tet_mesh *tet_mesh,
                                       struct node *node,
                                       u_int node_matindex);

int output_voronoi_element_dx(struct tet_mesh *tet_mesh,
                              struct node *node,
                              char *outfile_name);

int check_normals(struct tet_mesh *tet_mesh);

#endif
