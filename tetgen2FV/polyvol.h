#ifndef POLYVOL_H
#define POLYVOL_H

double polyvol_1(struct tet_mesh *tet_mesh, struct node *node);
double polyvol_by_material_1(struct tet_mesh *tet_mesh,
                           struct node *node,
                           u_int node_matindex);
double polyvol_2(struct tet_mesh *tet_mesh, struct node *node);
double polyvol_by_material_2(struct tet_mesh *tet_mesh,
                           struct node *node,
                           u_int node_matindex);
double meshvol(struct tet_mesh *tet_mesh);

#endif
