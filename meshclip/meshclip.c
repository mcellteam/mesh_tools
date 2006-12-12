#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "meshclip.h"
#include "vector.h"

#define NUM_DICES 20
#define EPSILON_1 1.0e-12
#define EPSILON_2 1.0e-9

extern FILE *yyin;
struct object *op;
struct polyhedron *input_mesh,*clipping_mesh;
struct volume *volume;
struct wall *wall_head;
char *infile;
char *input_mesh_filename;
char *clipping_mesh_filename;
char fully_outside_filename[64];
char fully_inside_filename[64];
char on_edge_filename[64];
int line_num;
int skip_freq;
int i,n_part;
int n_ligand_types = 0;
double vol_infinity,pos,pos1,pos2,part_delta,x_delta,y_delta,z_delta;
double *dblp;

int motion_map[6];
int wall_count[27] = {	
  3, 2, 3,  2, 1, 2,  3, 2, 3,	
  2, 1, 2,  1, 0, 1,  2, 1, 2,	
  3, 2, 3,  2, 1, 2,  3, 2, 3	
};	
int wall_map[27][3] = {	
  BOT, FRNT, LFT,  BOT, FRNT, -1,  BOT, FRNT, RT,	
  BOT, LFT, -1,  BOT, -1, -1,  BOT, RT, -1,	
  BOT, BCK, LFT,  BOT, BCK, -1,  BOT, BCK, RT,	

  FRNT, LFT, -1,  FRNT, -1, -1,  FRNT, RT, -1,	
  LFT, -1, -1,  -1, -1, -1,  RT, -1, -1,	
  BCK, LFT, -1,  BCK, -1, -1,  BCK, RT, -1,	
	
  TP, FRNT, LFT,  TP, FRNT, -1,  TP, FRNT, RT,	
  TP, LFT, -1,  TP, -1, -1,  TP, RT, -1,	
  TP, BCK, LFT,  TP, BCK, -1,  TP, BCK, RT	
};	


/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{

	if (argc<3) {
      	  fprintf(stderr,
            "Usage: %s input_mesh_filename clipping_mesh_filename\n",argv[0]);
	  exit(1);
	}

	input_mesh_filename=argv[1];
	clipping_mesh_filename=argv[2];
        sprintf(fully_outside_filename,"fully_outside.m");
        sprintf(fully_inside_filename,"fully_inside.m");
        sprintf(on_edge_filename,"on_edge.m");

        if ((volume=(struct volume *)malloc(sizeof(struct volume)))==NULL) {
          fprintf(stderr,"meshclip: cannot store volume data\n");
          exit(1);
        }      
        if ((volume->x_partitions=(double *)malloc
            ((4+NUM_DICES)*sizeof(double)))==NULL) {
          fprintf(stderr,"meshclip: cannot store volume partitions\n");
          exit(1);
        }      
        if ((volume->y_partitions=(double *)malloc
            ((4+NUM_DICES)*sizeof(double)))==NULL) {
          fprintf(stderr,"meshclip: cannot store volume partitions\n");
          exit(1);
        }      
        if ((volume->z_partitions=(double *)malloc
            ((4+NUM_DICES)*sizeof(double)))==NULL) {
          fprintf(stderr,"meshclip: cannot store volume partitions\n");
          exit(1);
        }      
        vol_infinity=sqrt(DBL_MAX)/4;
        volume->n_x_partitions=4+NUM_DICES;
        volume->n_y_partitions=4+NUM_DICES;
        volume->n_z_partitions=4+NUM_DICES;
        volume->x_walls=NULL;
        volume->y_walls=NULL;
        volume->z_walls=NULL;
        volume->n_subvol=0;
        volume->n_x_subvol=0;
        volume->n_y_subvol=0;
        volume->n_z_subvol=0;
        volume->subvol=NULL;
  
        wall_head = NULL;

	if ((yyin=fopen(input_mesh_filename,"r"))==NULL) {
	  fprintf(stderr,"meshclip: error opening file: %s\n",
            input_mesh_filename);
	  exit(1);
	} 
	fflush(stdout);
        infile=input_mesh_filename;
        line_num=0;
	if (yyparse()) {
	  fprintf(stderr,"meshclip: error parsing file: %s\n",
            input_mesh_filename);
	  exit(1);
	} 
	fclose(yyin);
        input_mesh=(struct polyhedron *)op->contents;

	if ((yyin=fopen(clipping_mesh_filename,"r"))==NULL) {
	  fprintf(stderr,"meshclip: error opening file: %s\n",
            clipping_mesh_filename);
	  exit(1);
	} 
	fflush(stdout);
        infile=clipping_mesh_filename;
        line_num=0;
	if (yyparse()) {
	  fprintf(stderr,"meshclip: error parsing file: %s\n",
            clipping_mesh_filename);
	  exit(1);
	} 
	fclose(yyin);
        clipping_mesh=(struct polyhedron *)op->contents;

        pos1=clipping_mesh->llf.x-(2*EPSILON_2);
        pos2=clipping_mesh->urb.x+(2*EPSILON_2);
        part_delta=(pos2-pos1)/(NUM_DICES+1);
        x_delta=part_delta;
        n_part=volume->n_x_partitions;
        dblp=volume->x_partitions;
        dblp[0]=-vol_infinity;
        dblp[n_part-1]=vol_infinity;
        for (i=1;i<n_part-1;i++) {
          dblp[i]=pos1+(i-1)*part_delta;
        }

        pos1=clipping_mesh->llf.y-(2*EPSILON_2);
        pos2=clipping_mesh->urb.y+(2*EPSILON_2);
        part_delta=(pos2-pos1)/(NUM_DICES+1);
        y_delta=part_delta;
        n_part=volume->n_y_partitions;
        dblp=volume->y_partitions;
        dblp[0]=-vol_infinity;
        dblp[n_part-1]=vol_infinity;
        for (i=1;i<n_part-1;i++) {
          dblp[i]=pos1+(i-1)*part_delta;
        }

        pos1=clipping_mesh->llf.z-(2*EPSILON_2);
        pos2=clipping_mesh->urb.z+(2*EPSILON_2);
        part_delta=(pos2-pos1)/(NUM_DICES+1);
        z_delta=part_delta;
        n_part=volume->n_z_partitions;
        dblp=volume->z_partitions;
        dblp[0]=-vol_infinity;
        dblp[n_part-1]=vol_infinity;
        for (i=1;i<n_part-1;i++) {
          dblp[i]=pos1+(i-1)*part_delta;
        }

	fprintf(stderr,"Partitioning volume...");
        if (partition_volume(volume)) {
	  fprintf(stderr,"meshclip: error partitioning volume\n");
	  exit(1);
        }
	fprintf(stderr,"Done\n");

	fprintf(stderr,"Initializing clipping geometry...");
        if (init_geom(clipping_mesh)) {
	  fprintf(stderr,"meshclip: error initializing geometry\n");
	  exit(1);
        }
	fprintf(stderr,"Done\n");

	fprintf(stderr,"Decomposing volume...");
        if (decompose_volume(volume,wall_head)) {
	  fprintf(stderr,"meshclip: error decomposing volume\n");
	  exit(1);
        }
	fprintf(stderr,"Done\n");

	fprintf(stderr,"Clipping...\n");
        if (clip_mesh(volume,input_mesh)) {
	  fprintf(stderr,"meshclip: error clipping mesh\n");
	  exit(1);
        }
	fprintf(stderr,"Done clipping\n");

	exit(0);
}


void swap_double(x,y)
double *x,*y;
{
  double temp;

  temp=*x;
  *x=*y;
  *y=temp;
}


sort_dbl_array(array,n)
  double *array;
  int n;
{
  int i,done;

  n=n-1;
  done=0;
  while (!done) {
    done=1;
    for (i=0;i<n;i++) {
      if (array[i] > array[i+1]) { 
        done=0;
        swap_double(&array[i],&array[i+1]);
      }
    }
  }
  return;
}


void cube_corners(p1,p2,corner)
struct vector3 *p1,*p2,*corner;
{
  double dx,dy,dz;

  dx=p2->x-p1->x;
  dy=p2->y-p1->y;
  dz=p2->z-p1->z;
  if (dx<0) {
    swap_double(&p1->x,&p2->x);
  }
  if (dy<0) {
    swap_double(&p1->y,&p2->y);
  }
  if (dz<0) {
    swap_double(&p1->z,&p2->z);
  }
  corner[0].x=p1->x;
  corner[0].y=p1->y;
  corner[0].z=p1->z;
  corner[1].x=p1->x;
  corner[1].y=p1->y;
  corner[1].z=p2->z;
  corner[2].x=p1->x;
  corner[2].y=p2->y;
  corner[2].z=p1->z;
  corner[3].x=p1->x;
  corner[3].y=p2->y;
  corner[3].z=p2->z;
  corner[4].x=p2->x;
  corner[4].y=p1->y;
  corner[4].z=p1->z;
  corner[5].x=p2->x;
  corner[5].y=p1->y;
  corner[5].z=p2->z;
  corner[6].x=p2->x;
  corner[6].y=p2->y;
  corner[6].z=p1->z;
  corner[7].x=p2->x;
  corner[7].y=p2->y;
  corner[7].z=p2->z;
}


void cube_face(corner,face,i)
  struct vector3 *corner,**face;
  int i;
{
      /* Build each face using right-hand rule */
      if (i==TP) {
        /* top face */
        face[0]=&corner[1];
        face[1]=&corner[5];
        face[2]=&corner[7];
        face[3]=&corner[3];
      }
      else if (i==BOT) {
        /* bottom face */
        face[0]=&corner[0];
        face[1]=&corner[2];
        face[2]=&corner[6];
        face[3]=&corner[4];
      }
      else if (i==FRNT) {
        /* front face */
        face[0]=&corner[0];
        face[1]=&corner[4];
        face[2]=&corner[5];
        face[3]=&corner[1];
      }
      else if (i==BCK) {
        /* back face */
        face[0]=&corner[2];
        face[1]=&corner[3];
        face[2]=&corner[7];
        face[3]=&corner[6];
      }
      else if (i==LFT) {
        /* left face */
        face[0]=&corner[0];
        face[1]=&corner[1];
        face[2]=&corner[3];
        face[3]=&corner[2];
      }
      else if (i==RT) {
        /* right face */
        face[0]=&corner[4];
        face[1]=&corner[6];
        face[2]=&corner[7];
        face[3]=&corner[5];
      }
      return;
}


void cube_faces(corner,face)
  struct vector3 *corner,*face[6][4];
{

  /* Build each face using right-hand rule */
  /* top face */
  face[0][0]=&corner[1];
  face[0][1]=&corner[5];
  face[0][2]=&corner[7];
  face[0][3]=&corner[3];
 
  /* bottom face */
  face[1][0]=&corner[0];
  face[1][1]=&corner[2];
  face[1][2]=&corner[6];
  face[1][3]=&corner[4];
 
  /* front face */
  face[2][0]=&corner[0];
  face[2][1]=&corner[4];
  face[2][2]=&corner[5];
  face[2][3]=&corner[1];
 
  /* back face */
  face[3][0]=&corner[2];
  face[3][1]=&corner[3];
  face[3][2]=&corner[7];
  face[3][3]=&corner[6];
 
  /* left face */
  face[4][0]=&corner[0];
  face[4][1]=&corner[1];
  face[4][2]=&corner[3];
  face[4][3]=&corner[2];
 
  /* right face */
  face[5][0]=&corner[4];
  face[5][1]=&corner[6];
  face[5][2]=&corner[7];
  face[5][3]=&corner[5];
}


struct wall *init_wall(face,face_vertex_normal,n_verts)
struct vector3 **face;
struct vector3 **face_vertex_normal;
int n_verts;
{
  struct wall *wp;
  struct vector3 v1,v2,nv1,nv2;
  double *dblp;
  double tmp;
  int i;

  if ((wp=(struct wall *)malloc(sizeof(struct wall)))==NULL) {
    return(NULL);
  }
  if ((dblp=(double *)malloc(n_verts*sizeof(double)))==NULL) {
    return(NULL);
  }
  wp->wall_type=NULL;
  wp->side=0;
  wp->n_vert=n_verts;
  wp->vert=face;
  wp->vert_normal=face_vertex_normal;
  wp->length=dblp;
  wp->next_wall=NULL;
  vectorize(wp->vert[0],wp->vert[1],&v1);
  vectorize(wp->vert[0],wp->vert[n_verts-1],&v2);
  nv1.x=v1.x;
  nv1.y=v1.y;
  nv1.z=v1.z;
  nv2.x=v2.x;
  nv2.y=v2.y;
  nv2.z=v2.z;
  normalize(&nv1);
  normalize(&nv2);
  cross_prod(&nv1,&nv2,&wp->normal);
  normalize(&wp->normal);
  wp->length[0]=vect_length(&v1);
  wp->length[n_verts-1]=vect_length(&v2);
  wp->length_first=wp->length[0];
  wp->length_last=wp->length[n_verts-1];
  wp->r_length_first=1.0/wp->length[0];
  wp->r_length_last=1.0/wp->length[n_verts-1];
  for (i=1;i<n_verts-1;i++) {
    vectorize(wp->vert[i],wp->vert[i+1],&v1);
    wp->length[i]=vect_length(&v1);
  }
  wp->d=wp->normal.x*wp->vert[0]->x+wp->normal.y*wp->vert[0]->y
    +wp->normal.z*wp->vert[0]->z;
  wp->projection=0;
  tmp = wp->normal.x;
  if (tmp*tmp < wp->normal.y*wp->normal.y) {
    wp->projection=1;
    tmp = wp->normal.y;
  }
  if (tmp*tmp < wp->normal.z*wp->normal.z) {
    wp->projection=2;
  }
  return(wp);
}


int partition_volume(volp)
  struct volume *volp;
{
  struct subvolume *subvolp;
  struct vector3 p1,p2;
  struct vector3 *verts,**face;
  struct wall_list *wlp;
  struct wall *wp;
  int i,j,k,l,m,n,nx,ny,nz,nxy,nx_parts,ny_parts,nz_parts,n_subvol;
 
  nx_parts=volp->n_x_partitions;
  ny_parts=volp->n_y_partitions;
  nz_parts=volp->n_z_partitions;
  p1.x=-vol_infinity;
  p1.y=-vol_infinity;
  p1.z=-vol_infinity;
  p2.x=vol_infinity;
  p2.y=vol_infinity;
  p2.z=vol_infinity;
  if (nx_parts==0) {
    if ((volp->x_partitions=(double *)malloc(2*sizeof(double)))==NULL) {
      fprintf(stderr,"meshclip: Cannot store volume data\n");
      return(1);
    }       
    volp->n_x_partitions=2;
    volp->x_partitions[0]=p1.x;
    volp->x_partitions[1]=p2.x;
  }
  if (ny_parts==0) {
    if ((volp->y_partitions=(double *)malloc(2*sizeof(double)))==NULL) {
      fprintf(stderr,"meshclip: Cannot store volume data\n");
      return(1);
    }       
    volp->n_y_partitions=2;
    volp->y_partitions[0]=p1.y;
    volp->y_partitions[1]=p2.y;
  }
  if (nz_parts==0) {
    if ((volp->z_partitions=(double *)malloc(2*sizeof(double)))==NULL) {
      fprintf(stderr,"meshclip: Cannot store volume data\n");
      return(1);
    }       
    volp->n_z_partitions=2;
    volp->z_partitions[0]=p1.z;
    volp->z_partitions[1]=p2.z;
  }
  cube_corners(&p1,&p2,volp->corner);
  nx_parts=volp->n_x_partitions;
  ny_parts=volp->n_y_partitions;
  nz_parts=volp->n_z_partitions;
  nx=nx_parts-1;
  ny=ny_parts-1;
  nz=nz_parts-1;
  nxy=nx*ny;
  n_subvol=nx*ny*nz;
  volp->n_x_subvol=nx;
  volp->n_y_subvol=ny;
  volp->n_z_subvol=nz;
  volp->n_subvol=n_subvol;

  if ((volp->x_walls=(struct wall **)malloc
       (nx_parts*sizeof(struct wall *)))==NULL) {
    fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
    return(1);
  }
  if ((volp->y_walls=(struct wall **)malloc
       (ny_parts*sizeof(struct wall *)))==NULL) {
    fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
    return(1);
  }
  if ((volp->z_walls=(struct wall **)malloc
       (nz_parts*sizeof(struct wall *)))==NULL) {
    fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
    return(1);
  }
  if ((subvolp=(struct subvolume *)malloc
       (n_subvol*sizeof(struct subvolume)))==NULL) {
    fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
    return(1);
  }
  volp->subvol=subvolp;
  for (i=0;i<n_subvol;i++) {
    subvolp[i].wall_list=NULL;
    for (j=0;j<6;j++) {
      subvolp[i].walls[j]=NULL;
    }
  }

  if ((verts=(struct vector3 *)malloc
       (4*sizeof(struct vector3)))==NULL) {
    fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
    return(1);
  }
  if ((face=(struct vector3 **)malloc
       (4*sizeof(struct vector3 *)))==NULL) {
    fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
    return(1);
  }
  for (i=0;i<4;i++) {
    face[i]=&verts[i];
  }
  for (i=0;i<nx_parts;i++) {
    verts[0].x=volp->x_partitions[i];
    verts[0].y=volp->y_partitions[0];
    verts[0].z=volp->z_partitions[0];

    verts[1].x=volp->x_partitions[i];
    verts[1].y=volp->y_partitions[0];
    verts[1].z=volp->z_partitions[nz];

    verts[2].x=volp->x_partitions[i];
    verts[2].y=volp->y_partitions[ny];
    verts[2].z=volp->z_partitions[nz];

    verts[3].x=volp->x_partitions[i];
    verts[3].y=volp->y_partitions[ny];
    verts[3].z=volp->z_partitions[0];
    if ((wp=init_wall(face,NULL,4))==NULL) {
      fprintf(stderr,"meshclip: cannot store wall data");
      return(1);
    }
    if ((wp->wall_type=(byte *)malloc
         ((1+n_ligand_types)*sizeof(byte)))==NULL) {
      fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
      return(1);
    }
    for (j=0;j<(1+n_ligand_types);j++) {
      wp->wall_type[j]=SUBVOL;
    }
    volp->x_walls[i]=wp;
  }

  if ((verts=(struct vector3 *)malloc
       (4*sizeof(struct vector3)))==NULL) {
    fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
    return(1);
  }
  if ((face=(struct vector3 **)malloc
       (4*sizeof(struct vector3 *)))==NULL) {
    fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
    return(1);
  }
  for (i=0;i<4;i++) {
    face[i]=&verts[i];
  }
  for (i=0;i<ny_parts;i++) {
    verts[0].x=volp->x_partitions[0];
    verts[0].y=volp->y_partitions[i];
    verts[0].z=volp->z_partitions[0];

    verts[1].x=volp->x_partitions[nx];
    verts[1].y=volp->y_partitions[i];
    verts[1].z=volp->z_partitions[0];

    verts[2].x=volp->x_partitions[nx];
    verts[2].y=volp->y_partitions[i];
    verts[2].z=volp->z_partitions[nz];

    verts[3].x=volp->x_partitions[0];
    verts[3].y=volp->y_partitions[i];
    verts[3].z=volp->z_partitions[nz];
    if ((wp=init_wall(face,NULL,4))==NULL) {
      fprintf(stderr,"meshclip: cannot store wall data");
      return(1);
    }
    if ((wp->wall_type=(byte *)malloc
         ((1+n_ligand_types)*sizeof(byte)))==NULL) {
      fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
      return(1);
    }
    for (j=0;j<(1+n_ligand_types);j++) {
      wp->wall_type[j]=SUBVOL;
    }
    volp->y_walls[i]=wp;
  }

  if ((verts=(struct vector3 *)malloc
       (4*sizeof(struct vector3)))==NULL) {
    fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
    return(1);
  }
  if ((face=(struct vector3 **)malloc
       (4*sizeof(struct vector3 *)))==NULL) {
    fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
    return(1);
  }
  for (i=0;i<4;i++) {
    face[i]=&verts[i];
  }
  for (i=0;i<nz_parts;i++) {
    verts[0].x=volp->x_partitions[0];
    verts[0].y=volp->y_partitions[0];
    verts[0].z=volp->z_partitions[i];

    verts[1].x=volp->x_partitions[0];
    verts[1].y=volp->y_partitions[ny];
    verts[1].z=volp->z_partitions[i];

    verts[2].x=volp->x_partitions[nx];
    verts[2].y=volp->y_partitions[ny];
    verts[2].z=volp->z_partitions[i];

    verts[3].x=volp->x_partitions[nx];
    verts[3].y=volp->y_partitions[0];
    verts[3].z=volp->z_partitions[i];
    if ((wp=init_wall(face,NULL,4))==NULL) {
      fprintf(stderr,"meshclip: cannot store wall data");
      return(1);
    }
    if ((wp->wall_type=(byte *)malloc
         ((1+n_ligand_types)*sizeof(byte)))==NULL) {
      fprintf(stderr,"meshclip: Cannot store subvolume partition data\n");
      return(1);
    }
    for (j=0;j<(1+n_ligand_types);j++) {
      wp->wall_type[j]=SUBVOL;
    }
    volp->z_walls[i]=wp;
  }

  l=0;
  for (k=0;k<nz;k++) {
    for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
        subvolp[l].x1=volp->x_partitions[i];
        subvolp[l].x2=volp->x_partitions[i+1];
        subvolp[l].y1=volp->y_partitions[j];
        subvolp[l].y2=volp->y_partitions[j+1];
        subvolp[l].z1=volp->z_partitions[k];
        subvolp[l].z2=volp->z_partitions[k+1];
        if (k<(nz-1)) {
          subvolp[l].walls[TP]=volp->z_walls[k+1];
        }
        else {
          subvolp[l].walls[TP]=NULL;
        }
        if (k>0) {
          subvolp[l].walls[BOT]=volp->z_walls[k];
        }
        else {
          subvolp[l].walls[BOT]=NULL;
        }
        if (j>0) {
          subvolp[l].walls[FRNT]=volp->y_walls[j];
        }
        else {
          subvolp[l].walls[FRNT]=NULL;
        }
        if (j<(ny-1)) {
        subvolp[l].walls[BCK]=volp->y_walls[j+1];
        }
        else {
          subvolp[l].walls[BCK]=NULL;
        }
        if (i>0) {
          subvolp[l].walls[LFT]=volp->x_walls[i];
        }
        else {
          subvolp[l].walls[LFT]=NULL;
        }
        if (i<(nx-1)) {
          subvolp[l].walls[RT]=volp->x_walls[i+1];
        }
        else {
          subvolp[l].walls[RT]=NULL;
        }
        subvolp[l].wall_list=NULL;

        l++;
      }
    }
  }

/* indexed by wall_index [0-5] for subvolume wall faces */
  motion_map[0]=nxy;
  motion_map[1]=-nxy;
  motion_map[2]=-nx;
  motion_map[3]=nx;
  motion_map[4]=-1;
  motion_map[5]=1;

  return(0);
}


int init_geom(php)
  struct polyhedron *php;
{
  struct polygon_list *plp;
  struct polygon *pop;
  struct wall *wp;
  int j;

  plp=php->polygon_list;
  while (plp!=NULL) {
    pop=plp->polygon;
    if ((wp=init_wall(pop->vertex_array,NULL,pop->n_verts))==NULL) {
      fprintf(stderr,"meshclip: cannot store wall data");
      return(1);
    }
    if ((wp->wall_type=(byte *)malloc
         ((1+n_ligand_types)*sizeof(byte)))==NULL) {
      fprintf(stderr,"meshclip: Cannot store wall data\n");
      return(1);
    }
    for (j=0;j<(1+n_ligand_types);j++) {
      wp->wall_type[j]=RFLCT;
    }
    wp->next_wall=wall_head;
    wall_head=wp;
    plp=plp->next;
  }

  return(0);
}


/* Find which pair of partitions a point lies between. */
/*   Return index of lower partition */
/*   or return -index if point is coincident with a partition. */ 
int find_range(u,u_range,n_u_range)
double u;
double *u_range;
int n_u_range;
{
  int n_lower,n_upper,n_mid;

  n_lower=0;
  n_upper=n_u_range-1;
  while (n_lower!=(n_upper-1)) {
    n_mid=(n_lower+n_upper)/2;
    if (u>u_range[n_mid]) {
      n_lower=n_mid;
    }
    else if (u<u_range[n_mid]) {
      n_upper=n_mid;
    }
    else {
      return(-n_mid);
    }
  }
  if (u==u_range[n_lower]) {
    return(-n_lower);
  }
  if (u==u_range[n_lower+1]) {
    return(-(n_lower+1));
  }

  return(n_lower);
}

/* find range of partitions spanned by a polygon */
find_polygon_range(umin_part,umax_part,vmin_part,vmax_part,
                   wmin_part,wmax_part,volp,wp)
  int *umin_part,*umax_part,*vmin_part,*vmax_part,*wmin_part,*wmax_part;
  struct volume *volp;
  struct wall *wp;
{
  double u,v,w,umin_bb,umax_bb,vmin_bb,vmax_bb,wmin_bb,wmax_bb;
  double pos,new_pos;
  int i; 

  /* form bounding box of polygon */ 
  umin_bb=DBL_MAX;
  umax_bb=-DBL_MAX;
  vmin_bb=DBL_MAX;
  vmax_bb=-DBL_MAX;
  wmin_bb=DBL_MAX;
  wmax_bb=-DBL_MAX;
  for (i=0;i<wp->n_vert;i++) {
    u=wp->vert[i]->x;
    v=wp->vert[i]->y;
    w=wp->vert[i]->z;
    if (u<umin_bb) {
      umin_bb=u;
    }
    if (u>umax_bb) {
      umax_bb=u;
    }
    if (v<vmin_bb) {
      vmin_bb=v;
    }
    if (v>vmax_bb) {
      vmax_bb=v;
    }
    if (w<wmin_bb) {
      wmin_bb=w;
    }
    if (w>wmax_bb) {
      wmax_bb=w;
    }
  }
    
  /* find partition range of polygon bounding box */
  *umin_part=find_range(umin_bb,volp->x_partitions,volp->n_x_partitions);
  *umax_part=find_range(umax_bb,volp->x_partitions,volp->n_x_partitions);
  *vmin_part=find_range(vmin_bb,volp->y_partitions,volp->n_y_partitions);
  *vmax_part=find_range(vmax_bb,volp->y_partitions,volp->n_y_partitions);
  *wmin_part=find_range(wmin_bb,volp->z_partitions,volp->n_z_partitions);
  *wmax_part=find_range(wmax_bb,volp->z_partitions,volp->n_z_partitions);

  return;
}


repartition_volume(volp)
  struct volume *volp;
{
  struct subvolume *subvolp;
  struct vector3 p1,p2;
  struct vector3 **verts;
  struct wall_list *wlp;
  struct wall *wp;
  int i,j,k,l,m,n,nx,ny,nz,nxy,nx_parts,ny_parts,nz_parts,n_subvol,adjacent;

  sort_dbl_array(volp->x_partitions);
  sort_dbl_array(volp->y_partitions);
  sort_dbl_array(volp->z_partitions);

  nx_parts=volp->n_x_partitions;
  ny_parts=volp->n_y_partitions;
  nz_parts=volp->n_z_partitions;
  nx=nx_parts-1;
  ny=ny_parts-1;
  nz=nz_parts-1;
  nxy=nx*ny;
  n_subvol=nx*ny*nz;
  volp->n_x_subvol=nx;
  volp->n_y_subvol=ny;
  volp->n_z_subvol=nz;
  volp->n_subvol=n_subvol;

  for (i=0;i<nx_parts;i++) {
    wp=volp->x_walls[i];
    verts=wp->vert;
    verts[0]->x=volp->x_partitions[i];
    verts[0]->y=volp->y_partitions[0];
    verts[0]->z=volp->z_partitions[0];
      
    verts[1]->x=volp->x_partitions[i];
    verts[1]->y=volp->y_partitions[0];
    verts[1]->z=volp->z_partitions[nz];
      
    verts[2]->x=volp->x_partitions[i];
    verts[2]->y=volp->y_partitions[ny];
    verts[2]->z=volp->z_partitions[nz];
      
    verts[3]->x=volp->x_partitions[i];
    verts[3]->y=volp->y_partitions[ny];
    verts[3]->z=volp->z_partitions[0];
    wp->d=wp->normal.x*wp->vert[0]->x+wp->normal.y*wp->vert[0]->y
      +wp->normal.z*wp->vert[0]->z;
  } 
  for (i=0;i<ny_parts;i++) {
    wp=volp->y_walls[i];
    verts=wp->vert;
    verts[0]->x=volp->x_partitions[0];
    verts[0]->y=volp->y_partitions[i];
    verts[0]->z=volp->z_partitions[0];
      
    verts[1]->x=volp->x_partitions[nx];
    verts[1]->y=volp->y_partitions[i];
    verts[1]->z=volp->z_partitions[0];
      
    verts[2]->x=volp->x_partitions[nx];
    verts[2]->y=volp->y_partitions[i];
    verts[2]->z=volp->z_partitions[nz];
      
    verts[3]->x=volp->x_partitions[0];
    verts[3]->y=volp->y_partitions[i];
    verts[3]->z=volp->z_partitions[nz];
    wp->d=wp->normal.x*wp->vert[0]->x+wp->normal.y*wp->vert[0]->y
      +wp->normal.z*wp->vert[0]->z;
  } 
  for (i=0;i<nz_parts;i++) {
    wp=volp->z_walls[i];
    verts=wp->vert;
    verts[0]->x=volp->x_partitions[0];
    verts[0]->y=volp->y_partitions[0];
    verts[0]->z=volp->z_partitions[i];
  
    verts[1]->x=volp->x_partitions[0];
    verts[1]->y=volp->y_partitions[ny];
    verts[1]->z=volp->z_partitions[i];
  
    verts[2]->x=volp->x_partitions[nx];
    verts[2]->y=volp->y_partitions[ny];
    verts[2]->z=volp->z_partitions[i];
  
    verts[3]->x=volp->x_partitions[nx];
    verts[3]->y=volp->y_partitions[0];
    verts[3]->z=volp->z_partitions[i];
    wp->d=wp->normal.x*wp->vert[0]->x+wp->normal.y*wp->vert[0]->y
      +wp->normal.z*wp->vert[0]->z;
  } 
  l=0;
  subvolp=volp->subvol;
  for (k=0;k<nz;k++) {
    for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
        subvolp[l].x1=volp->x_partitions[i];
        subvolp[l].x2=volp->x_partitions[i+1];
        subvolp[l].y1=volp->y_partitions[j];
        subvolp[l].y2=volp->y_partitions[j+1];
        subvolp[l].z1=volp->z_partitions[k];
        subvolp[l].z2=volp->z_partitions[k+1];
        if (k<(nz-1)) {
          subvolp[l].walls[TP]=volp->z_walls[k+1];
        }
        else {
          subvolp[l].walls[TP]=NULL;
        }
        if (k>0) {
          subvolp[l].walls[BOT]=volp->z_walls[k];
        }
        else {
          subvolp[l].walls[BOT]=NULL;
        }
        if (j>0) {
          subvolp[l].walls[FRNT]=volp->y_walls[j];
        }
        else {
          subvolp[l].walls[FRNT]=NULL;
        }
        if (j<(ny-1)) {
        subvolp[l].walls[BCK]=volp->y_walls[j+1];
        }
        else {
          subvolp[l].walls[BCK]=NULL;
        }
        if (i>0) {
          subvolp[l].walls[LFT]=volp->x_walls[i];
        }
        else {
          subvolp[l].walls[LFT]=NULL;
        }
        if (i<(nx-1)) {
          subvolp[l].walls[RT]=volp->x_walls[i+1];
        }
        else {
          subvolp[l].walls[RT]=NULL;
        }
        l++;
      }
    }
  }
  
  return;
}


/* determine if a 2-D point lies within the bounds of a 2-D clipping window */
int clip_point_2D(umin,umax,vmin,vmax,u1,v1)
  double umin,umax,vmin,vmax,u1,v1;
{

  if (u1>=umin) {
    if (u1<=umax) {
      if (v1>=vmin) {
        if (v1<=vmax) {
          return(1);
        }
      }
    }
  }
  return(0);
}


/* determine if a 3-D point lies within the bounds of a 3-D clipping box */
int clip_point_3D(umin,umax,vmin,vmax,wmin,wmax,u1,v1,w1)
  double umin,umax,vmin,vmax,wmin,wmax,u1,v1,w1;
{

  if (clip_point_2D(umin,umax,vmin,vmax,u1,v1)) {
    if (clip_point_2D(umin,umax,wmin,wmax,u1,w1)) {
      if (clip_point_2D(vmin,vmax,wmin,wmax,v1,w1)) {
        return(1);
      }
    }
  }
  return(0);
}


int clip_line(p,q,r1,r2) 
  double p,q,*r1,*r2;
{
  double r;

  if (p>0.0) {
    r=q/p;
    if (r<(*r1)) {
      return(0);
    }
    if (r<(*r2)) {
      *r2=r;
      return(1);
    }
  }
  else if (p<0.0) {
    r=q/p;
    if (r>(*r2)) {
      return(0);
    }
    if (r>(*r1)) {
      *r1=r;
      return(1);
    }
  }
  else if (q<0.0) {
    return(0);
  }
  return(1);
}


int clip_polygon(subvol,verts,nverts)
  struct subvolume *subvol;
  struct vector3 **verts;
  int nverts;
{
  double x1,x2,y1,y2,xw,yw,u1,u2,v1,v2,du,dv,b,m,t;
  double umin,umax,vmin,vmax,umin_bb,umax_bb,vmin_bb,vmax_bb;
  double r1,r2;
  int i,intersect;

/* intersect xy projections of subvol and polygon */
  umin=subvol->x1;
  umax=subvol->x2;
  vmin=subvol->y1;
  vmax=subvol->y2;
  intersect=0;
  i=0;
  while (intersect!=1 && i<nverts) {
    if (i==(nverts-1)) {
      u1=verts[i]->x;
      u2=verts[0]->x;
      v1=verts[i]->y;
      v2=verts[0]->y;
    } 
    else {
      u1=verts[i]->x;
      u2=verts[i+1]->x;
      v1=verts[i]->y;
      v2=verts[i+1]->y;
    }
    du=u2-u1;
    dv=v2-v1;

    if (clip_point_2D(umin,umax,vmin,vmax,u1,v1)) {
      intersect=1;
    }
    else if (clip_point_2D(umin,umax,vmin,vmax,u2,v2)) {
      intersect=1;
    }
    else if (du==0 && dv==0 && clip_point_2D(umin,umax,vmin,vmax,u1,v1)) {
      intersect=1;
    }
    else {
      r1=0.0;
      r2=1.0;
      if (clip_line(-du,u1-umin,&r1,&r2)) {
        if (clip_line(du,umax-u1,&r1,&r2)) {
          if (clip_line(-dv,v1-vmin,&r1,&r2)) {
            if (clip_line(dv,vmax-v1,&r1,&r2)) {
              intersect=1;
            }
          }
        }
      }
    }
    i++;
  } 
  if (intersect!=1) {
  /* Check to see if the bounding box of projected polygon completely */
  /*   surrounds the bounding box of the projected subvolume */
    /* form bounding box of projected polygon */
    umin_bb=DBL_MAX;
    umax_bb=-DBL_MAX;
    vmin_bb=DBL_MAX;
    vmax_bb=-DBL_MAX;
    for (i=0;i<nverts;i++) {
      u1=verts[i]->x;
      v1=verts[i]->y;
      if (u1<umin_bb) {
        umin_bb=u1;
      }
      if (u1>umax_bb) {
        umax_bb=u1;
      }
      if (v1<vmin_bb) {
        vmin_bb=v1;
      }
      if (v1>vmax_bb) {
        vmax_bb=v1;
      }
    }
    if (clip_point_2D(umin_bb,umax_bb,vmin_bb,vmax_bb,umin,vmin)) {
      if (clip_point_2D(umin_bb,umax_bb,vmin_bb,vmax_bb,umin,vmax)) {
        if (clip_point_2D(umin_bb,umax_bb,vmin_bb,vmax_bb,umax,vmax)) {
          if (clip_point_2D(umin_bb,umax_bb,vmin_bb,vmax_bb,umax,vmin)) {
            intersect=1;
          }
        }
      }
    }
    if (intersect!=1) {
      return(0);
    }
  }

/* intersect xz projections of subvol and polygon */
  umin=subvol->x1;
  umax=subvol->x2;
  vmin=subvol->z1;
  vmax=subvol->z2;
  intersect=0;
  i=0;
  while (intersect!=1 && i<nverts) {
    if (i==(nverts-1)) {
      u1=verts[i]->x;
      u2=verts[0]->x;
      v1=verts[i]->z;
      v2=verts[0]->z;
    } 
    else {
      u1=verts[i]->x;
      u2=verts[i+1]->x;
      v1=verts[i]->z;
      v2=verts[i+1]->z;
    }
    du=u2-u1;
    dv=v2-v1;
    if (clip_point_2D(umin,umax,vmin,vmax,u1,v1)) {
      intersect=1;
    }
    else if (clip_point_2D(umin,umax,vmin,vmax,u2,v2)) {
      intersect=1;
    }
    else if (du==0 && dv==0 && clip_point_2D(umin,umax,vmin,vmax,u1,v1)) {
      intersect=1;
    }
    else {
      r1=0;
      r2=1;
      if (clip_line(-du,u1-umin,&r1,&r2)) {
        if (clip_line(du,umax-u1,&r1,&r2)) {
          if (clip_line(-dv,v1-vmin,&r1,&r2)) {
            if (clip_line(dv,vmax-v1,&r1,&r2)) {
              intersect=1;
            }
          }
        }
      }
    }
    i++;
  } 
  if (intersect!=1) {
  /* Check to see if the bounding box of projected polygon completely */
  /*   surrounds the bounding box of the projected subvolume */
    /* form bounding box of projected polygon */
    umin_bb=DBL_MAX;
    umax_bb=-DBL_MAX;
    vmin_bb=DBL_MAX;
    vmax_bb=-DBL_MAX;
    for (i=0;i<nverts;i++) {
      u1=verts[i]->x;
      v1=verts[i]->z;
      if (u1<umin_bb) {
        umin_bb=u1;
      }
      if (u1>umax_bb) {
        umax_bb=u1;
      }
      if (v1<vmin_bb) {
        vmin_bb=v1;
      }
      if (v1>vmax_bb) {
        vmax_bb=v1;
      }
    }
    if (clip_point_2D(umin_bb,umax_bb,vmin_bb,vmax_bb,umin,vmin)) {
      if (clip_point_2D(umin_bb,umax_bb,vmin_bb,vmax_bb,umin,vmax)) {
        if (clip_point_2D(umin_bb,umax_bb,vmin_bb,vmax_bb,umax,vmax)) {
          if (clip_point_2D(umin_bb,umax_bb,vmin_bb,vmax_bb,umax,vmin)) {
            intersect=1;
          }
        }
      }
    }
    if (intersect!=1) {
      return(0);
    }
  }

/* intersect yz projections of subvol and polygon */
  umin=subvol->y1;
  umax=subvol->y2;
  vmin=subvol->z1;
  vmax=subvol->z2;
  intersect=0;
  i=0;
  while (intersect!=1 && i<nverts) {
    if (i==(nverts-1)) {
      u1=verts[i]->y;
      u2=verts[0]->y;
      v1=verts[i]->z;
      v2=verts[0]->z;
    } 
    else {
      u1=verts[i]->y;
      u2=verts[i+1]->y;
      v1=verts[i]->z;
      v2=verts[i+1]->z;
    }
    du=u2-u1;
    dv=v2-v1;
    if (clip_point_2D(umin,umax,vmin,vmax,u1,v1)) {
      intersect=1;
    }
    else if (clip_point_2D(umin,umax,vmin,vmax,u2,v2)) {
      intersect=1;
    }
    else if (du==0 && dv==0 && clip_point_2D(umin,umax,vmin,vmax,u1,v1)) {
      intersect=1;
    }
    else {
      r1=0;
      r2=1;
      if (clip_line(-du,u1-umin,&r1,&r2)) {
        if (clip_line(du,umax-u1,&r1,&r2)) {
          if (clip_line(-dv,v1-vmin,&r1,&r2)) {
            if (clip_line(dv,vmax-v1,&r1,&r2)) {
              intersect=1;
            }
          }
        }
      }
    }
    i++;
  } 
  if (intersect!=1) {
  /* Check to see if the bounding box of projected polygon completely */
  /*   surrounds the bounding box of the projected subvolume */
    /* form bounding box of projected polygon */
    umin_bb=DBL_MAX;
    umax_bb=-DBL_MAX;
    vmin_bb=DBL_MAX;
    vmax_bb=-DBL_MAX;
    for (i=0;i<nverts;i++) {
      u1=verts[i]->y;
      v1=verts[i]->z;
      if (u1<umin_bb) {
        umin_bb=u1;
      }
      if (u1>umax_bb) {
        umax_bb=u1;
      }
      if (v1<vmin_bb) {
        vmin_bb=v1;
      }
      if (v1>vmax_bb) {
        vmax_bb=v1;
      }
    }
    if (clip_point_2D(umin_bb,umax_bb,vmin_bb,vmax_bb,umin,vmin)) {
      if (clip_point_2D(umin_bb,umax_bb,vmin_bb,vmax_bb,umin,vmax)) {
        if (clip_point_2D(umin_bb,umax_bb,vmin_bb,vmax_bb,umax,vmax)) {
          if (clip_point_2D(umin_bb,umax_bb,vmin_bb,vmax_bb,umax,vmin)) {
            intersect=1;
          }
        }
      }
    }
    if (intersect!=1) {
      return(0);
    }
  }
  
  return(1);

}


/* Decompose object walls into spatial partitions */
/*   Return 0 on success */
/*   Return 1 on malloc failure */
int decompose_volume(volp,wp)
struct volume *volp;
struct wall *wp;
{
  struct subvolume *subvolp;
  struct wall *wp_start;
  struct wall_list *wlp,*sv_wlp;
  struct vector3 llf_vert,urb_vert;
  double u,v,w;
  double pos,new_pos;
  int umin_part,umax_part,vmin_part,vmax_part,wmin_part,wmax_part;
  int i,j,k,l,m,transp,llf_subvol,urb_subvol,reposition; 
  int nx,nxy;

  wp_start=wp;
  reposition=0;
  while (wp!=NULL) {
    /* find partition range of polygon */
    find_polygon_range(&umin_part,&umax_part,&vmin_part,&vmax_part,
                       &wmin_part,&wmax_part,volp,wp);
    if (umin_part<0) {
      pos=volp->x_partitions[-umin_part];
      new_pos=pos-(2.0*EPSILON_2);
      printf("meshclip warning: to avoid coincidence, X partition at %1.15g repositioned to %1.15g\n",pos,new_pos);
      volp->x_partitions[-umin_part]=new_pos;
      reposition=1;
    }
    if (umax_part<0) {
      pos=volp->x_partitions[-umax_part];
      new_pos=pos+(2.0*EPSILON_2);
      printf("meshclip warning: to avoid coincidence, X partition at %1.15g repositioned to %1.15g\n",pos,new_pos);
      volp->x_partitions[-umax_part]=new_pos;
      reposition=1;
    }
    if (vmin_part<0) {
      pos=volp->y_partitions[-vmin_part];
      new_pos=pos-(2.0*EPSILON_2);
      printf("meshclip warning: to avoid coincidence, Y partition at %1.15g repositioned to %1.15g\n",pos,new_pos);
      volp->y_partitions[-vmin_part]=new_pos;
      reposition=1;
    }
    if (vmax_part<0) {
      pos=volp->y_partitions[-vmax_part];
      new_pos=pos+(2.0*EPSILON_2);
      printf("meshclip warning: to avoid coincidence, Y partition at %1.15g repositioned to %1.15g\n",pos,new_pos);
      volp->y_partitions[-vmax_part]=new_pos;
      reposition=1;
    }
    if (umin_part<0) {
      pos=volp->x_partitions[-umin_part];
      new_pos=pos-(2.0*EPSILON_2);
      printf("meshclip warning: to avoid coincidence, X partition at %1.15g repositioned to %1.15g\n",pos,new_pos);
      volp->x_partitions[-umin_part]=new_pos;
      reposition=1;
    }
    if (wmax_part<0) {
      pos=volp->z_partitions[-wmax_part];
      new_pos=pos+(2.0*EPSILON_2);
      printf("meshclip warning: to avoid coincidence, Z partition at %1.15g repositioned to %1.15g\n",pos,new_pos);
      volp->z_partitions[-wmax_part]=new_pos;
      reposition=1;
    }

    if (reposition) {
      repartition_volume(volp);
      printf("meshclip: Restarting volume decomposition process...\n");
      reposition=0;
      wp=wp_start;
    }
    else {
      wp=wp->next_wall;
    }
  }

  nx=volp->n_x_subvol;
  nxy=nx*volp->n_y_subvol;
  wp=wp_start;
  while (wp!=NULL) {
    /* find partition range of polygon */
    find_polygon_range(&umin_part,&umax_part,&vmin_part,&vmax_part,
                       &wmin_part,&wmax_part,volp,wp);
    /* find all subvolumes which this polygon intersects */
    for (k=wmin_part;k<=wmax_part;k++) {
      for (j=vmin_part;j<=vmax_part;j++) {
        for (i=umin_part;i<=umax_part;i++) {
          l=(k*nxy)+(j*nx)+i;
          subvolp=&volp->subvol[l];

          /* if polygon intersects this subvolume */
          if (clip_polygon(subvolp,wp->vert,wp->n_vert)) {
            if ((wlp=(struct wall_list *)malloc
	         (sizeof(struct wall_list)))==NULL) {
              return(1);
            }
            wlp->wall=wp;
            transp=0;
            for (m=0;m<1+n_ligand_types;m++) {
              transp = (transp || (wp->wall_type[m]==TRANSP));
            }
            /* put transparent walls at head of list */
            if (transp) {
              wlp->next=subvolp->wall_list;
              subvolp->wall_list=wlp;
            }
            /* otherwise put wall at end of list */
            else {
              sv_wlp=subvolp->wall_list;
              if (sv_wlp!=NULL) {
	        while (sv_wlp->next!=NULL) {
	          sv_wlp=sv_wlp->next;
	        }
	        sv_wlp->next=wlp;
	        wlp->next=NULL;
              }
              else {
	        wlp->next=subvolp->wall_list;
	        subvolp->wall_list=wlp;
              }
            }
          }

        }
      }
    }
    wp=wp->next_wall;
  }

  return(0);
}

/* Find which subvolume a point lies within. */
/*   Return index of subvolume */
/*   or return -index if point is coincident with a subvolume wall. */ 
/*   For speed, first check to see if point lies within a */
/*   subvolume specified by "guess".  Skip this check if guess=-1. */
int find_subvol(volp,p,guess)
struct volume *volp;
struct vector3 *p;
int guess;
{
struct subvolume *subvolp;
int i,subvol,n_x_part,n_y_part,n_z_part,n_x_subvol,n_y_subvol,n_z_subvol;
int x_part,y_part,z_part;
double u,v,w,umin,umax,vmin,vmax,wmin,wmax;

  u=p->x;
  v=p->y;
  w=p->z;
  
  subvolp=volp->subvol;
  if (guess!=-1) {
    umin=subvolp[guess].x1;
    umax=subvolp[guess].x2;
    vmin=subvolp[guess].y1;
    vmax=subvolp[guess].y2;
    wmin=subvolp[guess].z1;
    wmax=subvolp[guess].z2;
    if (clip_point_2D(umin,umax,vmin,vmax,u,v)) {
      if (clip_point_2D(umin,umax,wmin,wmax,u,w)) {
        if (clip_point_2D(vmin,vmax,wmin,wmax,v,w)) {
          if (u==umin || u==umax || v==vmin || v==vmax || w==wmin || w==wmax) {
            return(-guess);
          }
          return(guess);
        } 
      } 
    } 
  }

  n_x_part=volp->n_x_partitions;
  n_y_part=volp->n_y_partitions;
  n_z_part=volp->n_z_partitions;
  n_x_subvol=volp->n_x_subvol;
  n_y_subvol=volp->n_y_subvol;
  n_z_subvol=volp->n_z_subvol;
  x_part=find_range(u,volp->x_partitions,n_x_part);
  y_part=find_range(v,volp->y_partitions,n_y_part);
  z_part=find_range(w,volp->z_partitions,n_z_part);
  subvol=(z_part*n_x_subvol*n_y_subvol)+(y_part*n_x_subvol)+x_part;
  if (x_part<0 || y_part<0 || z_part<0) {
    return(-subvol);
  }
  return(subvol);
}


int ray_trace (point,ray,subvol)
  struct vector3 *point,*ray;
  int subvol;
{
  struct subvolume *subvolp;
  struct wall *wp,**wpp,*wall_hit,*coincident_wall,*prev_coincident_wall;
  struct wall_list *wlp;
  struct vector3 curr_point;
  struct vector3 *vert;
  struct wall *prev_wall;
  double dx,dy,dz;
  double x,y,z,nx,ny,nz,t,t_min,t_inv,fuzz,a,b,c,d,ah,bh,ch;
  double wx,wy,wz,a1,a2,a3,b1,b2,b3,n_dot_delta,dot,ray_p;
  double det,prev_det,p1,p2;
  double x1,x2,y1,y2,z1,z2,pos_x,pos_y,pos_z;
  int region_x,region_y,region_z,n_walls,region_index,wall_index;
  int j,k,unbind;
  int wall_type,wall_hit_type,prev_wall_type,wall_hit_index;
  int it,coincident_t;
  int hit_status;
  byte lig_type;
  byte projection,n_vert,inside;

  /* ray trace the motion of a point */
  /* move until delta is used up */
/*
  fuzz=0.999999;
*/
  hit_status=NO_COLLISION;
  fuzz=1.0;
  lig_type=0;
  prev_wall=NULL;
  prev_coincident_wall=NULL;
  curr_point.x=point->x;
  curr_point.y=point->y;
  curr_point.z=point->z;
  dx=ray->x;
  dy=ray->y;
  dz=ray->z;
  it=0;
  while (dx!=0.0 || dy!=0.0 || dz!=0.0) {
    it++;
    x=curr_point.x;
    y=curr_point.y;
    z=curr_point.z;
    t_min=1000.0;
    wall_hit=NULL;
    coincident_wall=NULL;

    /* intersect motion with appropriate subset of subvolume walls */
    if (volume->n_subvol>1) {
      subvolp=&volume->subvol[subvol];
      wpp=subvolp->walls;
      x1=subvolp->x1;
      x2=subvolp->x2;
      y1=subvolp->y1;
      y2=subvolp->y2;
      z1=subvolp->z1;
      z2=subvolp->z2;
      pos_x=x+dx;
      pos_y=y+dy;
      pos_z=z+dz;
      if (pos_x>=x2) {
        region_x=2;
      }
      else if (pos_x<=x1) {
        region_x=0;
      }
      else {
        region_x=1;
      }

      if (pos_y>=y2) {
        region_y=6;
      }
      else if (pos_y<=y1) {
        region_y=0;
      }
      else {
        region_y=3;
      }

      if (pos_z>=z2) {
        region_z=18;
      }
      else if (pos_z<=z1) {
        region_z=0;
      }
      else {
        region_z=9;
      }
      region_index=region_x+region_y+region_z;

      n_walls=wall_count[region_index];
      for (j=0;j<n_walls;j++) {
        wall_index=wall_map[region_index][j];
        wp=wpp[wall_index];
        if (wp!=NULL) {
          if (wp!=prev_wall) {
            /* compute t value for collision with subvolume wall */
	    a=wp->normal.x;
	    b=wp->normal.y;
	    c=wp->normal.z;
	    d=wp->d;
	    n_dot_delta=a*dx+b*dy+c*dz;
	    if (n_dot_delta==0) {
	      t=1000;
	    }
	    else {
	      t=(d-a*x-b*y-c*z)/n_dot_delta;
	    }
      
	    if (t>=0.0 && t<=1.0 && t<=t_min) {
	      wall_hit=wp;
	      wall_hit_type=wp->wall_type[lig_type];
              wall_hit_index=wall_index;
	      t_min=t;
	      ah=a;
	      bh=b;
	      ch=c;
	      dot=n_dot_delta;
            }
          }
        }
      }
    }
    

    /* intersect motion with each wall in wall_list*/
    wlp=volume->subvol[subvol].wall_list;
    while (wlp!=NULL) {
      wp=wlp->wall;
      if (wp!=prev_wall) {
	/* compute t value */
	a=wp->normal.x;
	b=wp->normal.y;
	c=wp->normal.z;
	d=wp->d;
	n_dot_delta=a*dx+b*dy+c*dz;
	if (n_dot_delta==0) {
	  t=1000;
	}
	else {
	  t=(d-a*x-b*y-c*z)/n_dot_delta;
	}

	coincident_t = (fabs(t_min-t) < EPSILON_1);
	if (t>=0.0 && t<=1.0 && (t<=t_min || coincident_t)) {
	  /* take 2D projection of wall and hit point */
	  /*   and determine if point is within the wall boundaries */
	  projection = wp->projection;
	  vert=wp->vert[0];
	  if (projection == 0) {
	    p1=y+t*dy;
	    p2=z+t*dz;
	    a1=vert->y-p1;
	    a2=vert->z-p2;
	  }
	  else if (projection == 1) {
	    p1=x+t*dx;
	    p2=z+t*dz;
	    a1=vert->x-p1;
	    a2=vert->z-p2;
	  }
	  else {
	    p1=x+t*dx;
	    p2=y+t*dy;
	    a1=vert->x-p1;
	    a2=vert->y-p2;
          }
	  j=1;
	  prev_det=0;
	  inside=1;
	  n_vert = wp->n_vert;
	  while (inside && j<n_vert) {
	    vert=wp->vert[j];
	    if (projection == 0) {
	      b1=vert->y-p1;
	      b2=vert->z-p2;
	    }
	    else if (projection == 1) {
	      b1=vert->x-p1;
	      b2=vert->z-p2;
	    }
	    else {
	      b1=vert->x-p1;
	      b2=vert->y-p2;
            }
	    det=a1*b2-a2*b1;
	    a1=b1;
	    a2=b2;
	    if (j>1) {
	      inside=(!((det<0 && prev_det>=0) || (det>=0 && prev_det<0)));
	    }
	    prev_det=det;
	    j++;
	  }
	  if (inside) {
	    vert=wp->vert[0];
	    if (projection == 0) {
	      b1=vert->y-p1;
	      b2=vert->z-p2;
	    }
	    else if (projection == 1) {
	      b1=vert->x-p1;
	      b2=vert->z-p2;
	    }
	    else {
	      b1=vert->x-p1;
	      b2=vert->y-p2;
            }
	    det=a1*b2-a2*b1;
	    inside=(!((det<0 && prev_det>=0) || (det>=0 && prev_det<0)));
	    prev_det=det;
	  }
	
	  if (inside) {

	    /* check t and hit loc values to see if we have a valid hit */
	    /* and find minimum valid t value between 0 and 1 */

	    /* if 2 walls coincide with one another */
	    /* give hit priority to transparent walls */
	    wall_type=wp->wall_type[lig_type];
	    if (wall_type!=TRANSP && (coincident_t || prev_wall_type!=TRANSP)) {
	      coincident_wall=wp;
	    }
	    if (wall_type==TRANSP) {
	      wall_hit=wp;
	      wall_hit_type=wall_type;
	      t_min=t;
	      ah=a;
	      bh=b;
	      ch=c;
	      dot=n_dot_delta;
	    }
	    /* register hit to non-transparent wall */
	    /* iff this wall does not coincide with transparent wall */
	    else if (!coincident_t) {
	      if (prev_wall==NULL) {
	        wall_hit=wp;
	        wall_hit_type=wall_type;
	        t_min=t;
	        ah=a;
	        bh=b;
	        ch=c;
	        dot=n_dot_delta;
	      }
	      else if (t > EPSILON_1 || prev_wall_type!=TRANSP) {
	        wall_hit=wp;
	        wall_hit_type=wall_type;
	        t_min=t;
	        ah=a;
	        bh=b;
	        ch=c;
	        dot=n_dot_delta;
	      }
	    }
          }
        }
      }
      wlp=wlp->next;
    }
    if (wall_hit==NULL) {
      /* no collision with wall */
      curr_point.x=x+dx;
      curr_point.y=y+dy;
      curr_point.z=z+dz;
      dx=0.0;
      dy=0.0;
      dz=0.0;
      hit_status=NO_COLLISION;
    }
    else {
      /* collision with wall */
      prev_wall=wall_hit;
      prev_wall_type=wall_hit_type;
      if (wall_hit_type==TRANSP) {
	/* wall is transparent */
	nx=x+t_min*dx;
	ny=y+t_min*dy;
	nz=z+t_min*dz;
	t_inv=1.0-t_min;
	curr_point.x=nx;
	curr_point.y=ny;
	curr_point.z=nz;
	/* to continue motion through surface: */
/*
	dx=dx*t_inv;
	dy=dy*t_inv;
	dz=dz*t_inv;
*/
	/* to end motion at point of collision: */
        dx=0.0;
        dy=0.0;
        dz=0.0;
	if (dot<0) {
          hit_status=ENTERING;
	}
	else {
          hit_status=LEAVING;
	}
      }
      else if (wall_hit_type==RFLCT) {
	/* reflect off wall */
	nx=x+t_min*dx;
	ny=y+t_min*dy;
	nz=z+t_min*dz;
	t_inv=1.0-t_min;
	curr_point.x=nx;
	curr_point.y=ny;
	curr_point.z=nz;
	/* to reflect direction of motion: */
/*
	ray_p=-2.0*dot*t_inv;
	dx=ah*ray_p+dx*t_inv;
	dy=bh*ray_p+dy*t_inv;
	dz=ch*ray_p+dz*t_inv;
*/
	/* to end motion at point of collision: */
        dx=0.0;
        dy=0.0;
        dz=0.0;
	if (dot<0) {
          hit_status=ENTERING;
	}
	else {
          hit_status=LEAVING;
	}
      }
      else if (wall_hit_type==SUBVOL) {
	/* wall is a subvolume boundary */
	nx=x+t_min*dx;
	ny=y+t_min*dy;
	nz=z+t_min*dz;
	t_inv=1.0-t_min;
	curr_point.x=nx;
	curr_point.y=ny;
	curr_point.z=nz;
	dx=dx*t_inv;
	dy=dy*t_inv;
	dz=dz*t_inv;
        subvol=subvol+motion_map[wall_hit_index];
        hit_status=NO_COLLISION;
      }
    }
  }
  
  return(hit_status);
}


int clip_mesh(volp,php)
  struct volume *volp;
  struct polyhedron *php;
{
  FILE *fully_outside_file,*fully_inside_file,*on_edge_file;
  struct polygon_list *plp;
  struct polygon *pop;
  struct vertex_list *vlp,*vl0,*vl1,*vl2;
  struct vector3 *vert;
  struct vector3 ray;
  double umin,umax,vmin,vmax,wmin,wmax,u1,v1,w1;
  int hit_status,polygon_status;
  int fully_outside_vertex_count;
  int fully_inside_vertex_count;
  int on_edge_vertex_count;
  int fully_outside_polygon_count;
  int fully_inside_polygon_count;
  int on_edge_polygon_count;
  int i,subvol,n_x_part,n_y_part,n_z_part,n_x_subvol,n_y_subvol,n_z_subvol;

  /* classify vertices */
  fprintf(stderr,"Classifying vertices...\n");
  n_x_part=volp->n_x_partitions;
  n_y_part=volp->n_y_partitions;
  n_z_part=volp->n_z_partitions;
  umin=volp->x_partitions[1];
  umax=volp->x_partitions[n_x_part-2];
  vmin=volp->y_partitions[1];
  vmax=volp->y_partitions[n_y_part-2];
  wmin=volp->z_partitions[1];
  wmax=volp->z_partitions[n_z_part-2];
/*
  printf("umin = %g  umax = %g\n",umin,umax);
  printf("vmin = %g  vmax = %g\n",vmin,vmax);
  printf("wmin = %g  wmax = %g\n",wmin,wmax);
*/
  ray.x=1*(volp->x_partitions[n_x_part-2]-volp->x_partitions[1]);
  ray.y=0;
  ray.z=0;
  vlp=php->unique_vertex;
  while (vlp!=NULL) {
    vert=vlp->vertex;
    u1=vert->x;
    v1=vert->y;
    w1=vert->z;
    if (clip_point_3D(umin,umax,vmin,vmax,wmin,wmax,u1,v1,w1)) {
      subvol=find_subvol(volp,vert,-1);
      hit_status=ray_trace(vert,&ray,subvol);
      if (hit_status==NO_COLLISION) {
        vlp->vertex_status=FULLY_OUTSIDE;
      }
      else if (hit_status==ENTERING) {
        vlp->vertex_status=FULLY_OUTSIDE;
      }
      else if (hit_status==LEAVING) {
        vlp->vertex_status=FULLY_INSIDE;
      }
    }
    else {
      vlp->vertex_status=FULLY_OUTSIDE;
    }
    vlp=vlp->next;
  }
  fprintf(stderr,"Done\n");

  /* classify polygons and vertices */
  fprintf(stderr,"Classifying polygons and vertices...");
  plp=php->polygon_list;
  while (plp!=NULL) {
    pop=plp->polygon;
    vl0=pop->vertex_list_array[0];
    vl1=pop->vertex_list_array[1];
    vl2=pop->vertex_list_array[2];
    polygon_status=vl0->vertex_status+vl1->vertex_status+vl2->vertex_status;
    if (polygon_status==0) {
      /* mark polygon and vertices as fully_outside */
      pop->polygon_status=FULLY_OUTSIDE;
      vl0->fully_outside_member=1;
      vl1->fully_outside_member=1;
      vl2->fully_outside_member=1;
    }
    else if (polygon_status==3) {
      /* mark polygon and vertices as fully_inside */
      pop->polygon_status=FULLY_INSIDE;
      vl0->fully_inside_member=1;
      vl1->fully_inside_member=1;
      vl2->fully_inside_member=1;
    }
    else {
      /* mark polygon and vertices as on_edge */
      pop->polygon_status=ON_EDGE;
      vl0->on_edge_member=1;
      vl1->on_edge_member=1;
      vl2->on_edge_member=1;
    }
    plp=plp->next;
  }
  fprintf(stderr,"Done\n");

  /* output vertices to appropriate files according to classification */
  fprintf(stderr,"Writing vertices...");
  if ((fully_outside_file=fopen(fully_outside_filename,"w"))==NULL) {
    fprintf(stderr,"meshclip: error opening file: %s\n",
      fully_outside_filename);
    return(1);
  } 
  if ((fully_inside_file=fopen(fully_inside_filename,"w"))==NULL) {
    fprintf(stderr,"meshclip: error opening file: %s\n",
      fully_inside_filename);
    return(1);
  } 
/*
  if ((on_edge_file=fopen(on_edge_filename,"w"))==NULL) {
    fprintf(stderr,"meshclip: error opening file: %s\n",
      on_edge_filename);
    return(1);
  } 
*/
  fully_outside_vertex_count=0;
  fully_inside_vertex_count=0;
  on_edge_vertex_count=0;
  vlp=php->unique_vertex;
  while (vlp!=NULL) {
    vert=vlp->vertex;
    if (vlp->fully_outside_member) {
      fully_outside_vertex_count++;
      vlp->fully_outside_index=fully_outside_vertex_count;
      fprintf(fully_outside_file,"Vertex %d %.15g %.15g %.15g\n",
        fully_outside_vertex_count,vert->x,vert->y,vert->z);
    }
    if (vlp->fully_inside_member || vlp->on_edge_member) {
      fully_inside_vertex_count++;
      vlp->fully_inside_index=fully_inside_vertex_count;
      fprintf(fully_inside_file,"Vertex %d %.15g %.15g %.15g\n",
        fully_inside_vertex_count,vert->x,vert->y,vert->z);
    }
/*
    if (vlp->on_edge_member) {
      on_edge_vertex_count++;
      vlp->on_edge_index=on_edge_vertex_count;
      fprintf(on_edge_file,"Vertex %d %.15g %.15g %.15g\n",
        on_edge_vertex_count,vert->x,vert->y,vert->z);
    }
*/
    vlp=vlp->next;
  }
  fprintf(stderr,"Done\n");

  /* output polygons to appropriate files according to classification */
  fprintf(stderr,"Writing polygons...");
  fully_outside_polygon_count=0;
  fully_inside_polygon_count=0;
  on_edge_polygon_count=0;
  plp=php->polygon_list;
  while (plp!=NULL) {
    pop=plp->polygon;
    vl0=pop->vertex_list_array[0];
    vl1=pop->vertex_list_array[1];
    vl2=pop->vertex_list_array[2];
    if (pop->polygon_status==FULLY_OUTSIDE) {
      fully_outside_polygon_count++;
      fprintf(fully_outside_file,"Face %d %d %d %d\n",
        fully_outside_polygon_count,vl0->fully_outside_index,
        vl1->fully_outside_index,vl2->fully_outside_index);
    }
    else if (pop->polygon_status==FULLY_INSIDE
             || pop->polygon_status==ON_EDGE) {
      fully_inside_polygon_count++;
      fprintf(fully_inside_file,"Face %d %d %d %d\n",
        fully_inside_polygon_count,vl0->fully_inside_index,
        vl1->fully_inside_index,vl2->fully_inside_index);
    }
/*
    else if (pop->polygon_status==ON_EDGE) {
      on_edge_polygon_count++;
      fprintf(on_edge_file,"Face %d %d %d %d\n",
        on_edge_polygon_count,vl0->on_edge_index,
        vl1->on_edge_index,vl2->on_edge_index);
    }
*/
    plp=plp->next;
  }
  fprintf(stderr,"Done\n");

  fclose(fully_outside_file);
  fclose(fully_inside_file);
/*
  fclose(on_edge_file);
*/

  return(0);
}
