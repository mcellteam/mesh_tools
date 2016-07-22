%{
#include <stdio.h> 
#include <string.h> 
#include <sys/types.h>
#include <math.h>
#include "strfunc.h"
#include "tetgen2FV.h"
#include "voronoi.h"
#include "output_dx.h"
#include "output_edsim.h"
#include "tetgenelemparse.h"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *tetgenelemin;
extern char *infile_1;
extern char *infile_2;
extern char *dx_outfile_name;
extern char *renderman_outfile_name;
extern char *edsim_outfile_name;
extern char *curr_file;
extern int line_num;
extern struct tet_mesh *tet_mesh;

FILE *outfile;
int ival;
double rval;
char *cval;
char *var;
char *a_str,*rem_str;
char c,fmt_str[128],time_str[128];
char err_msg[128];
struct object *op,*curr_op;
struct meta_object *mop,*curr_mop,*wmop;
struct object_list *olp;
struct polyhedron *php;
struct polygon *pop;
struct polygon_list *plp,*polygon_head,*polygon_tail;
struct voronoi_volume_element *vvep;
struct voronoi_facet_list *vflp;
struct voronoi_facet *vfp;
struct node *np;
struct node **node_array;
struct node_list *nlp,*node_head,*node_tail;
struct tet *tetp;
struct tet_list *tlp,*tet_head,*tet_tail;
struct material_list *mlp,*found_mat;
struct node_material_info *nmip;
struct nurbs *nrbp;
struct double_list *dlp_head,*dlp;
struct ctlpt_list *clp_head,*clp;
struct ctlpt *cpp;
u_int node_count;
u_int child_nodes;
u_int tet_count;
u_int edge_count;
int node_index;
int polygon_count;
int output_count;
int multimat_tet_count;
int multimat_max;
int mat_found;
int data_offset;
double x,y,z;
float nx,ny,nz;
float fi;
u_int n1,n2,n3;
int node_1,node_2,node_3,node_4;
u_int i,j,k;


%}

%union {
int tok;
char *str;
double dbl;
struct object *obj;
} 

%{
  #include "tetgenelemlex.c"
%}


%token <tok> DIMENSION EDGE_SEGMENTS MESH3D POINTS SURFACE_ELEMENTS
%token <tok> VOLUME_ELEMENTS NEWLINE REAL INTEGER 
%type <dbl> int_arg real_arg num_arg 

%right '='
%left '+' '-'
%left '*' '/'
%left UNARYMINUS


%%


tetgen_elem_format: 
{
  fprintf(stderr,"tetgen2FV: parsing tetgen element file %s...\n",curr_file);
  tet_count=0;
  multimat_tet_count=0;
  multimat_max=0;
  tlp=NULL;
  tet_head=NULL;
  tet_tail=NULL;
}
	tetgen_elem_header
        element_list
{ 

  tet_mesh->ntets=tet_count;
  tet_mesh->tet_head=tet_head;

  fprintf(stderr,"          ... %d tet elements found\n",tet_count);

  /* split all interface nodes into material-specific children */
  child_nodes=0;
  for (i=0; i<tet_mesh->nnodes; i++)
  {
    np=tet_mesh->nodes[i];
    if (np->node_type==INTERFACE)
    {
      tet_mesh->interface_nodes++;
      if ((nmip=(struct node_material_info *)malloc
         (np->material_count*sizeof(struct node_material_info)))==NULL)
      {
        tetgenelemerror("Out of memory while storing node material info");
        return(1);
      }
      np->node_material=nmip;
      mlp=np->material_head;
      for (j=0;j<np->material_count;j++)
      {
        nmip[j].matnum=mlp->matnum;
        nmip[j].child_index=tet_mesh->nnodes+child_nodes;
        nmip[j].volume=0;
        nmip[j].nedges=0;
        nmip[j].ntets=0;
        nmip[j].nfacets=0;
        nmip[j].edge_head=NULL;
        child_nodes++;
        mlp=mlp->next;
      }
    }
  }

  fprintf(stderr,"          ... %d interface nodes found\n",
    tet_mesh->interface_nodes);

  if (construct_voronoi_mesh(tet_mesh))
  {
    tetgenelemerror("Fatal error while constructing voronoi mesh");
    return(1);
  }

/*
  if (output_dx_full(tet_mesh,dx_outfile_name))
*/
  if (output_dx_tet_mesh(tet_mesh,dx_outfile_name))
  {
    tetgenelemerror("Fatal error while writing DX output file");
    return(1);
  }

  if (output_renderman_full(tet_mesh,renderman_outfile_name))
  {
    tetgenelemerror("Fatal error while writing RenderMan output file");
    return(1);
  }

  if (output_edsim(tet_mesh,edsim_outfile_name))
  {
    tetgenelemerror("Fatal error while writing EDSIM output file");
    return(1);
  }
 
  edge_count=tet_mesh->nedges;

  fprintf(stderr,"\nvolume mesh:  %d nodes,  %d edges,  %d tets\n",node_count,edge_count,tet_count);

};


newline_list: NEWLINE
	| newline_list NEWLINE
;


tetgen_elem_header: int_arg int_arg int_arg newline_list
;


element_list: element
	| element_list element
;


element:	int_arg int_arg int_arg int_arg int_arg newline_list
{
  node_1=$<dbl>2;
  node_2=$<dbl>3;
  node_3=$<dbl>4;
  node_4=$<dbl>5;

  if ((tetp=(struct tet *)malloc(sizeof(struct tet)))==NULL)
  {
    tetgenelemerror("Out of memory while storing tet");
    return(1);
  }
  if ((tlp=(struct tet_list *)malloc(sizeof(struct tet_list)))==NULL)
  {
    tetgenelemerror("Out of memory while storing tet list");
    return(1);
  }

  tetp->tet_index=tet_count;
  tet_count++;
  tetp->node_index[0]=node_1-1;
  tetp->node_index[1]=node_2-1;
  tetp->node_index[2]=node_3-1;
  tetp->node_index[3]=node_4-1;
  tetp->cent.x=0;
  tetp->cent.y=0;
  tetp->cent.z=0;
  tetp->matnum=0;
  tetp->srfnum=-1;
  tetp->material_count=0;
  tetp->r=0;
  tetp->delaunay=1;
  tetp->tet_type=INTERIOR;
/*
  tetp->nboundaryfaces=0;
  tetp->boundary_head=NULL;
*/

  tlp->tet=tetp;
  if (tet_tail==NULL)
  {
    tet_tail=tlp;
  }
  tet_tail->next=tlp;
  tlp->next=NULL;
  tet_tail=tlp;
  if (tet_head==NULL)
  {
    tet_head=tlp;
  }
}
	| int_arg int_arg int_arg int_arg int_arg int_arg newline_list
{
  node_1=$<dbl>2;
  node_2=$<dbl>3;
  node_3=$<dbl>4;
  node_4=$<dbl>5;

  if ((tetp=(struct tet *)malloc(sizeof(struct tet)))==NULL)
  {
    tetgenelemerror("Out of memory while storing tet");
    return(1);
  }
  if ((tlp=(struct tet_list *)malloc(sizeof(struct tet_list)))==NULL)
  {
    tetgenelemerror("Out of memory while storing tet list");
    return(1);
  }

  tetp->tet_index=tet_count;
  tet_count++;
  tetp->node_index[0]=node_1-1;
  tetp->node_index[1]=node_2-1;
  tetp->node_index[2]=node_3-1;
  tetp->node_index[3]=node_4-1;
  tetp->cent.x=0;
  tetp->cent.y=0;
  tetp->cent.z=0;
  tetp->matnum=$<dbl>6;
  tetp->srfnum=-1;
  tetp->material_count=1;
  tetp->r=0;
  tetp->delaunay=1;
  tetp->tet_type=INTERIOR;
/*
  tetp->nboundaryfaces=0;
  tetp->boundary_head=NULL;
*/

  tlp->tet=tetp;
  if (tet_tail==NULL)
  {
    tet_tail=tlp;
  }
  tet_tail->next=tlp;
  tlp->next=NULL;
  tet_tail=tlp;
  if (tet_head==NULL)
  {
    tet_head=tlp;
  }

  /* add material property to each node in tet */
  for (i=0; i<4; i++)
  {
    np=tet_mesh->nodes[tetp->node_index[i]];

    /* look for matnum in material list for this node */
    found_mat=NULL;
    for (mlp=np->material_head; mlp!=NULL; mlp=mlp->next)
    {
      if (mlp->matnum==tetp->matnum)
      {
        found_mat=mlp;
      }
    }

    /* add material if not found */
    if (found_mat==NULL)
    {
      if ((mlp=(struct material_list *)malloc
         (sizeof(struct material_list)))==NULL)
      {
        tetgenelemerror("Out of memory while storing node material");
        return(1);
      }
      mlp->matnum=tetp->matnum;
      mlp->next=np->material_head;
      np->material_head=mlp;
      np->material_count++; 
      if (np->material_count>1)
      {
        np->node_type=INTERFACE;
      }
    }
    
  }

};


int_arg: INTEGER {$$=(double)ival;}
;


real_arg: REAL {$$=rval;}
;


num_arg: INTEGER {$$=(double)ival;}
	| REAL {$$=rval;}
;


%%


#undef tetgenelemwrap

int tetgenelemwrap()
{
        return(1);
} 


int tetgenelemerror(char *s)
{
	fprintf(stderr,"tetgen2FV: error on line: %d of file: %s  %s\n",
	        line_num,curr_file,s);
	fflush(stderr);
	return(1);
}

