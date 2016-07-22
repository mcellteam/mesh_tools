#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "strfunc.h"
#include "tetgennodeparse.h"
#include "tetgenfaceparse.h"
#include "tetgenelemparse.h"
#include "tetgen2FV.h"

extern FILE *tetgennodein;
extern FILE *tetgenfacein;
extern FILE *tetgenelemin;
char *infile_1;
char *infile_2;
char *dx_outfile_name;
char *renderman_outfile_name;
char *edsim_outfile_name;
char *curr_file;
int line_num=0;
byte input_mode;
struct tet_mesh *tet_mesh;

/* Begin main */
int main(argc,argv)
  int argc; 
  char *argv[];  
{
  char *argstr;
  char *tetgen_nodefile_name;
  char *tetgen_facefile_name;
  char *tetgen_elemfile_name;

	if (argc<3) {
      	  fprintf(stderr,"\nRead TetGen volume mesh file and convert to dx volume mesh format.\n ");
      	  fprintf(stderr,"  Output is written to dx_outfile_name.\n\n ");
      	  fprintf(stderr,"  Usage: %s vol_file_name dx_outfile_name\n\n",argv[0]);
	  exit(1);
	}
        else {
          argstr=argv[1];
          if (strcmp(argstr,"-h")==0) {
      	    fprintf(stderr,"\nRead TetGen volume mesh file and convert to dx volume mesh format.\n ");
      	    fprintf(stderr,"  Output is written to dx_outfile_name.\n\n ");
      	    fprintf(stderr,"  Usage: %s vol_file_name dx_outfile_name\n\n",argv[0]);
	    exit(1);
          }
        }

        if ((tet_mesh=(struct tet_mesh *)malloc(sizeof(struct tet_mesh)))==NULL)
        {
          fprintf(stderr,"Out of memory while storing tet_mesh");
          return(1);
        }
        tet_mesh->node_head=NULL;
        tet_mesh->nodes=NULL;
        tet_mesh->tet_head=NULL;
        tet_mesh->edge_hashtab=NULL;
        tet_mesh->voronoi_nodes=NULL;
        tet_mesh->ntets=0;
        tet_mesh->nnodes=0;
        tet_mesh->nvoronoinodes=0;
        tet_mesh->world_boundary_nodes=0;
        tet_mesh->interface_nodes=0;
        tet_mesh->nedges=0;
        tet_mesh->nboundarynodes=0;
        tet_mesh->nboundaryedges=0;
        tet_mesh->nboundaryfaces=0;
        tet_mesh->voronoi_boundary_nedges_tot=0;
        tet_mesh->voronoi_facet_nedges_tot=0;
        tet_mesh->voronoi_cap_nedges_tot=0;
        tet_mesh->voronoi_nfacets_tot=0;
        tet_mesh->voronoi_facet_ntriangles_tot=0;
        tet_mesh->voronoi_cap_ntriangles_tot=0;
        tet_mesh->nmaterials=0;
        tet_mesh->hashsize=0;
        tet_mesh->hashmask=0;

        input_mode=TETGEN_MODE;

	infile_1=argv[1];
	dx_outfile_name=argv[2];
	renderman_outfile_name="tet_mesh.rib";
	edsim_outfile_name="foo.edsim";

        if (input_mode==TETGEN_MODE)
        {
          tetgen_nodefile_name=my_strcat(infile_1,".node");
          tetgen_facefile_name=my_strcat(infile_1,".face");
          tetgen_elemfile_name=my_strcat(infile_1,".ele");

	  if ((tetgennodein=fopen(tetgen_nodefile_name,"r"))==NULL) {
	    fprintf(stderr,"tetgen2FV: error opening tetgen node file: %s\n",
              tetgen_nodefile_name);
	    exit(1);
	  } 
          curr_file=tetgen_nodefile_name;
	  fflush(stdout);
	  if (tetgennodeparse()) {
	    fprintf(stderr,"tetgen2FV: error parsing tetgen node file %s\n",
              curr_file);
	    exit(1);
	  } 
	  fclose(tetgennodein);

	  if ((tetgenfacein=fopen(tetgen_facefile_name,"r"))==NULL) {
	    fprintf(stderr,"tetgen2FV: error opening tetgen face file: %s\n",
              tetgen_facefile_name);
	    exit(1);
	  } 
          curr_file=tetgen_facefile_name;
	  fflush(stdout);
	  if (tetgenfaceparse()) {
	    fprintf(stderr,"tetgen2FV: error parsing tetgen face file %s\n",
              curr_file);
	    exit(1);
	  } 
	  fclose(tetgenfacein);

	  if ((tetgenelemin=fopen(tetgen_elemfile_name,"r"))==NULL) {
	    fprintf(stderr,"tetgen2FV: error opening tetgen element file: %s\n",
              tetgen_elemfile_name);
	    exit(1);
	  } 
          curr_file=tetgen_elemfile_name;
	  fflush(stdout);
	  if (tetgenelemparse()) {
	    fprintf(stderr,"tetgen2FV: error parsing tetgen element file %s\n",
              curr_file);
	    exit(1);
	  } 
	  fclose(tetgenelemin);

        }

	exit(0);
}
