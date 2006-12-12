#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "synu2mesh.h"

FILE *infile;
char *in_filename;
char *curr_filename;
int line_num=0;
int skip_freq;
struct vector3 translate;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{
  struct vector3 *vertices,*normals;
  struct polygon *polygons;
  int *edges;
  int nvertices,nedges,npolygons,flags,first_edge,last_edge,n_verts,edge_index;
  int i,j;
  char *argstr;
  char strbuf[256];

        if (argc>=2) {
          argstr=argv[1];
          if (strcmp(argstr,"-h")==0) {
            fprintf(stderr,"\nRead SYNU polygon mesh input and convert to mesh output.\n");
            fprintf(stderr,"  Read from stdin if synu_file_name is absent.\n");
            fprintf(stderr,"  Output is written to stdout.\n\n");
            fprintf(stderr,"  Usage: %s [-h] [synu_file_name]\n\n",argv[0]);
            fflush(stdout);
            exit(1);
          }
          else {
            in_filename=argv[1];
            if ((infile=fopen(in_filename,"r"))==NULL) {
              fprintf(stderr,"synu2mesh: error opening file: %s\n",in_filename);
              fflush(stdout);
              exit(1);
            }
            curr_filename=in_filename;
          }
        }
        else {
          infile=stdin;
          curr_filename="stdin";
        }

	fgets(strbuf,sizeof(strbuf),infile);
        fscanf(infile,"%d",&nvertices);
        fscanf(infile,"%d",&nedges);
        fscanf(infile,"%d",&npolygons);
        fscanf(infile,"%d",&flags);
        if (flags == 800 || flags == 0) {
          if ((vertices=(struct vector3 *)malloc(nvertices*sizeof(struct vector3)))==NULL) {
	    fprintf(stderr,"synu2mesh: error cannot store vertices from %s\n",curr_filename);
            exit(1);
          }
          if (flags == 800) {
            if ((normals=(struct vector3 *)malloc
              (nvertices*sizeof(struct vector3)))==NULL) {
	      fprintf(stderr,"synu2mesh: error cannot store normals from %s\n",
                curr_filename);
              exit(1);
            }
          }
          if ((edges=(int *)malloc(nedges*sizeof(int)))==NULL) {
	    fprintf(stderr,"synu2mesh: error cannot store edges from %s\n",curr_filename);
            exit(1);
          }
          if ((polygons=(struct polygon *)malloc(npolygons*sizeof(struct polygon)))==NULL) {
	    fprintf(stderr,"synu2mesh: error cannot store polygons from %s\n",curr_filename);
            exit(1);
          }
          for (i=0;i<nvertices;i++) {
            fscanf(infile,"%lf %lf %lf",&vertices[i].x,&vertices[i].y,&vertices[i].z);
            if (flags == 800) {
              fscanf(infile,"%lf %lf %lf",
                &normals[i].x,&normals[i].y,&normals[i].z);
            }
            fprintf(stdout,"Vertex %d %.15g %.15g %.15g\n",
              i+1,vertices[i].x,vertices[i].y,vertices[i].z);
          }
          for (i=0;i<nedges;i++) {
            fscanf(infile,"%d",&edges[i]);
          }
          first_edge=0;
          for (i=0;i<npolygons;i++) {
            fscanf(infile,"%d",&last_edge);
            n_verts=last_edge-first_edge+1;
            polygons[i].n_verts=n_verts;
            if ((polygons[i].vertex_index=(int *)malloc(n_verts*sizeof(int)))==NULL) {
	      fprintf(stderr,"synu2mesh: error cannot store polygons from %s\n",curr_filename);
              exit(1);
            }
            edge_index=first_edge;
            fprintf(stdout,"Face %d",i+1);
            for (j=0;j<n_verts;j++) {
              polygons[i].vertex_index[j]=edges[edge_index++];
              fprintf(stdout," %d",polygons[i].vertex_index[j]+1);
            }
            fprintf(stdout,"\n");
            first_edge=edge_index; 
          }
        }
        fprintf(stderr,"\npolygon mesh:  %d vertices & %d polygons\n",
          nvertices,npolygons);

	fclose(infile);

	exit(0);
}
