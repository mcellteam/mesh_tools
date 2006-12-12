#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "dx2mesh_2.h"

extern FILE *yyin;
char *infile;
char *curr_file;
int line_num=0;
int char_num=0;
int skip_freq;
int positions_item_count,positions_data_format,positions_data_offset;
int connections_item_count,connections_data_format,connections_data_offset;
struct vector3 translate;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{
  struct vertex_list *vertex_head,*vertex_tail,*vlp;
  struct vector3 *vecp;
  struct polygon_list *plp,*polygon_head,*polygon_tail;
  struct polygon *pop;
  float v1,v2,v3;
  int vertex_count,vertex_index,vert_1,vert_2,vert_3,polygon_count;
  int i;
  char *argstr;

        if (argc>=2) {
          argstr=argv[1];
          if (strcmp(argstr,"-h")==0) {
            fprintf(stderr,"\nRead DX input and convert to mesh output.\n");
            fprintf(stderr,"  Read from stdin if dx_file_name is absent.\n");
            fprintf(stderr,"  Output is written to stdout.\n\n");
            fprintf(stderr,"  Usage: %s [-h] [dx_file_name]\n\n",argv[0]);
            fflush(stdout);
            exit(1);
          }
          else {
            infile=argv[1];
            if ((yyin=fopen(infile,"r"))==NULL) {
              fprintf(stderr,"dx2mesh_2: error opening file: %s\n",infile);
              fflush(stdout);
              exit(1);
            }
            curr_file=infile;
          }
        }
        else {
          yyin=stdin;
          curr_file="stdin";
        }

	if (yyparse()) {
	  fprintf(stderr,"dx2mesh_2: error parsing mesh file %s\n",curr_file);
	  exit(1);
	} 
	fclose(yyin);

        if ((yyin=fopen(infile,"r"))==NULL) {
          fprintf(stderr,"dx2mesh_2: error opening file: %s\n",infile);
          fflush(stdout);
          exit(1);
        }

        vertex_head=NULL;
        vertex_tail=NULL;
        vertex_count=0;
        fseek(yyin,positions_data_offset,SEEK_SET);
        for (i=0;i<positions_item_count;i++) {
          if ((vecp=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
            fprintf(stderr,"dx2mesh_2: cannot store vertex");
            exit(1);
          }
          fread(&v1,sizeof(float),1,yyin);
          fread(&v2,sizeof(float),1,yyin);
          fread(&v3,sizeof(float),1,yyin);
          vecp->x=v1;
          vecp->y=v2;
          vecp->z=v3;
          if ((vlp=(struct vertex_list *)malloc
              (sizeof(struct vertex_list)))==NULL) {
            fprintf(stderr,"dx2mesh_2: cannot store vertex list");
            exit(1);
          }
          vlp->vertex_count=vertex_count++;
          vertex_index=vertex_count;
          vlp->vertex_index=vertex_index;
          vlp->vertex=vecp;
          vlp->normal=NULL;
          if (vertex_tail==NULL) {
            vertex_tail=vlp;
          }
          vertex_tail->next=vlp;
          vlp->next=NULL;
          vertex_tail=vlp;
          if (vertex_head==NULL) {
            vertex_head=vlp;
          }
        }

        polygon_head=NULL;
        polygon_tail=NULL;
        for (i=0;i<connections_item_count;i++) {
          fread(&vert_1,sizeof(int),1,yyin);
          fread(&vert_2,sizeof(int),1,yyin);
          fread(&vert_3,sizeof(int),1,yyin);
          if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
            fprintf(stderr,"dx2mesh_2: cannot store polygon");
            exit(1);
          }
          if ((plp=(struct polygon_list *)malloc
              (sizeof(struct polygon_list)))==NULL) {
            fprintf(stderr,"dx2mesh_2: cannot store polygon list");
            exit(1);
          }
          pop->n_verts=3;
          pop->vertex_index[0]=vert_1+1;
          pop->vertex_index[1]=vert_2+1;
          pop->vertex_index[2]=vert_3+1;
          plp->polygon=pop;
          if (polygon_tail==NULL) {
            polygon_tail=plp;
          }
          polygon_tail->next=plp;
          plp->next=NULL;
          polygon_tail=plp;
          if (polygon_head==NULL) {
            polygon_head=plp;
          }

        }

        vlp=vertex_head;
        while (vlp!=NULL) {
          printf("Vertex %d  %.15g %.15g %.15g\n",vlp->vertex_index,
            vlp->vertex->x,vlp->vertex->y,vlp->vertex->z);
          vlp=vlp->next;   
        }
        plp=polygon_head;
        polygon_count=0;
        while (plp!=NULL) {
          polygon_count++;
          printf("Face %d  %d %d %d\n",polygon_count,
            plp->polygon->vertex_index[0],plp->polygon->vertex_index[1],
            plp->polygon->vertex_index[2]);
          plp=plp->next;
        }
        fprintf(stderr,"\npolygon mesh: %d vertices & %d polygons\n",
          vertex_count,polygon_count);

	exit(0);
}
