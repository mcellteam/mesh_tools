#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mesh2ribwireframe.h"

extern FILE *yyin;
struct sym_table **main_sym_table;
char *infile;
char *curr_file;
int line_num;
int skip_freq;
double radius;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{
  char *argstr;

        
	if (argc>=2) {
	  argstr=argv[1];
          if (strcmp(argstr,"-h")==0) {
      	    fprintf(stderr,"\nRead mesh input and convert to rib wireframe output.\n");
      	    fprintf(stderr,"  Read from stdin if mesh_file_name is absent.\n");
      	    fprintf(stderr,"  Output is written to stdout.\n\n");
      	    fprintf(stderr,"  Usage: %s [-h] [mesh_file_name]\n\n",argv[0]);
	    fflush(stdout);
	    exit(1);
          }
          else {
	    infile=argv[1];
	    if ((yyin=fopen(infile,"r"))==NULL) {
	      fprintf(stderr,"mesh2ribwireframe: error opening file: %s\n",infile);
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

        main_sym_table=init_symtab(HASHSIZE);
        radius=0.04;

	if (yyparse()) {
	  fprintf(stderr,"mesh2ribwireframe: error parsing file: %s\n",curr_file);
	  exit(1);
	} 
	fclose(yyin);

	exit(0);
}

struct sym_table **init_symtab(size)
int size;
{
  struct sym_table **symtab;
  int i;
  symtab=(struct sym_table **)malloc(size*sizeof(struct sym_table *));
  for (i=0;i<size;symtab[i++]=NULL);
  return(symtab);
}

unsigned hash(sym)
  char *sym;
{
  unsigned hashval;
  
  for (hashval=0; *sym!='\0';sym++) {
    hashval=*sym+31*hashval;
  }
  return(hashval%HASHSIZE);
}
  
struct sym_table *retrieve_sym(sym,sym_type,hashtab)
  char *sym;
  unsigned int sym_type;
  struct sym_table **hashtab;
{
  struct sym_table *sp;
  
  for (sp=hashtab[hash(sym)]; sp!=NULL; sp=sp->next) {
    if (strcmp(sym,sp->name)==0 && sp->sym_type==sym_type) { 
      return(sp);
    }
  }
  return(NULL);
}

struct sym_table *store_sym(sym,sym_type,hashtab)
  char *sym;
  unsigned int sym_type;
  struct sym_table **hashtab;
{
  struct sym_table *sp;
  unsigned hashval;
  void *vp;
  struct edge *edgep;
  int i;

  vp=NULL;
  /* try to find sym in table */
  if ((sp=retrieve_sym(sym,sym_type,hashtab))==NULL) {  /* sym not found */
    if ((sp=(struct sym_table *)malloc(sizeof(struct sym_table)))==NULL) {
      return(NULL);
    }
    sp->name=sym;
    sp->sym_type=sym_type;
    hashval=hash(sym);
    sp->next=hashtab[hashval];
    hashtab[hashval]=sp;
    switch (sym_type) {
    case EDGE:
      if ((vp=(void *)malloc(sizeof(struct edge)))==NULL) {
        return(NULL);
      }
      edgep=(struct edge *)vp;
      edgep->sym=sp;
      edgep->vert_1=-1;
      edgep->vert_2=-1;
      break;
    }
    sp->value=vp;
  }
  else {                        /*sym found*/
    free((void *)sym);
  }
  return(sp);
}


