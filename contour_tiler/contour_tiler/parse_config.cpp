/*****************************************************************************/
/*****************************************************************************/
/**                                                                         **/
/** This SHASTRA software is not in the Public Domain. It is distributed on **/
/** a person to person basis, solely for educational use and permission is  **/
/** NOT granted for its transfer to anyone or for its use in any commercial **/
/** product.  There is NO warranty on the available software and neither    **/
/** Purdue University nor the Applied Algebra and Geometry group directed   **/
/** by C.  Bajaj accept responsibility for the consequences of its use.     **/
/**                                                                         **/
/*****************************************************************************/
/*****************************************************************************/

/* ---------------------------------------------------------------------
parse_config.c

 Parse the input file to read input variables
 
------------------------------------------------------------------------ */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <memory.h>

#include "common.h"
#include "myutil.h"
#include "parse_config.h"

#define BUFSIZE 10000
#define BUFSIZE10 1000

static char **Ptr;
static int Index;
extern int MATH_FUNC;
extern int FROM_CONTOURS;
extern int CONTOURS_ONLY;
extern int FILTER_TYPE;
extern int DO_APPROX_ONLY;
extern float FIRST_Z;
extern double MERGE_DISTANCE_SQUARE;
extern float Scale_z;
extern int MESH_FILE;
extern int SAME_LEVEL; /* if 1: only tile between the same level */
extern int DIFF_LEVEL_FILES; /* if 1: put diff. level into diff. files */
extern int NON_VALID; /* for visualization only, the triangulation may not be valid */

extern int STATIC_CHQ;

void my_clear_parse_config()
{
	Ptr = NULL;
	Index =0;
}

char * read_parse_file(char *filename)
{
    FILE *fi;
    int size;
    char *str;
	
    if((fi=fopen(filename,"r"))==NULL)
	{   
		printf("Could not open default configuration %s\n",filename); 
		exit(1); 
    }
    
	str=(char *)mymalloc(BUFSIZE*sizeof(*str));
    
	if((size=fread(str,1,BUFSIZE,fi))<=5) 
	{
        fprintf(stderr,"Error in reading %s\n",filename);
		exit(1);
    }
    
	fclose(fi);
    
	if(size==BUFSIZE)
	{
		str[size-1]=0;
		fprintf(stderr,"%s too big, only %d bytes are read\n",filename,BUFSIZE);
    }
    
	else str[size]=0;
    return str;
}

int parse_config(char *str)
{
    int i,j,k,flag,first;
    /* ignore all # up to end of file */
    Ptr=(char **)mymalloc(BUFSIZE*sizeof(Ptr[0]));
    flag=k=j=Index=0;
    first=1;
    
	for(i=0; str[i]!=0; i++) 
	{
        if(str[i]!='\n')
		{
			/* not end of line */
            if(flag==0) 
			{
                if(str[i]=='#')flag=1; /* begin comment to EOL */
                else 
				{
                    if(str[i]>' ') 
					{ 
						/* printable char */
                        if(first) 
						{
							first=0; 
							Ptr[Index++]=(char *)&(str[i]);
							if(Index>=BUFSIZE10) 
							{
								fprintf(stderr,"wrong config file\n");
								return(0);
							}
						}
					}
                    else 
					{
						first = 1; str[i]=0; 
					}
                }
            }
        }
        else 
		{
			/* EOL */
			first=1;
            k=flag=0;
			str[i]=0;
        }
    }
    return 1;
}

int read_database(char *filename, Interpo_struct *pstru, Name_struct *nstru)
{
    int i;
    char *str;
    str=read_parse_file(filename);
    if(str==NULL)return 0;
    parse_config(str);
    i=do_database(pstru, nstru);
    free(str);
    return i;
}

int tcp_read_database(char *str, Interpo_struct *pstru, Name_struct *nstru)
{
    int i;
    parse_config(str);
    i=do_database(pstru, nstru);
    return i;
}

int do_database(Interpo_struct *pstru, Name_struct *nstru)
{
    int j;
    memset(pstru,0,sizeof(*pstru));
    j=0;

    while(j<Index)
    {
        if(!strcmp(Ptr[j],"PREFIX:"))
        {
			if(++j>=Index)break;
			malloc_copy(&(nstru->prefix), Ptr[j]);
		}
        else if(!strcmp(Ptr[j],"SUFFIX:"))
        {
			if(++j>=Index)break;
			malloc_copy(&(nstru->suffix), Ptr[j]);
		}
        else if(!strcmp(Ptr[j],"OUTPUT_FILE_NAME:"))
        {
			if(++j>=Index)break;
			malloc_copy(&(nstru->output_file_name), Ptr[j]);
		}
        else if(!strcmp(Ptr[j],"TYPE:"))
        {
			if(++j>=Index)break;
			if(!strcmp(Ptr[j],"CTCTF"))pstru->type = CTCTF;
			else if(!strcmp(Ptr[j],"CTRAS"))pstru->type = CTRAS;
			else if(!strcmp(Ptr[j],"CTPGM"))pstru->type = CTPGM;
			else if(!strcmp(Ptr[j],"CTVOL"))pstru->type = CTVOL;
			else if(!strcmp(Ptr[j],"CTSLC"))pstru->type = CTSLC;
			else {printf("unknown type [%s], use default, CTPGM\n",Ptr[j]);}
        }
		/*
        else if(!strcmp(Ptr[j],"METHOD:"))
        {
		if(++j>=Index)break;
		if(!strcmp(Ptr[j],"MARCHING_CUBES"))pstru->method = MARCHING_CUBES;
		else if(!strcmp(Ptr[j],"TRILINEAR_MARCHING_CUBES"))
		pstru->method = TRILINEAR_MARCHING_CUBES;
		else if(!strcmp(Ptr[j],"MARCHING_TETRAHEDRA"))
		pstru->method = MARCHING_TETRAHEDRA;
		else if(!strcmp(Ptr[j],"SPIDER_WEB"))pstru->method = SPIDER_WEB;
		else {printf("unknown method [%s], use default, MARCHING_CUBES\n",
		Ptr[j]);}
        }
		*/
        else if(!strcmp(Ptr[j],"EPSILON:"))
		{
			if(++j>=Index)break;
			sscanf(Ptr[j],"%f",&(pstru->epsilon));
			/*	Do not use atof(), it does not work on some machine */
		}
        else if(!strcmp(Ptr[j],"DETZ:"))
		{
			if(++j>=Index)break;
			sscanf(Ptr[j],"%f",&pstru->detz);
		}
        else if(!strcmp(Ptr[j],"SCALEZ:"))
		{
			if(++j>=Index)break;
			sscanf(Ptr[j],"%f",&Scale_z);
		}
        else if(!strcmp(Ptr[j],"FIRST_Z:"))
		{
			if(++j>=Index)break;
			sscanf(Ptr[j],"%f",&FIRST_Z);
		}
        else if(!strcmp(Ptr[j],"MERGE_DISTANCE_SQUARE:"))
		{
			if(++j>=Index)break;
			sscanf(Ptr[j],"%lf",&MERGE_DISTANCE_SQUARE);
		}
        else if(!strcmp(Ptr[j],"THRES:"))
		{
			if(++j>=Index)break;
			sscanf(Ptr[j],"%f",&(pstru->thres));
		}
        else if(!strcmp(Ptr[j],"SLICE_RANGE:"))
		{
			if(++j>=Index)break;
			pstru->beg=atoi(Ptr[j]);
			if(++j>=Index)break;
			pstru->end=atoi(Ptr[j]);
		}
        else if(!strcmp(Ptr[j],"X_RANGE:"))
		{
			if(++j>=Index)break;
			pstru->x1=atoi(Ptr[j]);
			if(++j>=Index)break;
			pstru->x2=atoi(Ptr[j]);
		}
        else if(!strcmp(Ptr[j],"Y_RANGE:"))
		{
			if(++j>=Index)break;
			pstru->y1=atoi(Ptr[j]);
			if(++j>=Index)break;
			pstru->y2=atoi(Ptr[j]);
		}
        else if(!strcmp(Ptr[j],"CUBE_SIZE:"))
		{
			if(++j>=Index)break;
			pstru->xcube_vox=atoi(Ptr[j]);
			if(pstru->xcube_vox<=0 || pstru->xcube_vox>10)
				fprintf(stderr,"CUBE_SIZE, X =%d should be between [1,32]\n",
				pstru->xcube_vox);
			if(++j>=Index)break;
			pstru->ycube_vox=atoi(Ptr[j]);
			if(pstru->ycube_vox<=0 || pstru->ycube_vox>10)
				fprintf(stderr,"CUBE_SIZE, Y =%d should be between [1,32]\n",
				pstru->ycube_vox);
			if(++j>=Index)break;
			pstru->zcube_vox=atoi(Ptr[j]);
			if(pstru->zcube_vox<=0 || pstru->zcube_vox>10)
				fprintf(stderr,"CUBE_SIZE, Z =%d should be between [1,10]\n",
				pstru->zcube_vox);
		}
        else if(!strcmp(Ptr[j],"MATH_FUNC:"))
		{
			if(++j>=Index)break;
			MATH_FUNC=atoi(Ptr[j]);
		}
        else if(!strcmp(Ptr[j],"FROM_CONTOURS:"))
		{
			FROM_CONTOURS=1;
		}
        else if(!strcmp(Ptr[j],"DO_APPROX_ONLY:"))
		{
			DO_APPROX_ONLY=1;
		}
        else if(!strcmp(Ptr[j],"CONTOURS_ONLY:"))
		{
			CONTOURS_ONLY=1;
		}
        else if(!strcmp(Ptr[j],"FILE_FOR_MESH_GEN:"))
		{
			MESH_FILE=1;
		}
        else if(!strcmp(Ptr[j],"NON_VALID:"))
		{
			NON_VALID=1;
		}
        else if(!strcmp(Ptr[j],"SAME_LEVEL:"))SAME_LEVEL=1;
        else if(!strcmp(Ptr[j],"DIFF_LEVEL_FILES:"))DIFF_LEVEL_FILES=1;
        else if(!strcmp(Ptr[j],"FILTER_TYPE:"))
		{
			if(++j>=Index)break;
			FILTER_TYPE=atoi(Ptr[j]);
		}
		else 
		{
			printf("unknown command %s in configuration file\n",Ptr[j]);
			exit(1);
		}
		++j;
    }
    free(Ptr);
    return(1);
}



