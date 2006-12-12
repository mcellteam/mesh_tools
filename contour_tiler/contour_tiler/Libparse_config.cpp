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

#include <OB/CORBA.h> //KLC ???
#include "..\server\libHead.h" //KLC
#include "..\server\Prism_impl.h" //KLC

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

extern int get_no_of_poly_triangles();
extern char poly_file_name[300];

extern int STATIC_CHQ;

void my_clear_libparse_config()
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

int my_do_database(Interpo_struct *pstru, Name_struct *nstru, TileConfig myConf)
{
    memset(pstru,0,sizeof(*pstru));

	malloc_copy(&(nstru->prefix), myConf.prefix);
	malloc_copy(&(nstru->suffix), myConf.suffix);
	malloc_copy(&(nstru->output_file_name), myConf.outputFileName);
	
    if(!strcmp(myConf.type,"CTRAS"))
		pstru->type = CTRAS;
	else if(!strcmp(myConf.type,"CTPGM"))
		pstru->type = CTPGM;
	else if(!strcmp(myConf.type,"CTVOL"))
		pstru->type = CTVOL;
	else if(!strcmp(myConf.type,"CTSLC"))
		pstru->type = CTSLC;
	else
		pstru->type =0;

	if (myConf.epsilon != 0) pstru->epsilon = myConf.epsilon;

	if (myConf.detz != 0) pstru->detz= myConf.detz;

	if (myConf.scalez != 0) Scale_z= myConf.scalez;

	if (myConf.firstz != 0) FIRST_Z = myConf.firstz;

	if (myConf.mergeDistanceSquare != 0) MERGE_DISTANCE_SQUARE = myConf.mergeDistanceSquare;

	if (myConf.thres != 0)  pstru->thres = myConf.thres;

	pstru->beg = myConf.sliceRange[0];
	pstru->end = myConf.sliceRange[1];
	pstru->x1 = myConf.xRange[0];
	pstru->x2 = myConf.xRange[1];
	pstru->y1 = myConf.yRange[0];
	pstru->y2 = myConf.yRange[1];
    pstru->xcube_vox = myConf.cubeSize[0];
	pstru->ycube_vox = myConf.cubeSize[1];
	pstru->zcube_vox = myConf.cubeSize[2];
	MATH_FUNC = myConf.mathFunc;
	DO_APPROX_ONLY = myConf.doApproxOnly;
	CONTOURS_ONLY = myConf.contoursOnly;
	MESH_FILE = myConf.fileForMeshGen;
	NON_VALID = myConf.nonValid;
	SAME_LEVEL = myConf.sameLevel;
	DIFF_LEVEL_FILES = myConf.diffLevelFiles;

	if (myConf.filterType != 0) FILTER_TYPE = myConf.filterType;
	
	if (myConf.fromCountours == 1) FROM_CONTOURS =1; //?? its default 1. so it WILL b 1 :-(.

printf("prefix = %s\n", nstru->prefix);
printf("suffix = %s\n", nstru->suffix);
printf("opfname = %s\n", nstru->output_file_name);
printf("type = %d\n", pstru->type);
printf("eplison = %f\n", pstru->epsilon);
printf("detz = %f\n", pstru->detz);
printf("scalez= %f\n", Scale_z);
printf("type = %d\n", pstru->type);
printf("firstz = %f\n", FIRST_Z);
printf("mergedist = %f\n", MERGE_DISTANCE_SQUARE);
printf("thres = %f\n", pstru->thres);
printf("begin = %d\n", pstru->beg);
printf("end = %d\n", pstru->end);
printf("type = %d\n", pstru->type);
printf("xrange = %d and = %d\n", pstru->x1, pstru->x2);
printf("yrange = %d and = %d\n", pstru->y1, pstru->y2);
printf("x = %d and y= %d and z=%d\n", pstru->xcube_vox, pstru->ycube_vox, pstru->zcube_vox);
printf("mathfunc = %d\n", MATH_FUNC);
printf("contoursOnly = %d\n", CONTOURS_ONLY);
printf("doapprox = %d\n", DO_APPROX_ONLY);
printf("mathfunc = %d\n", MATH_FUNC);
printf("meshfile= %d\n", MESH_FILE);
printf("nonValid = %d\n", NON_VALID);
printf("samelevel = %d\n", SAME_LEVEL);
printf("diff = %d\n", DIFF_LEVEL_FILES);
printf("filterT = %d\n", FILTER_TYPE);
printf("fromC = %d\n", FROM_CONTOURS);

	return(1);
}

int my_read_database(Interpo_struct *pstru, Name_struct *nstru, TileConfig myConf)
{
    int i;
   
    i=my_do_database(pstru, nstru, myConf);
    
    return i;
}

HexOutput* my_write_output()
{
	//open the poly file n read from it....
	int i, j, triCounter=0, faceCounter=0, vertCounter=0;
	HexOutput* rest;
	int *mfaces;
	float* mverts;
	FILE* fp;
	float vx, vy, vz;
	
	triCounter = get_no_of_poly_triangles();

	rest  = (HexOutput*) (malloc (sizeof(HexOutput)));
	mfaces = (int*)(malloc((sizeof(int)) * (triCounter *3))); //thats the MAX that can b used
	mverts = (float*)(malloc((sizeof(float)) * (triCounter *9)));

	fp = fopen(poly_file_name, "r");
	if (fp == 0) 
	{
		printf("err in opening the poly file to read data..%s\n", poly_file_name);
		exit(0);
	}

	while (1) 
	{
		j = fscanf(fp, "%d", &i);
		if (j == EOF)	break;
		
		j = fscanf(fp, "%f %f %f", &vx, &vy, &vz);
		if (j == EOF)	break;
		mverts[vertCounter++] = vx;	mverts[vertCounter++] = vy;	mverts[vertCounter++] = vz;
		mfaces[faceCounter++] = (vertCounter-1)/3;

		j = fscanf(fp, "%f %f %f", &vx, &vy, &vz);
		if (j == EOF)	break;
		mverts[vertCounter++] = vx;	mverts[vertCounter++] = vy;	mverts[vertCounter++] = vz;
		mfaces[faceCounter++] = (vertCounter-1)/3;

		j = fscanf(fp, "%f %f %f", &vx, &vy, &vz);
		if (j == EOF)	break;
		mverts[vertCounter++] = vx;	mverts[vertCounter++] = vy;	mverts[vertCounter++] = vz;
		mfaces[faceCounter++] = (vertCounter-1)/3;
	}

	fclose(fp);

	rest->faces = mfaces;
	rest->verts = mverts;
	rest->numFaces = faceCounter;
	rest->numVerts = vertCounter;
	rest->hexes = NULL;
	rest->numHexes =0;
	printf("tris written =%d and verts are= %d\n", faceCounter, vertCounter);
	fflush(stdout);

	return rest;
}