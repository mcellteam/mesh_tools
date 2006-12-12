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

/*-----------------------------------------------------------------------
interp_read.c

 major functions:
 1. write polygons to files. 
 2. read slices and call functions to do work.
 3. print out the timing data. 
-----------------------------------------------------------------------*/


#include <stdio.h>
// #include <io.h>
#include <stdlib.h>
//#include <unistd.h> //KLC
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
/* #include "ct.h" */

#include "common.h" 
#include "contour.h" 
#include "gen_data.h" 
#include "read_slice.h" 
#include "image2cont2D.h" 
#include "myutil.h" 
#include "contour_read.h" 

extern int DEBUG;
extern int MATH_FUNC;
extern int FROM_CONTOURS;
extern int STORE_1ST_PASS;
extern Interpo_struct Istru;
extern Name_struct Nstru;
extern int COUNT_TRIANGLE_ONLY;
extern int NO_OUTPUT_FILE;
extern int DIFF_LEVEL_FILES; /* if 1: put diff. level into diff. files */
extern int CALC_VOLUME_ONLY;
extern int NON_VALID;

double Gvolume_ary[10]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
VertType *volumeVertAry=NULL;
int volumeVertAryInd;

/* In order to port to the Intel Paragon, the following file 
are indicated by file descriptor, so they can be opened by 
gopen() */

FILE *Tiled_fd=0, *Untiled_fd=0, *Debug_fd=0;
char poly_file_name[300];

// These are used when different files are required for different levels
#define MAX_LEVEL 20
FILE *TiledFdAry[MAX_LEVEL], *UntiledFdAry[MAX_LEVEL];

FILE *Mesh_fp=NULL, *Cpts_fp=NULL, *CptsFirst_fp=NULL;
float valuefunction(int, int, int);
void slicefunction(int, int);

static char Gfilename[200];

extern int STATIC_CHQ;

void my_clear_contour_read()
{
	strcpy(Gfilename, " ");
	init_find_contours_from_images();
}

void close_files()
{
    int i;

    if(Tiled_fd != 0)
	{
		fclose(Tiled_fd);  //klc
        Tiled_fd=0;
    }
    
	if(Untiled_fd !=0)
	{
		fclose(Untiled_fd); //klc
        Untiled_fd=0;
    }
    
	if(Debug_fd !=0)
	{
		fclose(Debug_fd); //klc
        Debug_fd=0;
    }
    
	for(i=0; i<MAX_LEVEL; i++)
	{
        if(TiledFdAry[i] != 0)	fclose(TiledFdAry[i]);
        TiledFdAry[i]=0;
        if(UntiledFdAry[i] != 0)	fclose(UntiledFdAry[i]);
        UntiledFdAry[i]=0;
    }
	
}

void file_open_error(char *s)
{
    fprintf(stderr,"open_files() Could not open %s for writing\n",s);
    exit(1);
}

void open_files(char *name)
{
    char name1[200],name2[200],name4[200];
    int oflag;
 
	oflag=O_CREAT | O_TRUNC | O_WRONLY;
    
	if(NO_OUTPUT_FILE)return;
    
	if(CALC_VOLUME_ONLY)return;
    
	strcpy(Gfilename,name);
    strcpy(name4,name);
    strcat(name4,"_1st.poly");
	
    if(DIFF_LEVEL_FILES)return;
    
	strcpy(name1,name);
    strcat(name1,".poly");
    strcpy(name2,name);
    strcat(name2,"_untiled.poly");
	
	strcpy(poly_file_name, name1);

#ifdef i860
    if((Tiled_fd=gopen(name1,oflag,M_LOG,0644))==-1) file_open_error(name1);
    if((Untiled_fd=gopen(name2,oflag,M_LOG,0644))==-1)file_open_error(name2);
#else 
    
	Tiled_fd = fopen(name1,"w");
	if (Tiled_fd == 0) file_open_error(name1);
	
	Untiled_fd = fopen(name1,"w");
	if (Untiled_fd == 0) file_open_error(name2);
    
	if(STORE_1ST_PASS)
	{
		if((Debug_fd = fopen(name4, "w"))==0) file_open_error(name1);
    }
#endif
}

FILE* openAnUntiled(int level)
{
    char name[230],exten[10];
    FILE* fd;
	
    sprintf(exten,"L%d",level);
    strcpy(name,Gfilename);
    strcat(name,exten);
    strcat(name,"_untiled.poly");
	
    fd = fopen(name, "w");
	if (fd == 0) file_open_error(name);
    return fd;
}

FILE* openATiled(int level)
{
    char name[230],exten[10];
    FILE* fd;
	
    sprintf(exten,"L%d",level);
    strcpy(name,Gfilename);
    strcat(name,exten);
    strcat(name,".poly");
	
    fd = fopen(name, "w");
	if (fd == 0) file_open_error(name);
    return fd;
}


void tile_all_from_contours(char *out_name, int beg_slice, int dimz, 
							double detz, char *prefix, char *suffix)
{
    int i;
    if(CALC_VOLUME_ONLY)
	{
		extern int SAME_LEVEL;
		SAME_LEVEL=1;
		NON_VALID=1;
		volumeVertAry=(VertType*)mymalloc(sizeof(VertType)*4096);
    }
	
	
#ifndef i860
    open_files(out_name);
#endif
	
	if(NON_VALID)
	{
		for(i=0; i<dimz; i++)
		{
			tile_from_contours(out_name, beg_slice+i, 2, detz, prefix, 
				suffix);
		}
		
		closed_all_tiling_files();
		
		if(CALC_VOLUME_ONLY)
		{
			printf("Level is the number of enclosing contours of a contour.\n");
			printf("For example, the outer contour is of level 0, and the nested contour\n");
			printf("is with level 1, The real volume is (sum of even levels) - (odd levels)\n");
			i=0;
			
			while(Gvolume_ary[i]>1)
			{
				printf("Level %d: volume = %f\n",i,Gvolume_ary[i]);
				i++;
			}
		}
		if(volumeVertAry)free(volumeVertAry);
		volumeVertAry = NULL;
		return;
	}
	
	/* make the file consistant.
	Reason: The tiling algorithm add points to break a contour segments
	into more segments. A contour is used twice. One for the upper section,
	the other for the bottom section. It could be added different vertices
	during the tiling of these two section. Thus the triangulation could
	be no longer valid (each edge is shared exactly twice). Non-valid
	triangulation is all right for visualization, but it can not be used
	in finite element analysis. So We need to make it valid.
	*/
	
    /* step1: Just make the new contours with add vertices. Each contour
	except the top and the bottom one are used twice, and has two
	new contours. Store them in /tmp/tiling.?
    */

    for(i=0; i<dimz; i++)
	{
		tile_from_contours_phase1(beg_slice+i, 2, prefix, suffix, 
			beg_slice, beg_slice+dimz);
    }
	
	printf("done phase one\n");
    
	/* step 2: take the union of the added vertices. */
    for(i=1; i<dimz; i++)
	{ /* read the /tmp/ and make it consistant */
		make_file_consistant(beg_slice+i);
		
    }
	
    for(i=0; i<dimz; i++)
	{
		//tile_from_contours(out_name, beg_slice+i, 2, detz, "/tmp/tiling.","");
		tile_from_contours(out_name, beg_slice+i, 2, detz, "tiling.","");
    }
#ifndef i860
    closed_all_tiling_files();
#endif
	system("rm tiling.*");
//	system("rm /tmp/tiling.*");
//	system("del tiling.*");
}



CTSlice init_volume_and_slice()
{
    CTVolume vol_in; /* a pointer to a structure */
    CTSlice	slice;
    if(MATH_FUNC==0) 
		vol_in = InitVolume(Istru.type, Nstru.prefix, Nstru.suffix,	Istru.beg, 
					Istru.end, Istru.detz);
    else 
		vol_in = myInitVolume16bits(Istru.type, Nstru.prefix, Nstru.suffix,	Istru.beg,
				Istru.end, Istru.detz);

    if(vol_in == NULL) 
	{
		fprintf(stderr,"InterpolateC1: Can not initialize input volume\n");
		fprintf(stderr,"PREFIX=%s SUFFIX=%s, from %d to %d\n",
			Nstru.prefix, Nstru.suffix,Istru.beg, Istru.end);
		exit(1);
    }
    pass_volume(vol_in);
	
	if(DEBUG)print_volume_structure(vol_in);
    if(MATH_FUNC==0) 
	{
		slice = myCTSliceRead(vol_in, Istru.beg, Istru.x1, Istru.y1,Istru.x2, Istru.y2);
		//slice = CTSliceRead(vol_in, Istru.beg, Istru.x1, Istru.y1,Istru.x2, Istru.y2); //KLC
	}
    else if(MATH_FUNC==4)/* special */
		slice = specialCTSliceRead(vol_in, Istru.beg, Istru.x1, Istru.y1,Istru.x2, 
				Istru.y2);
	
    else slice = myCTSliceRead(vol_in, Istru.beg, Istru.x1, Istru.y1,Istru.x2, Istru.y2);

    if(slice==NULL)
	{
		fprintf(stderr,"Read error in slice #%d, quit here\n",Istru.beg);
		exit(1);
    }
    return slice;
}

void init_find_contours_from_images()
{
    int dimx,dimy,dimz;
    static double vunit[3], orig[3];
    static CTSlice	slice=NULL;
	
    dimz=Istru.end - Istru.beg+1;
    slice=init_volume_and_slice();
    dimx=CTSliceWidth(slice)/Istru.xcube_vox;
    dimy=CTSliceHeight(slice)/Istru.ycube_vox;
    CTSliceFree(slice);
    slice=NULL;
    vunit[0]=Istru.xcube_vox;
    vunit[1]=Istru.ycube_vox;
    vunit[2]=Istru.detz*Istru.zcube_vox;
    orig[0]=0.0;
    orig[1]=0.0;
    orig[2]=vunit[2]*Istru.beg;
	
    initialization(Istru.beg, dimx, dimy, dimz, orig, vunit, Istru.thres);
    allocate_basic_memory();
    allocate_image_memory(dimx*dimy);
}

void free_find_contours_from_images()
{
    free_image_memory();
    free_basic_memory();
    free_pos_ary();
}

int  find_contours_from_images() 
{
    int k, dimz;	/* positions of the current  slice */
    
#ifndef i860
    init_find_contours_from_images();
#endif
	
    dimz=Istru.end - Istru.beg+1;
    
	if(COUNT_TRIANGLE_ONLY)
	{
		int tri_num=0;
		for(k=0; k<dimz-1; k++) 
		{
			tri_num+=count_triangles_one_slice(k+Istru.beg,
                valuefunction,slicefunction);
		}
		fprintf(stderr,"\nSlice %d-%d, Total Marching cubes tri. #: %d\n",Istru.beg, 
			Istru.end, tri_num);
		printf("\nSlice %d-%d, Total Marching cubes tri. #: %d\n",Istru.beg,
			Istru.end, tri_num);
    }	
    else
	{
		for(k=0; k<dimz; k++) 
		{
			printf(" --- reading & processing Slice %d\n", k+Istru.beg);
			
			read_and_approximate(k+Istru.beg, valuefunction,slicefunction);
			save_contours_pts_format(Nstru.output_file_name,k+Istru.beg);
		}
    }
	
#ifndef i860
    free_find_contours_from_images();
#endif
    return 1;
}



