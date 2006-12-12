#include <stdio.h>
#include <stdlib.h> //klc
#include <math.h>
#include <memory.h>
#include <malloc.h>
#include <string.h>
#include <sys/types.h>
//#include <unistd.h> //KLC

#include "common.h"
#include "myutil.h"
#include "tile_hash.h"
#include "tile_util.h"
#include "math_util.h"
#include "cross.h"
#include "group.h"
#include "contour_read.h"
#include "tiling_main.h"
#include "decompose.h"
#include "branch_util.h"
#include "branch.h"


extern  int SILENT;
static short **New_group_ary[2], *New_group_ind_ary[2];
static VertType *Lines_ary[2];
static int Group_num[2];
static int Untiled_triangle_number=0,Convex_triangle_number=0;
extern TileType **Best_tiling_table[2];
extern double Mid_zvalue;

extern int STATIC_CHQ;

void my_clear_branch()
{
	New_group_ary[0] = NULL;		New_group_ary[1] = NULL; 
	New_group_ind_ary[0] = NULL;	New_group_ind_ary[1] = NULL;
	Lines_ary[0] = NULL;			Lines_ary[1] = NULL;
	Group_num[0] = Group_num[1] = 0;
	Untiled_triangle_number=0;		Convex_triangle_number=0;	
}


int get_untiled_triangle_number()
{
    return Untiled_triangle_number;
}


int store_unused_contour(VertType **untiled_vert_ary, int *complete_index)
/* complete_index: the real
   return (complete_index+1)&0xffffe */
{
  #define ARY_SIZE 1000
  int topdown,group,sum,i,isum;
  int ind0,ind1,k,k1;
  VertType *vert_ary[ARY_SIZE],*temp_ary;
  int vert_ary_num[ARY_SIZE], vert_ary_index=0;
  char complete_flag[ARY_SIZE];

  for(topdown=0; topdown<2; topdown++) 
  {
    VertType *vp; 
    vp=Lines_ary[topdown];
  
	for(group=0; group<Group_num[topdown]; group++)
	{
		int temp=New_group_ind_ary[topdown][group];
		short *group_tab=New_group_ary[topdown][group];
		TileType *best_tiling_pointer;

		best_tiling_pointer = Best_tiling_table[topdown][group];
		vert_ary_num[vert_ary_index]=0;
		
		for(ind0=0; ind0<temp; ind0++)
		{
			if(!best_tiling_pointer[ind0].used)vert_ary_num[vert_ary_index]++; 
		}

		if(temp==vert_ary_num[vert_ary_index])
		{
			/*  completely unsed */
			complete_flag[vert_ary_index]=1;
			vert_ary_num[vert_ary_index]++;
			temp_ary=(VertType*) mymalloc (sizeof(VertType)*(temp+2));
			
			for(ind0=0; ind0<temp; ind0++)
			{
				if(topdown)		k=group_tab[ind0];
				else		k=group_tab[temp-1-ind0];
				memcpy(&(temp_ary[ind0]), &(vp[k]),sizeof(VertType));
			}

			memcpy(&(temp_ary[temp]), &(temp_ary[0]),sizeof(VertType));

			if(vert_ary_index>=ARY_SIZE-1)
			{
				fprintf(stderr,"vert_ary_index=%d overflow\n",vert_ary_index);
				exit(1);
		    }
			vert_ary[vert_ary_index++]=temp_ary;
		}
		else if(vert_ary_num[vert_ary_index]>0)
		{
			int temp_ind=0;
			complete_flag[vert_ary_index]=0;
			vert_ary_num[vert_ary_index]*=2;
			temp_ary=(VertType*) mymalloc(sizeof(VertType)*
			(vert_ary_num[vert_ary_index]+2)); /* each line has 2 V. */
			for(ind0=0; ind0<temp; ind0++)
			{
				if(best_tiling_pointer[ind0].used ==0)
				{ 
					/* not used yet */
					k=group_tab[ind0];
					ind1=ind0+1;
					if(ind1>=temp) ind1=0;
					k1=group_tab[ind1];

					if(topdown)
					{ 
						memcpy(&(temp_ary[temp_ind++]),
						&(Lines_ary[topdown][k]),sizeof(VertType));
						memcpy(&(temp_ary[temp_ind++]),
						&(Lines_ary[topdown][k1]), sizeof(VertType));
					}
					else
					{ 
						/* reverse direction */
						memcpy(&(temp_ary[temp_ind++]),
						&(Lines_ary[topdown][k1]), sizeof(VertType));
						memcpy(&(temp_ary[temp_ind++]),
						&(Lines_ary[topdown][k]),sizeof(VertType));
					}
				}
		    }

			if(temp_ind!=vert_ary_num[vert_ary_index] && !SILENT)
				printf("***** Error, temp_ind=%d, should be %d\n",
				temp_ind, vert_ary_num[vert_ary_index]);
			if(vert_ary_index>=ARY_SIZE-1)
			{
				fprintf(stderr,"store_unused_contour(), mem. overflow\n");
				exit(1);
			}
			vert_ary[vert_ary_index++]=temp_ary;
		}
	}
  } //end of the two for loops....
  
  isum=0;

  for(i=0; i<vert_ary_index; i++)	
	  isum+=vert_ary_num[i];

  *untiled_vert_ary=(VertType*) mycalloc(sizeof(VertType)*(isum+4));
  sum=0;
  
  for(i=0; i<vert_ary_index; i++)
  {
    if(complete_flag[i])
	{
		memcpy(&((*untiled_vert_ary)[sum]), vert_ary[i], sizeof(VertType)*vert_ary_num[i]);
		sum+=vert_ary_num[i];
    }
  }

  if(sum&1)sum++;
  *complete_index=sum;
  for(i=0; i<vert_ary_index; i++)
  {
	if(vert_ary_num[i]>0)
	{
		if(!complete_flag[i])
		{
			memcpy(&((*untiled_vert_ary)[sum]), vert_ary[i], sizeof(VertType)*vert_ary_num[i]);
			sum+=vert_ary_num[i];

			if(sum>isum+1)
			{
				printf("*** Error! over: sum=%d ,size=%d\n",sum,isum);
				exit(1);
			}
		}
	}
  } //end of for loop...

  for(i=0; i<vert_ary_index; i++)	free(vert_ary[i]);

  return sum;
}

int store_boundary_tiling(VertType **untiled_vert_ary, int index)
/* store the definitely boundary in one file, and the possible boundary
   at the end */
{
    #define BAK_ARY_SIZE 2000
    VertType vert_ary[BAK_ARY_SIZE],*temp_ary;
    int ind0,ind1,k,k1,sum;
    int tile_dir;
    int group0,group1;
    int vert_ind=0;

    while((pop_boundary_array(0,&group0,&ind0,&tile_dir, &group1,&ind1))
         || (pop_boundary_array(1,&group0,&ind0,&tile_dir, &group1,&ind1)))
	{
		/* check if this tiling is done */
		if(is_in_tile_hash_table(New_group_ary[0][group0][ind0],New_group_ary[1][group1][ind1]))
		{
			printf("*****Error G%d %d G%d %d dir=%d is in hash\n",
			group0,ind0,group1,ind1,tile_dir);
		}
		else 
		{
			if(tile_dir!=3)
			{
				k =New_group_ary[0][group0][ind0];
				k1=New_group_ary[1][group1][ind1];
				if(vert_ind>=BAK_ARY_SIZE-2)
				{
					fprintf(stderr,"Stack over, vert_ind=%d, store_b.._tiling()\n", vert_ind);
					exit(1);
				}
				
				if(vert_ind>=BAK_ARY_SIZE-4)
				{
					fprintf(stderr,"store_boundary_tiling(), mem. over\n");
					exit(1);
				}
				
				if(tile_dir==2)
				{
					memcpy(&(vert_ary[vert_ind++]),&(Lines_ary[1][k1]),sizeof(VertType));
					memcpy(&(vert_ary[vert_ind++]),&(Lines_ary[0][k]),sizeof(VertType));
				}
				else 
				{
					memcpy(&(vert_ary[vert_ind++]),&(Lines_ary[0][k]),sizeof(VertType));
					memcpy(&(vert_ary[vert_ind++]),&(Lines_ary[1][k1]),sizeof(VertType));
				}
		    }
		    else 
				fprintf(stderr,"G%d %d G%d %d dir=%d is in hash\n",
						group0,ind0,group1,ind1,tile_dir);
		}
    } //end of while loop...

    sum=index+vert_ind;
    temp_ary=(VertType*)mymalloc(sizeof(VertType)*(sum+4));
    memcpy(temp_ary, *untiled_vert_ary, sizeof(VertType)*index);
    memcpy(&temp_ary[index], vert_ary, sizeof(VertType)*vert_ind);
    free(*untiled_vert_ary);
    *untiled_vert_ary=temp_ary;
    return sum;
}

int get_convex_triangle_number()
{
    return Convex_triangle_number;
}

void  draw_convex_triangles(PolygonStruct **ret_poly_ary)
{
    int i;
    PolygonStruct *poly;
    int cnt=0;

    i=0;
    while (ret_poly_ary[i])
	{
		poly=ret_poly_ary[i];

		draw_tiled_triangle_no_mesh(&(poly->vertices[0]), 
			&(poly->vertices[1]),&(poly->vertices[2]));

		cnt++;
		i++;
    }
    Convex_triangle_number+=cnt;
}


void draw_untiled_contours(VertType *untiled_vert_ary,  
	short *untiled_group_ary, short num)
{
	extern int DIFF_LEVEL_FILES;
	extern int Current_level;
	extern FILE* UntiledFdAry[];
	extern FILE* Untiled_fd;
	extern int NO_OUTPUT_FILE;

    int ind,j,len;

    #define BUF_SIZE 40960

    char buf[BUF_SIZE];
    FILE* fd;
    if(NO_OUTPUT_FILE)return ;


    if(DIFF_LEVEL_FILES)
	{
		if(UntiledFdAry[Current_level] !=0) 
		{
			UntiledFdAry[Current_level]=openAnUntiled(Current_level); 
		}
		fd=UntiledFdAry[Current_level];
    }
    else 
	{
		if(Untiled_fd != 0)	return; /* no output desired */
		fd=Untiled_fd;
    }

    sprintf(buf,"%d\n",num);
    len=strlen(buf);

    for(j=0; j<num; j++)
	{
		if(Current_level &1)
		{ // need to reverse the sequence
			ind=untiled_group_ary[num-1-j];
		}
		else ind=untiled_group_ary[j];

		sprintf(&buf[len],"%lf %lf %lf %lf\n",untiled_vert_ary[ind][0],
			untiled_vert_ary[ind][1],untiled_vert_ary[ind][2],Mid_zvalue);
		len=strlen(buf);
		if(len>=BUF_SIZE -100)
		{
		   fprintf(stderr,"ERR draw_untiled_contours, buf overflow\n");
		   exit(1);
		}
    } //end of for loop...

    //j=write(fd,buf,len); //klc
	j = fwrite(buf, sizeof(char) , len , fd); //klc
    if(j!=len)
	{
		fprintf(stderr,"*** Error write, draw_untiled_.() should %d, but %d\n",	len,j);
		exit(1);
    }
}
/*
int triangulate_with_center(VertType *vert_ary, short *vert_index_ary,
	int poly_num, VertType *vt)

{
    int i,ind,ind1,j,num;

    extern Tile_triangle_num;
    extern FILE *Ofp;
    vert_index_ary[poly_num]=vert_index_ary[0];
    for(j=0; j<poly_num; j++){
	ind=vert_index_ary[j];
	ind1=vert_index_ary[j+1];
	draw_triangle(Ofp,&untiled_vert_ary[ind],&untiled_vert_ary[ind],vt);
	Tile_triangle_num++;
	}
    }
    
}
*/

/* triangulate the polygon vs index=start */
PolygonStruct ** breakAllIntoTris(int start, VertType *vert_ary,short *	group_ary, short poly_num)
{
    PolygonStruct **poly_ary;
    PolygonStruct *poly;
    int i,size,k1,kk;

    size=poly_num-2;
    poly_ary=(PolygonStruct**)mycalloc((size+1)*sizeof(*poly_ary));
    
	for(i=0; i<size; i++)
	{
		int ind_ary[3];
		poly=(PolygonStruct *)mycalloc(sizeof(*poly));
		poly_ary[i]=poly;
		poly->numpts=3;
		poly->normals=NULL;
		poly->vertices=(VertType *)mymalloc(sizeof(VertType)*3);
		/* The 3 triangles are start-i-2, start-i-1 and start */
		for(kk=0; kk<2; kk++)
		{
			ind_ary[kk] = (start -i -2 +kk +poly_num)%poly_num;
		}
		
		ind_ary[2] = start;
        
		for(kk=0; kk<3; kk++)
		{
			k1=group_ary[ind_ary[kk]];
            memcpy(&(poly->vertices[kk]), &(vert_ary[k1]), sizeof(VertType));
        }
    }
    poly_ary[i]=NULL;
    return poly_ary;
}

/* check it is a triangle or not */
PolygonStruct ** check_is_a_triangle(VertType *vert_ary,
	short *	group_ary, short poly_num)
{
    PolygonStruct *poly,**poly_ary;
    int cnt=0,kk,k1;
    double zval;

    if(poly_num!=3)	return NULL;

    zval=vert_ary[group_ary[0]][2];
    
	if(zval!=vert_ary[group_ary[1]][2])		cnt=1;
    else if(zval!=vert_ary[group_ary[2]][2])	cnt=1;

    if(cnt)
	{ 
		/* a legal triangle */
		if(!SILENT)printf("*** Warning, An untiled triangle is legal\n");
		poly_ary=(PolygonStruct**)mycalloc(2*sizeof(*poly_ary));
		poly=(PolygonStruct *)mymalloc(sizeof(*poly));
		poly_ary[0]=poly;
		poly->numpts=3;
		poly->normals=NULL;
		poly->vertices=(VertType *)mymalloc(sizeof(VertType)*poly_num);
		
		for(kk=0; kk<poly_num; kk++)
		{
			k1=group_ary[kk];
				memcpy(&(poly->vertices[kk]), &(vert_ary[k1]), sizeof(VertType));
		}
    }
    else
		return NULL;

    poly_ary[1]=NULL;
    
	return poly_ary;
}


/* check it is a quad or not */
PolygonStruct ** check_is_2_triangles(VertType *vert_ary, short *group_ary, short poly_num)
{
    int i, cnt=0,start;
    int flag_ary[4];
    double dist_ary[2];
    double zval;

    if(poly_num!=4)return NULL;

    zval=vert_ary[group_ary[0]][2];
    
	for(i=0; i<poly_num; i++)
	{
		if(vert_ary[group_ary[1]][2] == zval)	flag_ary[i]=0;
		else flag_ary[i]=1;
    }

    for(i=0; i<poly_num; i++)	cnt+=flag_ary[i];

    if(cnt ==4 || cnt==0)return NULL; /* They all on the same plane */
    
    for(i=0; i<2; i++)
	{
		if(flag_ary[i]==flag_ary[i+2]) dist_ary[i]=1e12;
		else 
		{
			int j1,j2;
			j1=group_ary[i];
			j2=group_ary[i+2];
			dist_ary[i]=line_distance(&vert_ary[j1], &vert_ary[j2]);
		}
    }
    
	if(dist_ary[0]<dist_ary[1])	start=0;
    else start=1;
    return breakAllIntoTris(start, vert_ary, group_ary, poly_num); 
}

int collect_untiled_contours(
            VertType *lines_ary0, VertType *lines_ary1,
            short ** new_group_ary0, short *new_group_ind_ary0, int group_num0,
            short ** new_group_ary1, short *new_group_ind_ary1, int group_num1,
	    VertType *untiled_vert_ary,
            short ** untiled_group_ary, 
	    short *untiled_group_ind_ary)
{
    int i,ary_index=0,complete_index=0,unused_index;
    int untiled_group_num;
    VertType vt;
    /* free memory */

    if(untiled_group_ary!=NULL) 
	{
		/* do not use free_group_ary() */
		int i=0;
		
		while(untiled_group_ary[i])
		{
			free(untiled_group_ary[i]);
			untiled_group_ary[i]=NULL;
		}
    }
    else
		printf("*** Err! untiled_group_ary=NULL\n");

    New_group_ary[0]=new_group_ary0;
    New_group_ary[1]=new_group_ary1;
    New_group_ind_ary[0]=new_group_ind_ary0;
    New_group_ind_ary[1]=new_group_ind_ary1;
    Lines_ary[0]=lines_ary0;
    Lines_ary[1]=lines_ary1;
    Group_num[0]=group_num0;
    Group_num[1]=group_num1;

    /* count how many line segments, predict the size of untiled */

    if(untiled_vert_ary)
	{
		free(untiled_vert_ary);
		untiled_vert_ary=NULL;
    }

    unused_index=store_unused_contour(&untiled_vert_ary, &complete_index);
    
	if(!SILENT)
		fprintf(stderr,"store_unsed_contour, complete_index=%d, unused_index=%d\n",
				complete_index,unused_index);
    
	ary_index=store_boundary_tiling(&untiled_vert_ary, unused_index);
    if(!SILENT)
		fprintf(stderr,"store_boundary_tiling, ary_index=%d\n",ary_index);
/*
{
FILE *ofp;
if((ofp=fopen("unorder.poly","w"))==NULL) {
        fprintf(stderr,"Could not open unorder.poly for writing\n");
        exit(1);

}
for(i=complete_index; i<ary_index; i=i+2){
    fprintf(ofp,"3\n%lf %lf %lf\n",untiled_vert_ary[i][0],
	untiled_vert_ary[i][1], untiled_vert_ary[i][2]);
    fprintf(ofp,"%lf %lf %lf\n",untiled_vert_ary[i+1][0],
	untiled_vert_ary[i+1][1], untiled_vert_ary[i+1][2]);
    fprintf(ofp,"%lf %lf %lf\n",
	(untiled_vert_ary[i+1][0]+untiled_vert_ary[i][0])/2.0,
	(untiled_vert_ary[i+1][1]+untiled_vert_ary[i][1])/2.0,
	(untiled_vert_ary[i+1][2]+untiled_vert_ary[i][2])/2.0);
}
close (ofp);
}

*/
    untiled_group_num=trace_untiled_contours(untiled_vert_ary, complete_index, unused_index,
											ary_index, untiled_group_ary, untiled_group_ind_ary);

    vt[2]=Mid_zvalue;
    
	for(i=0; i<untiled_group_num; i++)
	{
		int result;
		result = process_an_untiled_region( untiled_vert_ary,untiled_group_ary[i],
											untiled_group_ind_ary[i]);

		if(!result)
		{ 
			/* can not be processed, save into a file */
			fprintf(stderr,"--- Error!, Fail in median axis, poly_num=%d\n",
			untiled_group_ind_ary[i]);
			draw_untiled_contours(untiled_vert_ary, untiled_group_ary[i],untiled_group_ind_ary[i]);
		}
    }

    free(untiled_vert_ary);
    untiled_vert_ary=NULL;
	
    return(untiled_group_num);
}


