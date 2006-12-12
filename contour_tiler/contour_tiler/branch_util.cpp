#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <malloc.h>
#include <string.h>
#include <sys/types.h>
//#include <unistd.h> //KLC

#include "common.h" //KLC as its already included in the decompose.h file....
#include "branch_util.h"
#include "myutil.h"
#include "math_util.h"
#include "branch.h"
#include "decompose.h" //KLC
//#include "cross.h" //KLC included for the free_polygons_ary(poly_ary) method....

extern  int SILENT;
extern double Mid_zvalue;
extern int MESH_FILE;

extern int STATIC_CHQ;

void my_clear_branch_util()
{
	//hopefully nothing to do here..
}

/* file format for .cpts file:
num               This is a polygon.
x y z
...
...
num mode        The line segment, mode = 0, 1: If this slice is used as
...			top (0), then the contour need to be broken into
...			several contours. 
*/

/* return: 2: it is a hole
1: it has two slices, and is OK
0: wrong data
*/

/*
static int verify_adjacent_vertices(VertType *vert_ary, short *group_ary,
short poly_num, int topdown)
{
int i,cnt=0,flag;
int error=0,cnt1=0;
for(i=0; i<poly_num; i++){
VertType *v=&(vert_ary[group_ary[i]]);
if((*v)[2]>Mid_zvalue)flag=1;
else flag=0;
if(flag==topdown)cnt++;
else {
cnt=0;
cnt1++;
}
if(cnt>1) error=1;
}
if(!cnt1)return 2;
if(error && cnt1){
printf("Error! verify_adjacent_vertices(num=%d, topdown=%d)\n",
poly_num,topdown);
return 0;
}
return 1;
}
*/


static int verify_adjacent_vertices(VertType *vert_ary, short *group_ary,
									short poly_num, int topdown)
{
    int i,cnt0=0,cnt1=0;

    for(i=0; i<poly_num; i++)
	{
		VertType *v=&(vert_ary[group_ary[i]]);
		if((*v)[2]>Mid_zvalue)cnt1++;
		else cnt0++;
    }
    
	if(!cnt1 || !cnt0)return 2;
    return 1;
}
int trace_medial_axis(VertType *line_ary[2], int line_ind,
					  VertType *ret_ary[10], int ret_ind_ary[10])
{
    int i,j,start=0,aryi=0;
    char *flag;
    double max,min;
    VertType *v,*v0;
    max = Mid_zvalue+0.0001;
    min = Mid_zvalue-0.0001;
    flag=(char*)mycalloc(line_ind);
	
    /* the array is counted twice for both direction, get rid of one */
    for(i=0; i<line_ind; i++)
	{
		if(!flag[i])
		{
			int cnt=0;
		
			for(j=i+1; j<line_ind; j++)
			{
				if(!flag[j])
				{
					if(line_distance(&(line_ary[0][i]),&(line_ary[1][j])) <1e-6 &&
						line_distance(&(line_ary[0][j]),&(line_ary[1][i])) <1e-6)
					{
						cnt++;
						flag[j]=1;
					}
				}
			}
			if(cnt!=1)printf("Error! trace_medial_axis(),i=%d, cnt=%d\n",
				i,cnt);
		}
    }
	
    while(1)
	{
		int ind_ary[5], ind, done=0,mode,mode_ary[5];
		/* find the first point which is on the merging contour*/
		while(start<line_ind)
		{
			if(!flag[start])
			{
				for(mode=0; mode<2; mode++)
				{
					double val=line_ary[mode][start][2];
					if(val<min || val > max) goto break17; /* don't use break */
				}
			}
			start++;
		}
break17:;
		if(start>=line_ind)break; /* all done */
		
		ret_ary[aryi] = (VertType*)mymalloc(sizeof(VertType)*100); 
		ret_ind_ary[aryi]=0;
		
		/* use slow but simple method to trace */
		v0=&(line_ary[mode][start]);
		memcpy(&ret_ary[aryi][ret_ind_ary[aryi]], v0, sizeof(VertType));
		ret_ind_ary[aryi]++;
		v=&(line_ary[1-mode][start]);
		flag[start]=1;

		while(!done)
		{
			ind=0;
			for(i=0; i<line_ind; i++)
			{
				if(!flag[i])
				{ /* not used yet */
					if(line_distance(v, &(line_ary[0][i]))<1e-6)
					{
						mode_ary[ind]=0;
						ind_ary[ind++]=i;
					}
					else if(line_distance(v, &(line_ary[1][i]))<1e-6)
					{
						mode_ary[ind]=1;
						ind_ary[ind++]=i;
					}
				}
			}
			if(ind==0)
			{ 
				/* all done */
				memcpy(&ret_ary[aryi][ret_ind_ary[aryi]], v, sizeof(VertType));
				ret_ind_ary[aryi]++;
				done=1;
				aryi++;
			}
			else if(ind==1)
			{ 
				/* no confusion */
				v0=&(line_ary[mode_ary[0]][ind_ary[0]]);
				memcpy(&ret_ary[aryi][ret_ind_ary[aryi]], v0, sizeof(VertType));
				ret_ind_ary[aryi]++;
				v=&(line_ary[1-mode_ary[0]][ind_ary[0]]);
				flag[ind_ary[0]]=1;
			}
			else
			{
				/* more than one choice */ 
				/* find the CCW one */
				printf("Error! trace_medial_axis(), ind=%d>1, NOT implemented\n",ind); 
				free(flag);
				return 0;
				/*
				VertType *u0,*u,w0,w1;
				double max_angle=-1e6;
				int maxi;
				w0[0]=(*v)[0]-(*v0)[0];
				w0[1]=(*v)[0]-(*v0)[1];
				normalize_vector(&w0);
				for(k=0; k<ind; k++){
				double angle;
				u0=&(line_ary[mode_ary[k]][ind_ary[k]]);
				u=&(line_ary[1-mode_ary[k]][ind_ary[k]]);
				w1[0]=(*u)[0]-(*u0)[0];
				w1[1]=(*u)[0]-(*u0)[1];
				normalize_vector(&w1);
				angle=calc_angle(&w0,&w1);
				if(max_angle<angle){
				max_angle=angle;
				maxi = k;
				}
				}
				k=ind_ary[maxi];
				v0=&(line_ary[mode_ary[maxi]][k]);
				memcpy(&ret_ary[aryi][ret_ind_ary[aryi]], v0, sizeof(VertType));
				ret_ind_ary[aryi]++;
				v=&(line_ary[1-mode_ary[maxi]][k]);
				flag[k]=1;
				*/
			}
		}
    }
    ret_ary[aryi] = NULL;
    free(flag);
    return aryi;
}

static int find_median_axis( PolygonStruct ** poly_ary, int topdown, 
							VertType **ret_ary)
{
/* collect all the line segments 1) in the middle, 
	2) one is on the merging contour, and the other is middle */
    
    int i;
    PolygonStruct *poly;
    int tri_cnt;
    VertType *line_ary[2];
    int line_ind=0;
    double max,min;
	VertType *temp_ary[10]; 
	int temp_ind_ary[10];
	int aryi;
	
    max = Mid_zvalue+0.0001;
    min = Mid_zvalue-0.0001;
	
    tri_cnt=0;
    while (poly_ary[tri_cnt])	tri_cnt++;

    for(i=0; i<2; i++)line_ary[i]=(VertType*)mymalloc(sizeof(VertType)*tri_cnt*3);
	
    i=0;
    
	while (poly_ary[i])
	{
		int pos[3],j;
        poly=poly_ary[i];
		
		for(j=0; j<3; j++)
		{
			if(poly->vertices[j][2] > max)pos[j]=1;
			else if(poly->vertices[j][2] < min)pos[j]=0;
			else pos[j]=2;
		}
		
		for(j=0; j<3; j++)
		{
			int j1=(j+1)%3;
			int res=0;
			if(pos[j]==2)
			{
				if(pos[j1]==2)res=1;
				else if(pos[j1]==topdown)res=1;
			}
			else if(pos[j]==topdown)
			{
				if(pos[j1]==2)res=1;
			}
			if(res)
			{
				memcpy(&(line_ary[0][line_ind]), &(poly->vertices[j]),sizeof(VertType));
				memcpy(&(line_ary[1][line_ind]), &(poly->vertices[j1]),sizeof(VertType));
				line_ind++;
			}
		}
		i++;
    }
	
	aryi=trace_medial_axis(line_ary,line_ind, temp_ary, temp_ind_ary);
	/* write the sequence into file */
	if(aryi!=1)
	{
		printf("Error! find_median_axis(), aryi=%d\n",aryi);
	}
	
	for(i=1; i<aryi; i++)free(temp_ary[i]);
	*ret_ary=temp_ary[0];
	
    for(i=0; i<2; i++)free(line_ary[i]);
    return temp_ind_ary[0];
}

void convert_middle_line(VertType *ret_ary, int ret_ary_ind, double val)
{
    int i;
    double max,min;
    max = Mid_zvalue+0.0001;
    min = Mid_zvalue-0.0001;
	
    /* need to convert all the middle point to the topdown slice */
    for(i=0; i<ret_ary_ind; i++)
	{
		double val0;
		val0 = ret_ary[i][2];
		if(val0<max && val0>min)ret_ary[i][2]=val;
    }
}

/* let all triangle lay down on the merging slice */
void convert_middle_tri_ary(PolygonStruct ** poly_ary, double val)
{
    PolygonStruct *poly;
    int i,j;
    double max,min;
    max = Mid_zvalue+0.0001;
    min = Mid_zvalue-0.0001;
    i=0;
    
	while(poly_ary[i])
	{
		poly = poly_ary[i];
		
		for(j=0; j<3; j++)
		{
			double val0 = poly->vertices[j][2]; 
			if(val0>min && val0<max) 
				poly->vertices[j][2]=val;
		}
		i++;
    }
}

/* need to append the current medial axis to .cpts file */
int append_median_axis(VertType *ary, int ary_ind, int topdown)
{
	extern char Out_name[];
	extern int Current_slice_num;
    char str[256];
    FILE *ofp;
    int i;
	
	
    if(topdown)		sprintf(str,"%s%d.cpts",Out_name, Current_slice_num);
    else	sprintf(str,"%s%d.cpts",Out_name, Current_slice_num-1);

    if((ofp=fopen(str,"a"))==NULL) 
	{
		/* append */
		fprintf(stderr,"Error! cannot open %s for appending\n",str);
		return 0;
    }

    fprintf(ofp,"%d %d\n",ary_ind,topdown);
    
	for(i=0; i<ary_ind; i++)
	{
		fprintf(ofp,"%f %f %f\n",ary[i][0],ary[i][1],ary[i][2]);
    }
    
	fclose(ofp);
    return 1;
}

/* If a CCW middle triangle has a top vertex, then project the
middle vertices to the bottom slice (smaller Z).

 If a CW middle triangle has a bottom vertex, ....
	..................top slice.
*/

/* Rule for reading .cpts files

 #define LINE_UP 0
 #define LINE_DW 1
 #define LINE    2
 #define RANDOM_POINT_UP 4
 #define RANDOM_POINT_DW 5
 #define RANDOM_POINT 6
 #define CONTOUR 8
 
  NOTE: the points might be repeated, so should get rid of repeated point
  The triangle is arranged in CW sense. So invert it
*/
int append_Steiner_lines(PolygonStruct **poly_ary, int isCCW)
{
	extern char Out_name[];
	extern int Current_slice_num;
	extern int Current_level;
	extern double Zvalues[3];
    char str[256];
    FILE *ofp[2];
    int i,j,count;
    char *flag;
    double max,min;
    VertType *line_ary[200][2];
    VertType *pt_ary[100],*pt1_ary[100];
    int line_ind=0, anchor_slice,pt_ind=0,pt1_ind=0;

    /* 0: bottom slice, 1: top slice (large Z), 2: middle slice */
    
    if(Current_level &1)	isCCW=1-isCCW;
    
	isCCW=1-isCCW; /* invert it */
    anchor_slice = isCCW;
	
    for(i=0; i<2; i++)
	{
		sprintf(str,"%s%d.cpts",Out_name, Current_slice_num-1+i);
		if((ofp[i]=fopen(str,"a"))==NULL) 
		{
			/* append */
			fprintf(stderr,"Error! cannot open %s for appending\n",str);
			return 0;
		}
    }
	
    max = Mid_zvalue+0.0001;
    min = Mid_zvalue-0.0001;
	
    i=0;

    while (poly_ary[i])
	{
		PolygonStruct *poly;
		VertType *p[3];
		int slice[3]; /* 0,1,2: bottom, top, middle */
		int cnt[3],line[2];
		
        poly=poly_ary[i];
		for(j=0; j<3; j++)
		{
			p[j] = &(poly->vertices[j]);
			cnt[j]=0;
		}
		
		for(j=0; j<3; j++)
		{
			if((*p[j])[2]>max)slice[j]=1;
			else if((*p[j])[2]<min)slice[j]=0;
			else slice[j]=2;
		}
		
		for(j=0; j<3; j++)cnt[slice[j]]++;
		
		if(!cnt[2])
		{
			printf("Error! no middle vertex\n");
		}
		
		line[0]=-1;
		if(cnt[anchor_slice]==1)
		{ /* has one anchor (top for CCW) vertex */
			for(j=0; j<3; j++)
			{
				if(slice[j]==anchor_slice)
				{
					line[0]=(j+1)%3;
					line[1]=(j+2)%3;
				}
			}
			/* need to save the line segment */
			if(line[0]>=0)
			{
				line_ary[line_ind][0]=p[line[0]];
				line_ary[line_ind][1]=p[line[1]];
				line_ind++;
			}
		}
		else if(cnt[anchor_slice]==2)
		{
			/* has one anchor (top for CCW) vertex */
			for(j=0; j<3; j++)
			{
				if(slice[j]==2) pt_ary[pt_ind++] = p[j];
			}
		}
		
		for(j=0; j<3; j++)
		{
			/* store all middle point */
			if(slice[j]==2) pt1_ary[pt1_ind++] = p[j];
		}
		i++;
    }
	
    /* write the middle points as reference points only */
    count=0;
    flag=(char*)mycalloc(pt1_ind);
    for(j=0; j<pt1_ind; j++)
	{
		int doit=1;
		for(i=0; (i<j) && doit; i++)
		{ 
			if(line_distance(pt1_ary[j], pt1_ary[i])<1e-10)doit=0;
		}
		
		if(doit)
		{
			flag[j]=1;
			count++;
		}
		else flag[j]=0;
    }
    if(count)
	{
		fprintf(ofp[anchor_slice],"%d 7\n",count); /* reference  */
		for(j=0; j<pt1_ind; j++)
		{
			if(flag[j])
			{
				fprintf(ofp[anchor_slice],"%f %f %f\n",
					(*pt1_ary[j])[0],(*pt1_ary[j])[1], (*pt1_ary[j])[2]);
			}
		}
    }
	
    /* write the middle points as reference points only */
    count=0;
    for(j=0; j<pt_ind; j++)
	{
		int doit=1;
		
		for(i=0; (i<j) && doit; i++)
		{ 
			if(line_distance(pt_ary[j], pt_ary[i])<1e-10)doit=0;
		}
		
		if(doit)
		{
			for(i=0; (i<pt1_ind) && doit; i++)
			{ 
				if(line_distance(pt_ary[j], pt1_ary[i])<1e-10)doit=0;
			}
		}

		if(doit)
		{
			flag[j]=1;
			count++;
		}
		else flag[j]=0;
    }

    if(count)
	{
		fprintf(ofp[anchor_slice],"%d 7\n",count); /* reference  */
		for(j=0; j<pt_ind; j++)
		{
			if(flag[j])
			{
				fprintf(ofp[anchor_slice],"%f %f %f\n",
					(*pt_ary[j])[0],(*pt_ary[j])[1], (*pt_ary[j])[2]);
			}
		}
    }
	
	
    /* write lines to the proper slice */
    for(j=0; j<line_ind; j++)
	{
		int doit=1;
		for(i=0; (i<j) && doit; i++)
		{ /* check whether it is repeated or not */
			if(line_distance(line_ary[i][0], line_ary[j][0])<1e-10 &&
				line_distance(line_ary[i][1], line_ary[j][1])<1e-10)doit=0;
			if(doit && line_distance(line_ary[i][1], line_ary[j][0])<1e-10 &&
				line_distance(line_ary[i][0], line_ary[j][1])<1e-10)doit=0;
		}
		
		if(doit)
		{
			fprintf(ofp[1-anchor_slice],"2 %d\n",isCCW);
			fprintf(ofp[1-anchor_slice],"%f %f %f\n",
				(*line_ary[i][0])[0],(*line_ary[i][0])[1],
				Zvalues[1-anchor_slice]);
			fprintf(ofp[1-anchor_slice],"%f %f %f\n",
				(*line_ary[i][1])[0],(*line_ary[i][1])[1],
				Zvalues[1-anchor_slice]);
		}
    }

    /* write points to the proper slice */
    for(j=0; j<pt_ind; j++)
	{
		int doit=1;
		for(i=0; (i<line_ind) && doit; i++)
		{ 
			if(line_distance(pt_ary[j], line_ary[i][0])<1e-10 || 
				line_distance(pt_ary[j], line_ary[i][1])<1e-10)doit=0;
		}
		
		for(i=0; (i<j) && doit; i++)
		{ 
			if(line_distance(pt_ary[j], pt_ary[i])<1e-10)doit=0;
		}
		
		if(doit)
		{
			fprintf(ofp[1-anchor_slice],"1 %d\n",isCCW+4);
			fprintf(ofp[1-anchor_slice],"%f %f %f\n",
				(*pt_ary[i])[0],(*pt_ary[i])[1],
				Zvalues[1-anchor_slice]);
		}
    }
	
    fclose(ofp[0]);
    fclose(ofp[1]);
	free(flag);
    return 1;
}


/* need to append the triangles to .tile file */
int append_triangle_ary(PolygonStruct **poly_ary, int mode)
{
    extern char Out_name[];
    extern int Current_slice_num;
    extern int Current_level;
    int i,j;
    extern FILE *Mesh_fp;
	
    if(Mesh_fp==NULL)return 0;
	/*
    char str[256];
    sprintf(str,"%s%d.tile",Out_name, Current_slice_num-1);
    if((ofp=fopen(str,"a"))==NULL) {
	fprintf(stderr,"Error! cannot open %s for appending\n",str);
	return 0;
    }
	*/
    i=0;

    while (poly_ary[i])
	{
		PolygonStruct *poly;
		VertType *p[3];
		
        poly=poly_ary[i];
		for(j=0; j<3; j++)	p[j] = &(poly->vertices[j]);
		
		if(Current_level &1)
		{ // need to reverse the sequence
			if(mode)draw_triangle_mode_FP(Mesh_fp,p[1],p[0],p[2]);
			else draw_triangle_FP(Mesh_fp,p[1],p[0],p[2]);
		}
		else 
		{
			if(mode)draw_triangle_mode_FP(Mesh_fp,p[0],p[1],p[2]);
			else draw_triangle_FP(Mesh_fp,p[0],p[1],p[2]);
		}
        i++;
    }
    return 1;
}


int process_an_untiled_region(VertType *vert_ary, short * group_ary, 
							  short poly_num)
{
    int isCCW=0,mode=0;
    PolygonStruct ** poly_ary;
	//   VertType *ret_ary; //KLC
	//   int ret_ary_ind; //KLC
	
    isCCW=is_contour_CCW(vert_ary, group_ary, poly_num);
    poly_ary=check_is_a_triangle(vert_ary, group_ary,poly_num);

    if(!poly_ary)
	{
		poly_ary=check_is_2_triangles(vert_ary, group_ary,poly_num);
		mode=1;
    }
	
    /* mode ==2: the medial axis is involved */
    if(!poly_ary)
	{
		poly_ary = decompose(vert_ary, group_ary,poly_num, Mid_zvalue);
		mode=2;
    }
	
    if(!poly_ary)return 0; /* cannot be processed */ 
	
    draw_convex_triangles(poly_ary);
    if(MESH_FILE && mode <2)append_triangle_ary(poly_ary,0);
	
    if(!MESH_FILE ||  mode <2)return 1;
    /* need to make sure that the merging contour has no two adjacent v. */
	/*
    int res=verify_adjacent_vertices(vert_ary, group_ary,poly_num, isCCW);
    if(!res){
	return 0;
    }
    if(res==2){ 
	append_triangle_ary(poly_ary, 1);
	free_polygons_ary(poly_ary);
	return 1;
    }
	
	 ret_ary_ind = find_median_axis( poly_ary, isCCW, &ret_ary);
	 convert_middle_line(ret_ary, ret_ary_ind, ret_ary[0][2]); 
	 append_median_axis(ret_ary,ret_ary_ind, isCCW);
	 free(ret_ary);
	*/
	
	/*    convert_middle_tri_ary(poly_ary, ret_ary[0][2]); */
    append_triangle_ary(poly_ary,1); /* normal */
    append_Steiner_lines(poly_ary,isCCW);
    free_polygons_ary(poly_ary);
    return 1;
}


