#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "common.h"
#include "math_util.h"
#include "scale.h"

static double Scale,ShiftX,ShiftY;
static double CenterX, CenterY;
static double BaseZ, DiffZ;
extern int STATIC_CHQ;

void my_clear_scale()
{
	Scale = ShiftX = ShiftY= 0;
	CenterX = CenterY=0;
	BaseZ = DiffZ = 0;
}

int scale_contour( VertType *lines_ary, short ** new_group_ary,
				  short *new_group_ind_ary, int group_num)
{
    int i,j,k,maxi=0;

    for(j=0; j<group_num; j++)
	{
		for(i=0; i<new_group_ind_ary[j]; i++)
		{
			k=new_group_ary[j][i];
			if(k>maxi)maxi=k;
		}
    }

    for(i=0; i<=maxi; i++)
	{
		lines_ary[i][0]=(lines_ary[i][0]-CenterX)*Scale+ShiftX;
		lines_ary[i][1]=(lines_ary[i][1]-CenterY)*Scale+ShiftY;
    }
    return 1;
}

int scale_back_triangle(VertType *p)
{
    double z_ratio;
    double x,y;
    
	if((*p)[2]>=BaseZ)return 0; 
    z_ratio=(BaseZ-(*p)[2])/DiffZ;
    x=((*p)[0]-ShiftX)/Scale+CenterX;
    y=((*p)[1]-ShiftY)/Scale+CenterY;
    
	if(z_ratio<0.999)
	{
		x=(x-(*p)[0])*z_ratio+(*p)[0];
		y=(y-(*p)[1])*z_ratio+(*p)[1];
    }
    
	(*p)[0]=x;
    (*p)[1]=y;
    return 1;
}




int find_bounding_box(double box[2][2], VertType *lines_ary,short ** new_group_ary, 
					  short *new_group_ind_ary, int group_num)
{
    double maxx=-1000000.0, minx=1000000.0;
    double maxy=-1000000.0, miny=1000000.0;
    int i,j,k;
    VertType *vp;

    for(j=0; j<group_num; j++)
	{
		for(i=0; i<new_group_ind_ary[j]; i++)
		{
			k=new_group_ary[j][i];
			vp=&lines_ary[k];
			if((*vp)[0]>maxx)maxx=(*vp)[0];
			if((*vp)[0]<minx)minx=(*vp)[0];
			if((*vp)[1]>maxy)maxy=(*vp)[1];
			if((*vp)[1]<miny)miny=(*vp)[1];
		}
    }

    box[0][0]=minx;
    box[0][1]=maxx;
    box[1][0]=miny;
    box[1][1]=maxy;
    
	printf("bounding box [%lf %lf], [%lf %lf]\n",box[0][0],box[0][1],
		box[1][0],box[1][1]);
    
	return 1;
}

int calculate_factor(double box1[2][2], double box2[2][2], double x_ratio, 
					 double y_ratio)
{
    FILE *ofp;
    double lx1,ly1,cx1,cy1;
    double lx2,ly2,cx2,cy2;
    double ratio1,ratio2;

    lx1=box1[0][1]-box1[0][0];
    ly1=box1[1][1]-box1[1][0];
    cx1=(box1[0][1]+box1[0][0])/2.0;
    cy1=(box1[1][1]+box1[1][0])/2.0;
    CenterX=cx1;
    CenterY=cy1;
    lx2=box2[0][1]-box2[0][0];
    ly2=box2[1][1]-box2[1][0];
    cx2=(box2[0][1]+box2[0][0])/2.0;
    cy2=(box2[1][1]+box2[1][0])/2.0;
    ratio1=lx2/lx1;
    ratio2=ly2/ly1;
    Scale=ratio1;
    
	if(Scale>ratio2)Scale=ratio2;
    
	ShiftX=cx2+ lx2*x_ratio;
    ShiftY=cy2+ ly2*y_ratio;
    printf("ALL_CORR, scale=%lf, shift=[%lf %lf]\n",Scale,ShiftX,ShiftY);
    printf("ALL_CORR, BaseZ=%lf, DiffZ=%lf\n",BaseZ,DiffZ);
	
    if((ofp=fopen("scale.coeff","w"))==NULL) 
	{
        fprintf(stderr,"can not open scale.coeff to store scale info\n");
        exit(1);
    }
    
	fprintf(ofp,"%lf %lf %lf\n", Scale,ShiftX,ShiftY);
    fprintf(ofp,"%lf %lf\n", CenterX, CenterY);
    fprintf(ofp,"%lf %lf\n", BaseZ, DiffZ);
    fclose(ofp);
    
	return 1;
}

int scale_shift( double xshift, double yshift,
				VertType *lines_ary0, VertType *lines_ary1,
				short ** new_group_ary0, short *new_group_ind_ary0, int group_num0,
				short ** new_group_ary1, short *new_group_ind_ary1, int group_num1)
{
    double box1[2][2], box2[2][2];

    DiffZ=lines_ary0[0][2]-BaseZ;
    BaseZ=lines_ary1[0][2];
    
	find_bounding_box(box1, lines_ary0, new_group_ary0, new_group_ind_ary0, group_num0);
    
	find_bounding_box(box2, lines_ary1, new_group_ary1, new_group_ind_ary1, group_num1);
    
	calculate_factor(box1,box2,xshift,yshift);
    
	scale_contour(lines_ary0, new_group_ary0, new_group_ind_ary0, group_num0);
    
	return 1;
}

