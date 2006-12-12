#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <memory.h>

#include "common.h"
#include "dec_type.h"
#include "math_util.h"
#include "decom_util.h"

extern int STATIC_CHQ;

void my_clear_decom_util()
{
	check_3_sharp_triangle_angle(NULL, NULL, NULL);
}

//KLC
/*
void free_polygons_ary(PolygonStruct **polygon_ary)
{
int i;
PolygonStruct *poly;
if(polygon_ary==NULL)return;
i=0;
while(polygon_ary[i]){
poly=polygon_ary[i];
free(poly->vertices);
poly->vertices=NULL;
if(poly->normals){
free(poly->normals);
poly->normals=NULL;
}
free(polygon_ary[i]);
polygon_ary[i]=NULL;
i++;
}
free(polygon_ary);
}
*/

void calc_vector_from_cir(VertType *ret_v, VertType *vert_ary, 
						  short *vert_ind_ary, int poly_num, int index)
{
    int j,ind0,ind1;
    VertType *p0,*p1;

    ind0=vert_ind_ary[index];
    j=index+1;
    if(j>=poly_num)j=0;
    ind1=vert_ind_ary[j];
    p0=&vert_ary[ind0];
    p1=&vert_ary[ind1];
    (*ret_v)[0]=(*p1)[0]-(*p0)[0];
    (*ret_v)[1]=(*p1)[1]-(*p0)[1];
}


int is_convex(VertType *vert_ary,short *vert_ind_ary,int poly_num)
// the vert_ind_ary[] should have a size of at least poly_num+1 
{
	int i,sign,pre_sign=0;
	VertType v1,v0;
	double res;
	
	vert_ind_ary[poly_num]=vert_ind_ary[0];
	calc_vector_from_cir(&v0,vert_ary,vert_ind_ary,poly_num,poly_num-1);
	
	for(i=0; i<poly_num; i++) 
	{
        calc_vector_from_cir(&v1,vert_ary,vert_ind_ary,poly_num,i);
        res=outer_product_2d(v0,v1);
        if(res > 0.1)	sign=1;
        else if(res < -0.1)	sign=2;
        else	sign=0; // netrual
        if(pre_sign==0 && sign)	pre_sign=sign;
        if(sign && pre_sign !=sign)	return 0;

        memcpy(v0,v1,sizeof(v0));
	}
	return 1;
}


void put_a_poly_into_poly_structure(PolygonStruct** poly_ary,
									VertType *vert_ary, int *poly_ary_index, int size,
									short * vert_ind_ary, int poly_num)
{
    PolygonStruct *poly;
    int i;
	
    if(poly_ary==NULL)
	{
        fprintf(stderr,"Internal Error, poly_ary=NULL\n");
        exit(1);
    }
    
	if(*poly_ary_index>=size)
	{
		fprintf(stderr,"Mem. overflow put_a_poly_into..(), enlarge MAX_TRI_NUM\n");
		
        exit(1);
    }
    
	poly=(PolygonStruct *)mymalloc(sizeof(*poly));
    poly_ary[(*poly_ary_index)++]=poly;
    poly->numpts=poly_num;
    poly->normals=NULL;
    poly->vertices=(VertType *)mymalloc(sizeof(VertType)*poly_num);
    
	for(i=0; i<poly_num; i++)
	{
        int ind;
        ind=vert_ind_ary[i];
        memcpy(&(poly->vertices[i]), &(vert_ary[ind]), sizeof(VertType));
    }
}

double dot_product_2d(VertType *v1, VertType *v2)
{
    return((*v1)[0]*(*v2)[0]+(*v1)[1]*(*v2)[1]); 
}
/*
int normalize_vector(VertType *v1)
{
double dist,x,y;
x=(*v1)[0];
y=(*v1)[1];
dist=sqrt(x*x+y*y);
if(dist<0.0001){
fprintf(stderr,"Warning! normalize_vector(), dist=%lf\n",dist);
return 0;
}
(*v1)[0]=x/dist;
(*v1)[1]=y/dist;
return 1;
}
*/
int  calc_norm_vector_from_2vp(VertType *ret_v, VertType *p0, VertType *p1)
{
    int res;
    (*ret_v)[0]=(*p1)[0]-(*p0)[0];
    (*ret_v)[1]=(*p1)[1]-(*p0)[1];
    res=normalize_vector(ret_v);
    
	return res;
}


int is_center_line_not_cross_poly(VertType *vt, 
								  VertType *vert_ary, short *vert_ind_ary, int poly_num)
{
    int i,j,j1,res;
    VertType *vt1;
    
	for(i=0; i<poly_num; i++)
	{
        vt1=&vert_ary[i];
        for(j=0; j<poly_num; j++)
		{
            j1=j+1;
            if(j1>=poly_num)j1=0;
            if(i!=j && i!=j1)
			{
                res=find_intersection(vt,vt1,
					&vert_ary[vert_ind_ary[j]],
					&vert_ary[vert_ind_ary[j1]],NULL,NULL,NULL);
                if(res)return 0;
            }
        }
    }
    return 1;
}


int check_3_sharp_triangle_angle(VertType *vt1, VertType *vt2, 
								 VertType *vt3)
{
    VertType v[3];
    int res,i,j;
    static double thres=0.984; /* 10 degree */
    static double thres1=-0.984; /* 10 degree */
    double dval;
    
	res=calc_norm_vector_from_2vp(&v[0], vt1, vt2);
    
	if(!res)return 1; /* fail */
    
	res=calc_norm_vector_from_2vp(&v[1], vt2, vt3);
    
	if(!res)return 1; /* fail */
    
	res=calc_norm_vector_from_2vp(&v[2], vt3, vt1);
    
	if(!res)return 1; /* fail */
    
	for(i=0; i<3; i++)
	{
		j=i+1; 
		if(j>=3)j=0;
		dval=dot_product_2d(&v[i],&v[j]);
		if(dval>thres || dval < thres1)
			return 1;
    }
    return 0;
}

int is_angle_not_sharp(VertType *vert_ary,
					   short *vert_ind_ary,int poly_num, VertType *vt)
{
    int i,j,res;
    
	for(i=0; i<poly_num; i++)
	{
        j=i+1;
        if(j>=poly_num)	j=0;
        res=check_3_sharp_triangle_angle(&vert_ary[vert_ind_ary[i]],
            &vert_ary[vert_ind_ary[j]], vt);
        if(res)return 0;
    }
    return 1;
}


void find_polygon_center(VertType *vert_ary, 
						 short *vert_ind_ary,int poly_num, VertType *ret_vt)
						 /* the vert_ind_ary[] should have a size of at least poly_num+1 */
{
    double x=0.0, y=0.0;
    int i;
	
    for(i=0; i<poly_num; i++) 
	{
        x+= vert_ary[vert_ind_ary[i]][0];
        y+= vert_ary[vert_ind_ary[i]][1];
    }
    
	x=x/poly_num;
    y=y/poly_num;
    (*ret_vt)[0]=x;
    (*ret_vt)[1]=y;
}


int find_quad_center(VertType *vert_ary, short *vert_ind_ary,
					 int poly_num, VertType *ret_vt)
{
    int res;
    double dist0=1e14,dist1=1e14;
    VertType *p0,*p1,*p2,*p3,pt1,pt2;
	
    if(poly_num!=4)
	{
		printf("Error! find_quad_center(poly_num=%d)\n",poly_num);
		return 0;
    }
    
	p0=&(vert_ary[vert_ind_ary[0]]);
    p1=&(vert_ary[vert_ind_ary[1]]);
    p2=&(vert_ary[vert_ind_ary[2]]);
    p3=&(vert_ary[vert_ind_ary[3]]);
	
    pt1[0]=((*p0)[0] + (*p2)[0])/2.0;
    pt1[1]=((*p0)[1] + (*p2)[1])/2.0;
    pt2[0]=((*p1)[0] + (*p3)[0])/2.0;
    pt2[1]=((*p1)[1] + (*p3)[1])/2.0;
    
	res=is_inside_contour(&pt1,vert_ary, vert_ind_ary, 4, NULL,NULL);
    if(res)dist0=line_distance_2d(p0,p2);
    
	res=is_inside_contour(&pt2,vert_ary, vert_ind_ary, 4, NULL,NULL);
    
	if(res)dist1=line_distance_2d(p1,p3);
    
	if(dist0 <1e10 && dist1<1e10)
	{
		(*ret_vt)[0]=(pt1[0]+pt2[0])/2.0;
		(*ret_vt)[1]=(pt1[1]+pt2[1])/2.0;
		return 1;
    }
    
	if(dist0>dist1)
	{
		(*ret_vt)[0]=pt2[0];
		(*ret_vt)[1]=pt2[1];
		return 1;
    }
    
	if(dist0<dist1)
	{
		(*ret_vt)[0]=pt1[0];
		(*ret_vt)[1]=pt1[1];
		return 1;
    }
    else return 0;
}

int is_line_cross_polygon(int ind0, int ind1,VertType *vert_ary, 
						  short *vert_ind_ary, int poly_num)
						  /* Is this line [ind0, ind1] crosses any other line segment? */
{
    int i,j,res;
    VertType *p0,*p1;
    p0=&vert_ary[vert_ind_ary[ind0]];
    p1=&vert_ary[vert_ind_ary[ind1]];
    
	for(i=0; i<poly_num; i++)
	{
		j=i+1;
		if(j>=poly_num)j=0;
		if(ind0!=i && ind0!=j && ind1!=i && ind1!=j)
		{
			/* not neighboring */
			res=find_intersection(p0,p1,&vert_ary[vert_ind_ary[i]],
				&vert_ary[vert_ind_ary[j]],NULL,NULL,NULL);
			if(res)return 1;
		}
    }
    return 0;
}

double calc_edge_distance_2d_of_a_poly(VertType *vert_ary, short *vert_ind_ary, 
									   int poly_num, double *dist_ary)
{
    int i,j;
    double dist;

    for(i=0; i<poly_num; i++)
	{
        j=i+1;
        if(j>=poly_num)j=0;
        dist = line_distance_2d(&vert_ary[vert_ind_ary[i]],
			&vert_ary[vert_ind_ary[j]]);
        dist_ary[i]=sqrt(dist);
    }
    
	dist=0.0;
    for(i=0; i<poly_num; i++)dist+=dist_ary[i];
    return dist;
}
void swap_int(int *a, int *b)
{
    int temp;
    temp=*b;
    *b=*a;
    *a=temp;
}


