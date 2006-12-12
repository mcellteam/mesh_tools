#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <memory.h> 

#include "common.h" 
#include "dec_type.h" 
#include "math_util.h" 
#include "decom_util.h"
#include "decompose.h"


#define MIN_DIST 0.02 
#define MIN_VECT 0.00005 
#define MAX_TRI_NUMBER 100000 
#define MAX_POLY_NUM 20 
static double Mid_zvalue; 
static VertType *Vert_ary_vor; 
static VertType *Vert_ary_extra; 
static int Vert_ary_extra_size;
static int Vert_ary_extra_index;
static PolygonStruct **Triangle_ary; 
static int Triangle_ary_index; 
static short Crossed_pair[2][100]; 
static int Crossed_pair_ind; 
static short *Poly_index_ary[MAX_POLY_NUM]; 
static short Poly_size_ary[MAX_POLY_NUM]; 
static int Num_of_poly; 
static int Security_cnt, Security_limit;

extern int STATIC_CHQ;

void my_clear_decompose()
{
	Mid_zvalue =0;
	Vert_ary_vor = NULL;
	Vert_ary_extra = NULL;
	Vert_ary_extra_size = 0;
	Vert_ary_extra_index =0;
	Triangle_ary = NULL;
	Triangle_ary_index =0;
	Crossed_pair_ind=0;	 //no need to init the Crossed_pair arr...
	//Poly_index_ary = NULL; // no need if Num_of_poly =0;
	Num_of_poly =0;
	Security_cnt = Security_limit=0;
}


void check_polygon_self_crossed(
								VertType *vert_ary, short *vert_ind_ary, int poly_num)
{
    int i,j,i1,j1,res;
    double ratio1,ratio2;
    VertType pt;

    for(i=2; i<poly_num; i++)
	{
        i1=(i+1)%poly_num;
        for(j=0; j<i-1; j++)
		{
            j1=(j+1)%poly_num;
            res= find_intersection( &vert_ary[vert_ind_ary[i]],
                &vert_ary[vert_ind_ary[i1]],&vert_ary[vert_ind_ary[j]],
				&vert_ary[vert_ind_ary[j1]],&pt,&ratio1,&ratio2);
        
			if(res==1)
			{
				if(!(ratio1>0.9 || ratio1<0.1 ||ratio2>0.9 || ratio2<0.1))
					printf("Warning check_polygon_self_crossed(), ratio1=%lf, ratio2=%lf\n",
                    ratio1,ratio2);
            }
			
			if(res==1 || res==5)
			{
				Crossed_pair[0][Crossed_pair_ind]=i;
				Crossed_pair[1][Crossed_pair_ind]=j;
				Crossed_pair_ind++;
			}
        }
    }
}


void add_a_triangle_to_array(PolygonStruct *poly)
{
    if(Triangle_ary_index>=MAX_TRI_NUMBER)
	{
		fprintf(stderr,"Mem. overflow breaking_poly..(), enlarge MAX_TRI_NUMBER\n");
        exit(1);
    }
    
	Triangle_ary[Triangle_ary_index++]=poly;
}

PolygonStruct * create_a_poly_from_points(VertType *p_ary[], int num)
{
    PolygonStruct *poly;
    int i;
	
    if(num<3)return NULL;
    
	poly=(PolygonStruct *)mymalloc(sizeof(*poly));
    poly->numpts=num;
    poly->normals=NULL;
    poly->vertices=(VertType *)mymalloc(sizeof(VertType)*num);
    
	for(i=0; i<num; i++) 
		memcpy(&(poly->vertices[i]), p_ary[i], sizeof(VertType));
    return poly;
}

void breaking_poly_by_center(short *vert_ind_ary, int poly_num, VertType *vt)
/* return the number of polygons */
{
    int i,j;
    PolygonStruct *poly;
    VertType *p_ary[3];
    
	for(i=0; i<poly_num; i++)
	{
        int ind0,ind1;
		
        ind0=vert_ind_ary[i];
		j=i+1;
		if(j>=poly_num)j=0;
        ind1=vert_ind_ary[j];
		if(ind0>=10000)	p_ary[0]=&(Vert_ary_extra[ind0-10000]);
		else	 p_ary[0]=&(Vert_ary_vor[ind0]);
		if(ind1>=10000)	p_ary[1]=&(Vert_ary_extra[ind1-10000]);
		else	p_ary[1]=&(Vert_ary_vor[ind1]);
		p_ary[2]=vt;
		
		poly=create_a_poly_from_points(p_ary,3);
		
		add_a_triangle_to_array(poly);		
    }
}

int find_shortest_crossed(short *vert_ind_ary, int poly_num, 
						  int *min_low_ind, int *min_high_ind)
{
    int i,j,k0,k1;
    int ind0[2],ind1[2];
    int poly_num1=poly_num-1;
    double min_dist=1e12,dist;
    int min0,min1;
    VertType *p1, *p0;
	
    for(i=0; i<Crossed_pair_ind; i++)
	{
		int flag=0;
		ind0[0]=Crossed_pair[0][i];
		ind1[0]=Crossed_pair[1][i];
		ind0[1]= ind0[0]+1;
		ind1[1]= ind1[0]+1;
	
		if(ind0[1]>=poly_num)ind0[1]=0;
		if(ind1[1]>=poly_num)ind1[1]=0;
		
		for(j=0; j<4; j++)
		{
			int temp;
			/* try all 4 lines */
			k0=ind0[j&1];
			k1=ind1[(j>>1)&1];
			if(k1>k0)temp=k1-k0;
			else temp=k0-k1;
		
			if(temp!=1 && temp!=poly_num1)
			{ /* can not be neighboring */
				p0=&Vert_ary_vor[vert_ind_ary[k0]];
				p1=&Vert_ary_vor[vert_ind_ary[k1]];
				if(fabs((*p0)[2]-(*p1)[2])>0.001)
				{
					flag=1;
					/* different z plane */
					dist=line_distance_2d(p0,p1);
			
					if(dist<min_dist)
					{
						min_dist=dist;
						min0=k0;
						min1=k1;
					}
				}
			}
		}
		if(!flag)return 0;
    }
    if(min0>min1)swap_int(&min0,&min1);
    *min_low_ind=min0;
    *min_high_ind=min1;
    if(min_dist<1e10)return 1;
    else return 0;
}

void calc_middle_point(VertType *pt, VertType *p1, VertType *p2)
{
	(*pt)[0]=((*p1)[0]+(*p2)[0])/2.0;
	(*pt)[1]=((*p1)[1]+(*p2)[1])/2.0;
	(*pt)[2]=Mid_zvalue;
}

int find_any_best_diagnoal(short *vert_ind_ary, int poly_num, 
						   int *min_low_ind, int *min_high_ind)
						   /* the diagonal must in different slice */
{
    double dist,min_dist=1e14;
    int i,j;
	
    double cir_dist, cir_dist2, dist_ary[MAX_TRI_NUMBER];
    cir_dist=calc_edge_distance_2d_of_a_poly(Vert_ary_vor, vert_ind_ary, 
		poly_num, dist_ary);
    cir_dist2=cir_dist/2.0;

    for(i=2; i<poly_num; i++)
	{
		for(j=i-2; j>=0; j--)
		{
			if(vert_ind_ary[i]<10000 && vert_ind_ary[j]<10000)
			{
				/* none is extra vertex */
				VertType *p1,*p2;
				int kk;
				double part_dist=0;
				p1=&(Vert_ary_vor[vert_ind_ary[j]]);
				p2=&(Vert_ary_vor[vert_ind_ary[i]]);
				dist=line_distance_2d(p1,p2);
				dist=sqrt(dist);
				/* put some weight to balanced partition */
				for(kk=j; kk<i; kk++)part_dist+=dist_ary[kk]; 
				if(part_dist>cir_dist2)part_dist=cir_dist - part_dist;
				dist=dist/part_dist;
			
				if(min_dist>dist)
				{
					/* is this line cross any other line ? */
					int res;
					res=is_line_cross_polygon(j,i,Vert_ary_vor, 
						vert_ind_ary, poly_num);
					if(!res)
					{
						/* the middle point should be inside polygon */ 
						VertType pt;
						calc_middle_point(&pt, p1, p2);
						/* check if the middle point is inside polygon? */
						if(is_inside_contour(&pt, Vert_ary_vor, 
							vert_ind_ary, (short) poly_num,NULL,NULL))  
						{
							min_dist=dist;
							*min_low_ind=j;
							*min_high_ind=i;
						}
					}
				}
			}
		}
    }
    if(min_dist<1e12)return 1;
    else return 0;
}

int find_best_diagnoal(short *vert_ind_ary, int poly_num, 
					   int *min_low_ind, int *min_high_ind)
					   /* the diagonal must in different slice */
{
    short *high_ary, *low_ary; /* store vertices of different plane */
    int high_index=0, low_index=0;
    double dist,min_dist=1e14;
    int i,j,k;
	
    double cir_dist, cir_dist2, dist_ary[MAX_TRI_NUMBER];
    cir_dist=calc_edge_distance_2d_of_a_poly(Vert_ary_vor, vert_ind_ary, 
		poly_num, dist_ary);
    cir_dist2=cir_dist/2.0;
    j=0;
    k=0;

    for(i=0; i<poly_num; i++)
	{
		if(Vert_ary_vor[vert_ind_ary[i]][2]>Mid_zvalue)j++;
		else k++;
    }
    
	high_ary=(short *)mymalloc(sizeof(short *)*(j+1));
    low_ary=(short *)mymalloc(sizeof(short *)*(k+1));
    
	for(i=0; i<poly_num; i++)
	{
		if(Vert_ary_vor[vert_ind_ary[i]][2]>Mid_zvalue)
			high_ary[high_index++]=i;
		else low_ary[low_index++]=i;
    }
    
	for(i=0; i<high_index; i++)
	{
		for(j=0; j<low_index; j++)
		{
			int ind0,ind1;
			ind0=high_ary[i];
			ind1=low_ary[j];
			if(ind0>ind1)swap_int(&ind0,&ind1);
			k=ind1-ind0;
			
			if(k!=1 && k!=(poly_num-1))
			{ /* not neighboring */
				VertType *p1,*p2;
				int kk;
				double part_dist=0;
				p1=&(Vert_ary_vor[vert_ind_ary[ind0]]);
				p2=&(Vert_ary_vor[vert_ind_ary[ind1]]);
				dist=line_distance_2d(p1,p2);
				/* put some weight to balanced partition */
				for(kk=ind0; kk<ind1; kk++)part_dist+=dist_ary[kk]; 
				if(part_dist>cir_dist2)part_dist=cir_dist - part_dist;
				dist=dist/part_dist;
				if(min_dist>dist)
				{
					/* is this line cross any other line ? */
					int res;
					res=is_line_cross_polygon(ind0,ind1,Vert_ary_vor, 
						vert_ind_ary, poly_num);
					if(!res)
					{
						/* the middle point should be inside polygon */ 
						VertType pt;
						calc_middle_point(&pt, p1, p2);
						/* check if the middle point is inside polygon? */
						if(is_inside_contour(&pt, Vert_ary_vor, 
							vert_ind_ary, (short) poly_num,NULL,NULL))  
						{
							min_dist=dist;
							*min_low_ind=ind0;
							*min_high_ind=ind1;
						}
					}
				}
			}
		}
    }
    free(high_ary);
    free(low_ary);
    if(min_dist<1e12)return 1;
    else return 0;
}

void decompose_sub(short *vert_ind_ary, int poly_num) 
/* It recursively breaks the polygons into triangles or samller 
polygons. The breaking line must between two vertices which has 
different Z values */
{
    int i,j,k,res;
    int min_high_ind,min_low_ind;
    
    if(poly_num==3)
	{ /* It is a basic triangle */
		put_a_poly_into_poly_structure(Triangle_ary, Vert_ary_vor,
			&Triangle_ary_index,
			MAX_TRI_NUMBER, vert_ind_ary, poly_num);
		return;
    }
	
    j=0;
    k=0;
    Crossed_pair_ind=0;
    check_polygon_self_crossed(Vert_ary_vor,vert_ind_ary,poly_num);
    
	if(Crossed_pair_ind>1)
	{
		fprintf(stderr,"Warning, crossed self %d times\n",Crossed_pair_ind);
		printf("Warning, crossed self %d times\n",Crossed_pair_ind);
    }
    
	if(Crossed_pair_ind>0)
	{ /* processed the crossed line seg first */
		res=find_shortest_crossed(vert_ind_ary, poly_num, 
			&min_low_ind, &min_high_ind);
    }
    else res=0;
	
    if(res)
	{ /* It can be broke down to two polygons */
		int num1,num2;
		short *vert_ind_ary1,*vert_ind_ary2;
		num1=min_high_ind-min_low_ind+1;
		num2=min_low_ind+poly_num-min_high_ind+1;
		vert_ind_ary1=(short*)mymalloc(sizeof(short)*(num1+1));
		for(i=0; i<num1; i++) vert_ind_ary1[i]=vert_ind_ary[min_low_ind+i];
		decompose_sub(vert_ind_ary1, num1); 
		free(vert_ind_ary1);
		vert_ind_ary2=(short*)mymalloc(sizeof(short)*(num2+1));
		for(i=0; i<num2; i++) 
			vert_ind_ary2[i]=vert_ind_ary[(min_high_ind+i)%poly_num];
		decompose_sub(vert_ind_ary2, num2); 
		free(vert_ind_ary2);
    }
    else 
	{ /* can not break furthermore, store it for return */
		if(poly_num==4)
		{ /* It is four sided */
			VertType vt;
			int res;
			res=find_quad_center(Vert_ary_vor, vert_ind_ary, poly_num, &vt); 
			vt[2] = Mid_zvalue;
			if(res)
			{
				breaking_poly_by_center(vert_ind_ary, poly_num, &vt);
				return;
			}
		}
		
		Poly_size_ary[Num_of_poly]=poly_num;
		Poly_index_ary[Num_of_poly]=
			(short *)mymalloc(sizeof(short)*(poly_num+1));
		memcpy(Poly_index_ary[Num_of_poly],
			vert_ind_ary,poly_num*sizeof(*vert_ind_ary));
		Num_of_poly++;
		return;
    }
}

int define_a_triangle(VertType *p1_ary[],int ind0, int ind1, int ind2)
{
    PolygonStruct *poly;
    VertType *p_ary[3];
    p_ary[0]=p1_ary[ind0];
    p_ary[1]=p1_ary[ind1];
    p_ary[2]=p1_ary[ind2];
    poly=create_a_poly_from_points(p_ary, 3);
    add_a_triangle_to_array(poly);
    return 1;
}

void breaking_convex_quad_into_two_tri(VertType *v0, VertType *v1,
									   VertType *v2, VertType *v3)
{
    VertType *p_ary[4];
    p_ary[0]=v0;
    p_ary[1]=v1;
    p_ary[2]=v2;
    p_ary[3]=v3;
	
    if(line_distance_2d(v0,v2)<line_distance_2d(v1,v3))
	{
		define_a_triangle(p_ary,0,1,2);
		define_a_triangle(p_ary,2,3,0);
    }
    else
	{
		define_a_triangle(p_ary,0,1,3);
		define_a_triangle(p_ary,1,2,3);
    }
}

int special_breaking_a_quad(short *vert_ind_ary, int poly_num)
{
    int i,j;
    VertType *p1_ary[6];
	
    for(i=0; i<poly_num; i++)
	{
		j=vert_ind_ary[i];
		if(j<10000) p1_ary[i]=&Vert_ary_vor[j];
		else p1_ary[i]=&Vert_ary_extra[j-10000];
    }
    
	breaking_convex_quad_into_two_tri(p1_ary[0],p1_ary[1],p1_ary[2],p1_ary[3]);
    breaking_convex_quad_into_two_tri(p1_ary[0],p1_ary[3],p1_ary[4],p1_ary[5]);
    return 1;
}

int special_breaking_a_triangle(short *vert_ind_ary, int poly_num)
{
    int i,j;
    VertType *p1_ary[6];
	
    for(i=0; i<poly_num; i++)
	{
		j=vert_ind_ary[i];
		if(j<10000) p1_ary[i]=&Vert_ary_vor[j];
		else p1_ary[i]=&Vert_ary_extra[j-10000];
    }
    
	if(poly_num==4)
	{
		define_a_triangle(p1_ary,0,1,2);
		define_a_triangle(p1_ary,2,3,0);
		return 1;
    }
    
	if(poly_num==5)
	{
		define_a_triangle(p1_ary,0,1,2);
		breaking_convex_quad_into_two_tri(p1_ary[0],p1_ary[2],p1_ary[3],p1_ary[4]);
    }
    
	if(poly_num==6)
	{
		define_a_triangle(p1_ary,0,1,5);
		define_a_triangle(p1_ary,1,2,3);
		define_a_triangle(p1_ary,3,4,5);
		define_a_triangle(p1_ary,1,3,5);
		return 1;
    }
    return 1;
}

int rearrange_quad_points(short *vert_ind_ary, int poly_num)
/* re-arrange the point sequence 
poly_num1==6, then [0], [3] are the middle plane points */
{
    int temp,i;
    short ind_ary[6];
    if(poly_num!=6)return 0;
    temp=0;
    
	while(temp<poly_num && (vert_ind_ary[temp]<10000 || 
		vert_ind_ary[(temp+3)%poly_num]<10000))temp++;
    if(temp>=poly_num)
	{
		return 0;
    }
    
	for(i=0; i<poly_num; i++)
	{
		ind_ary[i]=vert_ind_ary[(temp+i)%poly_num];
    }
    
	for(i=0; i<poly_num; i++) vert_ind_ary[i]=ind_ary[i];
    return 1;
}


int rearrange_triangle_points(short *vert_ind_ary, int poly_num)
/* re-arrange the point sequence 
if poly_num1==4, then vert_ind_ary[0] is the middle plane point
if poly_num1==5, then [0] and [2] are the middle plane points
if poly_num1==6, then [0], [2] and [4] are the middle plane points */
{
    int temp,i;
    short ind_ary[6];
    
	if(poly_num==6)
	{
		if(vert_ind_ary[0]<10000 || vert_ind_ary[2]<10000 
			|| vert_ind_ary[4]<10000)
		{
			fprintf(stderr,"Internal Error in rearrange_triangle 6...\n");
			return 0;
		}
		return 1;
    }

    if(poly_num==4)
	{
		temp=0;
		while(vert_ind_ary[temp]<10000)temp++;
    }
    else if(poly_num==5)
	{
		temp=0;
		while(temp<poly_num && (vert_ind_ary[temp]<10000 || 
			vert_ind_ary[(temp+2)%poly_num]<10000))temp++;
		if(temp>=poly_num)
		{
			fprintf(stderr,"Internal error, rearrange_tri, 5\n");
			return 0;
		}
    }

    for(i=0; i<poly_num; i++)
	{
		ind_ary[i]=vert_ind_ary[(temp+i)%poly_num];
    }
    
	for(i=0; i<poly_num; i++) vert_ind_ary[i]=ind_ary[i];
    return 1;
}

int special_case_triangle(short *vert_ind_ary, int poly_num)
{
    int res;
    rearrange_triangle_points(vert_ind_ary,poly_num);
    res=special_breaking_a_triangle(vert_ind_ary, poly_num);
    return res;
}

int special_case_quad(short *vert_ind_ary, int poly_num)
{
    int res;
    res=rearrange_quad_points(vert_ind_ary,poly_num);
    if(!res)return 0;
    res=special_breaking_a_quad(vert_ind_ary, poly_num);
    return res;
}

int triangulate_a_poly_by_center(short *vert_ind_ary1, int poly_num1 )
{
    int i,j,res;
    VertType ret_pt;
    int poly_num=0;
    short vert_ind_ary[MAX_TRI_NUMBER];
	
	/*KLC: perhaps the vert_ind_ary is not getting initialised correctly...*/
	for (i=0; i<MAX_TRI_NUMBER; i++)
		vert_ind_ary[i]=0;
	/*KLC*/
    
	for(i=0; i<poly_num1; i++)
	{
		if(vert_ind_ary1[i]<10000)
			vert_ind_ary[poly_num++]=vert_ind_ary1[i];
    }
    
	if(poly_num==3 && poly_num1>3)
	{
		special_case_triangle(vert_ind_ary1, poly_num1);
		return 1;
    }
	
    if(poly_num==4 && poly_num1==6 && is_convex(Vert_ary_vor,vert_ind_ary,poly_num))
	{
		res=special_case_quad(vert_ind_ary1, poly_num1);
		if(res)return 1;
    }
	j=is_convex(Vert_ary_vor,vert_ind_ary,poly_num); 
	
	if(poly_num==4 && poly_num1==6 && is_convex(Vert_ary_vor,vert_ind_ary,poly_num))
	{
		res=special_case_quad(vert_ind_ary1, poly_num1);
		if(res)return 1;
    }
    j=is_convex(Vert_ary_vor,vert_ind_ary,poly_num);
	
    
	if(j)
	{ /* no reflex point, is a convex */
		find_polygon_center(Vert_ary_vor,vert_ind_ary,poly_num, &ret_pt);
		ret_pt[2]=Mid_zvalue;
		breaking_poly_by_center(vert_ind_ary1, poly_num1, &ret_pt);
    }
    else 
	{
		int res;
		find_polygon_center(Vert_ary_vor,vert_ind_ary, poly_num, &ret_pt);
		ret_pt[2]=Mid_zvalue;
		res=is_inside_contour(&ret_pt,Vert_ary_vor, 
			vert_ind_ary, (short) poly_num, NULL,NULL);
	
		if(!res)return 0;
		
		res=is_center_line_not_cross_poly(&ret_pt,Vert_ary_vor,
			vert_ind_ary, poly_num); 
		
		if(!res)return 0;
		
		res = is_angle_not_sharp(Vert_ary_vor, 
			vert_ind_ary, poly_num, &ret_pt);
		
		if(!res)return 0;
		
		breaking_poly_by_center(vert_ind_ary1, poly_num1, &ret_pt);
    }
    return 1;
}

int find_median_axis_sub(short *vert_ind_ary, int poly_num)
{
    int min_low_ind, min_high_ind;
    int res,i,ret_val=1;
    short *vert_ind_ary1, *mapping_ary;
    int poly_num1=0;
    int temp_index;
	
    if(Security_cnt>Security_limit)return 0; /* some problem */
	
    vert_ind_ary1=(short*)mymalloc(sizeof(short)*(poly_num+1));
    mapping_ary=(short*)mymalloc(sizeof(short)*(poly_num+1));
	
    for(i=0; i<poly_num; i++)
	{
		if(vert_ind_ary[i]<10000)
		{
			vert_ind_ary1[poly_num1]=vert_ind_ary[i];
			mapping_ary[poly_num1]=i;
			poly_num1++;
		}
    }
	
	if(Vert_ary_extra_index)
	{ // not the first time, try to triangulate 
		res=triangulate_a_poly_by_center(vert_ind_ary,  poly_num);
		if(res)
		{
			ret_val=1;
			goto end_label;
		}
    }
	
    res=find_any_best_diagnoal(vert_ind_ary1, poly_num1, 
		&min_low_ind, &min_high_ind);
    if(!res)
	{
		fprintf(stderr,"Err, fail in find_median_axis_sub()\n"); 
		ret_val=0;
		goto end_label;
    }
	
    if(res)
	{ /* It can be broken down to two polygons */
		int num1,num2;
		short *temp_ary;
		min_low_ind=mapping_ary[min_low_ind];
		min_high_ind=mapping_ary[min_high_ind];
		if(Vert_ary_extra_index>=Vert_ary_extra_size)
		{
			printf("Memory overwriting, index=%d, size=%d\n", 
				Vert_ary_extra_index, Vert_ary_extra_size); 
		}
		
		calc_middle_point(&Vert_ary_extra[Vert_ary_extra_index++], 
			&Vert_ary_vor[vert_ind_ary[min_low_ind]], 
			&Vert_ary_vor[vert_ind_ary[min_high_ind]]);  
		
		num1=min_high_ind-min_low_ind+1;
		num2=min_low_ind+poly_num-min_high_ind+1;
		temp_ary=(short*)mymalloc(sizeof(short)*(num1+2));
		
		for(i=0; i<num1; i++) temp_ary[i]=vert_ind_ary[min_low_ind+i];
		
		temp_index=Vert_ary_extra_index+9999; 
		temp_ary[num1]=temp_index;
		Security_cnt++;
		res=find_median_axis_sub(temp_ary, num1+1); 
		
		if(!res)
		{
			free(temp_ary);
			ret_val=0;
			goto end_label;
		}
		
		free(temp_ary);
		temp_ary=(short*)mymalloc(sizeof(short)*(num2+2));
		
		for(i=0; i<num2; i++) 
			temp_ary[i]=vert_ind_ary[(min_high_ind+i)%poly_num];
		temp_ary[num2]=temp_index;
		Security_cnt++;
		res=find_median_axis_sub(temp_ary, num2+1); 
		
		if(!res)
		{
			free(temp_ary);
			ret_val=0;
			goto end_label;
		}
		free(temp_ary);
    }
    ret_val=1;
	
end_label:
    free(vert_ind_ary1);
    free(mapping_ary);
    return ret_val;
}

int find_median_axis(short *vert_ind_ary, int poly_num) 
{
    int res;
    Vert_ary_extra_size=poly_num*2+2;
    Vert_ary_extra=(VertType*)mymalloc(sizeof(VertType)*Vert_ary_extra_size);
    Vert_ary_extra_index=0;
    
    Security_cnt=0;
    Security_limit=poly_num*2;
    res=find_median_axis_sub(vert_ind_ary, poly_num);
    free(Vert_ary_extra);
    Vert_ary_extra=NULL;
    return res;
}



PolygonStruct ** decompose( 
						   VertType *vert_ary, short *vert_ind_ary, 
						   int poly_num, double mid_z)
						   /* decompose the polygon into triangles and smaller polygons */
						   /* the vert_ind_ary[] should have a size of at least poly_num+1 */
{
	/*    PolygonStruct **ret_poly_ary; */
    int i,res;
	/*    VertType ret_pt; 
    short reflex_ary[MAX_TRI_NUMBER];
	*/
    int error=0;
	
    Num_of_poly=0;
    Mid_zvalue=mid_z;
    Vert_ary_vor=vert_ary;
    i=0;
	
    Triangle_ary=(PolygonStruct**)mymalloc((sizeof *Triangle_ary)*MAX_TRI_NUMBER);
    Triangle_ary_index=0; 
    
	if(poly_num>3)decompose_sub(vert_ind_ary, poly_num); 
    else 
	{
		Poly_size_ary[Num_of_poly]=poly_num;
		Poly_index_ary[Num_of_poly]=
			(short *)mymalloc(sizeof(short)*(poly_num+1));
		memcpy(Poly_index_ary[Num_of_poly],
			vert_ind_ary,poly_num*sizeof(*vert_ind_ary));
		Num_of_poly++;
    }
    
	i=0;
    for(i=0; (i<Num_of_poly) && error==0; i++)
	{
//		if(Poly_index_ary[i]) {
			res=triangulate_a_poly_by_center(Poly_index_ary[i], Poly_size_ary[i]);
			if(!res) 
			{ /* find its median axis */
				int res1;
				res1=find_median_axis(Poly_index_ary[i], Poly_size_ary[i]);
				if(!res1)error=1;
			}
//		}
    }

    for(i=0; i<Num_of_poly; i++)
	{
		free(Poly_index_ary[i]);
		Poly_index_ary[i]=NULL;
    }
    
	if(error)
	{
		free_polygons_ary(Triangle_ary);
		
		printf("---- Error! decompose(), Fail to process this polygon\n");
        printf("---- poly_num=%d, 1st pt=%lf %lf %lf\n",poly_num,
			vert_ary[0][0],vert_ary[0][1],vert_ary[0][2]);
		
		return NULL;
    }
    return Triangle_ary;
}

void free_polygons_ary(PolygonStruct **polygon_ary)
{
    int i;
    PolygonStruct *poly;
    if(polygon_ary==NULL)return;
    i=0;

    while(polygon_ary[i])
	{
        poly=polygon_ary[i];
		if(poly->vertices)
			free(poly->vertices);
		poly->vertices=NULL;
    
		if(poly->normals)
		{
			free(poly->normals);
			poly->normals=NULL;
		}
		
		free(polygon_ary[i]);
        polygon_ary[i]=NULL;
        i++;
    }
    free(polygon_ary);
}

