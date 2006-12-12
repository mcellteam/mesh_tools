#include <stdio.h>
// #include <io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
//#include <unistd.h> //KLC
#include <memory.h>

#include "string.h"
#include "common.h"
#include "math_util.h"
#ifndef PI
#define PI 3.14159
#endif

/* return the area enclosed by the contour */
double calculate_area(VertType *lines_ary, short *group_ary,
					  short group_size)
{
    double area=0;
    int i,j;
    VertType *p0, *p1;

    for(i=0; i<group_size; i++)
	{
        p0=&(lines_ary[group_ary[i]]);
        j=(i+1)%group_size;
        p1=&(lines_ary[group_ary[j]]);
        area+=((*p0)[0] * (*p1)[1] - (*p0)[1] * (*p1)[0]);
    }
    
	return area/2.0;
}


int normalize_vector(VertType *v1)
{
    double dist,x,y;
    x=(*v1)[0];
    y=(*v1)[1];
    dist=sqrt(x*x+y*y);
    
	if(dist<1e-6)
	{
		printf("Error! normalize_vector(), length^2=%f\n",dist);
		return 0;
    }
    
	(*v1)[0]=x/dist;
    (*v1)[1]=y/dist;
    
	return 1;
}

void bubble_sort(double *ratio_ary, int index)
{
    int i,j;
    double temp;
    
	for(i=1; i<index; i++)
	{
		for(j=0; j<i; j++)
		{
			if(ratio_ary[i]<ratio_ary[j])
			{
				/* swap */
				temp=ratio_ary[i];
				ratio_ary[i]=ratio_ary[j];
				ratio_ary[j]=temp;
			}
		}
    }
}

double line_distance(VertType *p1, VertType *p2)
{
    double dx,dy,dz,dist;
    dx=(*p2)[0]-(*p1)[0];
    dy=(*p2)[1]-(*p1)[1];
    dz=(*p2)[2]-(*p1)[2];
    
	dist=dx*dx+dy*dy+dz*dz;
    
	return dist;
}


double line_distance_2d(VertType *p1, VertType *p2)
{
    double dx,dy,dist;
    dx=(*p2)[0]-(*p1)[0];
    dy=(*p2)[1]-(*p1)[1];
    
	dist=dx*dx+dy*dy;
    
	return dist;
}

double outer_product_verttype(VertType *v1, VertType *v2)
{
    return (*v1)[0]*(*v2)[1]-(*v1)[1]*(*v2)[0];
}

double inner_product_verttype(VertType *v1, VertType *v2)
{
    return (*v1)[0]*(*v2)[0]+(*v1)[1]*(*v2)[1];
}

double outer_product_2d(double *v1, double *v2)
{
    return v1[0]*v2[1]-v1[1]*v2[0];
}

int calc_line_equation(VertType *pt0, VertType *pt1, double *a,
					   double *b, double *c)
					   /* f() >0 if it is on the left side */
{
    double dx,dy;
    dx=(*pt1)[0]-(*pt0)[0];
    dy=(*pt1)[1]-(*pt0)[1];

    if(fabs(dx)<=0.000001 && fabs(dy)<=0.000001)
	{
		fprintf(stderr,"I. Error! same pts in calc_line_equation(%lf, %lf)\n",
			(*pt1)[0],(*pt1)[1]);
		return 0;
    }
    *a=-dy;
    *b=dx;
    *c=(*pt0)[0]*dy - (*pt0)[1]*dx;
    
	return 1;
}


int is_inside_triangle(double *pt, double vertice[3][2])
/* This program check if the pt is inside the triangle of vertice[][] */
/* All four points need to be in the sampe plane */
/* return: 1: yes, 0: no */
{
    double v[3][2],res;
    int i,j,cnt=0;
    
	for(i=0; i<3; i++)
	{
        for(j=0; j<2; j++)v[i][j]=vertice[i][j]-pt[j];
    }
    
	res=outer_product_2d(v[0],v[1]);
    if(res<=-0.000001)cnt++; 
	else if(res>=0.000001)cnt--; 

    res=outer_product_2d(v[1],v[2]);

    if(res<=-0.000001)cnt++; 
    else if(res>=0.000001)cnt--; 
    
	res=outer_product_2d(v[2],v[0]);
    
	if(res<=-0.000001)cnt++; 
    else if(res>=0.000001)cnt--; 
	
    if(cnt==3 || cnt == -3)return 1;
    else return 0;
}


int is_pt_on_line(VertType *pt0,VertType *pt1, VertType *pt2)
/* is the pt0 on the line seg pt1-pt2 */
{
    double dx1,dy1,dx2,dy2;
    
	dx1=(*pt0)[0]-(*pt1)[0];
    dx2=(*pt0)[0]-(*pt2)[0];
    
	/* dx1 and dx2 must have different sign */
    if(dx1*dx2>0.000001)return 0;

    dy1=(*pt0)[1]-(*pt1)[1];
    dy2=(*pt0)[1]-(*pt2)[1];

    /* dy1 and dy2 must have different sign */
    if(dy1*dy2>0.000001)return 0;
    
	/* must at the same line */
    dx1=dy1*dx2-dx1*dy2;
    if(dx1<0.0001 && dx1>-0.0001)return 1;
    else return 0;
}

int check_bounding_boxes_intersection(VertType *pt0, VertType *pt1,
									  VertType *pt2, VertType *pt3, int dim)
									  /* find the bounding box to quickly check intersection */
									  /* dim: 001: x, 010: y, 100: z */
{
    double max_x0, min_x0;
    double max_x1, min_x1;
    double x0,x1,x2,x3;
    int i,ind=0;
	
    for(i=0; i<3; i++)
	{
		if(dim&1)
		{
			x0=(*pt0)[ind];
			x1=(*pt1)[ind];
			x2=(*pt2)[ind];
			x3=(*pt3)[ind];
			
			if(x0>x1)
			{
				max_x0=x0; min_x0=x1; 
			}
			else
			{
				max_x0=x1; min_x0=x0; 
			}
			
			if(x2>x3)
			{
				max_x1=x2; min_x1=x3; 
			}
			else
			{
				max_x1=x3; min_x1=x2; 
			}
			
			if(min_x0 > max_x1+0.00001)return 0;
			if(min_x1 > max_x0+0.00001)return 0;
		}
		ind++;
		dim=dim>>1;
    }
    return 1;
}

int find_intersection(VertType *pt0, VertType *pt1, 
					  VertType *pt2, VertType *pt3, VertType *ret_pt, 
					  double *ret_ratio, double *ret_ratio2)
					  /* If they intersect at the end point, it is considered as not intersected */
					  /* the return ratio is respective to the pt0, pt1 */
{
    double x,y,x0,y0,x1,y1,x2,y2,x3,y3,dx0,dy0,dx2,dy2;
    double c1,c0,delta,ratio1,ratio2;
    double length1, length2;
    int cnt,res;
	
	
    res=check_bounding_boxes_intersection(pt0,pt1,pt2,pt3,3);
    if(!res)return 0;
	
    x0=(*pt0)[0];
    y0=(*pt0)[1];
    x1=(*pt1)[0];
    y1=(*pt1)[1];
    x2=(*pt2)[0];
    y2=(*pt2)[1];
    x3=(*pt3)[0];
    y3=(*pt3)[1];
	
    length1=line_distance_2d(pt0,pt1);

    if(length1<0.00000001)
	{
		if(line_distance_2d(pt2,pt3)<0.00000001)return 0; /* same point */
		/* see if they overlap */
		if(line_distance_2d(pt0,pt2)<0.00000001)return 3; /* intersect one end */
		if(line_distance_2d(pt0,pt3)<0.00000001)return 3; /* intersect one end */
		if(is_pt_on_line(pt0,pt2,pt3))return 1; /* yes intersection */
		return 0;
    }
    
	length2=line_distance_2d(pt2,pt3);
    if(length2<0.00000001)
	{
		/* see if they overlap */
		if(line_distance_2d(pt2,pt0)<0.00000001)return 3; /* intersect one end */
		if(line_distance_2d(pt2,pt1)<0.00000001)return 3; /* intersect one end */
		if(is_pt_on_line(pt2,pt0,pt1))return 1; /* yes intersection */
		return 0;
    }
	
    dx0=x1-x0;
    dy0=y1-y0;
    dx2=x3-x2;
    dy2=y3-y2;
    delta=(dy0*dx2) -(dx0*dy2);
    
	if(delta<0.0000001 && delta >-0.0000001)
	{ 
		/* parallel,  or overlapped  */
		
		if(fabs(x3-x0)>0.01)
		{
			dx2=x3-x0;
			dy2=y3-y0;
		}
		else 
		{
			dx2=x3-x1;
			dy2=y3-y1;
		}
		
		delta=(dy0*dx2) -(dx0*dy2);
		delta/=(fabs(dx0)+fabs(dy0)+fabs(dx2)+fabs(dy2)); /* normalize it */
		if(fabs(delta)>0.0000001) return 0; /* parallel */ 
		/* overlapped */
		{
			double max,min,len;
			
			if(fabs(x1-x0)>fabs(y1-y0))
			{
				max=x0;
				min=x0;
				if(x1>max)max=x1;
				if(x1<min)min=x1;
				if(x2>max)max=x2;
				if(x2<min)min=x2;
				if(x3>max)max=x3;
				if(x3<min)min=x3;
				len=(fabs(x1-x0)+fabs(x3-x2))-(max-min);
				ratio1=len/fabs(x1-x0);
				ratio2=len/fabs(x3-x2);
			}
			else 
			{
				max=y0;
				min=y0;
				if(y1>max)max=y1;
				if(y1<min)min=y1;
				if(y2>max)max=y2;
				if(y2<min)min=y2;
				if(y3>max)max=y3;
				if(y3<min)min=y3;
				len=(fabs(y1-y0)+fabs(y3-y2))-(max-min);
				ratio1=len/fabs(y1-y0);
				ratio2=len/fabs(y3-y2);
			}
			
			if(ret_ratio) *ret_ratio=ratio1;
			if(ret_ratio2) *ret_ratio2=ratio2;
			return 5; /* overlapping */
		}
    }
    c0=x0*dy0 - y0*dx0;
    c1=x2*dy2 - y2*dx2;
    x=(dx2*c0-dx0*c1)/delta;
    y=(dy2*c0-dy0*c1)/delta;
	
    if(x1!=x0) ratio1=(x-x0)/(x1-x0);
    else ratio1=(y-y0)/(y1-y0);
    
	if(ret_ratio!=NULL)*ret_ratio=ratio1;
    
	if(ratio1<-0.00000001 || ratio1 >1.00000001) return 0;
	
    if(x3!=x2) ratio2=(x-x2)/(x3-x2);
    else  ratio2=(y-y2)/(y3-y2);
    
	if(ret_ratio2!=NULL)*ret_ratio2=ratio2;

    if(ratio2<-0.00000001 || ratio2 >1.00000001)
	{
		return 0;
    }
    
	cnt=0;
    
	if(ratio1<=0.00000001 || ratio1>=0.9999999) cnt++;
    
	if(ratio2<=0.00000001 || ratio2>=0.9999999)cnt++;
	
    if(ret_pt)
	{
		/*  need to return the result */
		(*ret_pt)[0]=x;
		(*ret_pt)[1]=y;
    }
	
    if(cnt==2)return 3; /* the intersection is right at the both end points */
    else if(cnt==1)return 2; /* one end is on the other line */
    
    return 1;
}

int solve_2_equations(double a0, double b0, double c0,
					  double a1, double b1, double c1, double *ret_t1, double *ret_t2)
{
    double delta,delta1,delta2;
    delta=(a0*b1) -(a1*b0);
    delta1=(c0*b1) -(c1*b0);
    delta2=(a0*c1) -(a1*c0);

    if(fabs(delta)<0.001)
	{
		/* parallel,  or overlapped  */
		if(fabs(delta1)>0.001)return 0; /* no solution   */
		
		if(fabs(delta2)>0.001)return 0; /* no solution   */
		
		return 2; /* infinite solution */
    }

    *ret_t1=delta1/delta;
    *ret_t2=delta2/delta;
    return 1;
}

int find_intersection_3D(VertType *pt0, VertType *pt1, 
						 VertType *pt2, VertType *pt3, VertType *ret_pt, 
						 double *ret_ratio)
						 /* If they intersect at the end point, it is considered as not intersected */
						 /* the return ratio is respective to the pt0, pt1 */
{
    double t1,t2,t1_ary[3],t2_ary[3];
    double a0,a1,a2,b0,b1,b2,c1,c2,c0;
    int i,res,index=0,cnt;
    
    if(line_distance(pt0,pt2)<0.000001 && 
		line_distance(pt1,pt3)<0.000001)return 3; /* overlapping */
    res=check_bounding_boxes_intersection(pt0,pt1,pt2,pt3,7);
    
	if(!res)return 0;
	
    a0=(*pt1)[0]-(*pt0)[0]; 
    a1=(*pt1)[1]-(*pt0)[1]; 
    a2=(*pt1)[2]-(*pt0)[2];
	
    b0=(*pt2)[0]-(*pt3)[0]; 
    b1=(*pt2)[1]-(*pt3)[1]; 
    b2=(*pt2)[2]-(*pt3)[2];
    /* solve the equation: pt0+v1*t1=pt2+v2*t2  */
    /* v1*t1 - v2*t2 = (pt2 - pt0)=0 */  
    c0=(*pt2)[0]-(*pt0)[0]; 
    c1=(*pt2)[1]-(*pt0)[1]; 
    c2=(*pt2)[2]-(*pt0)[2]; 
	
    /* 3 equations: a0* t1 + b0* t2 = c0 */
    res=solve_2_equations(a0,b0,c0,a1,b1,c1,&t1,&t2);
    if(res==0)return 0; /* no solution */
    
	if(res==1)
	{
		t1_ary[index]=t1; 
		t2_ary[index++]=t2;
    }
    
	res=solve_2_equations(a2,b2,c2,a1,b1,c1,&t1,&t2);
    if(res==0)return 0; /* no solution */
    
	if(res==1)
	{
		t1_ary[index]=t1; 
		t2_ary[index++]=t2;
    }
    
	res=solve_2_equations(a2,b2,c2,a0,b0,c0,&t1,&t2);
    if(res==0)return 0; /* no solution */
    
	if(res==1)
	{
		t1_ary[index]=t1; 
		t2_ary[index++]=t2;
    }
    
	if(index==0)
	{
		/* two lines are overlapping */
		/*	printf("overlapping or int. one pt\n"); */
		return 3; /* intersect at one end, or overlapped */
    }

    for(i=1; i<index; i++)
	{
		double diff;
	
		diff = t1_ary[i] -t1_ary[i-1];
		if(fabs(diff)>0.001)return 0;

		diff = t2_ary[i] -t2_ary[i-1];
		if(fabs(diff)>0.001)return 0;
    }
    
	if(t1<= -0.9999999999)return 0;
    if(t2<= -0.9999999999)return 0;
    if(t1>= 1.00000000001)return 0;
    if(t2>= 1.00000000001)return 0;
    cnt=0;
    if(t1>=0.999999 || t1<=0.000001)cnt++;
    if(t2>=0.999999 || t2<=0.000001)cnt++;
    if(cnt==2)return 4; /* intersect at both ends */
    /* printf("has intersection, t1=%lf t2=%lf\n",t1_ary[0],t2_ary[0]); */
    return 1;
}


int is_inside_contour(VertType *pt, VertType *lines_ary,
					  short *group_tab, short ary_ind, double *dist, int *ret_ind)
{
/* test the intersection of the line from (0,0) to pt with the 
contour segment. The chance that this line goes through right at
	a verice point is little. */
	
    VertType pt1,int_pt;
    double min_dist=100000000.0,d,dx,dy,ratio,ratio2;
    int i,j,cnt=0,res,min_ind=-1;
	
    /* find bounding box of the contour, so we can determine a point outside */
    double min_x=1e14,min_y=1e14;
    VertType *p;

    for(i=0; i<ary_ind; i++)
	{
		p=&(lines_ary[group_tab[i]]);
		if((*p)[0]<min_x)min_x=(*p)[0];
		if((*p)[1]<min_y)min_y=(*p)[1];
    }
	
    pt1[1]=min_y-0.87654;
    /* use the slope=1.123456, so the return distance is reliable */
    pt1[0]=(pt1[1]-(*pt)[1])*1.123456+(*pt)[0];
    
	for(i=0; i<ary_ind; i++)
	{
		j=i+1;
		if(j>=ary_ind)j=0;
		res=find_intersection(&(lines_ary[group_tab[i]]),
			&(lines_ary[group_tab[j]]),pt, &pt1, &int_pt,&ratio,&ratio2);
		if(res==3)
		{
			/* intersect at both end point */
			if(ret_ind)
			{
				if(ratio>0.5) *ret_ind=j;
				else *ret_ind=i;
			}
			
			if(dist!=NULL)*dist=0.0;
			return 3;
		}
		else if(res && line_distance_2d(pt,&int_pt)<0.00000001)
		{
			if(dist!=NULL)*dist=0.0;
			if(ret_ind)
			{
				if(ratio<0.01)*ret_ind=i; /* intersect at i */
				else  if(ratio>0.99)*ret_ind=j; /* intersect at j */
												/*
												else {
												printf("*** Error Intersect not end pt, ratio=%lf\n",ratio);
												}
				*/
			}
			return 2; /* right on the contour */
		}
		if(res)
		{
			dx=int_pt[0]-(*pt)[0];
			dy=int_pt[1]-(*pt)[1];
			d=dx*dx + dy*dy;
		
			if(d<min_dist)
			{ 
				min_dist=d;
				min_ind=i;
			}
			cnt++;
		}
    }
    if(dist!=NULL)*dist=min_dist;
    if(ret_ind!=NULL)*ret_ind=min_ind;
	/*
    if(min_dist<0.000000001)return 2;  right on the contour */
    return (cnt&1);
}


int is_contour_CCW(VertType *lines_ary,
				   short *group_tab, short ary_ind)
{
    double area;
    area = calculate_area(lines_ary, group_tab, ary_ind); 
    if(area>=0)return 1;
    else return 0;
}

/* calculate the angle between two vector, must be normalized */
/* reutn an angle between -PI and PI */
double calc_angle(VertType *u, VertType *v)
{
    double inner, outer,angle;
    
	outer=outer_product_verttype(u,v);
    inner=inner_product_verttype(u,v);
    angle=asin(outer);

    if(inner<0)
	{
		if(angle>0)angle=PI-angle;
		else angle=-PI-angle;
    }

    return angle;    
}


void draw_small_square(FILE *ofp, double x, double y, double z)
{
    fprintf(ofp,"4\n");
    fprintf(ofp,"%lf %lf %lf\n",x-0.1,y-0.1,z);
    fprintf(ofp,"%lf %lf %lf\n",x-0.1,y+0.1,z);
    fprintf(ofp,"%lf %lf %lf\n",x+0.1,y+0.1,z);
    fprintf(ofp,"%lf %lf %lf\n",x+0.1,y-0.1,z);
}
void draw_square(FILE *ofp, double x, double y, double z)
{
    fprintf(ofp,"4\n");
    fprintf(ofp,"%lf %lf %lf\n",x-0.3,y-0.3,z);
    fprintf(ofp,"%lf %lf %lf\n",x-0.3,y+0.3,z);
    fprintf(ofp,"%lf %lf %lf\n",x+0.3,y+0.3,z);
    fprintf(ofp,"%lf %lf %lf\n",x+0.3,y-0.3,z);
}
void draw_grid(FILE *ofp, double x, double y, double z)
{
    fprintf(ofp,"4\n");
    fprintf(ofp,"%lf %lf %lf\n",x-0.22,y-0.22,z);
    fprintf(ofp,"%lf %lf %lf\n",x-0.22,y+0.22,z);
    fprintf(ofp,"%lf %lf %lf\n",x+0.22,y+0.22,z);
    fprintf(ofp,"%lf %lf %lf\n",x+0.22,y-0.22,z);
}

void draw_line_seg(FILE *ofp, double x1, double y1, double z1,
				   double x2, double y2, double z2)
{
    fprintf(ofp,"3\n");
    fprintf(ofp,"%lf %lf %lf\n",x1,y1,z1);
    fprintf(ofp,"%lf %lf %lf\n",x2,y2,z2);
    fprintf(ofp,"%lf %lf %lf\n",(x1+x2)/2.0,(y1+y2)/2.0,(z1+z2)/2.0);
}
void draw_vert2vert(FILE *ofp, VertType *p1,VertType *p2)
{
    fprintf(ofp,"3\n");
    fprintf(ofp,"%lf %lf %lf\n",(*p1)[0],(*p1)[1],(*p1)[2]);
    fprintf(ofp,"%lf %lf %lf\n",(*p2)[0],(*p2)[1],(*p2)[2]);
    fprintf(ofp,"%lf %lf %lf\n",((*p1)[0]+(*p2)[0])*0.5,((*p1)[1]+(*p2)[1])*0.5,
        ((*p1)[2]+(*p2)[2])*0.5);
}

void draw_triangle(FILE* ofd, VertType *p1,  VertType *p2, VertType *p3)
{
	
/* In order to port to the Intel Paragon Supercomputer, use 
	file descriptor instead of file pointer */
    int len,j;
    char buf[1024];
    VertType v1,v2,v3;
	
    if(ofd<0)return;
    
	memcpy(&v1,p1,sizeof(v1));
    memcpy(&v2,p2,sizeof(v2));
    memcpy(&v3,p3,sizeof(v3));
	/*
	extern int ALL_CORRESPOND;
    if(ALL_CORRESPOND){
	scale_back_triangle(&v1);
	scale_back_triangle(&v2);
	scale_back_triangle(&v3);
    }
	*/
    sprintf(buf,"3\n%.15g %.15g %.15g\n%.15g %.15g %.15g\n%.15g %.15g %.15g\n",
		v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2]);
    len=strlen(buf);
    //j=write(ofd, buf, len); //klc
	j = fwrite(buf, sizeof(char) , len , ofd); //klc
    
	if(len!=j)
	{
		fprintf(stderr,"*** Error, draw_triangle(), should write %d, but %d\n",len,j);
		exit(1);
    }
}

void draw_triangle_FP(FILE *ofp, VertType *p1,  VertType *p2, VertType *p3)
{
    if(!ofp)return;
    fprintf(ofp,"3\n");
    fprintf(ofp,"%lf %lf %lf\n", (*p1)[0], (*p1)[1], (*p1)[2]);
    fprintf(ofp,"%lf %lf %lf\n", (*p2)[0], (*p2)[1], (*p2)[2]);
    fprintf(ofp,"%lf %lf %lf\n", (*p3)[0], (*p3)[1], (*p3)[2]);
}

void draw_triangle_mode_FP(FILE *ofp, VertType *p1,  VertType *p2, VertType *p3)
{
    if(!ofp)return;
    fprintf(ofp,"3 1\n");
    fprintf(ofp,"%lf %lf %lf\n", (*p1)[0], (*p1)[1], (*p1)[2]);
    fprintf(ofp,"%lf %lf %lf\n", (*p2)[0], (*p2)[1], (*p2)[2]);
    fprintf(ofp,"%lf %lf %lf\n", (*p3)[0], (*p3)[1], (*p3)[2]);
}


