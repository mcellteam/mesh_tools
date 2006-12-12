#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#include "common.h"
#include "myutil.h"
#include "filter.h"
#include "mymarch.h"
#include "math_util.h"
#include "group.h"
#include "approximate.h"
#include "correspond.h"
#include "image2cont2D.h"

#define VERT_TRI_BLOCK_SIZE (4096*3) 
#define VERT_BLOCK_SIZE 2048 


extern int FILTER_TYPE,DEBUG;
static VertType *Loc_lines_ary;
static int Loc_lines_ary_size;
static int Loc_lines_ind ;

extern int  Current_slice_num;
/* put the following variables as globle to save calling time */
static short **Group_ary,*Group_ind_ary;
static short **New_group_ary,*New_group_ind_ary;
static int Group_num;
extern int Base_slice_num;

extern double Tolerance;
static double Thres; 
int Begin_slice_num;
static float *Fdata;
static int Dimx,Dimy,Dimz,Detx=1,Dety=1,Detz=1;
static float *Xpos,*Ypos,*Zpos;
static double Xunit=1.0,Yunit=1.0, Zunit=3.0;
extern double Current_zvalue;

extern int STATIC_CHQ;

void my_clear_image2cont2D()
{
	Loc_lines_ary = NULL;
	Loc_lines_ary_size = 0;
	Loc_lines_ind =0;
	Group_ary = NULL;	Group_ind_ary=NULL;
	Group_num=0;
	Thres =0.0;
	Fdata = NULL;
	Dimx = Dimy = Dimz =0;
	Detx=Dety=Detz=1;
	Xpos = NULL;	Ypos = NULL;	Zpos = NULL;	
	Xunit=1.0; 	Yunit=1.0; Zunit=3.0;

	count_triangles_one_slice(STATIC_CHQ, NULL, NULL);
}


void save_contours_pts_format(char *s, int slice_num)
{
    char name[100];
    FILE *ofp;
    short *group_tab, group_ind,group_num;
    VertType *lines_ary;
    int i,j,group;
    sprintf(name,"%s%d.pts",s,slice_num);

    if((ofp=fopen(name,"w"))==NULL) 
	{
        fprintf(stderr,"save_contours_pts_format() Could not open %s for writing\n",name);
        exit(1);
    }
    
	group_num=Group_num;
    lines_ary=Loc_lines_ary;
    
	for(group=0; group<group_num; group++)
	{
        group_tab=New_group_ary[group];
        group_ind=New_group_ind_ary[group];
        fprintf(ofp, "%d\n",group_ind);
    
		for(i=0; i<group_ind; i++)
		{
            j=group_tab[i];
            fprintf(ofp,"%lf %lf %lf\n",lines_ary[j][0],
                lines_ary[j][1],lines_ary[j][2]);
        }
    }

    fclose(ofp);
    fprintf(stderr,"Saved contour file %s\n",name);
	
}

void save_marching_contours_pts_format(char *s, int slice_num)
{
    char name[100];
    FILE *ofp;
    short *group_tab, group_ind,group_num;
    VertType *lines_ary;
    int i,j,group;
    sprintf(name,"%s%d.pts",s,slice_num);
    
	if((ofp=fopen(name,"w"))==NULL) 
	{
        fprintf(stderr,"save_contours_pts_format() Could not open %s for writing\n",name);
        exit(1);
    }
    
	group_num=Group_num;
    lines_ary=Loc_lines_ary;
    
	for(group=0; group<group_num; group++)
	{
        group_tab=Group_ary[group];
        group_ind=Group_ind_ary[group];
        fprintf(ofp, "%d\n",group_ind);
    
		for(i=0; i<group_ind; i++)
		{
            j=group_tab[i];
            fprintf(ofp,"%lf %lf %lf\n",lines_ary[j][0],
                lines_ary[j][1],lines_ary[j][2]);
        }
    }

    fclose(ofp);
    fprintf(stderr,"Saved contour file %s\n",name);
	
}

void set_boundary_zero(float *fary, int dimx, int dimy)
{
    int i,j,ind,ind1,area;
    area=dimx*dimy;
    
	for(i=0; i<dimx; i++)fary[i]=0; 
    
	for(i=area-dimx; i<area; i++)fary[i]=0; 
    
	ind=dimx-1;
    ind1=0;
    
	for(j=0; j<dimy; j++)
	{
		fary[ind]=0;
		fary[ind1]=0;
		ind+=dimx;
		ind1+=dimx;
    }
}

void get_a_slice(float *fary, int slice_num, float (*valfun)(int, int, int),
				 void (*slicefun)(int, int))
{
    int i,j,ind=0;

    if (slicefun != NULL) (*slicefun)(slice_num, 1);

    for(j=0; j<Dimy; j++)
	{
		for(i=0; i<Dimx; i++) 
		{
			float val;
			val=(float)((*valfun)(i,j,slice_num));
			fary[ind++]=val;
		}
    }

    if (slicefun != NULL) (*slicefun)(slice_num, 0);
    
	if(FILTER_TYPE==1)
	{
		fprintf(stderr,"Applied a 3*3 averaging low pass filter\n");
		lowpass_filter33(fary,Dimx,Dimy); /* apply median filter */
		fprintf(stderr,"Finished filtering\n");
    }
    else if(FILTER_TYPE==2)
	{
		fprintf(stderr,"Applied a 3*3 median filter\n");
		median_filter33(fary,Dimx,Dimy); /* apply median filter */
		fprintf(stderr,"Finished median filtering\n");
    }
    set_boundary_zero(fary,Dimx,Dimy);
}

static void local_alloc_Lines_ary_if_empty()
{
    if(Loc_lines_ary==NULL) 
	{
        Loc_lines_ary_size=VERT_BLOCK_SIZE;
        Loc_lines_ary=(VertType *)mycalloc(sizeof(*Loc_lines_ary)*
			Loc_lines_ary_size);
    }
}

static int local_put_lines_on_map(int line_num, 
								  double lines[4][SPACE_VERT_DIM])
{
    int i,line,temp;
	
    local_alloc_Lines_ary_if_empty();
    temp=Loc_lines_ind;

	/*
	if(DEBUG)fprintf(stderr,"after local_alloc_Lines_ary_if_empty()\n");
	*/
    
	for(line=0; line<line_num; line++) 
	{
        int k,k1;
        k=line*2;
        for(i=k; i<k+2; i++)
		{
			/* leave extra 1000 space for breaking contour */
            if(Loc_lines_ary_size<=Loc_lines_ind+1000)
			{/* enlarge the buffer */
                int temp;
                VertType *ary;
                temp=Loc_lines_ind;
                Loc_lines_ary_size+=VERT_BLOCK_SIZE;
                ary=Loc_lines_ary;
                Loc_lines_ary=(VertType *)
					mycalloc(sizeof(*Loc_lines_ary)*Loc_lines_ary_size);
                memcpy(Loc_lines_ary,ary,sizeof(*ary)*temp);
                free(ary);
            }
			k1=Loc_lines_ind;
			memcpy(&Loc_lines_ary[k1],lines[i],
				sizeof(double)*SPACE_VERT_DIM);
            Loc_lines_ind++;
        }
    }
    return temp;
}


void shift_group_sequence(short *group_tab, int tab_size, int origin)
{
    int i,j;
    short ary[4000];

    if(origin>=tab_size)
		fprintf(stderr,"I. Error, tab_size=%d but origin is %d (big)\n",
		tab_size,origin);
    for(i=0; i<tab_size; i++)
	{
		j=i-origin;
		if(j<0)j=j+tab_size;
		ary[j]=group_tab[i];
    }
    memcpy(group_tab, ary, sizeof(short)*tab_size);
    group_tab[tab_size]=group_tab[0];
}


void reverse_group_direction(short *group_tab, int ind)
{
    int i,j,ind1,k;
    short temp;
    j=ind>>1;
    ind1=ind-1;
    
	for(i=0; i<j; i++)
	{
		k=ind1-i;
		temp=group_tab[i];
		group_tab[i]= group_tab[k];
		group_tab[k]= temp;
    }
    group_tab[ind]=group_tab[0]; /* make it to be circular */
}



static void adjust_contour_origin(VertType *lines_ary, short **group_ary,
								  short *group_ind_ary, int group_num)
								  /* pick up the left most point as the start point, This is required to 
								  produce the exactly same sequence while it is processed at top and
								  bottom. */
{
    int i,j,k,min_ind;
    double x,y, min_x, min_y;

    for(i=0; i<group_num; i++)
	{
		min_x=10000000.0;
		min_y=10000000.0;
		
		for(k=0; k<group_ind_ary[i]; k++)
		{
			j=group_ary[i][k]; /* get the first point */
			x=lines_ary[j][0];
			y=lines_ary[j][1];
			
			if(x<min_x+0.0001)
			{
				if(x>min_x-0.0001)
				{
					/* consider to be equal */
					if(y<=min_y)
					{
						min_y=y;
						min_ind=k;
					}
				} /* a new min, pt */
				else
				{
					min_x=x;
					min_y=y;
					min_ind=k;
				}
			}
		}
		/* printf("shift sequence to be %d (%lf, %lf), size=%d\n",
		min_ind,min_x,min_y,group_ind_ary[i]);
		*/
		if(min_ind) shift_group_sequence(group_ary[i],group_ind_ary[i],min_ind);
    }
}

static void adjust_contour_directions( VertType *lines_ary, short **group_ary,
									  short *group_ind_ary, int group_num)
{
    int i,j,j1,res;
    double x0,y0,x1,y1,val;
    int ix0,iy0,reverse=0;
    VertType pt;
    signed char is_inside_positive;
	
	
    for(i=0; i<group_num; i++)
	{
		j=group_ary[i][0]; /* get the first point */
		x0=lines_ary[j][0];
		y0=lines_ary[j][1];
		j1=group_ary [i][1]; /* get the second point */
		x1=lines_ary[j1][0];
		y1=lines_ary[j1][1];
		ix0=(int)x0;
		iy0=(int)y0;
		pt[0]=ix0;
		pt[1]=iy0;
		
		res=is_inside_contour(&pt, lines_ary, group_ary[i], group_ind_ary[i],
			NULL,NULL);
		val=Fdata[iy0*Dimx+ix0]-Thres;
		if(res) 
		{ /* It is inside the contour */
			if(val>0)is_inside_positive=1;
			else is_inside_positive=0;
		}
		else 
		{
			if(val>0)is_inside_positive=0;
			else is_inside_positive=1;
		}
		printf("Polarity of group %d is %d\n",
			i,is_inside_positive);
		
		/* Warning: the following code does not work */
		/* make the right side intensity > 0 */
		reverse=0;
		if((double)ix0==x0)
		{ /* It is on X grid line */
			if(x1>x0 && val <0) reverse=1; /* (ix0,iy0) on  right side */
			else if(x1<x0 && val >0) reverse=1; /* (ix0,iy0) on left side */
		}
		else if((double)iy0==y0)
		{ /* It is on Y grid line */
			if(y1>y0 && val >0) reverse=1; /* (ix0,iy0) on  left side */
			else if(y1<y0 && val <0) reverse=1; /* (ix0,iy0) on right side */
		}
		else fprintf(stderr,"I. Error, pt=(%lf %lf)\n",x0,y0);
		
		if(reverse) reverse_group_direction(group_ary[i],group_ind_ary[i]);
    }
}


void free_pos_ary()
{
    free(Xpos); Xpos=NULL;
    free(Ypos); Ypos=NULL;
    free(Zpos); Zpos=NULL;
}

int initialization(int begin_slice_num,int dimx, int dimy, int dimz, 
				   double *orig, double *voxel_unit, double thres )
{
    int i;
    Thres=thres;
    Dimx=dimx; Dimy=dimy; Dimz=dimz;

    if(voxel_unit!=NULL) 
	{ 
		Xunit=voxel_unit[0]; 
		Yunit=voxel_unit[1]; 
		Zunit=voxel_unit[2]; 
    } 
    
	Begin_slice_num=begin_slice_num;
    Xpos=(float*)mycalloc(Dimx*sizeof(float));
    Ypos=(float*)mycalloc(Dimy*sizeof(float));
    Zpos=(float*)mycalloc(Dimz*sizeof(float));
    for(i=0; i<Dimx; i++)Xpos[i]=(float)(orig[0]+Xunit*i);
    for(i=0; i<Dimy; i++)Ypos[i]=(float)(orig[1]+Yunit*i);
    for(i=0; i<Dimz; i++)Zpos[i]=(float)(orig[2]+Zunit*i);
    init_marching_cubes(thres, Xpos, Ypos, Zpos,
        dimx, dimy, dimz, Detx, Dety, Detz);
    return 1;
}


void allocate_image_memory(int area)
{
    Fdata=(float*)mycalloc(area*sizeof(float));
    Group_ary=(short **)mycalloc(MAX_GROUP_NUM*sizeof(*Group_ary));
    Group_ind_ary=(short *)mycalloc(MAX_GROUP_NUM*sizeof(*Group_ind_ary));
    New_group_ary=(short **)mycalloc(MAX_GROUP_NUM*sizeof(*New_group_ary));
    New_group_ind_ary=(short *)mycalloc(MAX_GROUP_NUM*sizeof(*New_group_ind_ary));
}
void free_image_memory()
{
    free(Fdata);
    Fdata=NULL;
    free(Group_ary);
    Group_ary=NULL;
    free(Group_ind_ary);
    Group_ind_ary=NULL;
}
/* got contour from one slice. */
int process_one_slice(float *ary)
{
    int i,j;
    int dimx_detx;
    Loc_lines_ind=0;
    dimx_detx=Dimx-Detx;
	
	if(DEBUG)fprintf(stderr," process_one_slice, Dimx %d,Dimy %d Detx %d Dety=%d\n",
		Dimx,Dimy,Detx,Dety);
    
	for(j=0; (j+Dety)<Dimy; j=j+Dety) 
	{
		for(i=0; i<dimx_detx; i=i+Detx)
		{
			double lines[4][SPACE_VERT_DIM];
			double box[4], vals[4];
			int line_num,temp;
			box[0]=Xpos[i]; box[1]=Xpos[i+1];
			box[2]=Ypos[j]; box[3]=Ypos[j+1];
			temp=j*Dimx+i;
			vals[0]=ary[temp];
			vals[1]=ary[temp+1];
			vals[2]=ary[temp+Dimx+1];
			vals[3]=ary[temp+Dimx];
			line_num=calc_lines(box,vals,lines); 
			if(line_num>0)
			{
				int i1;
				if(DEBUG)fprintf(stderr," j=%d i=%d, calc_lines()=%d\n",j,i,line_num);
			
				for(i1=(line_num<<1)-1; i1>=0; i1--)
				{
					lines[i1][2]=Current_zvalue;
					lines[i1][2]=Current_zvalue;
				}
				
				if(DEBUG)fprintf(stderr," before  local_put_lines_on_map()\n");
				local_put_lines_on_map(line_num,lines);
			}
		}
    }
	if(DEBUG)fprintf(stderr,"before trace_contours\n");
    Group_num=trace_contours(Loc_lines_ary,Loc_lines_ind,
		Group_ary, Group_ind_ary);
		/*
		adjust_contour_directions(Loc_lines_ary,Group_ary,
        Group_ind_ary, Group_num);
	*/
	
    /* adjust_contour_origin must after adjust_contour_direction(), if put
	in front, it wont work */
	if(DEBUG)fprintf(stderr,"before adjust_contour_origin\n");
    adjust_contour_origin(Loc_lines_ary,Group_ary,
        Group_ind_ary, Group_num);
	
	
		/*
		printf("%d groups\n",Group_num);
		for(i=0; i<Group_num; i++){
		short *temp_ary;
		int ind,j,k;
		temp_ary=Group_ary[i];
		ind=Group_ind_ary[i];
		printf(" %d pts\n",ind);
		for(k=0; k<ind; k++){
        j=temp_ary[k];
        printf(" %f %f %f\n",Loc_lines_ary[j][0],Loc_lines_ary[j][1],
		Current_zvalue);
		}
		}
	*/
    return 1;
}




int read_and_approximate(int k,
						 float (*valfun)(int, int, int),
						 void (*slicefun)(int, int))
{
    int i,j,res;
    double toler_ary[MAX_GROUP_NUM];
    char flag_ary[MAX_GROUP_NUM],error;
	
    Current_zvalue=Zpos[k-Begin_slice_num];
    Current_slice_num=Base_slice_num+k;
	/* fprintf(stderr,"read_and_approximate(), k=%d, Zval=%lf, slice#=%d\n",
    k, Current_zvalue,Current_slice_num);
	*/
	if(DEBUG)fprintf(stderr,"before get_a_slice()\n");
    get_a_slice(Fdata,k-Begin_slice_num,valfun,slicefun);
	if(DEBUG)fprintf(stderr,"before process_one_slice()\n");
    process_one_slice(Fdata);
	
    /* approximate contours, the Group_ary & Group_ind_ary will
	   be changed */
	   /* do both slices, because the previous might be a little different
	   from this one */
	
	if(DEBUG)fprintf(stderr,"before approximate_contours()\n");
	/* save_marching_contours_pts_format("test", k); */ 

	if(Tolerance<=0.001)
	{ // no approximation, simply copy
		for(i=0; i<Group_num; i++)
		{
			New_group_ind_ary[i]=Group_ind_ary[i];
			New_group_ary[i]=(short *)mymalloc(sizeof(short)*Group_ind_ary[i]+4);
			for(j=0; j<Group_ind_ary[i]; j++)New_group_ary[i][j]=Group_ary[i][j];
		}
		return 1;
	}
	
    for(i=0; i<Group_num; i++)	toler_ary[i]=Tolerance;

    while(1)
	{
		approximate_contours(toler_ary,
			Loc_lines_ary,Loc_lines_ind,
			Group_ary, Group_ind_ary, Group_num,
			New_group_ary, New_group_ind_ary);
	
		for(i=0; i<Group_num; i++)	flag_ary[i]=0;

		for(i=1; i<Group_num; i++)
		{
			for(j=0; j<i; j++)
			{
				res=check_if_two_contours_intersect(Loc_lines_ary,Loc_lines_ary,
					New_group_ary[i], New_group_ind_ary[i],New_group_ary[j],
					New_group_ind_ary[j],0, NULL, NULL, NULL);
				if(res)
				{ /* intersected */
					flag_ary[i]=1;
					flag_ary[j]=1;
					printf("Warning! G%d intersect G%d\n",i,j); 
				}
			}
		}
		error=0;

		for(i=0; i<Group_num; i++)
		{
			if(flag_ary[i])
			{
				toler_ary[i]/=2.0; /* tighten the tolerance */
				error=1;
				printf("Warning! G%d tolerance reduced to %lf\n",i,toler_ary[i]); 
			}
		}

		if(!error)return 1; /* no intersection */
		
		for(i=0; i<Group_num; i++)
		{
			free(New_group_ary[i]);
			New_group_ary[i]=NULL;
		}
		
    }
}


int count_triangles_one_slice(int k,
							  float (*valfun)(int, int, int),
							  void (*slicefun)(int, int))
{
    int i,j;
    int tri_num=0;
    static float *fdata[2]={NULL,NULL};

    if(fdata[0]==NULL)
	{
		fdata[0]=(float*)mymalloc(sizeof(float)*Dimx*Dimy);
		fdata[1]=(float*)mymalloc(sizeof(float)*Dimx*Dimy);
    }
    
	Current_slice_num=Base_slice_num+k;
    get_a_slice(fdata[0],k-Begin_slice_num,valfun,slicefun);
    get_a_slice(fdata[1],k-Begin_slice_num+1,valfun,slicefun);
	
    for(j=0; j<Dimy-1; j++) 
	{
        for(i=0; i<Dimx-1; i++)
		{
            tri_num+=count_triangles(fdata[0],fdata[1],i,j);
        }
	}

    printf("slice %d, tri_num=%d\n",k,tri_num);
    return tri_num;
}
