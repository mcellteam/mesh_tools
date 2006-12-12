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

/*
This program sort 2-D points into grids, then remove the redundant points.
It group the line segments into contours.
*/

#define  DEBUG  1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <memory.h>

#include "common.h"
#include "myutil.h"
#include "group.h"

#define VERT_BLOCK_SIZE 4096
#define LINE_BLOCK_SIZE 2048
#define CUBE_PER_DIM 64 
#define CUBE_PER_PLANE (CUBE_PER_DIM*CUBE_PER_DIM)
#define NEIGHBOR_NUM 50
#define GROUP_OFFSET 10000 


extern int SILENT;
static VertType *Vert_ary;
static int Processed_mark;
static int *Processed;
static int Vert_ind=0;
static int *Pindex;
static int *Matrix;
static int Total_vertice_num=0;
static double Bound_box[PLANE_VERT_DIM][2],Det[PLANE_VERT_DIM],
Det_large[PLANE_VERT_DIM];
static double Epsilon,Epsilon2,EpsilonSquare;

extern int STATIC_CHQ;

void my_clear_group()
{
	Vert_ary = NULL;
	Processed_mark =0;
	Processed = NULL;
	Vert_ind =0;
	Pindex =NULL;
	Matrix = NULL;
	Total_vertice_num =0;
	Epsilon = Epsilon2 = EpsilonSquare=0.0;
	// Hopefully no need to re-init the Bound_box Det and Det_large Ds...
}

void find_bound_box()
{
	int   i,j;
	double min[PLANE_VERT_DIM], max[PLANE_VERT_DIM];
	
	for (i = 0; i < PLANE_VERT_DIM; i++)
	{
		min[i] = 1000000.0;
		max[i] = -1000000.0;
	}
	
	for (i = 0; i < Vert_ind; i++)
	{
		for (j = 0; j < PLANE_VERT_DIM; j++)
		{
			if (Vert_ary[i][j] <= min[j]) min[j] = Vert_ary[i][j];
			if (Vert_ary[i][j] >= max[j]) max[j] = Vert_ary[i][j];
		}
	}

	for (i = 0; i < PLANE_VERT_DIM ; i++)
	{
		double val;
		val=(max[i]-min[i])*0.01;
		if(val<(1.1*Epsilon))val=1.1*Epsilon; 
		Bound_box[i][0]=min[i]-val;
		Bound_box[i][1]=max[i]+val;
	}
	/*
	fprintf(stderr,"Bounding [%f %f] [%f %f]\n",Bound_box[0][0],
	Bound_box[0][1],Bound_box[1][0],Bound_box[1][1]);
	*/
}

/*---------------------------------------------------------------------------
ReorderData -- reorder the data according to a value
----------------------------------------------------------------------------*/
void
ReorderData(int begin,int end,int xyz, double value,int *end1, int *begin1)
/* begin,end;                the begin and end indeces of the range
to be ordered                           */
/* xyz;                       order according to x,y,z ccordinate if
xyz == 0, 1, 2, respectively            */
/* value;                     threshold value                         */
/* *end1,*begin1;              partition ends                          */
{
	int    j,k;
	if(begin>end)
	{ /* no data inside */
		*end1 = end;
		*begin1 = begin;
		return;
	}
	
	j = begin;
	k = end;
	while(k>=j)
	{ /* Note: must use >= not just > */
		while((j <= end) && (Vert_ary[Pindex[j]][xyz] <= value))j++;
		while(((k >= begin) && Vert_ary[Pindex[k]][xyz] > value))k--;
		if(k>j)
		{
			int temp;
			temp=Pindex[j];
			Pindex[j]=Pindex[k];
			Pindex[k]=temp;
		}
	}
	*end1 = j - 1;
	*begin1 = j;
}


/*---------------------------------------------------------------------------
put_into_subcubes(begin,end,a1,a2,b1,b2)
put the data in Vert_ary into subcubes.
----------------------------------------------------------------------------*/
void
put_into_subcubes(int begin, int end, double a1, 
				  double a2, double b1, double b2)
{
	int   i,j,dir,begin1,end1;
	double midval;
	
	if ( (a2 - a1 <= Det_large[0]) && (b2 - b1 <= Det_large[1]))
	{
		int temp;
		i = (int)((a1 - Bound_box[0][0])/Det[0] + 0.5);
		j = (int)((b1 - Bound_box[1][0])/Det[1] + 0.5);
		
		temp=(j*CUBE_PER_DIM + i)<<1;
		Matrix[temp] = begin;
		Matrix[temp + 1] = end;
		return;
	}
	/* partition direction                         */
	
	if (a2 - a1 >= Det_large[0] ) dir = 0;
	else if (b2 - b1 >= Det_large[1] ) dir = 1;
	
	/* subdivide the set of vertice in x-direction  */
	if (dir == 0) 
	{
		midval=(a1+a2)/2.0;
		ReorderData(begin,end,0,midval,&end1,&begin1);
		put_into_subcubes(begin,end1,a1,midval,b1,b2);
		put_into_subcubes(begin1,end,midval,a2,b1,b2);
	}
	else if (dir == 1) 
	{
		/* subdivide the set of vertice in y-direction  */
		midval=(b1+b2)/2.0;
		ReorderData(begin,end,1,midval,&end1,&begin1);
		put_into_subcubes(begin,end1,a1,a2,b1,midval);
		put_into_subcubes(begin1,end,a1,a2,midval,b2);
	}
}

void verify_sorted_data()
{
    int i,j,ind,i1,k1;
    int begin,end;
    int total_pt=0;
    double dval;
    
	for(j=0; j<CUBE_PER_DIM; j++)
		for(i=0; i<CUBE_PER_DIM; i++)
		{
			ind=(j*CUBE_PER_DIM+i)*2;
			begin=Matrix[ind];
			end=Matrix[ind+1];
			for(i1=begin; i1<=end; i1++)
			{
				total_pt++;
				k1=Pindex[i1];
				if(Processed[k1]==-1)Processed[k1]=ind;		
				else printf("Error! %d processed, should be in %d, did in %d\n",
					k1,Processed[k1],ind);
				dval=Det[0]*i+Bound_box[0][0]-0.000001;
				if(Vert_ary[k1][0]<dval || Vert_ary[k1][0]>dval+Det[0]+0.000002)
					printf("Error X, (%d %d) val=%f, [%f %f)\n",i,j,
					Vert_ary[k1][0],dval, dval+Det[0]);
				dval=Det[1]*j+Bound_box[1][0]-0.000001;
				if(Vert_ary[k1][1]<dval || Vert_ary[k1][1]>dval+Det[1]+0.000002)
					printf("Error Y, (%d %d) val=%f, [%f %f)\n",i,j,
					Vert_ary[k1][1],dval, dval+Det[1]);
			}
		}

		if(total_pt!=Vert_ind)
		{
			printf("Error! only %d of %d vertices verified\n", total_pt,Vert_ind);
			for(i=0; i<Vert_ind; i++)
			{
				j=Pindex[i];
				if(Processed[j]==-1)printf("%d ",i);
			}
			printf("\n");
		}
		for(i=0; i<Vert_ind; i++)Processed[i]=-1;
}


int  Neighbors(double *point,int *neiborindex)
{
    int i,beginx,beginy,endx,endy;
    double vx,vy,valx,valy;
    int   ii,jj,begin,end,owmany;
	
    /* 1 determine the start and end index        */
    /* 1.1 determint the center interval index    */
    valx = point[0] - Bound_box[0][0];
    valy = point[1] - Bound_box[1][0];
    beginx = (int)((valx - Epsilon2)/Det[0]);
    beginy = (int)((valy - Epsilon2)/Det[1]);
    endx = (int)((valx + Epsilon2)/Det[0]);
    endy = (int)((valy + Epsilon2)/Det[1]);
	
	
    if(beginx<0)beginx=0;
    if(endx>=CUBE_PER_DIM)endx = CUBE_PER_DIM - 1;
    if(beginy<0)beginy=0;
    if(endy>=CUBE_PER_DIM)endy = CUBE_PER_DIM - 1;
	
    owmany = 0;
    for (jj = beginy; jj <= endy; jj++)
	{
		for (ii = beginx; ii <= endx; ii++)
		{
            int temp;
            temp = (CUBE_PER_DIM*jj + ii)<<1;
            begin = Matrix[temp];
            end = Matrix[temp +1];
            for (i = begin; i <= end ; i++) 
			{
				int pt;
				pt=Pindex[i];
				if(Processed[pt]<Processed_mark)
				{
					vx = point[0] - Vert_ary[pt][0];
					vy = point[1] - Vert_ary[pt][1];
					vx = vx*vx + vy*vy;
					if(vx<EpsilonSquare) 
					{
						neiborindex[owmany] = pt;
						owmany++;
						if (owmany >= NEIGHBOR_NUM)
						{
							if(!SILENT)printf("*** Error, %d vertare within %f\n",
								NEIGHBOR_NUM,Epsilon);
							return 0;
							
						}
					}
				}
            }
		}
    }
    return owmany;
}

/* remove_redundant_vertice() should not be used */
void remove_redundant_vertice()
{
    int total_pt=0,ni;
    int k,owmany,cur_pt;
    int neighbor_ary[NEIGHBOR_NUM];
    for(cur_pt=0; cur_pt<Vert_ind; cur_pt++)
	{
		if(Processed[cur_pt]<Processed_mark)
		{ 
			/* not processed yet */
			total_pt++;
			Processed[cur_pt]=cur_pt; /* point to itself */
			owmany=Neighbors((double *)(&(Vert_ary[cur_pt])), neighbor_ary);
			if(owmany>0)
			{
				for(k=0; k<owmany; k++) 
				{
					ni=neighbor_ary[k];
					Processed[ni]=cur_pt; /* this vertice is merged */ 
				}
			}
			/********************************* */
			if(owmany!=1 && !SILENT)
				fprintf(stderr,"Error, cur_pt=%d, owmany=%d, [%f %f]\n",
				cur_pt,owmany,Vert_ary[cur_pt][0],Vert_ary[cur_pt][1]);
			
		}
    }
	/*
    fprintf(stderr,"Merge %d vertices into %d vertices\n",Vert_ind,total_pt);
	*/
    Total_vertice_num=total_pt;
}

/*
void label_group(int offset, short *vert_group_ary, VertType *vert_ary, 
short **group_ary, short *group_ind_ary, int group_num)
{
int group;
int j,k,owmany,cur_pt;
double point[PLANE_VERT_DIM];
int neighbor_ary[NEIGHBOR_NUM];
for(group=0; group<group_num; group++){
for(cur_pt=0; cur_pt<group_ind_ary[group]; cur_pt++){
int i1;
i1=group_ary[group][cur_pt];
memcpy(point,&(vert_ary[i1]),sizeof(point));
owmany=Neighbors(point, neighbor_ary);
if(owmany>0){
int ni;
for(k=0; k<owmany; k++) {
ni=neighbor_ary[k];
vert_group_ary[ni]=group+offset;
}
}
}
}
}
*/

int tri_region_growing(char *slice_group[2], int start_index,
					   short *vert_group_ary)
{
    int group,cnt=0,vert;
    int i,j,owmany;
    int neighbor_ary[NEIGHBOR_NUM];
    int pichist[MAX_TRI_NUM],index=0; /* use region growing to walk */
	
    /* store the 3 vertice of this triangles into stack */
    pichist[index++]=start_index++;

    while(index>=0)
	{
		int j3;
		cnt++;
		j=pichist[--index];
		j3=(j/3)*3;
		Processed[j]=Processed_mark;
		for(vert=j3; vert<j3+3; vert++)
		{
			/* pick the other two vertices */
			if(Processed[vert]<Processed_mark)
			{
				/* the other vert. */
				group=vert_group_ary[j];
				if(group>=0)
				{
					if(group<GROUP_OFFSET) slice_group[0][group]=1;
					else slice_group[1][group-GROUP_OFFSET]=1;
				}
				
				pichist[index++]=vert;
				Processed[vert]=Processed_mark;
				/* also put in its corres. vertice */	
				owmany=Neighbors((double *)(&(Vert_ary[vert])), neighbor_ary);
				for(i=0; i<owmany; i++)
				{
					j=neighbor_ary[i];
					pichist[index++]=j;
					Processed[j]=Processed_mark;
				}
			}
		}
    }
    return -1; /* no correspondence */
}

int trace_contours(VertType *vert_ary, 
				   int vert_ind, short **group_ary, short *group_ind_ary)
				   /* input line_ary[] & vert_ary[], the output is in group_ary[][] */
{
    int j,i,k;
    int *table_ary;
    int ind, group_ind=0;
	
	
    Vert_ind=vert_ind;
    Vert_ary=vert_ary;
    
	if(Vert_ind<=2 )
	{
		if(!SILENT)printf("No line segments, or only 1 line seg.\n");
		return 0;
    }

    Processed_mark=0;
    find_bound_box();
    /* define the distance to merge vertices */
    Epsilon=0.00001;
    Epsilon2=Epsilon*2;
    EpsilonSquare=Epsilon*Epsilon;
    
	for(i=0; i<PLANE_VERT_DIM; i++)
	{
		Det[i]=(Bound_box[i][1]-Bound_box[i][0])/CUBE_PER_DIM;
		Det_large[i]=Det[i]*1.1;
    }
	
    Processed=(int *)mycalloc(sizeof(*Processed)*Vert_ind);
    Pindex=(int*)mycalloc(sizeof(*Pindex)*Vert_ind);
    Matrix=(int*)mycalloc(sizeof(*Matrix)*2*CUBE_PER_PLANE);
    
	for(i=0; i<Vert_ind; i++)
	{
		Pindex[i]=i;
		Processed[i]=-999999;
    }
    
	put_into_subcubes(0,Vert_ind-1,Bound_box[0][0],Bound_box[0][1],
		Bound_box[1][0],Bound_box[1][1]);
    /* verify_sorted_data(); */
    remove_redundant_vertice();
    free(Pindex);
    free(Matrix);
    Pindex=NULL;
    Matrix=NULL;
	
	
    /* The contour is supposed to be closed, so every 2 pts merge into one */
    table_ary=(int *)mycalloc(sizeof(int)*Vert_ind*2);
    
	for(i=Vert_ind*2-1; i>=0;  i--)table_ary[i]=-1;
    
	for(i=0; i<Vert_ind; i++)
	{
		j=Processed[i];
		k=j<<1;
		if(table_ary[k]==-1)table_ary[k]=i;
		else if(table_ary[k+1]==-1)table_ary[k+1]=i;
		else if(!SILENT)printf("Internal Error 1, k=%d\n",k);
    }
    
	ind=0;
    
	while(ind<Vert_ind)
	{
		short store_ary[2000],store_ind, *short_ary;
		int p2,old_pt; 
		int not_done=1,next_pt;
	
		while(table_ary[ind<<1]==-1 && ind<Vert_ind)ind++;
		
		if(ind>=Vert_ind)break;
		
		next_pt=ind;
		store_ind=0;
		
		while(not_done)
		{
			if(store_ind>=2000-1)
			{
				fprintf(stderr, "Error! store_ary[] overflow in trace_contours\n"); 
				exit(1);
			}
			
			if(next_pt&1)p2=next_pt-1; /* the next point is the other point */
				else p2=next_pt+1;
			store_ary[store_ind++]=next_pt;
			old_pt=next_pt;
			k=(Processed[p2])<<1;

			if(table_ary[k]!=-1 && table_ary[k]!=p2)
			{
				next_pt=table_ary[k];
			}
			else if(table_ary[k+1]!=-1 && table_ary[k+1]!=p2) 
			{
				next_pt=table_ary[k+1];
			}
			else 
			{
				not_done=0;
			}
			
			k=(Processed[old_pt])<<1;
			table_ary[k]=-1;
			table_ary[k+1]=-1;
		}
		/* if this contour has only 4 points, then delete it */
		if(store_ind>4)
		{
			short_ary = (short *)mycalloc( (store_ind+1)*sizeof(*short_ary));
			if(group_ary[group_ind] !=NULL)free( group_ary[group_ind]); 
			group_ary[group_ind] = short_ary; 
			memcpy(short_ary,store_ary, store_ind*sizeof(*store_ary));
			group_ind_ary[group_ind]=store_ind;
			group_ind++;
		
			if(group_ind>=MAX_GROUP_NUM-1)
			{
				fprintf(stderr,"Error! group overflow, group # =%d\n",group_ind);
				exit(1);
			}
		}
    }
    for(i=Vert_ind*2-1; i>=0;  i--)
		if(table_ary[i]!=-1)printf("Error! table[%d]=%d\n",i,table_ary[i]);
		free(table_ary); table_ary=NULL;
		free(Processed); Processed=NULL;
		return group_ind;
}


int trace_untiled_contours(VertType *vert_ary, int complete_ind, 
						   int unused_ind, int vert_ind, 
						   short **group_ary, short *group_ind_ary)
						   /* input line_ary[] & vert_ary[], the output is in group_ary[][] */
						   /* for vertice in [0]-[complete_ind-1], they are completed isolated contours */
{
    int j,i;
    char *flag_ary,*flag_ary_bak;
    int ind, group_ind=0;
    double  x,y,z,x1,y1,z1;
    short *short_ary;
    int is_printed=0,sharing_pt=0,ind_bak;
    int nei_ind2cnt=0;
	
	
    if(vert_ind==0)return 0;
    if(vert_ind<=2)
	{
		if(!SILENT)printf("I. Error! in trace_untiled_contours. too few line seg.\n");
		return 0;
    }
	
    /* malloc for extra one byte */
    flag_ary=(char *)mycalloc(sizeof(*flag_ary)*(vert_ind+1));
    flag_ary_bak=(char *)mycalloc(sizeof(*flag_ary_bak)*(vert_ind+1));
    memset(flag_ary,0,sizeof(*flag_ary)*(vert_ind+1));
	
    /* process the complete contour */
    ind=0;
    while(ind<complete_ind-2)
	{
        int start_ind=ind;
        ind++;
        while(memcmp((char*)(&vert_ary[start_ind]),(char*)(&vert_ary[ind]),
			sizeof(VertType))) {
			ind++;
			if(ind>complete_ind-2)
				break;
		}

		if(ind>=complete_ind)
		{
			if(!SILENT)fprintf(stderr,"I. Err, in trace_unt.. ind=%d > %d\n",
				ind, complete_ind);
			exit(1);
		}
		
		j=ind-start_ind;
        short_ary = (short *)mycalloc((j+4)*sizeof(*short_ary));
        if(group_ary[group_ind] !=NULL)
		{
			free( group_ary[group_ind]); 
			group_ary[group_ind]=NULL;
		}
        
		group_ary[group_ind] = short_ary; 
		
		for(i=start_ind; i<ind; i++)short_ary[i-start_ind]=i;
		
		short_ary[j]=short_ary[0];
		group_ind_ary[group_ind]=j;
		group_ind++;
		ind++;
    }
	
    ind=complete_ind;
    
	while(ind<vert_ind)
	{
#define ARY_SIZE2 4000
		short store_ary[ARY_SIZE2],store_ind;
		int p2,old_pt; 
		int not_done=1,next_pt,start_pt;
		int neighbor[5],nei_ind;
		/* skip the used points */

		while(flag_ary[ind] && ind<vert_ind)ind+=2;
		
		if(ind>=vert_ind)break;
		
		memcpy(flag_ary_bak,flag_ary,vert_ind);
		next_pt=ind;
		start_pt=ind;
		store_ind=0;
		
		while(not_done)
		{
			if(store_ind>=ARY_SIZE2-1)
			{
				printf("Error! mem. over trace_untiled_contours\n"); 
				exit(1);
			}
			
			p2=next_pt+1;
			store_ary[store_ind++]=next_pt;
			old_pt=next_pt;
			/* mark the line to be used */
			if(next_pt!=start_pt)flag_ary[next_pt]=1; /* not mark start pt */
			flag_ary[p2]=1;
			
			/* find the connected point */
			x=vert_ary[p2][0]-0.00001;
			y=vert_ary[p2][1]-0.00001;
			z=vert_ary[p2][2]-0.00001;
			x1=x+0.00002;
			y1=y+0.00002;
			z1=z+0.00002;
			nei_ind=0;
			
			for(i=complete_ind; i<vert_ind; i=i+2)
			{
				/* look for same point */
				if(flag_ary[i]==0 && /* not used yet */
					vert_ary[i][0]>x && vert_ary[i][0] < x1 &&
					vert_ary[i][1]>y && vert_ary[i][1] < y1
					&& vert_ary[i][2]>z && vert_ary[i][2] < z1){
					neighbor[nei_ind++]=i;
				}
			}
			
			for(i=0; i<nei_ind; i++)
			{
				if(neighbor[i]==start_pt)not_done=0; /* done */
			}
			
			if(!not_done)break;
			
			if(nei_ind>2)
				if(!SILENT)
					fprintf(stderr,"Err, nei_ind=%d, should at most 2\n",nei_ind);
				if(nei_ind==2 && !nei_ind2cnt)
				{
					if(!SILENT)
						fprintf(stderr,"nei_ind=2, (%d %d)\n",neighbor[0],neighbor[1]);
					store_ind=0; /* ignore this group */
					not_done=0;
					/* restart from the sharing point to resolve where to go  */
					ind_bak=ind;
					ind=neighbor[0];
					sharing_pt=1;
					memcpy(flag_ary,flag_ary_bak,vert_ind);
					nei_ind2cnt++; /*prevent inifite loop */
				}
				else if(nei_ind==1){ /* keep going */
					next_pt=neighbor[0];
				}
				else 
				{
					/* error ignore this group */
					if(!SILENT)printf("I. Error! nei_ind=%d\n",nei_ind);
					if(!SILENT)printf("x,y,z=[%lf %lf %lf]\n",vert_ary[p2][0],
						vert_ary[p2][1],vert_ary[p2][2]);
					flag_ary[start_pt]=1;
					if(!is_printed && !SILENT)
					{
						for(i=0; i<vert_ind; i++)
						{
							printf("[%d], vertice=%lf %lf %lf ",i,
								vert_ary[i][0], vert_ary[i][1], vert_ary[i][2]);
							if(flag_ary[i]) printf("USED");
							if(i>=vert_ind)printf("BACKUP\n");
							else printf("\n");
						}
						is_printed=1;
					}
					store_ind=0; /* ignore this group */
					not_done=0;
				}
		}

        if(store_ind>=3)
		{
            flag_ary[start_pt]=1;
            short_ary = (short *)mycalloc( (store_ind+4)*sizeof(*short_ary));
            if(group_ary[group_ind] !=NULL)
			{
				free( group_ary[group_ind]); 
				group_ary[group_ind]=NULL;
			} 
            
			group_ary[group_ind] = short_ary; 
            memcpy(short_ary,store_ary, store_ind*sizeof(*store_ary));
			short_ary[store_ind]=short_ary[0];
			group_ind_ary[group_ind]=store_ind;
			group_ind++;
			
			if(group_ind>=MAX_GROUP_NUM-1)
			{
				fprintf(stderr,"Error! group overf, group # =%d\n",group_ind);
				exit(1);
			}
			
			if(sharing_pt)ind=ind_bak;
			sharing_pt=0;
			nei_ind2cnt=0;
		}
        else
		{
			if(store_ind)
			{
				if(!SILENT)fprintf(stderr,"I. error, store_ind=%d\n",store_ind);
				flag_ary[start_pt]=1; /* ignore it */
			}
		}
    }

    if(vert_ind && !SILENT)printf("untiled vert_ind=%d \n",vert_ind);
	
	/*
    printf("\n print all stack \n");
    for(i=0; i<vert_ind; i++){
	printf("[%d], vertice=%lf %lf %lf ",i,
	vert_ary[i][0], vert_ary[i][1], vert_ary[i][2]);
	if(flag_ary[i]) printf("USED ");
	if(i>=vert_ind)printf("BACKUP\n");
	else printf("\n");
    }
	*/

    free(flag_ary);
    free(flag_ary_bak);
    return group_ind;
}

