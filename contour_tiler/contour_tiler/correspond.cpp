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



#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <malloc.h>

#include "common.h"
#include "contour.h"
#include "myutil.h"
#include "math_util.h"
#include "correspond.h"

extern int DEBUG,SILENT;
extern short Hierachy[2][MAX_GROUP_NUM];
static VertType *Lines_ary[2];
static short **New_group_ary[2],*New_group_ind_ary[2];
static int Group_num[2];
extern short *Group_relation_table;


extern short *Group_connection_table;
extern signed char *Is_inside_positive[2];
extern TileType **Best_tiling_table[2];

extern short *Crossed_group[2];
extern short *Crossed_index[2];
extern int Crossed_num;
extern int NUMERICALLY_STABLE;

extern int STATIC_CHQ;

void my_clear_correspond()
{
	Lines_ary[0] = NULL;			Lines_ary[1] = NULL; 
	New_group_ary[0] = NULL;		New_group_ary[1] = NULL;
	New_group_ind_ary[0] = NULL;	New_group_ind_ary[1] = NULL;
	Group_num[0] = 0;				Group_num[1] = 0;
}


void print_group_relation_table(int group_num0, int group_num1)
{
    int i,j;

    printf("group_relation_table\n ");
    
	for(j=0; j<group_num0; j++) 
	{
        if(!(j%5))printf(" ");
        if(!(j%10))printf(" ");
        printf("%d",j%10);
    }
    
	printf("\n");
    
	for(j=0; j<group_num1; j++)
	{
        printf("%d",j%10);
        
		for(i=0; i<group_num0; i++) 
		{
            if(!(i%5))printf(" ");
            if(!(i%10))printf(" ");
            printf("%d",Group_relation_table[j*group_num0+i]);
        }
        
		printf("\n");
    }
}


void build_polarity_table()
{
	
    int i,topdown,k;
    int not_done;

    for(topdown=0; topdown<2; topdown++)
	{
		for(i=Group_num[topdown]; i>=0; i--) Is_inside_positive[topdown][i]=-1;
    }
    
	for(topdown=0; topdown<2; topdown++)
	{
		not_done=1;
		while(not_done)
		{
			not_done=0;
			for(i=Group_num[topdown]-1; i>=0; i--)
			{
				if(Is_inside_positive[topdown][i]==-1)
				{
					not_done=1;
					k=Hierachy[topdown][i];
					
					if(k==-1) /* the top level */
						Is_inside_positive[topdown][i]=1; /* positive */
					else 
					{ 
						/* someone's child */
						if(Is_inside_positive[topdown][k]>=0)
						{ /* parent done */
							Is_inside_positive[topdown][i]=
								1- Is_inside_positive[topdown][k];
						}
					}
				}
			}
		}
    }
}

void build_hierachy() /* the parent-child relation of contours on a slice */
{
    int topdown,res;
    int j,k,j1;
    VertType *pt;
    double dist;

    for(topdown=0; topdown<2; topdown++)
	{
		for(j=0; j<Group_num[topdown]; j++)
		{
			double min_dist=1000000.0;
			int nec_ind=-1;
			k=New_group_ary[topdown][j][0];
			pt=&(Lines_ary[topdown][k]);
		
			for(j1=0; j1<Group_num[topdown]; j1++)
			{
				if(j!=j1)
				{
					res=is_inside_contour(pt, Lines_ary[topdown], 
						New_group_ary[topdown][j1], New_group_ind_ary[topdown][j1],
						&dist,NULL);
					if(res==1)
					{
						if(dist<min_dist)
						{
							min_dist=dist;
							nec_ind=j1;
						}
					}
				}
				
			}
			Hierachy[topdown][j]=nec_ind;
		}
    }
}

int check_if_two_contours_intersect(VertType *lines_ary0, VertType *lines_ary1,
									short *group_tab0, short num0, short *group_tab1, int num1, 
									int index_ary_size, short *index_ary0, short *index_ary1, int *ret_num)
{
    int i,j,k1,res;
    int num=0, intersect=0;
    double ratio1,ratio2;
    VertType *p0,*p1,*p2,*p3;

    for(j=0; j<num0; j++)
	{
		p0=&(lines_ary0[group_tab0[j]]);
		p1=&(lines_ary0[group_tab0[j+1]]);
	
		for(i=0; i<num1; i++)
		{
			k1=group_tab1[i];
			p2=&(lines_ary1[group_tab1[i]]);
			p3=&(lines_ary1[group_tab1[i+1]]);
			res=find_intersection(p0, p1, p2, p3, NULL,&ratio1,&ratio2);
		
			if( res == 5 && ratio1 > 0.0001 && ratio2>0.0001)
			{
				/* overlapped */
				if(num<index_ary_size)
				{
					if(index_ary0) index_ary0[num]=j;
					if(index_ary1) index_ary1[num]=i;
					num++;
				}
			}
			else if(res==1 ) 
			{
				/* 1: crossed, 5: overlapped */
				if(ratio1<=0.9999 && ratio1>=0.0001 && 
					ratio2<=0.9999 && ratio2>=0.0001){
					if(num<index_ary_size)
					{
						if(index_ary0)index_ary0[num]=j;
						if(index_ary1)index_ary1[num]=i;
						num++;
					}
				}
			}
			if(res)intersect=1;
		}
    }
    if(ret_num)*ret_num=num;
    return intersect;
}

void check_intersection_all_contours()
{
    int i,j,k,j1,num,res;
    short index_ary[2][4000];
    int index_ary_size=4000;
    char hist[2][MAX_LINES_NUM];
	
	
#define TEMP_BUF_SIZE 4000
    
	if(NUMERICALLY_STABLE)
	{
		for(i=0; i<2; i++)
		{
			Crossed_group[i]=(short*)mymalloc(sizeof(short)*TEMP_BUF_SIZE);
			Crossed_index[i]=(short*)mymalloc(sizeof(short)*TEMP_BUF_SIZE);
		}
    }
	
	
    for(j=0; j<Group_num[0]; j++)
	{
		for(j1=0; j1<Group_num[1]; j1++)
		{
			res=check_if_two_contours_intersect(Lines_ary[0], Lines_ary[1],
				New_group_ary[0][j], New_group_ind_ary[0][j],New_group_ary[1][j1], 
				New_group_ind_ary[1][j1],index_ary_size, index_ary[0], 
				index_ary[1],&num);

			k=j1*Group_num[0]+j;
			
			if(res)		Group_relation_table[k]=1;
			else	Group_relation_table[k]=0;

			if(NUMERICALLY_STABLE)
			{
				memset(hist[0],0,MAX_LINES_NUM);
				memset(hist[1],0,MAX_LINES_NUM);
			
				for(i=0; i<num; i++)
				{
					int k1,k2;
					for(k1=0; k1<2; k1++)
					{
						k2=index_ary[k1][i];
						if(hist[k1][k2])
							fprintf(stderr,"G[%d][%d %d], index=%d, intersected >1, with %d\n", 
							k1,j,j1,k2,index_ary[1-k1][i]);
						hist[k1][k2]++;
					}
				}

				for(i=0; i<num; i++)
				{
					Crossed_group[0][Crossed_num+i]=j;
					Crossed_group[1][Crossed_num+i]=j1;
					Crossed_index[0][Crossed_num+i]=index_ary[0][i];
					Crossed_index[1][Crossed_num+i]=index_ary[1][i];
				}
				Crossed_num+=num;
			}
		}
    }

    if(NUMERICALLY_STABLE)
	{
		short *temp;
		for(i=0; i<2; i++)
		{
			temp=(short*)mymalloc(sizeof(short)*(Crossed_num+1));
			memcpy(temp, Crossed_group[i], sizeof(short)*Crossed_num);
			free(Crossed_group[i]);
			Crossed_group[i]=temp;
			
			temp=(short*)mymalloc(sizeof(short)*(Crossed_num+1));
			memcpy(temp, Crossed_index[i], sizeof(short)*Crossed_num);
			free(Crossed_index[i]);
			Crossed_index[i]=temp;
		}
    }
}

void check_projection_table()
{
    int topdown;
    int i,j,i2;
    TileType *best_tiling_tab;
    signed char por1,por0,por2,pre_por=-1;
    int ibak,jbak;
	
    for(topdown=0; topdown<2; topdown++)
	{
		for(j=0; j<Group_num[topdown]; j++)
		{
			best_tiling_tab = Best_tiling_table[topdown][j];
			for(i=0; i<New_group_ind_ary[topdown][j]; i++)
			{
				int i1;
				i1=i+1;
				if(i1>=New_group_ind_ary[topdown][j])i1=0;
				por0=best_tiling_tab[i].nec_polarity;
				por1=best_tiling_tab[i1].nec_polarity;
			
				if(por0!=2 && por1!=2 && por0!=por1)
				{
					if(!SILENT)
					{
						printf("*** Warning in check_projection_table(), topdown=%d, g%d p[%d]=%d != p[%d]=%d\n",
							topdown,j,i,por0,i1,por1);
						printf("It is caused by numerical instability. dist0=%lf, dist1=%lf\n",
							best_tiling_tab[i].dist,best_tiling_tab[i1].dist);
					}
					if(pre_por==0 || pre_por==1)por2=pre_por;
					else 
					{
						i2=i1+1;
						if(i2>=New_group_ind_ary[topdown][j])i2=0;
						por2=best_tiling_tab[i2].nec_polarity;
					}
					
					if(por2<2)
					{
						if(por2==por0)
						{ /* por1 has problem */
							ibak=i1;
							jbak=i;
						}
						else 
						{
							ibak=i;
							jbak=i1;
						}
						if(best_tiling_tab[ibak].dist < 0.0001)
						{
							if(!SILENT)printf("This warning remedied. Don't worry the result\n");
							best_tiling_tab[ibak].nec_polarity=
								best_tiling_tab[jbak].nec_polarity;
							best_tiling_tab[ibak].nec_group=
								best_tiling_tab[jbak].nec_group;
						}
					}
				}
				pre_por=por0;
			}
		}
    }
}



void build_projection_table()
{
    int topdown,topdown1;
    int i,j,res,j1;
    int min_group,ret_ind;
    double min_dist,d;
    VertType *pt;
    TileType *best_tiling_tab;

	
    for(topdown=0; topdown<2; topdown++)
	{
		topdown1=1-topdown;
		for(j=0; j<Group_num[topdown]; j++)
		{
			best_tiling_tab = Best_tiling_table[topdown][j];
			for(i=New_group_ind_ary[topdown][j]-1; i>=0; i--)
			{
				min_group=-1;
				min_dist=100000000;
				pt=&(Lines_ary[topdown][New_group_ary[topdown][j][i]]);
				for(j1=Group_num[topdown1]-1; j1>=0; j1--)
				{
				/*
				if(topdown==0 && j==8  && i==2 && j1==10)
				debug_pt();
					*/
					res=is_inside_contour(pt, Lines_ary[topdown1], 
						New_group_ary[topdown1][j1], 
						New_group_ind_ary[topdown1][j1],
						&d,&ret_ind);
					if(res==3 || res==2 || d<0.00000001)
					{ 
						/* this point is right on the contour */
						min_group=j1;
						j1=-1; /* quit the loop */
					}
					else if(res==1)
					{
						if(d<min_dist)
						{
							min_dist=d;
							min_group=j1;
						}
					}
                }

				best_tiling_tab[i].nec_group=min_group;
				if(res==3)
				{ /* overlapping vertex */
					best_tiling_tab[i].group=min_group;
					best_tiling_tab[i].ind=ret_ind;
					best_tiling_tab[i].nec_polarity=2;
					best_tiling_tab[i].dist=0;
				}
				else if(min_group>=0)
				{
					best_tiling_tab[i].nec_polarity=
						Is_inside_positive[1-topdown][min_group];
				}
				else /* the slice boundary */
					best_tiling_tab[i].nec_polarity=0;
			}
			memcpy(&(best_tiling_tab[New_group_ind_ary[topdown][j]]),
				&(best_tiling_tab[0]), sizeof(*best_tiling_tab));
		}
    }
}

int have_a_negative_projection_vertex(int topdown, int group)
{
    int i,sign;
    TileType *best_tiling_tab;
    
	sign=Is_inside_positive[topdown][group];
    best_tiling_tab = Best_tiling_table[topdown][group];

    for(i=New_group_ind_ary[topdown][group]-1; i>=0; i--)
	{
		if(best_tiling_tab[i].nec_polarity==sign)return 1;
    }
    return 0;
}

int have_a_positive_projection_vertex(int topdown, int group)
{
    int i,sign;
    TileType *best_tiling_tab;
    
	best_tiling_tab = Best_tiling_table[topdown][group];
    sign=Is_inside_positive[topdown][group];
    
	for(i=New_group_ind_ary[topdown][group]-1; i>=0; i--)
	{
		if(best_tiling_tab[i].nec_polarity!=sign &&
			best_tiling_tab[i].nec_polarity!=2)return 1;
    }
    return 0;
}

int have_a_vertex_nec_equal_to_another_contour(int topdown, int g0, int g1)
{
    int i;
    unsigned char g1_ch;
    TileType *best_tiling_tab;
    best_tiling_tab = Best_tiling_table[topdown][g0];
	
    if(g1>=0)g1_ch=g1;
    else g1_ch=-1; /* for comparsion */
	
    for(i=New_group_ind_ary[topdown][g0]-1; i>=0; i--)
	{
		if(best_tiling_tab[i].nec_group==g1_ch && 
			best_tiling_tab[i].nec_polarity!=2 )return 1;
    }
    return 0;
}

int check_condition2(int g0, int g1)
{
    int nec_g1,nec_g0;

    if(Is_inside_positive[0][g0]==
		Is_inside_positive[1][g1])return 0; /* same sign */

    if(!have_a_negative_projection_vertex(0,g0))return 0;

    if(!have_a_negative_projection_vertex(1,g1))return 0;

    nec_g0=Hierachy[0][g0];
    nec_g1=Hierachy[1][g1];
    
	if(!have_a_vertex_nec_equal_to_another_contour(0,g0,nec_g1))return 0;

    if(!have_a_vertex_nec_equal_to_another_contour(1,g1,nec_g0))return 0;

    return 1;
}

int is_the_nec_of_different_slice(int topdown,int g0,int g1)
/* if g1=NEC(g0') */
{
    int i,ind,ind0,topdown1;
    topdown1=1-topdown;
    
	if(topdown)	ind0=g0*Group_num[0]+g1;
    else	 ind0=g1*Group_num[0]+g0;

    if(!topdown && Group_relation_table[ind0]!=2)return 0;

    if(topdown && Group_relation_table[ind0]!=3)return 0;

    for(i=0; i<Group_num[topdown1]; i++)
	{
		if(topdown)	ind=g0*Group_num[0]+i;
		else	 ind=i*Group_num[0]+g0;

		/* because Hierachy[] is dependent on topdown, so separate it */
		if(topdown==1 && Group_relation_table[ind]==3 && ind!=ind0)
		{
			/* check it has an even smaller one contain g0' */
			if(Hierachy[topdown1][i]==g1)return 0;
		}

		if(topdown==0 && Group_relation_table[ind]==2 && ind!=ind0)
		{
			/* check it has an even smaller one contain g0' */
			if(Hierachy[topdown1][i]==g1)return 0;
		}
    }
    return 1;
}

int check_condition3(int topdown, int g0, int g1)
/* g0 of slice(topdown) is inside g1 */
{
    int nec_g0;
    int topdown1=1-topdown;

    if(Is_inside_positive[topdown][g0]!=Is_inside_positive[topdown1][g1])
		return 0; /* different sign */
    
	nec_g0=Hierachy[topdown][g0];
    
	if(!have_a_vertex_nec_equal_to_another_contour(topdown1,g1,nec_g0))return 0;
    
	if(!is_the_nec_of_different_slice(topdown,g0,g1))return 0;
    
	return 1;
}


void    make_all_group_counter_clockwise()
{
    int i,j;
    
	for(i=0; i<Group_num[0]; i++)
	{
		j=make_contour_counter_clockwise(Is_inside_positive[0][i],Lines_ary[0], 
				New_group_ary[0][i],New_group_ind_ary[0][i]);
		if(j&&DEBUG)
			printf("top slice, g%d, seq. reversed, is_positive=%d\n",i,
				Is_inside_positive[0][i]);
    }

    for(i=0; i<Group_num[1]; i++)
	{
		j=make_contour_counter_clockwise(Is_inside_positive[1][i],Lines_ary[1], 
			New_group_ary[1][i],New_group_ind_ary[1][i]);
		if(j&&DEBUG)
			printf("bot slice, g%d, seq. reversed, is_positive=%d\n",i,
			Is_inside_positive[1][i]);
    }
}


void build_contours_relation(VertType *lines_ary0, VertType *lines_ary1,
							 short ** group_ary0, short *group_ind_ary0, int group_num0,
							 short ** group_ary1, short *group_ind_ary1, int group_num1)
{
    int j,k,j1,ind;
    VertType *pt,*pt1;
    New_group_ary[0]=group_ary0;
    New_group_ary[1]=group_ary1;
    New_group_ind_ary[0]=group_ind_ary0;
    New_group_ind_ary[1]=group_ind_ary1;
    Lines_ary[0]=lines_ary0;
    Lines_ary[1]=lines_ary1;
    Group_num[0]=group_num0;
    Group_num[1]=group_num1;
	
    build_hierachy();
    build_polarity_table();
    make_all_group_counter_clockwise();
    check_intersection_all_contours();
    build_projection_table();

    for(j=0; j<Group_num[0]; j++)
	{
		for(j1=0; j1<Group_num[1]; j1++)
		{
			ind=j1*Group_num[0]+j;
			if(!Group_relation_table[ind])
			{
				k=New_group_ary[0][j][0];
				pt=&(Lines_ary[0][k]);
				k=New_group_ary[1][j1][0];
				pt1=&(Lines_ary[1][k]);
			
				if(is_inside_contour(pt, Lines_ary[1], 
					New_group_ary[1][j1], New_group_ind_ary[1][j1],
					NULL,NULL)) Group_relation_table[ind]=2; /* a < b */
				else if(is_inside_contour(pt1, Lines_ary[0], 
					New_group_ary[0][j], New_group_ind_ary[0][j],
					NULL,NULL)) Group_relation_table[ind]=3; /* a > b */
				else Group_relation_table[ind]=0;
			}
			else Group_connection_table[ind]=1;
		}
    }
    if(DEBUG)print_group_relation_table(Group_num[0], Group_num[1]);
	
    for(j=0; j<Group_num[0]; j++)
	{
		for(j1=0; j1<Group_num[1]; j1++)
		{
			ind=j1*Group_num[0]+j;
			if(!Group_connection_table[ind])
			{
				if(Group_relation_table[ind]==0) /* disjoint */
					Group_connection_table[ind]=check_condition2(j,j1);
				else if(Group_relation_table[ind]==2) /* a < b */
					Group_connection_table[ind]=check_condition3(0,j,j1);
				else if(Group_relation_table[ind]==3) /* a > b */
					Group_connection_table[ind]=check_condition3(1,j1,j);
			}
		}
    }
}

void build_contours_levels(VertType *lines_ary[2], 
						   short ** group_ary[2], short *group_ind_ary[2], int group_num[2],
						   short *level_ary[2])
{
    int j,k,j1; //,ind; //KLC
    //VertType *pt, *pt1; //KLC

    for(j=0; j<2; j++)
	{
		New_group_ary[j]=group_ary[j];
		New_group_ind_ary[j]=group_ind_ary[j];
		Lines_ary[j]=lines_ary[j];
		Group_num[j]=group_num[j];
    }
	
    build_hierachy();

    for(j=0; j<2; j++)
	{
		for(k=0; k<group_num[j]; k++)
		{
			int cnt=0;
			j1=k;
			while( Hierachy[j][j1]>=0)
			{
				cnt++;
				j1=Hierachy[j][j1]; /* its parents */
			}

			level_ary[j][k]=cnt;
		}
    }
}

