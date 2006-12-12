#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include "common.h"
#include "myutil.h"
#include "math_util.h"
#include "tile_hash.h"
#include "qsort.h"
#include "legal_region.h"
#include "tile_util.h"

#define UNTILED_ARY_SIZE 2000
extern int SECURITY_CHECK,SILENT;
extern int DEBUG;

extern int Dimx,Dimy,Dimz;

static int SAVE_TILE=0;
static BoundaryType **Untiled_ary[2];
static int Untiled_ary_index[2]={0,0};

static VertType *Lines_ary[2];
static short **New_group_ary[2],*New_group_ind_ary[2];
extern TileType **Best_tiling_table[2];
static int Group_num[2];
extern short *Group_connection_table;
extern int STATIC_CHQ;

void my_clear_tile_util()
{
	SAVE_TILE=0;
	Untiled_ary[0] = NULL;			Untiled_ary[1] = NULL; 
	Untiled_ary_index[0] = NULL;	Untiled_ary_index[1] = NULL; 
	Lines_ary[0] = NULL;			Lines_ary[1] = NULL;
	New_group_ary[0] = NULL;		New_group_ary[1] = NULL; 
	New_group_ind_ary[0] = NULL;	New_group_ind_ary[1] = NULL; 
	Group_num[0] =0;				Group_num[1] =0; 

	save_best_tile(NULL, NULL, NULL,NULL, STATIC_CHQ,NULL, NULL);
	
}

void free_group_ary(char **p)
{
    int i=0;
    if(p==NULL)return;

    while(p[i])
	{
        free(p[i]);
        p[i]=NULL;
		i++;
    }
}


void pass_tiling_parameters(VertType *lines_ary0, VertType *lines_ary1,
							short ** new_group_ary0, short *new_group_ind_ary0, int group_num0,
							short ** new_group_ary1, short *new_group_ind_ary1, int group_num1)
{
    New_group_ary[0]=new_group_ary0;
    New_group_ary[1]=new_group_ary1;
    New_group_ind_ary[0]=new_group_ind_ary0;
    New_group_ind_ary[1]=new_group_ind_ary1;
    Lines_ary[0]=lines_ary0;
    Lines_ary[1]=lines_ary1;
    Group_num[0]=group_num0;
    Group_num[1]=group_num1;
    Untiled_ary_index[0]=0;
    Untiled_ary_index[1]=0;
}


int push_boundary_array(int mode,short group0, short ind0, char tile_dir,
						short group1, short ind1)
{
    int i,empty_ind=-1;
    BoundaryType *temp;
    BoundaryType **untiled_ary;
    
	if(Untiled_ary[mode]==NULL)
	{
		Untiled_ary[mode]=
			(BoundaryType **)mycalloc(sizeof(BoundaryType *)*UNTILED_ARY_SIZE);
    }
	
    untiled_ary=Untiled_ary[mode];
    
	for(i=0; i<Untiled_ary_index[mode]; i++)
	{
		if(untiled_ary[i])
		{
			if(group0==untiled_ary[i]->group0 && 
				group1==untiled_ary[i]->group1 && 
				ind0==untiled_ary[i]->ind0 && 
				ind1==untiled_ary[i]->ind1) 
			{ /* it is already inside */ 
				if(tile_dir==3)return 0;
				if(tile_dir+untiled_ary[i]->tile_dir==3)
				{
					/* It comes with opposite direction, it is done */
					/* delete this line from stack */
					free(untiled_ary[i]);
					untiled_ary[i]=NULL; /* marked as skipped */
					add_an_entrance_to_tile_hash(
						New_group_ary[0][group0][ind0],
						New_group_ary[1][group1][ind1]);
					return 0;
				}
				
				if(tile_dir!=untiled_ary[i]->tile_dir &&
					untiled_ary[i]->tile_dir!=3)
					
					printf("**** I. Error, push_b._array (%d %d %d, %d %d %d),(%d %d)\n",
					group0, ind0, tile_dir,group1, ind1, tile_dir,untiled_ary[i]->tile_dir,
					untiled_ary[i]->tile_dir); 
				
				untiled_ary[i]->tile_dir=tile_dir;
				return 0;
			}
		}
		else /* untiled_ary[i]==NULL */
			if(empty_ind==-1)empty_ind=i; 
    }

    if(empty_ind<0)
	{
		empty_ind=Untiled_ary_index[mode];
		Untiled_ary_index[mode]++;
    }
    
	if(Untiled_ary_index[mode]>= UNTILED_ARY_SIZE)
	{
		fprintf(stderr,"overflow in push_boundary_array()\n"); 
		exit(1);
    }
    
	temp = (BoundaryType*)mycalloc(sizeof(BoundaryType));
    temp->group0=group0;
    temp->group1=group1;
    temp->ind0=ind0;
    temp->ind1=ind1;
    temp->tile_dir=tile_dir;
    untiled_ary[empty_ind]=temp;
    return 1;
}

int pop_boundary_array(int mode, int *group0, int *ind0, int *tile_dir,
					   int *group1,int *ind1)
{
    BoundaryType *temp;
    BoundaryType **untiled_ary;
	
    untiled_ary=Untiled_ary[mode];
    
	while(Untiled_ary_index[mode]>0)
	{
		Untiled_ary_index[mode]--;
		temp=untiled_ary[Untiled_ary_index[mode]];
	
		if(temp)
		{
			*group0=temp->group0;
			*group1=temp->group1;
			*ind0=temp->ind0;
			*ind1=temp->ind1;
			*tile_dir=temp->tile_dir;
			free(temp);
			untiled_ary[Untiled_ary_index[mode]]=NULL;
			return 1;
		}
		untiled_ary[Untiled_ary_index[mode]]=NULL;
    }
    return 0;
}

int find_tiling_direction_from_array(int mode, int group0, int group1, 
									 int ind0, int ind1)
{
    int j,res;
    BoundaryType *temp;
    BoundaryType **untiled_ary;
	
    untiled_ary=Untiled_ary[mode];
    j=0;
    
	while(j<Untiled_ary_index[mode])
	{
		temp=untiled_ary[j];
		if((int)temp>1)
		{
			if( group0==temp->group0 &&	group1==temp->group1 &&	ind0  ==temp->ind0 &&
					ind1  ==temp->ind1)
			{
				res=temp->tile_dir;
				free(temp);
				untiled_ary[j]=NULL;
				return res;
			}
		}
		j++;
    }
    return 3;
}


void security_check_untiled_stack(int mode)
{
    int i;
    BoundaryType **untiled_ary;
	
    if(Untiled_ary[mode]==NULL)
	{
	/*
	fprintf(stderr,"I. Err. security..(), Untiled_ary[%d]=NULL\n", mode);
		*/
		return;
    }
    
	if(Untiled_ary_index[mode]!=0)
		printf("*** Error. security..(), Untiled_ary_index[%d]=%d\n",
		mode,Untiled_ary_index[mode]);
    untiled_ary=Untiled_ary[mode];
    
	for(i=0; i<UNTILED_ARY_SIZE; i++)
	{
		if(untiled_ary[i])
		{
			printf("Error, free_untiled_stack(%d), un.._ary[%d]=%X\n",
				mode,i,untiled_ary[i]);
			if((int)untiled_ary[i]>1)free(untiled_ary[i]);
		}
    }
}

int is_pt_in_LS(VertType *p1, TileType* tile)
/* return 0,1,2: for contour, LS, RS */
{
    double x,y,val[2];
    int i,cnt=0;
    x=(*p1)[0];
    y=(*p1)[1];
	
    cnt=0;
    
	for(i=0; i<2; i++)
	{
		val[i] = tile->a[i]* x + tile->b[i]* y + tile->c[i];
		if(val[i] >-0.01)cnt++;
    }
	
    if(tile->and_or==1 )
	{ /* AND */
		if(cnt==2)
		{
			if(val[0]<=0.01 || val[1]<=0.01)return 0; /* on contour */
			else return 1; /* LS */
		}
		else return 2;     /* RS */
    }
    else if(tile->and_or==0 )
	{ /* OR */
		if(!cnt)return 2; /* RS */
		else if(val[0]<=0.01 && val[0]>=-0.01)return 0; /* on contour */
		else if(val[1]<=0.01 && val[1]>=-0.01)return 0; /* on contour */
		return 1; /* LS */
    }
    return 1;
}

int is_triangle_legal(VertType *p0, VertType *p1, VertType *p2, 
					  TileType *tile0, TileType *tile1)
					  /* p0,p1 must on the same contour, and be ordered 
					  check if *pt is on correct side of (p0,p1) */
					  /* tile1 is for *p1 */
{
    int in_LS;
    double val;
	
    /* security check */
	if(SECURITY_CHECK)
	{
		int res;
		double a,b,c,da,db,dc;
		
		if((*p0)[2]!=(*p1)[2])
		{
			printf("*** Error! is_triangle_legal(), p0[2]!=p1[2] (%lf != %lf)\n",
				(*p0)[2],(*p1)[2]);
			return 0;
		}
		
		res=calc_line_equation(p0, p1,&a,&b,&c);
		if(!res) 
			fprintf(stderr,"Same point: in is_triangle_legal\n"); 
		da=tile1->a[0]-a;
		db=tile1->b[0]-b;
		dc=tile1->c[0]-c;
		val=da*da+db*db+dc*dc;
		if(val>0.01)
			fprintf(stderr,"val=%lf in is_triangle_legal()\n",val);
	}
	
    val = tile1->a[0]* (*p2)[0] + tile1->b[0]* (*p2)[1] + tile1->c[0];
    if(val>0.0001)in_LS=1;
    else if(val<-0.0001)in_LS=2;
    else return 1;
   
	if(tile0->nec_polarity==0)
	{ 
		/* negative nec, should not in RS */
		if(in_LS==2)return 0;
    }
    else if(tile0->nec_polarity==1)
	{
		/* positive nec, should not in LS */
		if(in_LS==1)return 0;
    }
    else if(tile1->nec_polarity==0)
	{
		/* negative nec, should not in RS */
		if(in_LS==2)return 0;
    }
    else if(tile1->nec_polarity==1)
	{
		/* positive nec, should not in LS */
		if(in_LS==1)return 0;
    }

    return 1;
}
int is_line_legal(int mode, VertType *p0, VertType *p1, 
				  TileType *tile0, TileType *tile1)
				  /* p0 is on top slice */
				  /* if mode==1, then check the tiling crossing only */
{
    TileType *tile2;
    int i,j,k,debug;
    int in_LS[2];
    if(mode)goto check_tile_only;
	/*
    if(debug_check(p0,131.0, 124.798) && debug_check(p1, 134.0, 122.906))
	debug_pt();
	*/
    /* return 0,1,2: for contour, LS, RS */
    in_LS[0]=is_pt_in_LS(p1, tile0);
    
	if(tile0->nec_polarity==0)
	{ 
		/* negative nec, should be not in RS */
		if(in_LS[0]==2)return 0;
    }
    else if(tile0->nec_polarity==1)
	{
		/* positive nec, should not be in LS */
		if(in_LS[0]==1)return 0;
    }
    else if(tile0->nec_polarity==2)
	{
		/* overlapping pt */
        int group,ind; /* find its overlapping point */
        group=tile0->group;
        ind=tile0->ind;
        
		if(group>=0) 
		{
            tile2=&(Best_tiling_table[1][group][ind]);
			in_LS[1]=is_pt_in_LS(p1, tile2);
			if(in_LS[0]==in_LS[1] && in_LS[0])return 0; /* both on RS, or LS */
		}
		else printf("***Error is_line_legal(0) group=%d\n"); 
    }
	
    in_LS[0]=is_pt_in_LS(p0, tile1);

    if(tile1->nec_polarity==0)
	{
		/* negative nec, should be not in RS */
		if(in_LS[0]==2)return 0;
    }
    else if(tile1->nec_polarity==1)
	{
		/* positive nec, should not be in LS */
		if(in_LS[0]==1)return 0;
    }
    else if(tile1->nec_polarity==2)
	{
		/* overlapping pt */
        int group,ind; /* find its overlapping point */
        group=tile1->group;
        ind=tile1->ind;
        
		if(group>=0) 
		{
            tile2=&(Best_tiling_table[0][group][ind]);
			in_LS[1]=is_pt_in_LS(p0, tile2);
			if(in_LS[0]==in_LS[1] && in_LS[0])return 0; /* both on RS, or LS */
		}
		else printf("***Error is_line_legal(1) group=%d\n"); 
    }
	
	
    for(i=0; i<2; i++)
	{
		short **group_ary,*group_ind_ary;
		VertType *lines_ary;
		group_ary=New_group_ary[i];
		group_ind_ary=New_group_ind_ary[i];
		lines_ary=Lines_ary[i];
		
		for(j=0; j<Group_num[i]; j++)
		{
			for(k=0; k<group_ind_ary[j]; k++)
			{
				int k1,j0,j1;
				k1=k+1;
			
				if(k1>=group_ind_ary[j])k1=0;
				j0=group_ary[j][k];
				j1=group_ary[j][k1];
				
				if(find_intersection(p0, p1, &(lines_ary[j0]),
					&(lines_ary[j1]),NULL,NULL,NULL)==1)return 0;
			}
		}
    }
check_tile_only:;

		debug=0;
		/*
		if(Untiled_ary_index[0] || Untiled_ary_index[1]){
		debug_check(p0,p1);
		debug=1;
		}
		*/
		for(mode=0; mode<2; mode++)
		{
			BoundaryType *temp;
			BoundaryType **untiled_ary;
			untiled_ary=Untiled_ary[mode];
			for(i=0; i<Untiled_ary_index[mode]; i++)
			{
				temp=untiled_ary[i];
				if(temp)
				{
					/* not empty */
					int res;
					VertType *p2,*p3;
					short group0,group1,ind0,ind1;
					group0=temp->group0;
					group1=temp->group1;
					ind0=temp->ind0;
					ind1=temp->ind1;
					p2=&Lines_ary[0][New_group_ary[0][group0][ind0]];
					p3=&Lines_ary[1][New_group_ary[1][group1][ind1]];
					/*
					if((int)(*p0)[0]==114 && (int)(*p0)[1]==103 &&
					(int)(*p1)[0]==117 && (int)(*p1)[1]==102 &&
					(int)(*p2)[0]==118 && (int)(*p2)[1]==103 &&
					(int)(*p3)[0]==116 && (int)(*p3)[1]==101) 
					debug_pt();
					if((int)(*p2)[0]==114 && (int)(*p2)[1]==103 &&
					(int)(*p3)[0]==117 && (int)(*p3)[1]==102 &&
					(int)(*p0)[0]==118 && (int)(*p0)[1]==103 &&
					(int)(*p1)[0]==116 && (int)(*p1)[1]==101) 
					debug_pt();
					*/
					res=find_intersection(p0, p1, p2, p3, NULL,NULL,NULL);
					
					if(res==1 || res==2) return 0; /* intersection */
					
					if(res==5)
					{ /* overlapping, check it in 3D */
						res=find_intersection_3D(p0, p1, p2, p3, NULL,NULL);
						if(res==1)return 0; /* 3 is overlapping line */
					}
					/* else if(res==3){}  intersect at both end points */
				}
			}
		}
				return 1;
}





void find_best_tiling_group_to_group(int topdown,
									 VertType *lines_ary0, VertType *lines_ary1,
									 short *new_group_ary0, short ary_ind0, int group0,
									 short *new_group_ary1, short ary_ind1, int group1,
									 TileType *best_tiling_tab0, TileType *best_tiling_tab1)
{
#define ARY_SIZE1 2000
    int i,j,k;
    int i1,j1;
    int result,not_shortest_cnt=0;
    double min_dist;
    double dist_ary[ARY_SIZE1];
    int dist_ind_ary[ARY_SIZE1];
    /* use the stupidest search */

    for(j=0; j<ary_ind0; j++)
	{
		min_dist=best_tiling_tab0[j].dist;
		if(min_dist>0.000001)
		{
			if(ary_ind1>=ARY_SIZE1)
			{
				fprintf(stderr,"dist_ary[2000] overflow, ary_ind1=%d\n",ary_ind1);
				printf("dist_ary[2000] overflow, ary_ind1=%d\n",ary_ind1);
				exit(1);
			}
			
			j1=new_group_ary0[j];
			
			for(i=0; i<ary_ind1; i++)
			{
				i1=new_group_ary1[i];
				dist_ary[i]=line_distance_2d(&lines_ary0[j1], &lines_ary1[i1]);
				dist_ind_ary[i]=i;
			}
			
			quick_sort(dist_ary, dist_ind_ary, ary_ind1);
			/*
			if(j==65 && group0==1 && group1==1)
			debug_pt();
			*/
			i=0;
			
			while(i<ary_ind1)
			{ /* try the shortest first */
				k=dist_ind_ary[i]; /* the index of shortest */
				i1=new_group_ary1[k];
				if(dist_ary[k]>=min_dist)break;
			
				if(dist_ary[k]<=0.000001) 
				{ /* overlapping point */
					if(!(best_tiling_tab0[j].and_or&2)) /* overlapping */
						if(!SILENT)printf("**** Warning find_best_tiling_group_to_group(), dist=%.10lf topdown=%d, %d %d\n",
							dist_ary[k], topdown,group0,j); 
						result=1;
				}
				else 
				{
					if(topdown)
						result=is_line_legal(0, &lines_ary1[i1], &lines_ary0[j1],
						&(best_tiling_tab1[k]),
						&(best_tiling_tab0[j]));
					else 
						result=is_line_legal(0, &lines_ary0[j1], &lines_ary1[i1],
						&(best_tiling_tab0[j]),
						&(best_tiling_tab1[k]));
				}
				if(result) 
				{
					min_dist=dist_ary[k];
					best_tiling_tab0[j].dist=(float) min_dist;
					best_tiling_tab0[j].group=group1;
					best_tiling_tab0[j].ind=k;
					if(i>0)not_shortest_cnt++;
					break;
				}
				i++;
			}
		}
    }
	/*
    printf("ind0=%d ind1=%d, %d pts are not shortest\n",
	ary_ind0,ary_ind1,not_shortest_cnt);
	*/
}

void save_best_tile(char *s, VertType *lines_ary0, VertType *lines_ary1,
					short * group_ary, short num, short ** other_group_ary, 
					TileType *tile_table)
{
	
    int i1,k1;
    static int cnt=0;
    char name[60];
    FILE *ofp;
    sprintf(name,"%s%d.poly",s,cnt++);

    if((ofp=fopen(name,"w"))==NULL) 
	{
        fprintf(stderr,"Could not open %s for writing\n",name);
        exit(1);
    }
    
	for(i1=0; i1<num; i1++)
	{
		int k2,ind,group;
		k1=group_ary[i1];
	
		if(tile_table[i1].group>=0)
		{
			group=tile_table[i1].group;
			ind=tile_table[i1].ind;
			k2=other_group_ary[group][ind];
			draw_line_seg(ofp,lines_ary0[k1][0],
				lines_ary0[k1][1], lines_ary0[k1][2],
				lines_ary1[k2][0], lines_ary1[k2][1], lines_ary1[k2][2]);
		}
    }
    fclose(ofp);
}


int calc_best_tiling_table(
						   VertType *lines_ary0, VertType *lines_ary1,
						   short ** new_group_ary0, short *new_group_ind_ary0, int group_num0,
						   short ** new_group_ary1, short *new_group_ind_ary1, int group_num1,
						   short *group_connection_table) 
{
    int i,j,k;
	
	
    make_group_legal_equation_tables(
		lines_ary0, lines_ary1,
		new_group_ary0, new_group_ind_ary0, group_num0,
		new_group_ary1, new_group_ind_ary1, group_num1,
		Best_tiling_table);
	
    for(j=0; j<group_num0; j++)
	{
		for(i=0; i<group_num1; i++)
		{
			k=i*group_num0+j;
		
			if(group_connection_table[k])
			{
				/* It is connected */
				find_best_tiling_group_to_group(0,
					lines_ary0, lines_ary1, 
					new_group_ary0[j], new_group_ind_ary0[j], j,
					new_group_ary1[i], new_group_ind_ary1[i], i,
					Best_tiling_table[0][j],Best_tiling_table[1][i]);
			}
		}
		if(SAVE_TILE) save_best_tile("btileup",lines_ary0,lines_ary1,
			new_group_ary0[j], new_group_ind_ary0[j], new_group_ary1, 
			Best_tiling_table[0][j]);
    }
	
    for(j=0; j<group_num1; j++)
	{
		for(i=0; i<group_num0; i++)
		{
			k=j*group_num0+i;
			if(group_connection_table[k]) 
			{
				/* It is connected */
				find_best_tiling_group_to_group(1,
					lines_ary1, lines_ary0, 
					new_group_ary1[j], new_group_ind_ary1[j], j,
					new_group_ary0[i], new_group_ind_ary0[i], i,
					Best_tiling_table[1][j],Best_tiling_table[0][i]);
			}
		}

		if(SAVE_TILE) save_best_tile("btiledw",lines_ary1,lines_ary0,
			new_group_ary1[j], new_group_ind_ary1[j], new_group_ary0, 
			Best_tiling_table[1][j]);
    }
    return 1;	    
}

int find_start_tile_point(int topdown, int group, int index,
						  int *group0, int *group1, int *ind0, int *ind1)
						  /* if mode==0; then only search for best pair */
{
    int pre_group=-100, pre_pre_group=-100;
    int pre_ind,cnt=0;
    int next_group=-1, next_ind;
	
    topdown=topdown&1; /* either 1 or 0 */
	
	while(1)
	{
		if(++cnt > 100)
		{
			printf("**** Error, loop in find_start_tile_point()\n");
			return 0;
		}
		pre_pre_group=pre_group;
		pre_group=next_group;
		pre_ind=next_ind;
		
		if(Best_tiling_table[topdown][group][index].used==1)
		{
			int ind;
			ind=index-1;
			if(ind<0)ind=New_group_ind_ary[topdown][group]-1;
			/* check if both end line seg are used */
			if(Best_tiling_table[topdown][group][ind].used==1) return 0;
		}
		
		next_group=Best_tiling_table[topdown][group][index].group;
		
		if(next_group<0)return 0;
		
		next_ind=Best_tiling_table[topdown][group][index].ind;
		
		if(next_ind<0)return 0;
		
		if(pre_pre_group==next_group) 
		{
			/* it goes back to the same group */
			if(topdown)
			{
				/* the current is on the second */
				*group0=next_group;
				*group1=pre_group;
				*ind0=next_ind;
				*ind1=pre_ind;
			}
			else
			{
				*group1=next_group;
				*group0=pre_group;
				*ind1=next_ind;
				*ind0=pre_ind;
			}
			return 1;
		}
		topdown=1-topdown;
		group=next_group;
		index=next_ind;;
	}
}


int find_best_start_tile_point(int group0, int ind0,
							   int *ret_group1, int *ret_ind1)
{
    int group1,ind1,group2,ind2;
	
    group1=Best_tiling_table[0][group0][ind0].group;
    ind1=Best_tiling_table[0][group0][ind0].ind;

    if(group1<0 || ind1<0)return 0;

    group2=Best_tiling_table[1][group1][ind1].group;
    ind2=Best_tiling_table[1][group1][ind1].ind;
    
	if(group0==group2 && ind0==ind2)
	{
		/* they are best pair */
		*ret_group1=group1;
		*ret_ind1=ind1;
		return 1;
    }
    return 0;
}



Linklist *record_break_contour(int topdown, int group, int ind0)
{
#define BREAK_ARY_SIZE 100
    extern int NUMERICALLY_STABLE; 
    int ind1,topdown1;
    int j,j1,k,k1,i,i1,cnt=0,kkk;
    int group1;
    double ratio_ary[BREAK_ARY_SIZE],pre_ratio;
    double ratio_ary1[BREAK_ARY_SIZE];
    Linklist *linklist;
	
    topdown1=1-topdown;
    ind1=ind0+1;
    if(ind1>=New_group_ind_ary[topdown][group])ind1=0;
    j =New_group_ary[topdown][group][ind0];
    j1=New_group_ary[topdown][group][ind1];

    for(group1=0; group1<Group_num[topdown1]; group1++)
	{
		if(topdown)	kkk=group*Group_num[0]+group1;
		else	 kkk=group1*Group_num[0]+group;
		
		/*
		if(Group_connection_table[kkk]){
		*/

		for(i=0; i<New_group_ind_ary[topdown1][group1]; i++)
		{
			int res;
			i1=i+1;
			if(i1>=New_group_ind_ary[topdown1][group1])i1=0;
			k=New_group_ary[topdown1][group1][i];
			k1=New_group_ary[topdown1][group1][i1];
			res=find_intersection(&(Lines_ary[topdown][j]),
				&(Lines_ary[topdown][j1]), &(Lines_ary[topdown1][k]),
				&(Lines_ary[topdown1][k1]), NULL, &(ratio_ary[cnt]),NULL);
		
			if(res==1 || res==2 || res==5)
			{ 
				cnt++; 
				if(cnt>=BREAK_ARY_SIZE)
				{
					fprintf(stderr,"over BREAK_ARY_SIZE\n");
					exit(1);
				}
			}
		}
		/*
		}
		*/
    }
    if(!cnt)return NULL;
    bubble_sort(ratio_ary,cnt); /* bubble sort is slow, but N is small */
								/* the break point is between the intersection point, and the intersection
	points themselves */
    /* get rid of the duplicated point */
    j=0;
    i=0;

    while(i<cnt && (ratio_ary[i]<=0.0001 || ratio_ary[i]>=0.9999) )i++;
    
	if(i==cnt)return NULL;
    
	pre_ratio=ratio_ary[i];
    ratio_ary1[j++]=pre_ratio;
    i++;
    
	for(; i<cnt; i++)
	{
	/*	if(ratio_ary[i]>(pre_ratio+0.001) && ratio_ary[i] < 0.999) { 
	if(ratio_ary[i] < 0.9999) { 
		*/
		/* different point */
		ratio_ary1[j++]=ratio_ary[i];
		pre_ratio=ratio_ary[i];
		/*
		}
		*/
    }
    cnt=j;
    if(!cnt)return NULL;
	
    if(NUMERICALLY_STABLE)
	{
		if(cnt<=1)return NULL;
		j=0;
		
		for(i=1; i<cnt; i++)
		{
			ratio_ary1[j++]=(ratio_ary1[i-1]+ratio_ary1[i])/2.0;
		}
		
		cnt=j;
    }
	
    linklist=(Linklist*)mycalloc(sizeof(Linklist));
    linklist->index=ind0;
    linklist->ary_ind=cnt;
    linklist->ratio_ary=(double*)mycalloc(sizeof(double)*cnt);
    
	for(i=0; i<cnt; i++)	linklist->ratio_ary[i]=ratio_ary1[i];
    
	linklist->next=NULL;
    
	return linklist;
}

int break_contour_if_necessary( VertType *lines_ary0, VertType *lines_ary1,
							   short ** new_group_ary0, short *new_group_ind_ary0, int group_num0,
							   short ** new_group_ary1, short *new_group_ind_ary1, int group_num1,
							   int *lines_ind0, int *lines_ind1, 
							   int lines_ary_size0, int lines_ary_size1)
							   
							   /* if a contour line segment has more than one intersection with other contour
							   in other slice, it needs to be broken up. So legal tiling become possible */
{
    Linklist **group_link_list[2];
    int i,i1,group;
    int topdown;
    Linklist *linklist,**pre_pointer,*nextlist;
    int local_lines_ind[2];
    int local_lines_ary_size[2];
	
    New_group_ary[0]=new_group_ary0;
    New_group_ary[1]=new_group_ary1;
    New_group_ind_ary[0]=new_group_ind_ary0;
    New_group_ind_ary[1]=new_group_ind_ary1;
    Lines_ary[0]=lines_ary0;
    Lines_ary[1]=lines_ary1;
    Group_num[0]=group_num0;
    Group_num[1]=group_num1;
    local_lines_ind[0]=*lines_ind0;
    local_lines_ind[1]=*lines_ind1;
    local_lines_ary_size[0]=lines_ary_size0;
    local_lines_ary_size[1]=lines_ary_size1;
	
    /* order the linklist, so it is as the contour order */
    for(topdown=0; topdown<2; topdown++)
	{
		group_link_list[topdown]=(Linklist**)
			mycalloc(sizeof(Linklist**)*Group_num[topdown]);
	
		for(group=0; group<Group_num[topdown]; group++)
		{
			pre_pointer=&(group_link_list[topdown][group]);
		
			for(i=0; i<New_group_ind_ary[topdown][group]; i++)
			{
				linklist=record_break_contour(topdown,group,i);
				if(linklist!=NULL)
				{
					*pre_pointer=linklist;
					pre_pointer=&(linklist->next);
				}
			}
		}
    }
    /* begin to break contour segment */
    for(topdown=0; topdown<2; topdown++)
	{
		for(group=0; group<Group_num[topdown]; group++)
		{
			linklist=group_link_list[topdown][group];
			
			if(linklist!=NULL)
			{
				/* some break in this group */
				/* first count how many break so we can allocate mem */
				int cnt=0,new_index,ii,num;
				short *temp_group_ary, *group_tab;
				VertType *lines_tab;
				nextlist=linklist;
				
				while(nextlist!=NULL)
				{
					cnt+=nextlist->ary_ind;
					nextlist=nextlist->next;
				}
				
				new_index=New_group_ind_ary[topdown][group]+cnt;
				if(DEBUG)
					printf("New_group_ind_ary[%d][%d] change from %d to %d\n",
					topdown, group, New_group_ind_ary[topdown][group], new_index);
				temp_group_ary=(short*)mycalloc(sizeof(short)*(new_index+1));
				ii=0;
				group_tab=New_group_ary[topdown][group];
				num=New_group_ind_ary[topdown][group];
				lines_tab=Lines_ary[topdown];
				
				for(i=0; i<New_group_ind_ary[topdown][group]; i++)
				{
					temp_group_ary[ii++]=group_tab[i];
					
					if(linklist!=NULL && i==linklist->index)
					{
						/* break it */
						double x0,y0,x1,y1,z,dx,dy;
						int jj0,jj1,i9;
						i1=i+1;
						if(i1>=num)i1=0;
						jj0=group_tab[i];
						jj1=group_tab[i1];
						x0=lines_tab[jj0][0];
						y0=lines_tab[jj0][1];
						z=lines_tab[jj0][2];
						x1=lines_tab[jj1][0];
						y1=lines_tab[jj1][1];
						dx=x1-x0;
						dy=y1-y0;
						
						for(i9=0; i9<linklist->ary_ind; i9++)
						{
							lines_tab[local_lines_ind[topdown]][0]=
								x0+(linklist->ratio_ary)[i9]*dx;
							lines_tab[local_lines_ind[topdown]][1]=
								y0+(linklist->ratio_ary)[i9]*dy;
							lines_tab[local_lines_ind[topdown]][2]=z;
							temp_group_ary[ii++]= local_lines_ind[topdown];
							local_lines_ind[topdown]++;
							if(local_lines_ind[topdown]>=
								local_lines_ary_size[topdown]) 
							{
								fprintf(stderr,"I. error 109, Lines_ary overflow\n");
								printf("I. error 109, Lines_ary overflow\n");
								exit(1);
							}
						}

						free(linklist->ratio_ary);
						linklist->ratio_ary=NULL;
						nextlist = linklist;
						linklist=linklist->next;
						free(nextlist);
						nextlist=NULL;
					}
				}
				if(ii!=new_index)
					printf("**** Error, should %d, only %d\n",
					new_index, ii);
				New_group_ind_ary[topdown][group]=new_index;
				free(New_group_ary[topdown][group]);
				New_group_ary[topdown][group]=temp_group_ary;
				New_group_ary[topdown][group][new_index]=
					New_group_ary[topdown][group][0];
			}
		}
		free(group_link_list[topdown]);
		group_link_list[topdown]=NULL;
    }
    *lines_ind0 = local_lines_ind[0];
    *lines_ind1 = local_lines_ind[1];
    return 1;
}






