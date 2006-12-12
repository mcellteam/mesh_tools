#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
#include <math.h>
#include <fcntl.h>
#include "ct/ct.h"
#include "common.h"
#include "math_util.h"
#include "tile_hash.h"
#include "tile_util.h"
#include "contour_read.h"
#include "myutil.h"
#include "tiling_main.h"
#define UNTILED_ARY_SIZE 2000

extern int DEBUG,SILENT;
extern int Dimx,Dimy,Dimz;
extern short *Group_connection_table;

extern FILE* Debug_fd;
extern FILE* Mesh_fp;
int Optimum_tiling_mode;

extern TileType **Best_tiling_table[2];
static int Tile_triangle_num; /* record the tiled number */
void close(int);

static int no_of_poly_triangles =0;

extern int STATIC_CHQ;

void my_clear_tiling()
{
	Tile_triangle_num =0;

	check_and_process_shortest_regular_triangle(NULL, NULL, NULL, STATIC_CHQ, STATIC_CHQ,
												STATIC_CHQ, STATIC_CHQ, STATIC_CHQ,
												STATIC_CHQ, NULL, NULL,NULL, NULL, NULL, NULL);
	do_partial_tiling(STATIC_CHQ, STATIC_CHQ, STATIC_CHQ, STATIC_CHQ, STATIC_CHQ, NULL,
					  NULL, NULL, STATIC_CHQ, NULL, STATIC_CHQ);
	christiansen(NULL, NULL, NULL, NULL, STATIC_CHQ, NULL, NULL, STATIC_CHQ);
	
}

void debug_pt()
{
	if(!SILENT)fprintf(stderr,"debug, in shortest()\n");
}


void draw_tiled_triangle_sub(VertType *p1, VertType *p2, VertType *p3, int mode)
{
	extern int DIFF_LEVEL_FILES;
	extern int Current_level;
	extern FILE* TiledFdAry[];
	extern FILE* Tiled_fd;
	extern int NO_OUTPUT_FILE;
	int i,j,k;
	FILE* fd;
	extern VertType *volumeVertAry;
	extern int volumeVertAryInd;

    if(NO_OUTPUT_FILE)return ;

    if(DIFF_LEVEL_FILES)
	{
		if(TiledFdAry[Current_level] !=0)
		{
            TiledFdAry[Current_level]=openATiled(Current_level);
		}
		fd=TiledFdAry[Current_level];
    }
    else
	{
		fd=Tiled_fd;
    }

    if(volumeVertAry!=NULL)
	{ 
		/* write the middle line for volume calculate */
		VertType *p[3];
		int cnt=0;
		p[0]=p1; p[1]=p2; p[2]=p3; 

		for(i=0; i<3; i++)
		{
			j=(i+1)%3;
			if( ((*p[i])[2]<= ((*p[j])[2] -1e-4)) || ((*p[i])[2]>= ((*p[j])[2] +1e-4)))
			{ 
				/* yes they are different */
				for(k=0; k<3; k++)
				{
					volumeVertAry[volumeVertAryInd][k] =((*p[i])[k]+(*p[j])[k])/2.0;
				}
				
				volumeVertAryInd++;
				cnt++;
			}
		}
		
		if(cnt!=2)printf("Error! save the middle section, cnt=%d\n",cnt);
    }
			
    if(fd == 0)return;

    if(Current_level &1)
	{ 
		/* need to reverse the sequence*/
		draw_triangle(fd,p2,p1,p3);
		if(Debug_fd !=0)draw_triangle(Debug_fd,p2,p1,p3);
		if(mode && Mesh_fp)draw_triangle_FP(Mesh_fp,p2,p1,p3);
		no_of_poly_triangles++;
    }
    else 
	{
		draw_triangle(fd,p1,p2,p3);
		if(Debug_fd !=0 )draw_triangle(Debug_fd,p1,p2,p3);
		if(mode && Mesh_fp)draw_triangle_FP(Mesh_fp,p1,p2,p3);
		no_of_poly_triangles++;
    }
}
	    



void draw_tiled_triangle(VertType *p1, VertType *p2, VertType *p3)
{
    draw_tiled_triangle_sub(p1, p2, p3,1);
}

void draw_tiled_triangle_no_mesh(VertType *p1, VertType *p2, VertType *p3)
{
    draw_tiled_triangle_sub(p1, p2, p3,0);
}

int check_and_process_optimum_triangle_pair(VertType *lines_ary0, VertType *lines_ary1,
											int group0, int group1, int ind0, int ind1,
											short group_num0, short group_num1,
											short *group_ary0, short *group_ary1,
											TileType *best_tiling_table0, 
											TileType *best_tiling_table1,int *ret_pt0,
											int *ret_pt1,int *ret_edge_ind0, 
											int *ret_edge_ind1)

{
    int next_ind0,next_ind1,edge_ind0,edge_ind1;
    int best_ind0,best_group0;
    int best_ind1,best_group1;

    if(is_in_tile_hash_table(group_ary0[ind0], group_ary1[ind1]))	return 0;

    next_ind0=ind0+1;
    next_ind1=ind1+1;
    
	if(next_ind0>=group_num0)	next_ind0=0;

    if(next_ind1>=group_num1)	next_ind1=0;

    if(is_in_tile_hash_table(group_ary0[next_ind0], group_ary1[next_ind1]))
		return 0;

    edge_ind0=ind0;
    edge_ind1=ind1;

    /* check if the top line segment is used */ 
    if(best_tiling_table0[edge_ind0].used) /* used */
	    return 0;
    
	if(best_tiling_table1[edge_ind1].used) /* used */
	    return 0;
    
	best_ind0=best_tiling_table0[next_ind0].ind;
    best_group0=best_tiling_table0[next_ind0].group;
    best_ind1=best_tiling_table1[next_ind1].ind;
    best_group1=best_tiling_table1[next_ind1].group;

    if(best_ind1==next_ind0 && best_group1==group0 && best_ind0==next_ind1 && best_group0==group1)
	{
		/* The line is absolutely legal */
		*ret_pt0=next_ind0;
		*ret_pt1=next_ind1;
		*ret_edge_ind0=edge_ind0;
		*ret_edge_ind1=edge_ind1;
		return 1;
    }
    return 0;
}
    

int check_and_process_optimum_triangle(int topdown,VertType *lines_ary0, 
									   VertType *lines_ary1, int dir, int group0, 
									   int group1, int ind0, int ind1, short group_num0,
									   short *group_ary0, short *group_ary1,
									   TileType *best_tiling_table0, 
									   TileType *best_tiling_table1,int *new_pt)
{
    int next_ind0,edge_ind;
    int done=0;
    int best_ind1,best_group1;

    if(dir==1)
	{
		next_ind0=ind0+1;
		if(next_ind0>=group_num0)next_ind0=0;
		edge_ind=ind0;
    }
    else 
	{
		next_ind0=ind0-1;
		if(next_ind0<0)next_ind0=group_num0-1;
		edge_ind=next_ind0;
    }

    /* check if the top line segment is used */ 
    if(best_tiling_table0[edge_ind].used) /* used */
	    return 0;
    
	best_ind1=best_tiling_table1[ind1].ind;
    best_group1=best_tiling_table1[ind1].group;

    if(best_ind1==next_ind0 && best_group1==group0) /* an optimum triangle */
		done=1;
    else
	{
		int best_next_ind0, best_next_group0;
		best_next_ind0=best_tiling_table0[next_ind0].ind;
		best_next_group0=best_tiling_table0[next_ind0].group;
		if(ind1==best_next_ind0 && best_next_group0 == group1)done=1; 
    }

    if(done)
	{
		int res;
		/* check if the new tiling segment used */

		if(topdown)
			res=is_in_tile_hash_table(group_ary1[ind1], group_ary0[next_ind0]);
		else 
			res=is_in_tile_hash_table(group_ary0[next_ind0], group_ary1[ind1]);		
		if(res)return 0;

		if(dir==1)
			res=is_triangle_legal( &(lines_ary0[group_ary0[ind0]]), 
			&(lines_ary0[group_ary0[next_ind0]]),
					&(lines_ary1[group_ary1[ind1]]), 
			&best_tiling_table0[ind0],&best_tiling_table0[next_ind0]);
		else 
			res=is_triangle_legal( &(lines_ary0[group_ary0[next_ind0]]), 
			&(lines_ary0[group_ary0[ind0]]),
					&(lines_ary1[group_ary1[ind1]]), 
			&best_tiling_table0[next_ind0],&best_tiling_table0[ind0]);
		if(!res)return 0;
		
		if(topdown)
			res=is_line_legal(1,&(lines_ary1[group_ary1[ind1]]), 
			&(lines_ary0[group_ary0[next_ind0]]), NULL,NULL);
		else 
			res=is_line_legal(1,&(lines_ary0[group_ary0[next_ind0]]), 
			&(lines_ary1[group_ary1[ind1]]), NULL,NULL);

		if(res)
		{
			*new_pt=next_ind0;
			return 1;
		}
    }
    return 0;
}
   

int check_and_process_regular_triangle(int topdown,	VertType *lines_ary0, 
									   VertType *lines_ary1,int dir, int group0, 
									   int group1, int ind0, int ind1, short group_num0,
									   short *group_ary0, short *group_ary1,
									   TileType *best_tiling_table0,
									   TileType *best_tiling_table1,double *dist,
									   int *edge_ret)
{
    int next_ind0,edge_ind,best_next_group0;
    int res,best_group1,best_group0;
    int p0,p1,p2;
    int different_group=0;

    if(dir==1)
	{
		next_ind0=ind0+1;
		if(next_ind0>=group_num0)next_ind0=0;
		edge_ind=ind0;
    }
    else 
	{
		next_ind0=ind0-1;
		if(next_ind0<0)next_ind0=group_num0-1;
		edge_ind=next_ind0;
    }

    best_group1=best_tiling_table1[ind1].group;
    best_group0=best_tiling_table0[ind0].group;
    /* check if the top line segment is used */ 
    
	if(best_tiling_table0[edge_ind].used)
	{
		/* used */
		*dist=100000000.0;
		return -1;
	}

    best_next_group0=best_tiling_table0[next_ind0].group;

    if(best_next_group0 != group1)	different_group=1;

    if(Optimum_tiling_mode<2 && different_group)
	{
		*dist=100000000.0;
		return -1;
    }

    p0=group_ary0[ind0];
    p1=group_ary0[next_ind0];
    p2=group_ary1[ind1];

    if(topdown)
	    res=is_in_tile_hash_table(group_ary1[ind1], group_ary0[next_ind0]);
    else 
	    res=is_in_tile_hash_table(group_ary0[next_ind0], group_ary1[ind1]);
    
	if(res)
	{
        *dist=100000000.0;
        return -1;
    }

    if(dir==1)
        res=is_triangle_legal( &(lines_ary0[group_ary0[ind0]]),
                &(lines_ary0[group_ary0[next_ind0]]),
                &(lines_ary1[group_ary1[ind1]]),
                &best_tiling_table0[ind0],&best_tiling_table0[next_ind0]);
    else
        res=is_triangle_legal( &(lines_ary0[group_ary0[next_ind0]]),
                &(lines_ary0[group_ary0[ind0]]),
                &(lines_ary1[group_ary1[ind1]]),
                &best_tiling_table0[next_ind0],&best_tiling_table0[ind0]);
    if(!res)
	{
        *dist=100000000.0;
        return -1;
    }

    if(topdown)
	res=is_line_legal(0,&(lines_ary1[p2]), &(lines_ary0[p1]),
	    &(best_tiling_table1[ind1]),&(best_tiling_table0[next_ind0]));
    else res=is_line_legal(0,&(lines_ary0[p1]), &(lines_ary1[p2]),
	    &(best_tiling_table0[next_ind0]),&(best_tiling_table1[ind1]));
    
	if(res)
	{
		*dist=line_distance_2d(&lines_ary0[p1], &lines_ary1[p2]);
		*edge_ret=edge_ind;
		return next_ind0;
    }
    *dist=10000000;
    return -1;
}

int check_and_process_shortest_regular_triangle(int *tile_dir, VertType *lines_ary0, 
												VertType *lines_ary1,int group0, 
												int group1, int ind0, int ind1,
												short group_num0, short group_num1,
												short *group_tab0, short *group_tab1,
												TileType *top_best_tiling_table, 
												TileType *bot_best_tiling_table,
												int *new_pt0, int *new_pt1)
{
    double min_dist=100000.0, dist;
    int index=-1, topdown,edge_ind,final_edge_ind;
    int i,j;
    short dir;
    static int dir_array[3]={0,1,-1};
    

    for(j=1; j<=2; j++)
	{
		if(*tile_dir &j)
		{
			i=check_and_process_regular_triangle(0,	lines_ary0,lines_ary1,dir_array[j],
				group0, group1, ind0, ind1, group_num0,group_tab0,group_tab1,
				top_best_tiling_table, bot_best_tiling_table,&dist,&edge_ind);

			if(dist<min_dist)
			{
				min_dist=dist;
				index=i;
				topdown=0;
				final_edge_ind=edge_ind;
				dir=j;
			}
		}
    }
	
    for(j=1; j<=2; j++)
	{
		if(*tile_dir &j)
		{
			i=check_and_process_regular_triangle(1,lines_ary1,lines_ary0,dir_array[j],
				group1, group0, ind1, ind0, group_num1,group_tab1,group_tab0,
				bot_best_tiling_table, top_best_tiling_table,&dist,&edge_ind);

			if(dist<min_dist)
			{
				min_dist=dist;
				index=i;
				topdown=1;
				final_edge_ind=edge_ind;
				dir=j;
			}
		}
    }

    if(index<0)	return 0;

    if(topdown)
	{
		/* bottom */
		*new_pt0=ind0;
		*new_pt1=index;
		*tile_dir=dir;
		return 1;
    }
    else 
	{
		*new_pt1=ind1;
		*new_pt0=index;
		*tile_dir=dir;
		return 1;
    }
    return 0;
}
	

int is_both_edge_used(TileType *best_tiling_tab, int ind0, short group_num)
{
    int ind1;
    if(best_tiling_tab[ind0].used==0)return 0;
    ind1=ind0-1;
    if(ind1<0)ind1=group_num-1;
    if(best_tiling_tab[ind1].used==0)return 0;
    return 1;
}

void do_best_tiling(int group0, int group1, int ind0, int ind1, VertType *lines_ary0,
					VertType *lines_ary1, short * group_tab0, short group_num0,
					short * group_tab1, short group_num1)
{
    int res,new_ind0,new_ind1,edge_ind0,edge_ind1,is_top=0;
    TileType *top_best_tiling_tab;
    TileType *bot_best_tiling_tab;
    VertType *p0,*p3,*p1,*p2;
    VertType *pp1,*pp2,*pp3;
    VertType *pp4,*pp5,*pp6;
   
    top_best_tiling_tab=Best_tiling_table[0][group0];
    bot_best_tiling_tab=Best_tiling_table[1][group1];
/*
if(group0==1  && ind0==26 && group1==1 && (ind1==33 || ind1==31)) 
debug_pt();
*/
    res=check_and_process_optimum_triangle_pair(lines_ary0,lines_ary1, group0, group1, 
			ind0, ind1, group_num0,group_num1, group_tab0,group_tab1,top_best_tiling_tab, 
			bot_best_tiling_tab,&new_ind0, &new_ind1,&edge_ind0,&edge_ind1);

    if(!res)return;
    /* arrange the triangle in counter clockwise (CCW)direction */

    /*        p0-------p1
	      |        |
	      p2-------p3
    */	
    
	p0=&(lines_ary0[group_tab0[ind0]]);
    p1=&(lines_ary0[group_tab0[new_ind0]]);
    p2=&(lines_ary1[group_tab1[ind1]]);
    p3=&(lines_ary1[group_tab1[new_ind1]]);
    
	if(line_distance_2d(p0,p3)>line_distance_2d(p1,p2))
	{ 
		/* p1-p2 */
		pp1=p0; pp2=p2; pp3=p1;
		pp4=p2; pp5=p3; pp6=p1;
		is_top=1;
    }
    else 
	{
		pp1=p0; pp2=p2; pp3=p3;
		pp4=p0; pp5=p3; pp6=p1;
		is_top=0;
    }

    draw_tiled_triangle(pp1,pp2,pp3);
    Tile_triangle_num+=2;
    draw_tiled_triangle(pp4,pp5,pp6);
    top_best_tiling_tab[edge_ind0].used=1; /* used */
    bot_best_tiling_tab[edge_ind1].used=1; /* used */
    
	push_boundary_array(Optimum_tiling_mode&1, (short)group0,(short)ind0,2,(short)group1,(short)ind1); /* left boundary */
    push_boundary_array(Optimum_tiling_mode&1, (short)group0,(short)new_ind0,1,(short)group1,(short)new_ind1); /* right boundary */
}

void do_overlapping_pt_tiling(int group0, int group1, int ind0, int ind1,
							  VertType *lines_ary0, VertType *lines_ary1,
							  short * group_tab0, short group_num0,short * group_tab1, 
							  short group_num1)
{
    int new_ind0,new_ind1;
    TileType *top_best_tiling_tab;
    TileType *bot_best_tiling_tab;
    VertType *p0,*p3,*p1,*p2;
    double dist[2];
    int cnt;
/*
if(group0==4  && group1==3) 
debug_pt();
*/

    /* check if this is done */
    if(is_in_tile_hash_table(group_tab0[ind0], group_tab1[ind1]))return;

    top_best_tiling_tab=Best_tiling_table[0][group0];
    bot_best_tiling_tab=Best_tiling_table[1][group1];

    new_ind0=ind0-1;
    new_ind1=ind1-1;

    if(new_ind0<0)new_ind0=group_num0-1;
    if(new_ind1<0)new_ind1=group_num1-1;
    p0=&(lines_ary0[group_tab0[ind0]]);
    p2=&(lines_ary1[group_tab1[ind1]]);
    p1=&(lines_ary0[group_tab0[new_ind0]]);
    p3=&(lines_ary1[group_tab1[new_ind1]]);

    /*        p1-------p0
	      |        |
	      p3-------p2
    */	

    cnt=0;
    
	if(!top_best_tiling_tab[new_ind0].used && !is_in_tile_hash_table(group_tab0[new_ind0], group_tab1[ind1]))
	{
		cnt++;
		dist[0]=line_distance_2d(p1,p2);
    }

    if(!bot_best_tiling_tab[new_ind1].used && !is_in_tile_hash_table(group_tab0[ind0], group_tab1[new_ind1]))
	{
		cnt++;
		dist[1]=line_distance_2d(p0,p3);
    }

    if(cnt==2)
	{
		push_boundary_array(Optimum_tiling_mode&1, (short)group0,(short)ind0,1,(short)group1,(short)ind1);

		Tile_triangle_num++;
		
		if(dist[0]<dist[1])
		{ 
			/* choose the upper triangle */
			draw_tiled_triangle(p0,p1,p2);
			push_boundary_array(Optimum_tiling_mode&1,(short)group0,(short)new_ind0,2,(short)group1,(short)ind1);
			top_best_tiling_tab[new_ind0].used=1; /* used */
		}
		else 
		{
			/* p0-p3 */
			draw_tiled_triangle(p0,p3,p2);
			push_boundary_array(Optimum_tiling_mode&1,(short)group0,(short)ind0,2,(short)group1,(short)new_ind1);
			bot_best_tiling_tab[new_ind1].used=1; /* used */
		}
    }

/*  It might happen
    else if(cnt==1) 
	fprintf(stderr,"Error do_overlapping_pt_tiling()\n");
*/

    new_ind0=ind0+1;
    new_ind1=ind1+1;
    if(new_ind0>=group_num0)	new_ind0=0;
    if(new_ind1>=group_num1)	new_ind1=0;

    p1=&(lines_ary0[group_tab0[new_ind0]]);
    p3=&(lines_ary1[group_tab1[new_ind1]]);

    /*        p0-------p1
	      |        |
	      p2-------p3
    */	

    cnt=0;
    
	if(!top_best_tiling_tab[ind0].used && !is_in_tile_hash_table(group_tab0[new_ind0], group_tab1[ind1]))
	{
		cnt++;
		dist[0]=line_distance_2d(p1,p2);
    }

    if(!bot_best_tiling_tab[ind1].used && !is_in_tile_hash_table(group_tab0[ind0], group_tab1[new_ind1]))
	{
		cnt++;
		dist[1]=line_distance_2d(p0,p3);
    }

    if(cnt==2)
	{
		push_boundary_array(Optimum_tiling_mode&1,(short)group0,(short)ind0,2,(short)group1,(short)ind1);

		Tile_triangle_num++;

		if(dist[0]<dist[1])
		{ /* choose the upper triangle */
			draw_tiled_triangle(p0,p2,p1);
			push_boundary_array(Optimum_tiling_mode&1,
				(short)group0,(short)new_ind0,1,(short)group1,(short)ind1);
			top_best_tiling_tab[ind0].used=1; /* used */
		}
		else 
		{
			/* p0-p3 */
			draw_tiled_triangle(p0,p2,p3);
			push_boundary_array(Optimum_tiling_mode&1,(short)group0,(short)ind0,1,(short)group1,(short)new_ind1);
			bot_best_tiling_tab[ind1].used=1; /* used */
		}
    }
/*
    else if(cnt==1)
	fprintf(stderr,"Error do_overlapping_pt_tiling()\n");
*/
}


int debug_check(VertType *p1, double x, double y)
{
    double xval,yval,zval;
    xval= (*p1)[0];
    yval= (*p1)[1];
    zval= (*p1)[2];
    if(xval<=x+0.01 && xval >= x-0.01 &&
       yval<=y+0.01 && yval >= y-0.01)
	return 1;
    else return 0;
}

void do_partial_tiling(int group0, int group1, int ind0, int ind1, int tile_dir,
					   VertType *lines_ary0, VertType *lines_ary1, short * group_tab0,
					   short group_num0, short * group_tab1, short group_num1)
{
    int i,new_ind,new_ind1,did_one,loop,old_ind0,old_ind1;
    int seq,ind0_bak, ind1_bak,tile_dir_bak,done_num_ary[2];
    int is_break;
    static int dir_array[3]={0,1,-1};
    /* tile_dir: 3 : undefined, 1: increase, 2: decrease */
    TileType *top_best_tiling_table;
    TileType *bot_best_tiling_table;
   
    top_best_tiling_table=Best_tiling_table[0][group0];
    bot_best_tiling_table=Best_tiling_table[1][group1];

    /* tiling can go from both side, so use it twice */
    ind0_bak=ind0;
    ind1_bak=ind1;

	/*
	if(group0==1  && ind0==26 && group1==1 && (ind1==33 || ind1==31)) 
	debug_pt();
	if(group0==8 && (ind0==1 || ind0==2) && group1==10 && ind1==1)
	debug_pt();
	*/
    
	i=is_in_tile_hash_table(group_tab0[ind0], group_tab1[ind1]);

    if(i)
	{
		if(tile_dir==3)	return;
		printf("**** Error, g0=%d ind0=%d, g1=%d, ind1=%d, tile_dir=%d in hash\n",group0,ind0,group1,ind1,tile_dir);
		return;
    }

    /* need to check whether both end are used */
    if(tile_dir==3)
	{
		if(is_both_edge_used(Best_tiling_table[0][group0], ind0, group_num0))
			return;
		if(is_both_edge_used(Best_tiling_table[1][group1], ind1, group_num1))
			return;
    }

    /* need to chek is this line legal */
    if(!Optimum_tiling_mode)
	{
		if(tile_dir!=3)
			printf("*** Error, tile_dir=%d\n",tile_dir);
		tile_dir=find_tiling_direction_from_array(0,group0,group1,ind0,ind1);
    }

    if(tile_dir==3)
	{
		i=is_line_legal(0, &lines_ary0[group_tab0[ind0]], &lines_ary1[group_tab1[ind1]],
			&(top_best_tiling_table[ind0]),	&(bot_best_tiling_table[ind1]));
		  if(!i)return;
    }


	tile_dir_bak=tile_dir;
	if(tile_dir==3)loop=2; /* try both direction */
	else loop=1;
	done_num_ary[0]=0;
	done_num_ary[1]=0;

  for(seq=0; seq<loop; seq++)
  {
    int not_done=1;
    is_break=0;
    did_one=0;

    /* If it is second round tiling, do not reverse its direction */
    if(seq==1)
	{ /* reverse the direction */
		ind0=ind0_bak;
		ind1=ind1_bak;
		tile_dir=3-tile_dir;
    }

    while(not_done)
	{
		int stored_tile_dir,mode;
		not_done=0;
		old_ind0=ind0;
		old_ind1=ind1;
		if(did_one)
		{
			mode=Optimum_tiling_mode &1;
			stored_tile_dir=find_tiling_direction_from_array( mode, group0,group1,
				ind0,ind1);

			if(stored_tile_dir==3)
			{
				mode=1-mode;
				stored_tile_dir=find_tiling_direction_from_array(
					mode, group0,group1,ind0,ind1);
			}
			if(stored_tile_dir+tile_dir==3)
			{
				/* in the opposite direction */
				/* push to cancel */
				/*
				if(push_boundary_array(mode,group0,ind0,tile_dir,group1,ind1))
					fprintf(stderr,"Err, (%d %d) (%d %d) is not canceled\n",
					group0,ind0,group1,ind1);
				*/
				/* mark it is done */
				add_an_entrance_to_tile_hash( group_tab0[ind0], group_tab1[ind1]);
				is_break=1;
				break;
			}
		}
		/*	
		if(group0==1  && ind0==25 && group1==1 && (ind1==33 || ind1==31)) 
		debug_pt();
		if(group0==8 && (ind0==1 || ind0==2) && group1==10 && ind1==1)
		debug_pt();
		*/
			if(Optimum_tiling_mode==0)
			{
				for(i=1; i<=2 && (!not_done); i++)
				{ 
					if(tile_dir & i)
					{
						/* check the triangle formed by the top line seg */
						if(check_and_process_optimum_triangle(0,lines_ary0,lines_ary1,
							dir_array[i], group0, group1, ind0, ind1, group_num0,
							group_tab0,group_tab1,top_best_tiling_table, 
							bot_best_tiling_table,&new_ind))
						{
							tile_dir=i; 
							ind0=new_ind;
							not_done=1;
							did_one++;
						}
					}
				}

			  for(i=1; i<=2 &&(!not_done); i++)
			  { 
				if(tile_dir & i)
				{
				  /* check the triangle formed by the bot line seg */
					  if(check_and_process_optimum_triangle(1,lines_ary1,lines_ary0,
						  dir_array[i], group1, group0, ind1, ind0, group_num1,
						  group_tab1,group_tab0,bot_best_tiling_table, top_best_tiling_table,
						  &new_ind))
					  {
						tile_dir=i; 
						ind1=new_ind;
						not_done=1;
						did_one++;
					  }
				}
			  } 
			} //end of if (Optimum_tiling_mode==0)condition....


		  if(!not_done && Optimum_tiling_mode>0 && 
			  check_and_process_shortest_regular_triangle( &tile_dir,lines_ary0,
			  lines_ary1, group0, group1, ind0, ind1, group_num0, group_num1,
			  group_tab0,group_tab1,top_best_tiling_table, bot_best_tiling_table,
			  &new_ind,&new_ind1))
		  {
			ind0=new_ind;
			ind1=new_ind1;
			not_done=1;
			did_one++;
		  }

      if(not_done)
	  {
		int is_top,edge_ind;
		VertType *p3,*p1,*p2;
		/* arrange the triangle in counter clockwise direction */

		/*      bottom triangle		top triangle    
				 p1			  p3	p2
				 o			  o-----o
				/ \			   \   /
			   /   \		    \ /
			  o-----o		     o
			  p2    p3		     p1
 		*/
		if(ind0==old_ind0 && ind1!=old_ind1)
		{ /* the bottom triangle */
			is_top=0;
			p1=&(lines_ary0[group_tab0[ind0]]);
			if(tile_dir==1)
			{
				p2=&(lines_ary1[group_tab1[old_ind1]]);
				p3=&(lines_ary1[group_tab1[ind1]]);
				edge_ind=old_ind1;
			}
			else 
			{
				p3=&(lines_ary1[group_tab1[old_ind1]]);
				p2=&(lines_ary1[group_tab1[ind1]]);
				edge_ind=ind1;
			}

			/* If tile_dir_bak!=3, i.e. it has tiled before */
			if(did_one==1 && tile_dir_bak!=3)
				add_an_entrance_to_tile_hash(group_tab0[ind0], group_tab1[old_ind1]);
			else if(did_one>1) /* mark this can not be used again */
				add_an_entrance_to_tile_hash(group_tab0[ind0], group_tab1[old_ind1]);
		}

		else if(ind0!=old_ind0 && ind1==old_ind1)
		{ /* the top triangle */
			is_top=1;
			p1=&(lines_ary1[group_tab1[ind1]]);
			
			if(tile_dir==1)
			{
				p2=&(lines_ary0[group_tab0[ind0]]);
				p3=&(lines_ary0[group_tab0[old_ind0]]);
				edge_ind=old_ind0;
			}
			else 
			{
				p3=&(lines_ary0[group_tab0[ind0]]);
				p2=&(lines_ary0[group_tab0[old_ind0]]);
				edge_ind=ind0;
			}

			/* If tile_dir_bak!=3, i.e. it has tiled before */
			if(did_one==1 && tile_dir_bak!=3)
				add_an_entrance_to_tile_hash(group_tab0[old_ind0], group_tab1[ind1]);
			else if(did_one>1) /* mark this can not be used again */
				add_an_entrance_to_tile_hash(group_tab0[old_ind0], group_tab1[ind1]);
		} 
		else
			fprintf(stderr,"*** partial_tiling(), g: %d %d, old %d %d, %d %d\n",
			group0,group1,old_ind0,old_ind1,ind0,ind1);
		
		draw_tiled_triangle(p1,p2,p3);
		Tile_triangle_num++;
		if(is_top)
		{
			top_best_tiling_table[edge_ind].used=1; /* used */
		}
		else 
		{ 
			bot_best_tiling_table[edge_ind].used=1; /* used */
		}
	  }
		
	  done_num_ary[seq]=did_one;
		}

    if(seq==1 && done_num_ary[0] && done_num_ary[1])
	/* both side is done */
		add_an_entrance_to_tile_hash(group_tab0[ind0_bak], group_tab1[ind1_bak]);

    if(ind0==ind0_bak && ind1==ind1_bak && did_one)	return; 
	/*  completed in one direction */
    
	/* always push the best tiling into the stack. It is required for collecting untiled data */
    if(!is_break) push_boundary_array(Optimum_tiling_mode&1,(short) group0,(short)ind0,(char)tile_dir,(short)group1,(short)ind1);
    
	if(tile_dir==3)return; /* It tried both directions, no tiling done */
  }
}

void triangulate_quad( VertType *lines_ary0, VertType *lines_ary1,int group0, int ind0, 
					  int new_ind0, int kk0, int kk1, int group1, int ind1,
					  int new_ind1, int kk2, int kk3,TileType *best_tiling_table0,
					  TileType *best_tiling_table1)
{
    int is_top=0;
    TileType *top_best_tiling_tab;
    TileType *bot_best_tiling_tab;
    VertType *p0,*p3,*p1,*p2;
    VertType *pp1,*pp2,*pp3;
    VertType *pp4,*pp5,*pp6;
    int legal1=0, legal2=0;
  

    top_best_tiling_tab=Best_tiling_table[0][group0];
    bot_best_tiling_tab=Best_tiling_table[1][group1];
    /* arrange the triangle in counter clockwise (CCW)direction */

    /*        p0-------p1
              |        |
              p2-------p3
    */
    p0=&(lines_ary0[kk0]);
    p1=&(lines_ary0[kk1]);
    p2=&(lines_ary1[kk2]);
    p3=&(lines_ary1[kk3]);

    if(line_distance_2d(p0,p3)>line_distance_2d(p1,p2))
	{
		/* p1-p2 */
        pp1=p0; pp2=p2; pp3=p1;
        pp4=p2; pp5=p3; pp6=p1;
        is_top=1;
    }
    else
	{
        pp1=p0; pp2=p2; pp3=p3;
        pp4=p0; pp5=p3; pp6=p1;
        is_top=0;
    }
    
	if(!is_in_tile_hash_table(kk0,kk2) && 
		is_line_legal(0,p0,p2,&(best_tiling_table0[ind0]),&(best_tiling_table1[ind1])))
	{
		legal1=1;
		draw_tiled_triangle(pp1,pp2,pp3);
		Tile_triangle_num++;
    }

    if(!is_in_tile_hash_table(kk1,kk3) && 
		is_line_legal(0,p1,p3,&(best_tiling_table0[new_ind0]),&(best_tiling_table1[new_ind1])))
	{
		legal2=1;
		draw_tiled_triangle(pp4,pp5,pp6);
		Tile_triangle_num++;
		if(!is_top)	 top_best_tiling_tab[ind0].used=1; /* used */
		else  bot_best_tiling_tab[ind1].used=1; /* used */
    }

    if(legal1 && legal2)
	{
		top_best_tiling_tab[ind0].used=1; /* used */
		bot_best_tiling_tab[ind1].used=1; /* used */

		push_boundary_array(Optimum_tiling_mode&1,(short)group0,(short)ind0,2,(short)group1,(short)ind1); /* left boundary */
		push_boundary_array(Optimum_tiling_mode&1,(short)group0,(short)new_ind0,1,(short)group1,(short)new_ind1); /* right boundary */
		if(is_top) /* mark p1-p2 is used */
			add_an_entrance_to_tile_hash(kk1,kk2);
        else  /* mark p0-p3 is used */
	    add_an_entrance_to_tile_hash(kk0,kk3);
	return;
    }

    if(legal1)
	{
		push_boundary_array(Optimum_tiling_mode&1,(short)group0,(short)ind0,2,(short)group1,(short)ind1); /* left boundary */
		if(is_top) {
			top_best_tiling_tab[ind0].used=1; /* used */
			push_boundary_array(Optimum_tiling_mode&1,
			(short)group0,(short)new_ind0,1,(short)group1,(short)ind1); /* right boundary */
		}
		else  
		{
			bot_best_tiling_tab[ind1].used=1; /* used */
			push_boundary_array(Optimum_tiling_mode&1,
			(short)group0,(short)ind0,1,(short)group1,(short)new_ind1); /* right boundary */
		}
    }
    else if(legal2)
	{
		push_boundary_array(Optimum_tiling_mode&1,(short)group0,(short)new_ind0,1,(short)group1,(short)new_ind1); /* left boundary */
		if(is_top) 
		{
			bot_best_tiling_tab[ind1].used=1; /* used */
			push_boundary_array(Optimum_tiling_mode&1,(short)group0,(short)new_ind0,2,(short)group1,(short)ind1); /* right boundary */
		}
		else
		{
			top_best_tiling_tab[ind0].used=1; /* used */
			push_boundary_array(Optimum_tiling_mode&1,(short)group0,(short)ind0,2,(short)group1,(short)new_ind1); /* right boundary */
		}
    }
}

void do_cross_tiling( VertType *lines_ary0, VertType *lines_ary1,
            short ** new_group_ary0, short *new_group_ind_ary0, int group_num0,
            short ** new_group_ary1, short *new_group_ind_ary1, int group_num1)
{
    extern short *Crossed_group[2], *Crossed_index[2];
    extern int Crossed_num;

    int group0,group1,k,k1,l1,l,num;
    int ind0,ind1,new_ind0,new_ind1;
    TileType *best_tiling_table0;
    TileType *best_tiling_table1;

    for(num=0; num<Crossed_num; num++)
	{    
		group0=Crossed_group[0][num];
		group1=Crossed_group[1][num];
		ind0=Crossed_index[0][num];
		ind1=Crossed_index[1][num];
		best_tiling_table0=Best_tiling_table[0][group0];
		best_tiling_table1=Best_tiling_table[1][group1];

		if(best_tiling_table0[ind0].used)
			printf("*** Warning do_cross tiling, top, group=%d, ind=%d, used\n", 
				group0,ind0);
		else if(best_tiling_table1[ind1].used)
			printf("*** Warning do_cross tiling, bottom, group=%d, ind=%d, used\n",
				group1,ind1);
		else 
		{
			new_ind0=ind0+1;
			if(new_ind0>=new_group_ind_ary0[group0])new_ind0=0;
			new_ind1=ind1+1;
			if(new_ind1>=new_group_ind_ary1[group1])new_ind1=0;
			k=new_group_ary0[group0][ind0];
			k1=new_group_ary0[group0][new_ind0];
			l=new_group_ary1[group1][ind1];
			l1=new_group_ary1[group1][new_ind1];
			triangulate_quad(lines_ary0, lines_ary1,
			group0,ind0,new_ind0,k,k1,group1,ind1,new_ind1,l,l1,
			best_tiling_table0, best_tiling_table1);
		}
    }
    return;
}
    
int do_tiling( VertType *lines_ary0, VertType *lines_ary1,
            short ** new_group_ary0, short *new_group_ind_ary0, int group_num0,
            short ** new_group_ary1, short *new_group_ind_ary1, int group_num1)
{
    int i;
    int group,res;
    int group0,group1, ind0,ind1;
    int tile_ind;
    extern int Current_slice_num;
    extern int NUMERICALLY_STABLE;

    Tile_triangle_num=0;
    pass_tiling_parameters(lines_ary0, lines_ary1, new_group_ary0, new_group_ind_ary0,
		group_num0, new_group_ary1, new_group_ind_ary1, group_num1);

	if(!SILENT)fprintf(stderr,"before calc_best_tiling_table\n");
		calc_best_tiling_table( lines_ary0, lines_ary1, new_group_ary0, 
			new_group_ind_ary0, group_num0, new_group_ary1, new_group_ind_ary1,
			group_num1, Group_connection_table);

		 /* check_projection_table(); *//* after the best_tiling_table is done */

		/* There are three passes in tiling:
		 pass 0: Only do the optimum tiling
		 pass 1: regular tiling
		 pass 2: regular & different group tiling 
		 pass 3: ? 
		*/

    if(!SILENT)fprintf(stderr,"begin to do tiling\n");
    Optimum_tiling_mode=0;

    if(NUMERICALLY_STABLE) do_cross_tiling(lines_ary0, lines_ary1,
            new_group_ary0, new_group_ind_ary0, group_num0,
            new_group_ary1, new_group_ind_ary1, group_num1);
    if(Debug_fd>0)
	{
		fclose(Debug_fd);
		Debug_fd=0;
    }

    for(group0=0; group0<group_num0; group0++)
	{
		for(ind0=0; ind0<new_group_ind_ary0[group0]; ind0++)
		{
			res=find_best_start_tile_point(group0,ind0, &group1,&ind1);

			if(res) 
			{ /* found a start point */
				do_best_tiling(group0, group1, ind0, ind1,
						lines_ary0, lines_ary1,	new_group_ary0[group0],
						new_group_ind_ary0[group0],new_group_ary1[group1], 
						new_group_ind_ary1[group1]);
			}
		}
    }

    /* process overlapping point */
    
	for(group0=0; group0<group_num0; group0++)
	{
		TileType *best_tiling_tab;
		best_tiling_tab=Best_tiling_table[0][group0];

		for(ind0=0; ind0<new_group_ind_ary0[group0]; ind0++)
		{
			if(best_tiling_tab[ind0].dist<0.000001)
			{ /* overlapping point */
				group1=best_tiling_tab[ind0].group;
				ind1=best_tiling_tab[ind0].ind;
				if(ind1>=0)
					do_overlapping_pt_tiling(group0, group1, ind0, ind1,lines_ary0, 
					lines_ary1,	new_group_ary0[group0], new_group_ind_ary0[group0],
					new_group_ary1[group1], new_group_ind_ary1[group1]);
			}
		}
    }

    for(group=0; group<group_num0; group++)
	{
		for(i=0; i<new_group_ind_ary0[group]; i++)
		{
			res=find_start_tile_point(0,group,i, &group0,&group1,&ind0,&ind1);
			if(res) 
			{ /* found a start point */
				do_partial_tiling(group0, group1, ind0, ind1,3,lines_ary0, lines_ary1,
					new_group_ary0[group0], new_group_ind_ary0[group0],
					new_group_ary1[group1], new_group_ind_ary1[group1]);
			}
		}
    }

    if(Debug_fd>0)
	{
		fclose(Debug_fd);
		Debug_fd=0;
    }

	/*
	if(!SILENT)fprintf(stderr,"Finish 1st round of tiling\n");
	*/
    
	for(Optimum_tiling_mode=1; Optimum_tiling_mode<3; Optimum_tiling_mode++)
	{
		while(pop_boundary_array(1-(Optimum_tiling_mode&1),&group0,&ind0,&tile_ind,
			&group1,&ind1))
		{
			do_partial_tiling(group0, group1, ind0, ind1,tile_ind,lines_ary0, 
				lines_ary1, new_group_ary0[group0], new_group_ind_ary0[group0],
				new_group_ary1[group1], new_group_ind_ary1[group1]);
		}
		/*
			if(!SILENT)fprintf(stderr,"Finish %drd round of tiling\n",Optimum_tiling_mode+1);
			
		*/
		security_check_untiled_stack(1-(Optimum_tiling_mode&1)); 
    }
    return Tile_triangle_num;
}


int christiansen( VertType *lines_ary0, VertType *lines_ary1, short ** group_ary0, 
				 short *group_ind_ary0, int group_num0,short ** group_ary1, 
				 short *group_ind_ary1, int group_num1)
{
    int i;
    int num0,num1, ind0,ind1;
    double minlen,len,len0,len1;
    int min_ind;
    VertType *pt1,*pt0;
    static int first=1;
    short *g_ary0, *g_ary1;

    g_ary0=group_ary0[0];
    g_ary1=group_ary1[0];

    if(group_num0!=1 || group_num1!=1)
	{
		fprintf(stderr,"It must be one-to-one group, exit\n");
		exit(1);
    }

    Tile_triangle_num=0;
    num0=group_ind_ary0[0];
    num1=group_ind_ary1[0];

    pt1=&(lines_ary0[g_ary0[0]]);
    minlen=line_distance_2d(pt1,&(lines_ary1[g_ary1[0]]));
    min_ind=0;
    
	for(i=1; i<num1; i++)
	{
		len=line_distance_2d(pt1,&(lines_ary1[g_ary1[i]]));
		if(len<minlen)
		{
			minlen=len;
			min_ind=i;
		}
    }

    ind0=0;
    ind1=min_ind;
    pt0=&(lines_ary0[g_ary0[ind0]]);
    pt1=&(lines_ary1[g_ary1[ind1]]);
    
	while(!(ind0==0 && ind1==min_ind && !first))
	{
		VertType *pt00,*pt11;
		int k0,k1;
		first=0;
		k0=(ind0+1)%num0;
		k1=(ind1+1)%num1;
		pt00=&(lines_ary0[g_ary0[k0]]);
		pt11=&(lines_ary1[g_ary1[k1]]);
		len0=line_distance_2d(pt1,pt00);
		len1=line_distance_2d(pt0,pt11);

		if(len0<len1)
		{
			draw_tiled_triangle(pt0,pt1,pt00);
			pt0=pt00;
			ind0=k0;
		}
		else 
		{
			draw_tiled_triangle(pt1,pt0,pt11);
			pt1=pt11;
			ind1=k1;
		}
    }
    return 1;
}

int get_no_of_poly_triangles()
{
	return no_of_poly_triangles;
}
