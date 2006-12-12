
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

/*---------------------------------------------------------------
It approximate a closed contour by line segments according to a
tolerance.
---------------------------------------------------------------- */


#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "common.h" 
#include "myutil.h" 
#include "cross.h" 
#include "approximate.h" 

static double Tolerance=1.0;
static VertType *Lines_ary;

extern int STATIC_CHQ;

void my_clear_approximate()
{
	Tolerance = 1.0;
	Lines_ary = NULL;	
}

void adaptive_approx(int start, int end, short *ind_ary,char *flag_ary)
{
    double dx,dy,a,b,c,d,dist; /* the line is ax+by+1.0=0 */
    double max=0.0;
    int i;
    int max_ind=-1;

    if(end-start<=1)return;

    dx=Lines_ary[ind_ary[start]][0]-Lines_ary[ind_ary[end]][0];
    dy=Lines_ary[ind_ary[start]][1]-Lines_ary[ind_ary[end]][1];
    a=dy;
    b=-dx;
    c=-(a*(Lines_ary[ind_ary[start]][0])+b*(Lines_ary[ind_ary[start]][1]));
    a/=c;
    b/=c;
    /* the distance is |ax+by+c|/(a^2 + b^2)^0.5 */
    d=sqrt(a*a+ b*b);
    
	for(i=start+1; i<end; i++)
	{
		dist=a*Lines_ary[ind_ary[i]][0]+b*Lines_ary[ind_ary[i]][1]+1.0;

		if(dist<0.0)	dist=-dist;
		dist=dist/d;
		if(dist>=max)
		{
			max=dist-0.0001; /* took the last one to avoid numerical unstability */
			max_ind=i;
		}
    }

    if(max<=Tolerance)	return;
    flag_ary[max_ind]=1;
    adaptive_approx(start,max_ind,ind_ary,flag_ary);
    adaptive_approx(max_ind,end,ind_ary,flag_ary);
}

int approximate_line(int num, short *short_group_ary, char *flag_ary)
{
    int num2;
    if(num<=3) return num;
    num2=num/2;
    flag_ary[0]=1;
    flag_ary[num2]=1;
    flag_ary[num-1]=1;
    adaptive_approx(0,num2,short_group_ary,flag_ary);
    adaptive_approx(num2,num-1,short_group_ary,flag_ary);
    return 1;
}

int approximate_contours( double *toler_ary,VertType *lines_ary, int lines_ind,
						short **group_ary, short *group_ind_ary, int group_num,
			            short **new_group_ary, short *new_group_ind_ary)
{
    int i,group,pt_num,cnt;
    short *short_ary,*short_group_ary;
    char flag_ary[MAX_LINES_NUM];
    int self_intersect;

    Lines_ary=lines_ary;

    for(group=0; group<group_num; group++)
	{
		int ind;
		Tolerance=toler_ary[group];
		
		if(new_group_ary[group]!=NULL)
		{
			free(new_group_ary[group]);
			new_group_ary[group]=NULL;
		}
		
		pt_num=group_ind_ary[group];
		short_group_ary=group_ary[group];
		
		if(short_group_ary[0]!=short_group_ary[pt_num-1]) 
		{ 
			/* the end pts are different */ 
			short_group_ary[pt_num++]=short_group_ary[0];
		}
		self_intersect=-1;
		while(self_intersect)
		{
			if(self_intersect>=0)		free(short_ary);
			cnt=0;
			while(cnt<=3 )
			{
				cnt=0;
				for(i=pt_num; i>=0; i--)flag_ary[i]=0;
				approximate_line(pt_num,short_group_ary,flag_ary);
				for(i=pt_num-1; i>=0; i--)if(flag_ary[i])cnt++;
				Tolerance/=2.0; /* do it again to avoid cnt<=3 */
			}
			
			/* pt_num--;*/ /* the end point is repeated, get rid of it  */
			short_ary=(short *)mymalloc(sizeof(*short_ary)*(cnt+2));
			new_group_ary[group]=short_ary;
			ind=0;

			for(i=0; i<(pt_num-1); i++)if(flag_ary[i])
			{
				short_ary[ind++]=short_group_ary[i];
			}

			new_group_ind_ary[group]=ind;
			short_ary[ind]=short_group_ary[0];
			self_intersect=is_polygon_self_crossed(lines_ary, short_ary, ind);
	
			if(self_intersect)
			printf("Warning! Group %d self-inter, reduce epsilon to %lf\n",

			group,Tolerance); 
		}

		printf("G%d reduce %d pts into %d lines\n",
	    group,group_ind_ary[group],new_group_ind_ary[group]);
    }
    return 1;
}
	


