#include <stdio.h>
#include <math.h>

#include "common.h"
#include "math_util.h"
#include "legal_region.h"

extern signed char *Is_inside_positive[2];
extern VertType *Lines_ary[2];

void  calc_group_legal_equations(VertType *lines_ary,
								 short * new_group_tab, short group_ind,
								 TileType *best_tiling_table)
{
    int i,i0,i1,j,res;
    double x,y;

    for(i=0; i<group_ind; i++)
	{
		j=new_group_tab[i];
		i0=i-1;
	
		if(i0<0)i0=group_ind-1;
		i1=i+1;
		
		if(i1>=group_ind)i1=0;
		
		res=calc_line_equation(&(lines_ary[new_group_tab[i0]]), &(lines_ary[j]),
			&(best_tiling_table[i].a[0]),
			&(best_tiling_table[i].b[0]), &(best_tiling_table[i].c[0]));
		
		if(!res)
			fprintf(stderr,"Same point: index %d, %d, group=%X\n",i,i1,new_group_tab);
		
		res=calc_line_equation(&(lines_ary[j]),&(lines_ary[new_group_tab[i1]]), 
			&(best_tiling_table[i].a[1]),
			&(best_tiling_table[i].b[1]), &(best_tiling_table[i].c[1]));
		
		if(!res)
			fprintf(stderr,"Same point: index %d, %d, group=%X\n",i,i1,new_group_tab);
		
		x=lines_ary[new_group_tab[i1]][0];
		y=lines_ary[new_group_tab[i1]][1];
		
		if(best_tiling_table[i].a[0]*x +best_tiling_table[i].b[0]*y+
			best_tiling_table[i].c[0]>0) /* inside the legal region */
			best_tiling_table[i].and_or=1; /* narrow down region */
		else best_tiling_table[i].and_or=0;
    }
}



void make_group_legal_equation_tables(
									  VertType *lines_ary0, VertType *lines_ary1,
									  short ** new_group_ary0, short *new_group_ind_ary0, int group_num0,
									  short ** new_group_ary1, short *new_group_ind_ary1, int group_num1,
									  TileType **best_tiling_table[2])
{
    int j;
    //char **flag_group_ary[2];// KLC
	
    for(j=0; j<group_num0; j++)
        calc_group_legal_equations(lines_ary0,new_group_ary0[j], new_group_ind_ary0[j],
		best_tiling_table[0][j]);
	
    for(j=0; j<group_num1; j++)
        calc_group_legal_equations(lines_ary1,new_group_ary1[j], new_group_ind_ary1[j],
		best_tiling_table[1][j]);
}


