#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include <malloc.h>

#include "common.h"
#include "myutil.h"
#include "math_util.h"
#include "tile_util.h"
#include "branch.h"
#include "contour_read.h"
#include "image2cont2D.h"
#include "tile_hash.h"
#include "tiling_main.h"
#include "scale.h"
#include "correspond.h"
#include "approximate.h"
#include "contour.h"

#define VERT_TRI_BLOCK_SIZE (4096*3) 
#define VERT_BLOCK_SIZE 2048 

extern int ALL_CORRESPOND,SILENT;
extern int SAVE_CONTOURS; 
extern int CHANGE_Z;
extern double TURBERANCE;
extern int DEBUG;
extern int CONTOURS_ONLY;
extern int DO_APPROX_ONLY;
extern int MESH_FILE;
extern double Tolerance;
extern float Scale_z;
extern FILE *Mesh_fp, *Cpts_fp, *CptsFirst_fp;
extern int SAME_LEVEL;
extern int Current_level;
extern int Beg_slice;
/* put the following variables as globle to save calling time */

VertType *Lines_ary[2];
int Lines_ary_size[2];
int Lines_ind[2] ;
static VertType *Untiled_vertice;
static short **New_group_ary[2],*New_group_ind_ary[2];
static short **Backup_group_ary[2],*Backup_group_ind_ary[2];
static short *GLevel[2];
TileType **Best_tiling_table[2];
int Group_num[2];
int Backup_group_num[2];
short *Group_connection_table;
short *Group_relation_table;
int Base_slice_num;
double Zunit=0;
short Hierachy[2][MAX_GROUP_NUM];

/* must use signed, some compiler treat char as unsigned */
signed char *Is_inside_positive[2];
extern int NUMERICALLY_STABLE;


int Dimz;
double Current_zvalue, Mid_zvalue, Zvalues[3]; /* 0: bottom (small Z) 1:top */
short ** Untiled_group_ary; 
short *Untiled_group_ind_ary;
int Untiled_group_num;
int Current_slice_num=-1;
int Total_triangle_number;

char Out_name[256];

short *Crossed_group[2]={NULL,NULL}, *Crossed_index[2]={NULL,NULL};
int Crossed_num=0;

/* The following is for volume calculation only */
extern double Gvolume_ary[10];
extern VertType *volumeVertAry;
extern int volumeVertAryInd;
extern int NON_VALID;
extern int CALC_VOLUME_ONLY;
double calculate_approximate_volume();
double calculate_accurate_volume();
double Local_zunit;



#ifndef i860
static int mynode(){return 0;}
#endif

extern int STATIC_CHQ;

void my_clear_contour()
{
	Untiled_vertice = NULL;
	New_group_ary[0] = NULL;		New_group_ary[1] = NULL;
	New_group_ind_ary[0] = NULL;	New_group_ind_ary[1] = NULL;
	Backup_group_ary[0] = NULL;		Backup_group_ary[1] = NULL;
	GLevel[0] = NULL;				GLevel[1] = NULL;						

	read_contours_from_1_file_sub(STATIC_CHQ, NULL, STATIC_CHQ);
	tile_from_contours(NULL, STATIC_CHQ, STATIC_CHQ, STATIC_CHQ, NULL, NULL);
	tile_from_contours_phase1(STATIC_CHQ, STATIC_CHQ, NULL, NULL, STATIC_CHQ, STATIC_CHQ);
}	


double calculate_approximate_volume()
{
    int i,j;
    double area[2],max,min,mid,volume;

    for(j=0; j<2; j++)
	{
        area[j]=0.0;
        for(i=0; i<Group_num[j]; i++)
		{
            area[j]+= calculate_area(Lines_ary[j], New_group_ary[j][i],
				New_group_ind_ary[j][i]);
        }
    }

    if(area[0]<area[1])
	{
        max=area[1];
        min=area[0];
    }
    else 
	{
        max=area[0];
        min=area[1];
    }
    
	mid=min+(max-min)*0.25;
    if(min>1e-6)	volume = ((max+min+ 4.0*mid)/6.0) * Local_zunit;
    else	volume = ((max+min+ 4.0*mid)/12.0) * Local_zunit;
    return volume;
}

static void alloc_Lines_ary_if_empty(int index)
{
    if(Lines_ary[index]==NULL) 
	{
        Lines_ary_size[index]=VERT_BLOCK_SIZE;
        Lines_ary[index]=(VertType *)mycalloc(sizeof(*Lines_ary[0])*
			Lines_ary_size[index]);
    }
}
static int put_lines_on_map(int index,int line_num,
							double lines[4][SPACE_VERT_DIM])
{
    int i,line,temp;
	
    alloc_Lines_ary_if_empty(index);
    temp=Lines_ind[index];
    
	for(line=0; line<line_num; line++) 
	{
        int k,k1;
        k=line*3;
        for(i=k; i<k+2; i++)
		{
            /* leave extra 1000 space for breaking contour */
            if(Lines_ary_size[index]<=Lines_ind[index]+1000)
			{
				/* enlarge the buffer */
                int temp;
                VertType *ary;
                temp=Lines_ind[index];
                Lines_ary_size[index]+=VERT_BLOCK_SIZE;
                ary=Lines_ary[index];
                Lines_ary[index]=(VertType *)
                    mycalloc(sizeof(*Lines_ary[0])*Lines_ary_size[index]);
                memcpy(Lines_ary[index],ary,sizeof(*ary)*temp);
                free(ary);
            }
            k1=Lines_ind[index];
            memcpy(&Lines_ary[index][k1],lines[i],
                sizeof(double)*SPACE_VERT_DIM);
				/*
				for(j1=0; j1<SPACE_VERT_DIM; j1++)
                Lines_ary[index][k1][j1]=lines[i][j1];
            */
            Lines_ind[index]++;
        }
    }
    return temp;
}


int make_contour_counter_clockwise(signed char is_positive,VertType *vert_ary,
								   short *vert_ind_ary, short poly_num)
{
    int res;
    VertType v0,v1;
    VertType *p0,*p1;
    VertType pt;
	
    p0=&(vert_ary[vert_ind_ary[0]]);
    p1=&(vert_ary[vert_ind_ary[1]]);
    v0[0]=(*p1)[0]-(*p0)[0];
    v0[1]=(*p1)[1]-(*p0)[1];
    v1[0]=-v0[1];
    v1[1]=v0[0];
    normalize_vector(&v1);
	
    /* make the pt as the middle point of p0 & p1 */
    pt[0]=((*p0)[0]+(*p1)[0])/2.0;
    pt[1]=((*p0)[1]+(*p1)[1])/2.0;
    pt[0]=pt[0]+v1[0]*0.03; /* shift to its left side a little bit */
    pt[1]=pt[1]+v1[1]*0.03;
	
    res=is_inside_contour(&pt, vert_ary, vert_ind_ary, poly_num, NULL, NULL);
    if((!res && is_positive) || (res && !is_positive))
	{
        reverse_group_direction(vert_ind_ary, poly_num);
		return 1;
    }
    return 0;
}

void save_a_group(FILE *ofp, VertType *lines_ary, 
				  short *group_tab, int group_size) 
{
    int i,j, j1,k1;
    for(i=0; i<group_size; i++)
	{
		j=group_tab[i];
		k1=i+1;
		if(k1>=group_size)	k1=0;
		j1=group_tab[k1];
		/* mark the first point */
		if(i==0) 
			draw_square(ofp,lines_ary[j][0], lines_ary[j][1], lines_ary[j][2]);
		else if(i==1)
			draw_grid(ofp,lines_ary[j][0], lines_ary[j][1], lines_ary[j][2]);
		else 
			draw_small_square(ofp,lines_ary[j][0], lines_ary[j][1], lines_ary[j][2]);

		draw_line_seg(ofp, lines_ary[j][0], lines_ary[j][1], lines_ary[j][2],
			lines_ary[j1][0], lines_ary[j1][1], lines_ary[j1][2]);
    }
}

void save_a_contour(FILE *ofp, VertType *lines_ary, 
					short *group_tab, int group_size) 
{
    int i,j;
    fprintf(ofp,"%d\n",group_size);
    
	for(i=0; i<group_size; i++)
	{
		j=group_tab[i];
		fprintf(ofp,"%f %f %f\n",lines_ary[j][0], lines_ary[j][1],lines_ary[j][2]);
    }
}



void save_contours(int topdown, char *s, int slice_num)
{
    char name[100];
    FILE *ofp;
    short **group_ary, *group_ind_ary,group_num;
    VertType *lines_ary;
    int group;
    sprintf(name,"%s%d_cont.poly",s,slice_num);

    if((ofp=fopen(name,"w"))==NULL) 
	{
        fprintf(stderr,"save_contours() Could not open %s for writing\n",name);
        exit(1);
    }
    
	group_num=Group_num[topdown];
    lines_ary=Lines_ary[topdown];
    group_ary=New_group_ary[topdown];
    group_ind_ary=New_group_ind_ary[topdown];
    
	for(group=0; group<group_num; group++)
	{
		save_a_group(ofp, lines_ary, group_ary[group], group_ind_ary[group]);
    }
    fclose(ofp);
}

void save_all_groups(int topdown, char *s, int slice_num)
{
    char name[100];
    FILE *ofp;
    short **group_ary, *group_ind_ary,group_num;
    VertType *lines_ary;
    int group;
    group_num=Group_num[topdown];
    lines_ary=Lines_ary[topdown];
    group_ary=New_group_ary[topdown];
    group_ind_ary=New_group_ind_ary[topdown];
    
	for(group=0; group<group_num; group++)
	{
		sprintf(name,"%s%dg%d.poly",s,slice_num,group);
		if((ofp=fopen(name,"w"))==NULL) 
		{
            fprintf(stderr,"save_all_groups() Could not open %s for writing\n",name);
            exit(1);
		}
		save_a_group(ofp, lines_ary, group_ary[group], group_ind_ary[group]);
		fclose(ofp);
    }
}


int print_group_connection_table(int group_num0, int group_num1)
{
    int i,j;
    int res=1;
	
    printf("\n ");

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
			printf("%d",Group_connection_table[j*group_num0+i]);
			if(!Group_connection_table[j*group_num0+i])res=0;
		}
		printf("\n");
    }
    return res;
}

void allocate_basic_memory()
{
    int i;
	/*
    if(Untiled_group_ary){
	fprintf(stderr,"allocate_basic_memory(), Untiled_group_ary=%X !=NULL\n",
	Untiled_group_ary);
	exit(1);
    }
	*/
    if(Untiled_group_ary)free (Untiled_group_ary);
    if(Untiled_group_ind_ary)free (Untiled_group_ind_ary);
    Untiled_group_ary=(short **)
		mycalloc(MAX_GROUP_NUM*sizeof(*Untiled_group_ary));
    Untiled_group_ind_ary= (short *)
		mycalloc(MAX_GROUP_NUM*sizeof(*Untiled_group_ind_ary));
	
    for(i=0; i<2; i++)
	{
		if(Is_inside_positive[i])
			free(Is_inside_positive[i]);
		
		Is_inside_positive[i]= (signed char *)mycalloc(MAX_GROUP_NUM*sizeof(char));
	
		if(New_group_ary[i])free(New_group_ary[i]);
		New_group_ary[i]=(short **)mycalloc(MAX_GROUP_NUM*sizeof(*New_group_ary));
		
		if(New_group_ind_ary[i])free(New_group_ind_ary[i]);
		New_group_ind_ary[i]=
			(short*)mycalloc(MAX_GROUP_NUM*sizeof(*New_group_ind_ary[0]));
    }
}

void allocate_backup_memory()
{
    int i;
    for(i=0; i<2; i++)
	{
		if(Backup_group_ary[i])
			free(Backup_group_ary[i]);
			
		Backup_group_ary[i]=(short **)mycalloc(MAX_GROUP_NUM*sizeof(*Backup_group_ary));

		if(Backup_group_ind_ary[i])
			free(Backup_group_ind_ary[i]);

		Backup_group_ind_ary[i]=(short*)mycalloc(MAX_GROUP_NUM*sizeof(*Backup_group_ind_ary[0]));
		
		if(GLevel[i])free(GLevel[i]);
			GLevel[i]= (short*)mycalloc(MAX_GROUP_NUM*sizeof(*GLevel[0]));
    }
}

void free_basic_memory()
{
    int i;
	
    if(Untiled_group_ary)free_group_ary((char **)Untiled_group_ary);
    if(Untiled_group_ary)free(Untiled_group_ary);
    Untiled_group_ary=NULL;
    if(Untiled_group_ind_ary)free(Untiled_group_ind_ary);
    Untiled_group_ind_ary=NULL;
	
    for(i=0; i<2; i++)
	{
		if(Is_inside_positive[i])
			free(Is_inside_positive[i]);
		
		Is_inside_positive[i]=NULL;

		if(New_group_ary[i])
			free(New_group_ary[i]);
		
		New_group_ary[i]=NULL;
		
		if(New_group_ind_ary[i])	
			free(New_group_ind_ary[i]);

		New_group_ind_ary[i]=NULL;
    }
}

void free_backup_memory()
{
    int  i;

    for(i=0; i<2; i++)
	{
		if(Backup_group_ary[i])
			free(Backup_group_ary[i]);
		Backup_group_ary[i]=NULL;

		if(Backup_group_ind_ary[i])
			free(Backup_group_ind_ary[i]);
		Backup_group_ind_ary[i]=NULL;

		if(GLevel[i])free(GLevel[i]);
		GLevel[i]=NULL;
    }
}

void free_individual_backup_memory()
{
    int  i,j;

    for(i=0; i<2; i++)
	{
		for(j=0; j<MAX_GROUP_NUM; j++)
		{
			if(Backup_group_ary[i][j])free(Backup_group_ary[i][j]);
			Backup_group_ary[i][j]=NULL;
		}
    }
}

void group_reset(int index)
{
    Lines_ind[index]=0;
}

void read_contours_from_1_file(int slice_num, char *prefix, char *sufix,int ind)
{
    char name[200];
//	if(strstr(prefix, "ImageSet"))
//		sprintf(name,"%s%04d%s",prefix, slice_num, sufix);
//	else
		sprintf(name,"%s%d%s",prefix, slice_num, sufix);
    read_contours_from_1_file_sub(slice_num, name, ind);
}

int read_contours_from_1_file_sub(int slice_num, char *name, int ind)
{
    FILE *fp;
    int i,j,num, gindex=0;
    short **group_tab, *group_ind_tab;
    double pt[4][SPACE_VERT_DIM] ;
    VertType v,v1;
    static double base_z=-1;
    static int base_z_index;
    extern double Read_file_time;
    int merged_num=0;
	
    if((fp=fopen(name,"r"))==NULL) 
	{
        fprintf(stderr,"read_contours_from_1_file() Could not open %s for reading\n",name);
        exit(1);
    }
    
	alloc_Lines_ary_if_empty(ind);
    group_tab=New_group_ary[ind];
    group_ind_tab=New_group_ind_ary[ind];
	
    group_reset(ind);
    
	while((fscanf(fp,"%d",&num)==1))
	{
		extern double MERGE_DISTANCE_SQUARE;
		merged_num=0;
		group_tab[gindex]=(short *)mymalloc(sizeof(*group_tab[0])*(num+3));
	
		for(i=0; i<num; i++)
		{
			if(fscanf(fp,"%lf %lf %lf",&v[0],&v[1],&v[2])!=3)
			{
				fprintf(stderr,"Corrupt file %s\n",name);
				exit(1);
			}
			
			if(i==0 || line_distance_2d(&v,&v1)>=MERGE_DISTANCE_SQUARE)
			{
				pt[0][0]=v[0];
				pt[0][1]=v[1];
				pt[0][2]=v[2]*Scale_z;
				if(CHANGE_Z)pt[0][2]=Zunit*slice_num;
				/*
				pt[0][0]+=((rand()&0xff)-128)/128.0 * TURBERANCE;
				pt[0][1]+=((rand()&0xff)-128)/128.0 * TURBERANCE;
				*/
				j=put_lines_on_map(ind,1, pt);
				group_tab[gindex][merged_num++]=j;
				memcpy(&v1,&v,sizeof(v));
			}
		}
        if(num<3)
		{ 
			/* Z as a reference
				   printf("Warning! contour has %d points, ignored, Z is useful\n",num); 
			/* make it looped  for easy access */ 
			free(group_tab[gindex]);
		}
		else
		{
			group_tab[gindex][merged_num]=group_tab[gindex][0]; 
			/* make core dump if illegal access */ 
			group_tab[gindex][merged_num+1]=20000; 
			group_ind_tab[gindex]=merged_num;
            if(merged_num!=num)
				printf("C%d %d pints merged to %d points\n",gindex,num,merged_num);
			gindex++;
		}
    }
    fclose(fp);
    group_tab[gindex]=NULL;
    Group_num[ind]=gindex;
	
    if(base_z==-1)
	{
		if(Zunit>0)
		{
			base_z=pt[0][2];
			base_z_index=slice_num;
		}
    }
    return 1;
}

int read_contours_from_2_files(int k, int ind0, int ind1,
							   char *prefix, char *sufix)
{
    fprintf(stderr,"N%d --- reading & processing Slice %d %d\n",
		mynode(),Current_slice_num-1,Current_slice_num); 
    printf("N%d --- reading & processing Slice %d %d\n",
		mynode(),Current_slice_num-1,Current_slice_num); 

    read_contours_from_1_file(Current_slice_num-1,prefix,sufix,ind0);
    read_contours_from_1_file(Current_slice_num,prefix,sufix,ind1);
    return 1;
}

void initialize_tile_table(TileType *tile_table, int num)
{
    int j;
	
    for(j=0; j<=num; j++)
	{
        tile_table[j].group=-1;
        tile_table[j].ind=-1;
        tile_table[j].dist=1000000.0;
        tile_table[j].used=0;
    }
}

void allocate_best_tiling_table(int ind0, int ind1)
{
    int i,k;
    int group_num0, group_num1;
    group_num0 = Group_num[ind0];
    group_num1 = Group_num[ind1];
    
	for(k=0; k<2; k++)
	{
		if(Best_tiling_table[k]!=NULL)
		{
			i=0;
			while(Best_tiling_table[k][i]!=NULL)
			{
				free(Best_tiling_table[k][i]);
				Best_tiling_table[k][i]=NULL;
				i++;
			}
		}

		free(Best_tiling_table[k]);
		Best_tiling_table[k]=NULL;
    }

    Best_tiling_table[0]=(TileType **) mycalloc(sizeof(TileType *)*(group_num0+1));
    Best_tiling_table[1]=(TileType **) mycalloc(sizeof(TileType *)*(group_num1+1));

		/* Although calloc() suppose to clear all memory, but do not trust it.
	set the following to be NULL. */
    Best_tiling_table[0][group_num0]=NULL;
    Best_tiling_table[1][group_num1]=NULL;
    
	for(i=0; i<group_num0; i++)
	{
        Best_tiling_table[0][i]= (TileType*)mycalloc(sizeof(TileType)*(New_group_ind_ary[ind0][i]+1));
        initialize_tile_table(Best_tiling_table[0][i], 	New_group_ind_ary[ind0][i]);
    }

    for(i=0; i<group_num1; i++)
	{
        Best_tiling_table[1][i]= (TileType*)mycalloc(sizeof(TileType)*(New_group_ind_ary[ind1][i]+1));
        initialize_tile_table(Best_tiling_table[1][i], New_group_ind_ary[ind1][i]);
    }
}

int do_major_work(char *out_name, int k, int ind0, int ind1)
{
    int jj,i,res;
    int vertice_num[2];
	
    volumeVertAryInd=0;
    for(jj=0; jj<2; jj++)
	{
		vertice_num[jj]=0;
		for(i=0; i<Group_num[jj]; i++)
			vertice_num[jj]+=New_group_ind_ary[jj][i];
    }
    if(!SILENT)
		printf("Top vertice #: %d, bot= %d, sum=%d\n",vertice_num[0],vertice_num[1],
			vertice_num[0]+vertice_num[1]);
    init_tile_hash_table(vertice_num[0]+vertice_num[1]);
	
    res=do_tiling(Lines_ary[ind0], Lines_ary[ind1],New_group_ary[ind0], 
			New_group_ind_ary[ind0], Group_num[ind0],New_group_ary[ind1], 
			New_group_ind_ary[ind1], Group_num[ind1]);
	
#ifdef i860
    if(numnodes()>1)request_work();
#endif
	if(!DEBUG)
		Untiled_group_num= collect_untiled_contours(Lines_ary[ind0], Lines_ary[ind1],
			New_group_ary[ind0], New_group_ind_ary[ind0], Group_num[ind0],
			New_group_ary[ind1], New_group_ind_ary[ind1], Group_num[ind1],
			Untiled_vertice,Untiled_group_ary, Untiled_group_ind_ary);
	/* security_check_untiled_stack(1); */
	
    if(CALC_VOLUME_ONLY)
	{
		double volume1, volume2;
		volume1 = calculate_approximate_volume();
		printf("level %d approx. vol=%f\n",Current_level,volume1);
		if(Untiled_group_num==0)
		{ /* All are tiled, so we can calculate
			the accurate volume */
			static double appvol=0.0, accvol=0.0;
			appvol+=volume1;
		
			volume2=calculate_accurate_volume();
			accvol+=volume2;
			printf("accurate volume=%f, app./acc = %f\n",
				volume2, volume1/volume2);
			printf("accumulate accurate volume=%f, app. =%f, ratio=%f\n",
				accvol, appvol, appvol/accvol);
			Gvolume_ary[Current_level]+= volume2;
		}
		else Gvolume_ary[Current_level]+= volume1;
    }
	
    free_tile_hash_table();
    if(k==1)
	{
		extern int WriteIndividualFileInPoly;
		if(SAVE_CONTOURS)save_contours(0,out_name,Current_slice_num-1);
		if(WriteIndividualFileInPoly)
		{
			save_all_groups(0, "slice", 0);
			save_all_groups(1, "slice", 1);
		}
    }

    if(SAVE_CONTOURS)save_contours(k&1,out_name,Current_slice_num);
	
    free_group_ary((char**)New_group_ary[ind0]);
    free_group_ary((char**)New_group_ary[ind1]);
	
    return res;
}



/* flag  0 : It is preprocessing, just want to break the contour
flag  1 : normal 
flag  2 : It is postprocessing, does not need to break the contour 
*/
void handle_correspondence(int ind0, int ind1, int flag)
{
    extern double XSHIFT;
    
	if(ALL_CORRESPOND)
	{
		scale_shift(XSHIFT,0.0,
			Lines_ary[ind0], Lines_ary[ind1],
			New_group_ary[ind0], New_group_ind_ary[ind0], Group_num[ind0],
			New_group_ary[ind1], New_group_ind_ary[ind1], Group_num[ind1]);
    }
	
    if(flag <2)
	{
		break_contour_if_necessary(
			Lines_ary[ind0], Lines_ary[ind1],
			New_group_ary[ind0], New_group_ind_ary[ind0], Group_num[ind0],
			New_group_ary[ind1], New_group_ind_ary[ind1], Group_num[ind1],
			&(Lines_ind[ind0]),&(Lines_ind[ind1]),
			Lines_ary_size[ind0],  Lines_ary_size[ind1]);
		if(!flag)return;
    }
	
	allocate_best_tiling_table(ind0,ind1);
	
	if(Group_connection_table!=NULL)	free(Group_connection_table);
	Group_connection_table=(short*)mycalloc(
						sizeof(*Group_connection_table)*Group_num[ind0]*Group_num[ind1]);
	
	if(Group_relation_table!=NULL)	free(Group_relation_table);

	Group_relation_table=(short*)mycalloc(
						sizeof(*Group_relation_table)*Group_num[ind0]*Group_num[ind1]);

	build_contours_relation(Lines_ary[ind0], Lines_ary[ind1],New_group_ary[ind0], 
				New_group_ind_ary[ind0], Group_num[ind0],New_group_ary[ind1],
				New_group_ind_ary[ind1], Group_num[ind1]);

	if(ALL_CORRESPOND)
	{
		int res;
		res=print_group_connection_table(Group_num[ind0],Group_num[ind1]);
		/*
		if(!res){
		fprintf(stderr,"The ALL_CORRESPOND flag is on, but fail to make it\n");
		exit(1);
		}
		*/
	}
}

void static_save_contours_pts_format(char *s, int slice_num)
{
    char name[100];
    FILE *ofp;
    short *group_tab, group_ind,group_num;
    VertType *lines_ary;
    int i,j,group;
    sprintf(name,"%s%d.pts",s,slice_num);
    
	if((ofp=fopen(name,"w"))==NULL) 
	{
        fprintf(stderr,"static_save_contours_pts_format() Could not open %s for writing\n",name);
        exit(1);
    }
    
	group_num=Group_num[0];
    lines_ary=Lines_ary[0];
    
	for(group=0; group<Group_num[0]; group++)
	{
        group_tab=New_group_ary[1][group];
        group_ind=New_group_ind_ary[1][group];
        fprintf(ofp, "%d\n",group_ind);
        for(i=0; i<group_ind; i++)
		{
            j=group_tab[i];
            fprintf(ofp,"%lf %lf %lf\n",Lines_ary[0][j][0],
                Lines_ary[0][j][1],Lines_ary[0][j][2]);
        }
    }

    fclose(ofp);
    fprintf(stderr,"Saved contour file %s\n",name);
}


int do_approx_only(char *outname, int base_slice_num, int dimz, 
				   double zunit, char *prefix, char *sufix) 
{
    fprintf(stderr,"do_approx_only() disabled\n");
    exit(1);
	/*
    int i, j, k, kk ,poly_num=0, res;
    int total_tiled_triange=0; 
	
	 Dimz=dimz;
	 Zunit=zunit;
	 if(Zunit<=0 && CHANGE_Z){
	 printf("Zunit=%lf, wrong\n",Zunit);
	 exit(1);
	 }
	 Base_slice_num=base_slice_num;
	 
	  allocate_basic_memory();
	  
	   for(k=0; k<Dimz+1; k++) {
	   group_reset(0);
	   read_contours_from_1_file(Base_slice_num+k,prefix,sufix,0);
	   approximate_contours(Tolerance, Lines_ary[0], Lines_ind[0],
	   New_group_ary[0], New_group_ind_ary[0], Group_num[0],
	   New_group_ary[1], New_group_ind_ary[1]);
	   static_save_contours_pts_format(outname, Base_slice_num+k);
	   }
	   return 1;
	*/
    return 1;
}

void save_a_slice(FILE *ofp, int updown)
/* save the contour information into a file prepared for mesh generation */
{
    int j,k;
    int group, group_num,group_ind;
    VertType *lines_ary;
    group_num=Group_num[updown];
    lines_ary=Lines_ary[updown];
    
	for(group=0; group<Group_num[updown]; group++)
	{
		short *group_tab=New_group_ary[updown][group];
		group_ind=New_group_ind_ary[updown][group];
		fprintf(ofp, "%d\n",group_ind);
	
		for(k=0; k<group_ind; k++)
		{
			j=group_tab[k];
			fprintf(ofp,"%lf %lf %lf\n",lines_ary[j][0],
				lines_ary[j][1],lines_ary[j][2]);
		}
    }
}


int find_next_match(int ind, VertType *vAry, char *flag_ary, int size)
{
    int i;
    VertType *p;
    double min[2],max[2];
    
	p = & (vAry[ind]);

    for(i=0; i<2; i++)
	{
		min[i]=(*p)[i]-1e-4;
		max[i]=(*p)[i]+1e-4;
    }
    
	for(i=0; i<size; i++)
	{
		if(!flag_ary[i]){
			p = & (vAry[i]);
			if((*p)[0]>min[0] && (*p)[0]<max[0] && (*p)[1]>min[1] && (*p)[1]<max[1])
			{ 
				/* the same point */
				return i;
			}
		}
    }
	/*    printf("Error! no match, find_next_match(ind=%d, .. size=%d\n",ind,size);
	*/
    return -1;
}

double calculate_accurate_volume()
{
    /* first need to make the group_ary */
    int cnt=0;
    short group_ary[2000];
    char  flag_ary[4000]; 
    int i,j,ind=0,not_done=1;
    double area[2],middle_area=0.0;
    for(i=0; i<volumeVertAryInd; i++)flag_ary[i]=0;
	
    for(j=0; j<2; j++)
	{
		area[j]=0.0;
        for(i=0; i<Group_num[j]; i++)
		{
            area[j]+= calculate_area(Lines_ary[j], New_group_ary[j][i],
						New_group_ind_ary[j][i]);
        }
    }
	
    while(not_done)
	{
		ind=0;
		j=0;
		while(j<volumeVertAryInd)
		{
			if(flag_ary[j]==0)	break; /* not used */
			j=j+2;
		}

		if(j>=volumeVertAryInd)	break;

		group_ary[ind++] = j;
		flag_ary[j]=1; /* It is used */
		
		while(j>=0)
		{
			j=(j/2)*4 + 1 -j; /* take its pair, e.g. 4-5, 5-4 */
			group_ary[ind++]=j;
			flag_ary[j]=1;
			j=find_next_match(j, volumeVertAry, flag_ary, volumeVertAryInd);
			flag_ary[j]=1;
		}
        middle_area+= calculate_area(volumeVertAry, group_ary, (short)ind);
		cnt++;
    }
    printf("%d middle section loop\n",cnt);
    return (area[0]+area[1]+middle_area*4)/6.0 * Local_zunit;
}


int open_mesh_files(int k)
{
    char str[200];
    char mode[10];
    if(!MESH_FILE)return 0;
    if(!Current_level)	sprintf(mode,"w");
    else	sprintf(mode,"a");

    sprintf(str,"%s%d.cpts",Out_name, Current_slice_num);

    if((Cpts_fp=fopen(str,mode))==NULL) 
	{
        fprintf(stderr,"Could not open %s for writing, mode=%s\n",str,mode);
        exit(1);
    }
    
	if(Current_slice_num-1 == Beg_slice)
	{
		/* the first */
        sprintf(str,"%s%d.cpts",Out_name, Current_slice_num-1);
        if((CptsFirst_fp=fopen(str,mode))==NULL) 
		{
            fprintf(stderr,"Could not open %s for writing, mode=%s\n",str,mode);
            exit(1);
		}
    }
    return 1;
}



int tile_from_contours(char *out_name, int base_slice_num, int dimz, 
					   double zunit, char *prefix, char *sufix) 
					   /* section : 0: first, 1: the rest */
{
    static int processed_num=0;
    int k, kk;
    int level; 
	/*    int total_tiled_triange=0; */
	
    if(zunit>1e-6)Zunit=zunit;
    Dimz=dimz;
    Base_slice_num=base_slice_num;
    strcpy(Out_name, out_name);
	
    if(!processed_num)
	{
		Total_triangle_number=0;
		allocate_basic_memory();
		if(SAME_LEVEL)allocate_backup_memory();
    }
	
    kk=1;
    
	for(k=1; k<Dimz; k++) 
	{
		extern int CHRISTIANSEN_METHOD;
		int ind0,ind1;
		int res;
		extern float FIRST_Z;
		ind1=k&1;
		ind0=1-ind1; /* top slice */
		Current_slice_num=Base_slice_num+k;
		read_contours_from_2_files(k, ind0, ind1,prefix,sufix);
		if(MESH_FILE)
		{
			/* prepare file for mesh generation */
			char str[200];
			sprintf(str,"%s%d.tile",Out_name, Current_slice_num-1);
			if((Mesh_fp=fopen(str,"w"))==NULL) 
			{
				fprintf(stderr,"Could not open %s for writing\n",str);
				exit(1);
			}
		}
		
		
		if(CHRISTIANSEN_METHOD)
		{
			res=christiansen(
				Lines_ary[ind0], Lines_ary[ind1],
				New_group_ary[ind0], New_group_ind_ary[ind0], Group_num[ind0],
				New_group_ary[ind1], New_group_ind_ary[ind1], Group_num[ind1]);
			close_files();
			fprintf(stderr,"Done\n");
			exit(1);
		}

		Current_zvalue=Lines_ary[ind1][0][2];
		
		if(Lines_ind[0] && Lines_ind[1])
		{
			Zvalues[0]=Lines_ary[0][0][2];
			Zvalues[1]=Lines_ary[1][0][2];
			Mid_zvalue=(Lines_ary[0][0][2]+Lines_ary[1][0][2])/2.0;
			Zvalues[2]=Mid_zvalue;
		}
		else if(Lines_ind[0])
		{
			if(Zunit>0.001) 
			{
				Mid_zvalue=Lines_ary[0][0][2]+(Zunit/2.0);
				FIRST_Z= (float)Zunit;
				Zvalues[0]=Lines_ary[0][0][2];
				Zvalues[1]=Zvalues[0]+Zunit;
				Zvalues[2]=Mid_zvalue;
			}
			else
			{
				if(FIRST_Z<1e-4)
				{
					printf("Error! 1st slice is empty, you need to specify a point, or use DETZ\n");
					exit(1);
				}
			}
		}
		else if(Lines_ind[1])
		{
			if(Zunit>0.001) 
			{
				Mid_zvalue=Lines_ary[1][0][2]-(Zunit/2.0);
				FIRST_Z= (float)Zunit;
				Zvalues[1]=Lines_ary[1][0][2];
				Zvalues[0]=Zvalues[1]-Zunit;
				Zvalues[2]=Mid_zvalue;
			}
			else 
			{
				if(FIRST_Z<1e-4)
				{
					printf("Error! 2nd slice is empty, you need to specify DETZ in .config file or use a point\n");
					exit(1);
				}
			}
		}
		
		
		Current_level=0;
		if(SAME_LEVEL)
		{ /* only tile between same level */
			level=copy_input_to_backup();
		}
		else	level=0; /* tread everythind as the same level */

		for(Current_level=0; Current_level<=level; Current_level++)
		{
			if(SAME_LEVEL)copy_backup_to_input(Current_level);
			
			if(Group_num[0] || Group_num[1])
			{
				int flag; //is_mesh; KLC
				if(NON_VALID)flag=1;
				else flag=2;
				if(Zunit==0.0) 
				{
					Local_zunit=Lines_ary[1][0][2]-Lines_ary[0][0][2];
					Zunit = Local_zunit;
				}
				else Local_zunit=Zunit;
				
				handle_correspondence(ind0, ind1,flag); 
				/* save_both need to be after handle_corres... (split contour) */
				open_mesh_files(k);
				if(Cpts_fp)
				{
					save_a_slice(Cpts_fp,1); /* save the top one */
					fclose(Cpts_fp);
					Cpts_fp=NULL;
				}
				
				if(CptsFirst_fp)
				{
					save_a_slice(CptsFirst_fp,0); /* bottom one*/
					fclose(CptsFirst_fp);
					CptsFirst_fp=NULL;
				}
				res=do_major_work(out_name,k,ind0,ind1);
				
				Total_triangle_number+=res;
				
				if(NUMERICALLY_STABLE)
				{
					int i;
					for(i=0; i<2; i++)
					{
						if(Crossed_group[i])	free(Crossed_group[i]); 
						Crossed_group[i]=NULL;
						if(Crossed_index[i])	free(Crossed_index[i]); 
						Crossed_index[i]=NULL;
					}
					Crossed_num=0;
				}
			}
		}
		/* need to free all individual backup memory. */
		/*	if(SAME_LEVEL)free_individual_backup_memory(); */ 
		kk=1-kk;
		if(Mesh_fp)fclose(Mesh_fp);
		Mesh_fp=NULL;
    }
    processed_num++;
    return 1;
}




int tile_from_contours_phase1(int base_slice_num, int dimz, 
							  char *prefix, char *sufix,int begnum, int endnum) 
{
    static int processed_num=0;
    int k, kk;
    int level; 
	
    Base_slice_num=base_slice_num;
	
    if(!processed_num)
	{
		allocate_basic_memory();
		allocate_backup_memory();
    }
	
    kk=1;
    
	for(k=1; k<dimz; k++) 
	{
		int flag=0; /* This is pre-processing */
		int ind0,ind1;
		ind1=k&1;
		ind0=1-ind1; /* top slice */
		Current_slice_num=Base_slice_num+k;
		read_contours_from_2_files(k, ind0, ind1,prefix,sufix);
		
		Current_level=0;
		if(SAME_LEVEL)
		{ /* only tile between same level */
			level=copy_input_to_backup();
		}
		else level=0; /* tread everythind as the same level */
		
		for(Current_level=0; Current_level<=level; Current_level++)
		{
			if(SAME_LEVEL)copy_backup_to_input(Current_level);
			
			if(Group_num[0] || Group_num[1])
			{
				handle_correspondence(ind0, ind1, flag);
				
				if(NUMERICALLY_STABLE)
				{
					int i;
					for(i=0; i<2; i++)
					{
						if(Crossed_group[i])	free(Crossed_group[i]); 
						Crossed_group[i]=NULL;
						if(Crossed_index[i])	free(Crossed_index[i]); 
						Crossed_index[i]=NULL;
					}
					Crossed_num=0;
				}
			}
			if(SAME_LEVEL)copy_one_level_to_backup(Current_level);
		}
		
		if(!SAME_LEVEL)
		{ /* only tile between same level */
			level=copy_input_to_backup();
		}
		
		/* need to save the backup into files */ 
		
		save_backup_files(base_slice_num, begnum, endnum);
		
		
		/* need to free all individual backup memory. */
		free_individual_backup_memory(); 
		kk=1-kk;
		
    }
    processed_num++;
    return 1;
}


int closed_all_tiling_files()
{
    int j,k;
	
    close_files();
    j=get_untiled_triangle_number();
    k=get_convex_triangle_number();
    printf("Tri num:  %d (%d + %d (convex) +%d (untiled))\n",
        Total_triangle_number+j+k,Total_triangle_number,k,j);

    free_basic_memory();
    if(SAME_LEVEL)	free_backup_memory();
	
    return Total_triangle_number+j+k;
}

int  save_backup_files(int num, int begnum, int endnum)
{
    char name[2][256];
    FILE *ofp;
    int i;
	/* //KLC
    if(num==begnum)
	sprintf(name[0],"/tmp/tiling.%d",num);
    else 
	sprintf(name[0],"/tmp/tiling.%d.1",num);
    if(num+1==endnum)
	sprintf(name[1],"/tmp/tiling.%d",num+1);
    else 
	sprintf(name[1],"/tmp/tiling.%d.0",num+1);
	*/
	
	if(num==begnum)
		sprintf(name[0],"tiling.%d",num);
    else 
		sprintf(name[0],"tiling.%d.1",num);

    if(num+1==endnum)
		sprintf(name[1],"tiling.%d",num+1);
    else 
		sprintf(name[1],"tiling.%d.0",num+1);
	
    for(i=0; i<2; i++)
	{
		int group;
		if((ofp=fopen(name[i],"w"))==NULL) 
		{
			fprintf(stderr,"Cannot open %s for writing temp file\n",name[i]);
			exit(1);
		}
		
		for(group=0; group<Backup_group_num[i]; group++)
		{
			if(Backup_group_ary[i][group]==NULL)
			{
				fprintf(stderr,"Error! Backup_group_ary[%d][%d]=NULL\n",i,group);
				exit(1);
			}
			
			save_a_contour(ofp, Lines_ary[i], Backup_group_ary[i][group], 
				Backup_group_ind_ary[i][group]);
		}

		if(group==0)
		{
			fprintf(ofp,"1\n0 0 %f\n",Lines_ary[i][0][2]);
		}
        
		fclose(ofp);
    }
    return 1;
}


int  copy_input_to_backup()
{
    int i,j,max=-1;
    
	for(i=0; i<2; i++)
	{
		for(j=0; j<Group_num[i]; j++)
		{
			Backup_group_ind_ary[i][j]=New_group_ind_ary[i][j];
			Backup_group_ary[i][j]=New_group_ary[i][j];
			New_group_ary[i][j]=NULL;
		}
		
		Backup_group_ind_ary[i][j]=0;
		Backup_group_num[i]=Group_num[i];
    }

    build_contours_levels(Lines_ary, Backup_group_ary,Backup_group_ind_ary, 
		Group_num, GLevel);
    
	for(i=0; i<2; i++)
	{
		for(j=0; j<Group_num[i]; j++)
		{
			if(GLevel[i][j]>max)max=GLevel[i][j];
		}
    }
    return max;
}
void copy_backup_to_input(int level)
{
    int i,j;

    for(i=0; i<2; i++)
	{
		int k=0;
		int cnt=0;
		for(j=0; j<Backup_group_num[i]; j++)
		{
			if(GLevel[i][j]==level)
			{
				if(New_group_ary[i][k])	free(New_group_ary[i][k]);
				New_group_ary[i][k]=Backup_group_ary[i][j];
				Backup_group_ary[i][j]=NULL; 
				New_group_ind_ary[i][k]=Backup_group_ind_ary[i][j];
				Backup_group_ind_ary[i][j]=0;
				cnt++;
				k++;
			}
		}
		Group_num[i]=cnt;
    }
}


void copy_one_level_to_backup(int level)
{
    int i,j;
    
	for(i=0; i<2; i++)
	{
		int k=0;
		int cnt=0;
		
		for(j=0; j<Backup_group_num[i]; j++)
		{
			if(GLevel[i][j]==level)
			{
				if(Backup_group_ary[i][j])
				{
					printf("Error! copy_one_level_to_backup(%d), j=%d, not NULL\n",level,j);
					exit(1);
				}
				
				Backup_group_ary[i][j]=New_group_ary[i][k];
				New_group_ary[i][k]=NULL; 
				Backup_group_ind_ary[i][j]=New_group_ind_ary[i][k];
				Backup_group_ind_ary[i][j]= New_group_ind_ary[i][k];
				New_group_ind_ary[i][k]=0;
				cnt++;
				k++;
			}
		}
    }
}

/* make the file consistant.
Reason: The tiling algorithm add points to break a contour segments
into more segments. A contour is used twice. One for the upper section,
the other for the bottom section. It could be added different vertices
during the tiling of these two section. Thus the triangulation could
be no longer valid (each edge is shared exactly twice). Non-valid
triangulation is all right for visualization, but it can not be used
in finite element analysis. So We need to make it valid.
*/

int make_file_consistant(int num)
{
    int i,group,index[2];
    char name[200];
    VertType *p[2],*v;
    FILE *ofp;
	
    short *mygroup_tab;
    char *flag_ary;
    int mygroup_size;
	
    short *group_tab[2], group_size[2];
    //sprintf(name,"/tmp/tiling.%d.0",num); //KLC
	sprintf(name,"tiling.%d.0",num);
    read_contours_from_1_file_sub(num, name, 0);
    //sprintf(name,"/tmp/tiling.%d.1",num); //KLC
	sprintf(name,"tiling.%d.1",num);
    read_contours_from_1_file_sub(num, name, 1);
	
    //sprintf(name,"/tmp/tiling.%d",num); //KLC
	sprintf(name,"tiling.%d",num);

    if((ofp=fopen(name,"w"))==NULL) 
	{
        fprintf(stderr,"make_file_consistant(%d) Could not open %s for writing\n",num,name);
        exit(1);
    }
	
    /* find the matching of starting point */
    if(Group_num[0]!=Group_num[1])
	{
		printf("Error! make_file_consistant(%d), groupnum=%d %d, not equal\n",
			num,Group_num[0], Group_num[1]);
		exit(1);
    }
	
    for(group=0; group<Group_num[0]; group++)
	{
		for(i=0; i<2; i++)
		{
			group_tab[i]=New_group_ary[i][group];
			group_size[i]=New_group_ind_ary[i][group];
			p[i] = &(Lines_ary[i][group_tab[i][0]]);
		}
		
		mygroup_tab=(short*)mymalloc((group_size[0]+group_size[1])*sizeof(short));
		flag_ary=(char*)mymalloc(group_size[0]+group_size[1]);
		mygroup_size=0;
		/* the starting points should be the same */
		if(line_distance(p[0],p[1])>1e-6)
		{
			printf("Error! the first point should be the same, S%d,G%d\n",num,group); 
			exit(1);
		}
		
		mygroup_size=0;
		mygroup_tab[mygroup_size]=group_tab[0][0];
		flag_ary[mygroup_size++]=0;
		index[0]=1;
		index[1]=1;
		
		while(index[0]<group_size[0] || index[1]<group_size[1])
		{
			for(i=0; i<2; i++)
				p[i] = &(Lines_ary[i][group_tab[i][index[i]]]);
			
			if(line_distance(p[0],p[1])<1e-6)
			{ /* same point */
				mygroup_tab[mygroup_size]=group_tab[0][index[0]];
				flag_ary[mygroup_size++]=0;
				index[0]++;
				index[1]++;
			}
			else 
			{
				double dist[2];
				int i1=mygroup_size-1;
				v = &(Lines_ary[flag_ary[i1]][mygroup_tab[i1]]);
			
				for(i=0; i<2; i++)
				{
					dist[i]=line_distance(v,p[i]);
				}
				
				if(dist[0]<dist[1])
				{
					mygroup_tab[mygroup_size]=group_tab[0][index[0]];
                    flag_ary[mygroup_size++]=0;
                    index[0]++;
				}
				else 
				{
					mygroup_tab[mygroup_size]=group_tab[1][index[1]];
                    flag_ary[mygroup_size++]=1;
                    index[1]++;
				}
			}
		}

		if(index[0]<group_size[0] || index[1]<group_size[1])
		{
			printf("Error! index[]=%d %d, size=%d %d\n",
				index[0],index[1],group_size[0],group_size[1]);
		}
		// save this group
		fprintf(ofp,"%d\n",mygroup_size);
		for(i=0; i<mygroup_size; i++){
			v = &(Lines_ary[flag_ary[i]][mygroup_tab[i]]);
			fprintf(ofp,"%f %f %f\n",(*v)[0], (*v)[1], (*v)[2]);
		}
		free (mygroup_tab);
		free(flag_ary);
    }
    fclose (ofp);
    return 1;
}

