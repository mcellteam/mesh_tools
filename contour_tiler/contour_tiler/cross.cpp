
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <malloc.h>

#include "common.h"
#include "math_util.h"
#include "myutil.h"
#include "cross.h"

#define MIN_DIST 0.02
#define MIN_VECT 0.00005

static VertType *Vert_ary_vor;
extern int SILENT;
extern int STATIC_CHQ;

void my_clear_cross()
{
	Vert_ary_vor = NULL;
}


void calc_vector_from_2vp(VertType *ret_v, VertType *p0, VertType *p1)
{
    (*ret_v)[0]=(*p1)[0]-(*p0)[0];
    (*ret_v)[1]=(*p1)[1]-(*p0)[1];
}

/*
void calc_vector_from_cir(VertType *ret_v, short *vert_ind_ary, int index)
{
int ind0,ind1;
VertType *p0=NULL,*p1=NULL;
ind0=vert_ind_ary[index];
ind1=vert_ind_ary[index+1];
p0=&Vert_ary_vor[ind0];
p1=&Vert_ary_vor[ind1];
(*p1)[0]-(*p0)[0];

 (*ret_v)[0]=(*p1)[0]-(*p0)[0]; 
 (*ret_v)[1]=(*p1)[1]-(*p0)[1];
 }
*/

/*KLC.... this function was commented...*/
/*int is_convex(short *vert_ind_ary,int poly_num, VertType *ret_vt)
// the vert_ind_ary[] should have a size of at least poly_num+1 
{
int i,sign,pre_sign=0;
VertType v1,v0;
double res,x=0,y=0;

  Vert_ary_vor = ret_vt; //KLC as the Vert_ary_vor is not being initialized anywhere...
  
   calc_vector_from_cir(&v0,vert_ind_ary,poly_num-1);
   for(i=0; i<poly_num; i++) {
   calc_vector_from_cir(&v1,vert_ind_ary,i);
   res=outer_product_2d(v0,v1);
   if(res > 0.1)sign=1;
   else if(res < -0.1)sign=2;
   else sign=0; // netrual 
   if(pre_sign==0 && sign)pre_sign=sign;
   if(sign && pre_sign !=sign)return 0;
   memcpy(v0,v1,sizeof(v0));
   }
   for(i=0; i<poly_num; i++) {
   x+= Vert_ary_vor[vert_ind_ary[i]][0];
   y+= Vert_ary_vor[vert_ind_ary[i]][1];
   }
   x=x/poly_num; 
   y=y/poly_num; 
   (*ret_vt)[0]=x;
   (*ret_vt)[1]=y;
   return 1;
   }
*/

/*KLC.... this function was commented...*/
/*void free_polygons_ary(PolygonStruct **polygon_ary)
{
int i;
PolygonStruct *poly;
if(polygon_ary==NULL)return;
i=0;
while(polygon_ary[i]){
poly=polygon_ary[i];
free(poly->vertices);
poly->vertices=NULL;
if(poly->normals){
free(poly->normals);
poly->normals=NULL;
}
free(polygon_ary[i]);
polygon_ary[i]=NULL;
i++;
}
free(polygon_ary);
}*/
/**/


int is_polygon_self_crossed(
							VertType *vert_ary, short *vert_ind_ary, int poly_num)
{
    int i,j,i1,j1,res;
    double ratio1,ratio2;
    VertType pt;

    for(i=2; i<(poly_num-1); i++)
	{
        i1=(i+1);
        for(j=0; j<=(i-2); j++)
		{
            j1=(j+1)%poly_num;
            res= find_intersection( &vert_ary[vert_ind_ary[i]],
                &vert_ary[vert_ind_ary[i1]],&vert_ary[vert_ind_ary[j]],
				&vert_ary[vert_ind_ary[j1]],&pt,&ratio1,&ratio2);
        
			if(res==1)
			{
				if(!(ratio1>0.9 || ratio1<0.1 ||ratio2>0.9 || ratio2<0.1))
					printf("Warning is_polygon_self_crossed(), ratio1=%lf, ratio2=%lf\n",
                    ratio1,ratio2);
            }
			
			if(res)
				return 1;
        }
    }
    return 0;
}

