#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "meshalyzer.h"

extern FILE *yyin;
char *infile;
int line_num;
int total_vertex_count,total_face_count;
struct object *world_obj;

/* Begin main */
int main(argc,argv)
  int argc; 
  char *argv[];  
{
  struct surface_mesh *smp;

	if (argc<2) {
      	  printf("Usage: %s in_file_name\n",argv[0]);
	  exit(1);
	}

	infile=argv[1];

	world_obj=NULL;
	if ((world_obj=(struct object *)malloc(sizeof(struct object)))==NULL) {
	  printf("meshalyzer: cannot store world object\n");
	  exit(1);
        }
        world_obj->next=NULL;
        world_obj->parent=NULL;
        world_obj->first_child=NULL;
        world_obj->last_child=NULL;
        world_obj->name="world metaobject";
        world_obj->object_type=META_OBJ;
        world_obj->contents=NULL;

	if ((yyin=fopen(infile,"r"))==NULL) {
	  printf("meshalyzer: error opening file: %s\n",infile);
	  exit(1);
	} 
	fflush(stdout);
	if (yyparse()) {
	  printf("meshalyzer: error parsing file: %s\n",infile);
	  exit(1);
	} 
	fclose(yyin);


        smp=(struct surface_mesh *)world_obj->first_child->contents;

/*
        fprintf(stderr,"Sort vertex array...");
        qsort((void *)smp->vertex_array,smp->n_verts,sizeof(smp->vertex_array),compare_verts);
        fprintf(stderr,"Done\n");

        fprintf(stderr,"Merging duplicate vertices...\n");
        merge_duplicate_verts(smp,EPS_D);
        fprintf(stderr,"Done\n");
*/

        analyze_surface_mesh(smp);
/*
        output_surface_mesh(smp);
*/

	exit(0);
}



int compare_verts(const void *p1, const void *p2)
{
  struct vertex *v1,*v2;

  v1 = * (struct vertex **)p1;
  v2 = * (struct vertex **)p2;

  if (v1->projection<v2->projection)
  {
    return(-1);
  }
  if (v1->projection>v2->projection)
  {
    return(1);
  }
  return(0);

}



/**********************************************************************
distinguishable -- reports whether two doubles are measurably different

Parameters
        a -- first double
        b -- second double
        eps -- fractional difference that we think is different

Returns
        1 if the numbers are different, 0 otherwise
**********************************************************************/

int distinguishable(double a,double b,double eps)
{
  double c;

  c=a-b;

  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;

  if (b<a) eps*=a;
  else eps*=b;
  return (c>eps);
}



void merge_duplicate_verts(struct surface_mesh *smp, double epsilon)
{
  struct vertex *vp1,*vp2;
  int running_index;
  int identical_verts,epsilon_verts;
  int i,j;

  identical_verts=0;
  epsilon_verts=0;
  running_index=0;
  for(i=0; i<smp->n_verts; i++)
  {
    vp1=smp->vertex_array[i];
    if (vp1->merged_index==-1)
    {
/*
      fprintf(stderr,"%d running_index = %d\n",i,running_index);
*/
      vp1->merged_index=running_index;
      smp->merged_verts++;
      for (j=i+1; j<smp->n_verts && !distinguishable(vp1->projection,smp->vertex_array[j]->projection,epsilon); j++)
      {
        vp2=smp->vertex_array[j];
        if ( !distinguishable(vp1->vertex->x,vp2->vertex->x,epsilon)
             && !distinguishable(vp1->vertex->y,vp2->vertex->y,epsilon)
             && !distinguishable(vp1->vertex->z,vp2->vertex->z,epsilon) )
        {
          vp2->merged_index=running_index;
          if (vp1->vertex->x==vp2->vertex->x
              && vp1->vertex->y==vp2->vertex->y
              && vp1->vertex->z==vp2->vertex->z)
          {
            identical_verts++;
          }
          else
          {
            epsilon_verts++;
          }
        }
      }
      running_index++;
    }
  }

  fprintf(stderr,"  identical_verts = %d\n",identical_verts);
  fprintf(stderr,"  epsilon_verts = %d\n",epsilon_verts);
  fprintf(stderr,"  unique_verts = %d\n",smp->merged_verts);

  return;
}



void output_nurbs(op)
  struct object *op;
{
  struct nurbs *nrbp;
  struct double_list *dlp;
  struct ctlpt *cpp;
  struct ctlpt_list *clp;
  double val,umin,umax,vmin,vmax;
  int i,upts,vpts,uorder,vorder;

  nrbp=(struct nurbs *)op->contents;

  printf("NuPatch\n");

  upts=nrbp->upts;
  uorder=nrbp->uorder;
  printf("\t%d %d\n",upts,uorder);
  printf("\t[ ");
  i=0;
  dlp=nrbp->uknots;
  while(dlp) {
    val=dlp->val;
    if(i==uorder-1) {
      umin=val;
    }
    if(i==upts) {
      umax=val;
    }
    i++;
    printf("%.9g ",val);
    dlp=dlp->next;
  }
  printf("]\n");
  printf("\t%.9g %.9g\n",umin,umax);

  vpts=nrbp->vpts;
  vorder=nrbp->vorder;
  printf("\t%d %d\n",vpts,vorder);
  printf("\t[ ");
  i=0;
  dlp=nrbp->vknots;
  while(dlp) {
    val=dlp->val;
    if(i==vorder-1) {
      vmin=val;
    }
    if(i==vpts) {
      vmax=val;
    }
    i++;
    printf("%.9g ",val);
    dlp=dlp->next;
  }
  printf("]\n");
  printf("\t%.9g %.9g\n",vmin,vmax);

  clp=nrbp->ctlpts;
  if (clp) {
    if (clp->ctlpt->rational) {
      printf("\t\"Pw\" [\n");
    }
    else {
      printf("\t\"P\" [\n");
    }
  }
  while(clp) {
    cpp=clp->ctlpt;
    printf("\t%.9g %.9g %.9g",cpp->x,cpp->y,cpp->z);
    if(cpp->rational) {
      printf(" %.9g\n",cpp->w);
    }
    else {
      printf("\n");
    }
    clp=clp->next;
  }
  printf("\t]\n");
  
  return;
}


double polygon_area(struct polygon *pop)
{
  double area;
  u_int i;
  struct vector3 *v1,*v2,*v3,va,vb,vx;

  area = 0.0;

 /* For each triangle in the polygon, compute the normal vector*/

  v1 = pop->vertex_array[0]->vertex;
  for (i=0; i<pop->n_verts-2; i++)
  {

    v2 = pop->vertex_array[i+1]->vertex;
    v3 = pop->vertex_array[i+2]->vertex;

    vectorize(v2,v3,&va);
    vectorize(v2,v1,&vb);

 /* Compute the cross product */
    cross_prod(&va,&vb,&vx);

 /* Add the magnitude of the normal vector to the sum.  */
    area += 0.5*vect_length(&vx);
  }

  return(area);
}


/*
  Calculate volume of a manifold surface mesh by summing over
  all triangles of the surface mesh, the quantity:

              |  x1   x2   x3  |
    (1/6) det |  y1   y2   y3  |
              |  z1   z2   z3  |

  where (xn,yn,zn) is vertex n of a surface triangle.  The surface normals of
  all triangles of the surface mesh must be consistently oriented.
  Volume will be positive if all normals face outward, and negative if
  all normals face inward.
*/
double mesh_vol(struct surface_mesh *smp)
{
  struct polygon_list *plp;
  struct polygon *pop;
  struct vector3 *v1,*v2,*v3;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,det;
  double volume;
  u_int i;

  volume = 0.0;

 /* Sum over each polygon in the surface mesh */
 /* Assume polygon normals are consistently oriented */
 for (plp=smp->polygon_list; plp!=NULL; plp=plp->next)
 {
   pop=plp->polygon;

   /* For each triangle in the polygon, compute the normal vector*/
   v3 = pop->vertex_array[pop->n_verts-1]->vertex;
   for (i=0; i<pop->n_verts-2; i++)
   {
     v1 = pop->vertex_array[i]->vertex;
     v2 = pop->vertex_array[i+1]->vertex;

     x1=v1->x;
     y1=v1->y;
     z1=v1->z;
     x2=v2->x;
     y2=v2->y;
     z2=v2->z;
     x3=v3->x;
     y3=v3->y;
     z3=v3->z;

     /* compute determinant of oriented triangle */
     det=x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1);

     volume+=det;
   }

 }

  volume=volume/6.0;

  return(volume);
}


/*
  Calculate the aspect ratio of a triangle.
    The apsect ratio is defined as the length of longest edge
    divided by the shortest altitude.
*/
double triangle_aspect_ratio(struct polygon *pop)
{
  struct vector3 *v1,*v2,*v3,va,vb,vc,*vbase,*vopp,valt;
  double la,lb,lc;
  double amin,lmax,dot,aspect_ratio;

  /* Make triangle edge vectors */
  v1=pop->vertex_array[0]->vertex;
  v2=pop->vertex_array[1]->vertex;
  v3=pop->vertex_array[2]->vertex;
  vectorize(v1,v2,&va);
  vectorize(v2,v3,&vb);
  vectorize(v3,v1,&vc);

  /* Find length of longest edge */
  lmax=-DBL_MAX;
  la=vect_length(&va);
  lb=vect_length(&vb);
  lc=vect_length(&vc);
  if (la>lmax)
  {
    lmax=la;
    vbase=&va;
    vectorize(v1,v3,&vc);
    vopp=&vc;
  }
  if (lb>lmax)
  {
    lmax=lb;
    vbase=&vb;
    vectorize(v2,v1,&va);
    vopp=&va;
  }
  if (lc>lmax)
  {
    lmax=lc;
    vbase=&vc;
    vectorize(v3,v2,&vb);
    vopp=&vb;
  }

  /* Find shortest altitude */
  normalize(vbase);
  dot=dot_prod(vbase,vopp);
  valt.x=vopp->x-(dot*vbase->x);
  valt.y=vopp->y-(dot*vbase->y);
  valt.z=vopp->z-(dot*vbase->z);
  amin=vect_length(&valt);

  aspect_ratio = lmax/amin;

/*
  printf("\n");
  printf("la = %g  lb = %g  lc = %g\n",la,lb,lc);
  printf("lmax = %g  amin = %g  dot = %g\n",lmax,amin,dot);
  printf("aspect ratio = %g\n",aspect_ratio);
  printf("\n");
*/
 
  return(aspect_ratio); 
  
}


void bin_value(double value, struct histogram *histp)
{
  u_int start,mid,end;

  start=0;
  end=histp->size-1;

  while (end-start>1)
  {
    mid=(start+end)/2;
    if (value >= histp->bin_lower_limit[start]
        && value <= histp->bin_lower_limit[mid])
    {
      end=mid;
    }
    else
    {
      start=mid;
    }
  }
  if (value >= histp->bin_lower_limit[start]
       && value <= histp->bin_lower_limit[end])
  {
    histp->bin_count[start]++;
  }
  else
  {
    histp->bin_count[end]++;
  }
  
  return;
}


int analyze_surface_mesh(struct surface_mesh *smp)
{
  struct polygon_list *plp;
  struct polygon *pop;
  struct histogram *aspect_ratio_histogram;
  struct histogram *area_histogram;
  double area_avg,area_sd,area_sum,area_sum2,area_min,area_max,bin_width;
  struct vector3 bb_llf,bb_urb;
  int count;
  int i;

  
  if ((aspect_ratio_histogram=(struct histogram *)malloc(sizeof(struct histogram)))==NULL) {
    printf("meshalyzer: out of memory while creating aspect ratio histogram\n");
    return(1);
  }
  aspect_ratio_histogram->size=16;
  if ((aspect_ratio_histogram->bin_lower_limit=(double *)malloc(aspect_ratio_histogram->size*sizeof(double)))==NULL) {
    printf("meshalyzer: out of memory while creating aspect ratio histogram\n");
    return(1);
  }
  if ((aspect_ratio_histogram->bin_count=(u_int *)malloc(aspect_ratio_histogram->size*sizeof(u_int)))==NULL) {
    printf("meshalyzer: out of memory while creating aspect ratio histogram\n");
    return(1);
  }
  aspect_ratio_histogram->bin_lower_limit[0]=2.0/sqrt(3.0);
  aspect_ratio_histogram->bin_lower_limit[1]=1.5;
  aspect_ratio_histogram->bin_lower_limit[2]=2.0;
  aspect_ratio_histogram->bin_lower_limit[3]=2.5;
  aspect_ratio_histogram->bin_lower_limit[4]=3.0;
  aspect_ratio_histogram->bin_lower_limit[5]=4.0;
  aspect_ratio_histogram->bin_lower_limit[6]=6.0;
  aspect_ratio_histogram->bin_lower_limit[7]=10.0;
  aspect_ratio_histogram->bin_lower_limit[8]=15.0;
  aspect_ratio_histogram->bin_lower_limit[9]=25.0;
  aspect_ratio_histogram->bin_lower_limit[10]=50.0;
  aspect_ratio_histogram->bin_lower_limit[11]=100.0;
  aspect_ratio_histogram->bin_lower_limit[12]=300.0;
  aspect_ratio_histogram->bin_lower_limit[13]=1000.0;
  aspect_ratio_histogram->bin_lower_limit[14]=10000.0;
  aspect_ratio_histogram->bin_lower_limit[15]=100000.0;
  for (i=0; i<aspect_ratio_histogram->size; i++)
  {
    aspect_ratio_histogram->bin_count[i]=0;
  }

  smp->volume=mesh_vol(smp);
  printf("Mesh volume: %g\n",smp->volume);

  count=0;
  area_avg=0.0;
  area_sd=0.0;
  area_sum=0.0;
  area_sum2=0.0;
  area_min=DBL_MAX;
  area_max=-DBL_MAX;
  bb_llf.x=DBL_MAX;
  bb_llf.y=DBL_MAX;
  bb_llf.z=DBL_MAX;
  bb_urb.x=-DBL_MAX;
  bb_urb.y=-DBL_MAX;
  bb_urb.z=-DBL_MAX;
  for (plp=smp->polygon_list; plp!=NULL; plp=plp->next)
  {
    count++;
    pop=plp->polygon;
    pop->area=polygon_area(pop);
    //area_avg=(area_avg*count/(count+1))+pop->area/(count+1);
    area_avg=(area_avg*(count-1)/count)+pop->area/count;
    area_sum+=pop->area;
    area_sum2=area_sum2+(pop->area*pop->area);
    if (pop->area < area_min) area_min=pop->area;
    if (pop->area > area_max) area_max=pop->area;
    pop->aspect_ratio=triangle_aspect_ratio(pop);
    bin_value(pop->aspect_ratio,aspect_ratio_histogram);
    for (i=0;i<pop->n_verts;i++)
    {
      if (pop->vertex_array[i]->vertex->x<bb_llf.x) bb_llf.x=pop->vertex_array[i]->vertex->x;
      if (pop->vertex_array[i]->vertex->y<bb_llf.y) bb_llf.y=pop->vertex_array[i]->vertex->y;
      if (pop->vertex_array[i]->vertex->z<bb_llf.z) bb_llf.z=pop->vertex_array[i]->vertex->z;
      if (pop->vertex_array[i]->vertex->x>bb_urb.x) bb_urb.x=pop->vertex_array[i]->vertex->x;
      if (pop->vertex_array[i]->vertex->y>bb_urb.y) bb_urb.y=pop->vertex_array[i]->vertex->y;
      if (pop->vertex_array[i]->vertex->z>bb_urb.z) bb_urb.z=pop->vertex_array[i]->vertex->z;
    }
  }
  area_sd=sqrt((area_sum2-(area_sum*area_sum/count))/(count-1));

  printf("Mesh bounding box: [ %.15g, %.15g, %.15g ] [ %.15g, %.15g, %.15g ]\n",bb_llf.x,bb_llf.y,bb_llf.z,bb_urb.x,bb_urb.y,bb_urb.z);

  if ((area_histogram=(struct histogram *)malloc(sizeof(struct histogram)))==NULL) {
    printf("meshalyzer: out of memory while creating aspect ratio histogram\n");
    return(1);
  }
  area_histogram->size=16;
  if ((area_histogram->bin_lower_limit=(double *)malloc(area_histogram->size*sizeof(double)))==NULL) {
    printf("meshalyzer: out of memory while creating aspect ratio histogram\n");
    return(1);
  }
  if ((area_histogram->bin_count=(u_int *)malloc(area_histogram->size*sizeof(u_int)))==NULL) {
    printf("meshalyzer: out of memory while creating aspect ratio histogram\n");
    return(1);
  }
  bin_width=(4.0*area_sd)/(area_histogram->size-1.0);
  printf("\bbin_width = %.12g\n",bin_width);
  printf("\barea_avg = %.12g\n",area_avg);
  printf("\barea_sum = %.12g\n",area_sum);
  printf("\bcount = %d\n",count);
  area_histogram->bin_lower_limit[0]=0;
  area_histogram->bin_count[0]=0;
  area_histogram->bin_lower_limit[1]=0;
  area_histogram->bin_count[1]=0;
  for (i=2; i<area_histogram->size; i++)
  {
    area_histogram->bin_lower_limit[i]=area_avg+(((i-1)-(area_histogram->size/2.0))*bin_width);
    area_histogram->bin_count[i]=0;
  }

  for (plp=smp->polygon_list; plp!=NULL; plp=plp->next)
  {
    pop=plp->polygon;
    bin_value(pop->area,area_histogram);
  }

  printf("\n");
  printf("Mesh quality statistics:\n\n");
  printf("  Avg triangle area = %g +-%g\n",area_avg,area_sd);
  printf("  Smallest area:       %10.5g   |  Largest area:       %10.5g\n\n",area_min,area_max);
  printf("  Triangle aspect ratio histogram:\n");
  printf("  %.5g - %g       : %9d    | %6g - %g         : %9d\n",aspect_ratio_histogram->bin_lower_limit[0],aspect_ratio_histogram->bin_lower_limit[1],aspect_ratio_histogram->bin_count[0],aspect_ratio_histogram->bin_lower_limit[8],aspect_ratio_histogram->bin_lower_limit[9],aspect_ratio_histogram->bin_count[8]);
  printf("  %6g - %g         : %9d    | %6g - %g         : %9d\n",aspect_ratio_histogram->bin_lower_limit[1],aspect_ratio_histogram->bin_lower_limit[2],aspect_ratio_histogram->bin_count[1],aspect_ratio_histogram->bin_lower_limit[9],aspect_ratio_histogram->bin_lower_limit[10],aspect_ratio_histogram->bin_count[9]);
  printf("  %6g - %g       : %9d    | %6g - %g        : %9d\n",aspect_ratio_histogram->bin_lower_limit[2],aspect_ratio_histogram->bin_lower_limit[3],aspect_ratio_histogram->bin_count[2],aspect_ratio_histogram->bin_lower_limit[10],aspect_ratio_histogram->bin_lower_limit[11],aspect_ratio_histogram->bin_count[10]);
  printf("  %6g - %g         : %9d    | %6g - %g        : %9d\n",aspect_ratio_histogram->bin_lower_limit[3],aspect_ratio_histogram->bin_lower_limit[4],aspect_ratio_histogram->bin_count[3],aspect_ratio_histogram->bin_lower_limit[11],aspect_ratio_histogram->bin_lower_limit[12],aspect_ratio_histogram->bin_count[11]);
  printf("  %6g - %g         : %9d    | %6g - %g       : %9d\n",aspect_ratio_histogram->bin_lower_limit[4],aspect_ratio_histogram->bin_lower_limit[5],aspect_ratio_histogram->bin_count[4],aspect_ratio_histogram->bin_lower_limit[12],aspect_ratio_histogram->bin_lower_limit[13],aspect_ratio_histogram->bin_count[12]);
  printf("  %6g - %g         : %9d    | %6g - %g      : %9d\n",aspect_ratio_histogram->bin_lower_limit[5],aspect_ratio_histogram->bin_lower_limit[6],aspect_ratio_histogram->bin_count[5],aspect_ratio_histogram->bin_lower_limit[13],aspect_ratio_histogram->bin_lower_limit[14],aspect_ratio_histogram->bin_count[13]);
  printf("  %6g - %g        : %9d    | %6g - %g     : %9d\n",aspect_ratio_histogram->bin_lower_limit[6],aspect_ratio_histogram->bin_lower_limit[7],aspect_ratio_histogram->bin_count[6],aspect_ratio_histogram->bin_lower_limit[14],aspect_ratio_histogram->bin_lower_limit[15],aspect_ratio_histogram->bin_count[14]);
  printf("  %6g - %g        : %9d    | %6g -            : %9d\n",aspect_ratio_histogram->bin_lower_limit[7],aspect_ratio_histogram->bin_lower_limit[8],aspect_ratio_histogram->bin_count[7],aspect_ratio_histogram->bin_lower_limit[15],aspect_ratio_histogram->bin_count[15]);
  printf("  (Aspect ratio is longest edge divided by shortest altitude)\n");
  printf("\n");

  printf("  Triangle area histogram:\n");
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",area_histogram->bin_lower_limit[0],area_histogram->bin_lower_limit[1],area_histogram->bin_count[0],area_histogram->bin_lower_limit[8],area_histogram->bin_lower_limit[9],area_histogram->bin_count[8]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",area_histogram->bin_lower_limit[1],area_histogram->bin_lower_limit[2],area_histogram->bin_count[1],area_histogram->bin_lower_limit[9],area_histogram->bin_lower_limit[10],area_histogram->bin_count[9]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",area_histogram->bin_lower_limit[2],area_histogram->bin_lower_limit[3],area_histogram->bin_count[2],area_histogram->bin_lower_limit[10],area_histogram->bin_lower_limit[11],area_histogram->bin_count[10]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",area_histogram->bin_lower_limit[3],area_histogram->bin_lower_limit[4],area_histogram->bin_count[3],area_histogram->bin_lower_limit[11],area_histogram->bin_lower_limit[12],area_histogram->bin_count[11]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",area_histogram->bin_lower_limit[4],area_histogram->bin_lower_limit[5],area_histogram->bin_count[4],area_histogram->bin_lower_limit[12],area_histogram->bin_lower_limit[13],area_histogram->bin_count[12]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",area_histogram->bin_lower_limit[5],area_histogram->bin_lower_limit[6],area_histogram->bin_count[5],area_histogram->bin_lower_limit[13],area_histogram->bin_lower_limit[14],area_histogram->bin_count[13]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",area_histogram->bin_lower_limit[6],area_histogram->bin_lower_limit[7],area_histogram->bin_count[6],area_histogram->bin_lower_limit[14],area_histogram->bin_lower_limit[15],area_histogram->bin_count[14]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",area_histogram->bin_lower_limit[7],area_histogram->bin_lower_limit[8],area_histogram->bin_count[7],area_histogram->bin_lower_limit[15],area_max,area_histogram->bin_count[15]);
  printf("\n");

  return(0);
}


void output_surface_mesh(struct surface_mesh *smp)
{
  struct vertex *vp;
  struct vector3 *vecp;
  struct polygon_list *plp;
  struct polygon *pop;
  int running_index;
  int total_face_count;
  int i;

  running_index=0;
  for (i=0; i<smp->n_verts; i++)
  {
    vp=smp->vertex_array[i];
    if (vp->merged_index>=running_index)
    {
      vecp=vp->vertex;
      printf("Vertex %d %.15g %.15g %.15g\n",running_index+1,vecp->x,vecp->y,vecp->z);
      running_index++;
    }
  }

  total_face_count=0;
  for (plp=smp->polygon_list; plp!=NULL; plp=plp->next)
  {
    pop=plp->polygon;

    total_face_count++;
    printf("Face %d",total_face_count);
    for (i=0; i<pop->n_verts; i++)
    {
      printf(" %d",pop->vertex_array[i]->merged_index+1);
    }
    printf("\n");
  }

  return;
}

