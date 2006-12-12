/*
MCell (tm) Version 2.08 5/18/1998

Copyright (C) 1997,1998, The Salk Institute & Cornell University.
MCell was written jointly by T.M. Bartol Jr. & J.R. Stiles,
with design input for Monte Carlo algorithms from E.E. Salpeter.

  Acknowledgements:
    T.J. Sejnowski for development input and support
    (NSF Grant IBN-9603611), and M.M. Salpeter for fostering
    quantitative experimental applications.  Additional support
    from NIH Grant K08NS01776 (J.R. Stiles).

MCell is a scientific simulation software tool distributed freely to
registered beta-test research laboratories, and must be obtained as a
machine-specific executable program from the authors.  Copying and/or
editing of MCell without the authors' consent is expressly forbidden.
MCell is provided as is and the authors assume no responsibility for
simulation results obtained by users.  Any material published from
MCell simulations must acknowledge the authors and granting
institutions.
*/

/* 3D vector routines */

#include <math.h>
#include "vector.h"

#define MY_PI 3.14159265358979323846

void mult_matrix(m1,m2,om,l,m,n)
double m1[][4],m2[][4],om[][4];
unsigned short l,m,n;
{
  double tm[4][4];
  unsigned short i,j,k;

  for (i=0;i<l;i++) {
    for (j=0;j<n;j++) {
      tm[i][j]=0;
      for (k=0;k<m;k++) {
	tm[i][j]=tm[i][j]+(m1[i][k])*(m2[k][j]);
      }
    }
  }
  for (i=0;i<l;i++) {
    for (j=0;j<n;j++) {
      om[i][j]=tm[i][j];
    }
  }
}


void normalize(v)
struct vector3 *v;
{
  double length;

  length=vect_length(v);
  v->x=v->x/length;
  v->y=v->y/length;
  v->z=v->z/length;
}


void init_matrix(im)
double im[][4];
{

  im[0][0]=1;
  im[0][1]=0;
  im[0][2]=0;
  im[0][3]=0;
  im[1][0]=0;
  im[1][1]=1;
  im[1][2]=0;
  im[1][3]=0;
  im[2][0]=0;
  im[2][1]=0;
  im[2][2]=1;
  im[2][3]=0;
  im[3][0]=0;
  im[3][1]=0;
  im[3][2]=0;
  im[3][3]=1;
}

void scale_matrix(im,om,scale)
double im[][4];
double om[][4];
struct vector3 *scale;
{
  double sc[4][4];
  unsigned short l,m,n;

  sc[0][0]=scale->x;
  sc[0][1]=0;
  sc[0][2]=0;
  sc[0][3]=0;
  sc[1][0]=0;
  sc[1][1]=scale->y;
  sc[1][2]=0;
  sc[1][3]=0;
  sc[2][0]=0;
  sc[2][1]=0;
  sc[2][2]=scale->z;
  sc[2][3]=0;
  sc[3][0]=0;
  sc[3][1]=0;
  sc[3][2]=0;
  sc[3][3]=1;

  l=4;
  m=4;
  n=4;
  mult_matrix(im,sc,om,l,m,n);
}

void translate_matrix(im,om,translate)
double im[][4];
double om[][4];
struct vector3 *translate;
{
  double tm[4][4];
  unsigned short l,m,n;

  tm[0][0]=1;
  tm[0][1]=0;
  tm[0][2]=0;
  tm[0][3]=0;
  tm[1][0]=0;
  tm[1][1]=1;
  tm[1][2]=0;
  tm[1][3]=0;
  tm[2][0]=0;
  tm[2][1]=0;
  tm[2][2]=1;
  tm[2][3]=0;
  tm[3][0]=translate->x;
  tm[3][1]=translate->y;
  tm[3][2]=translate->z;
  tm[3][3]=1;

  l=4;
  m=4;
  n=4;
  mult_matrix(im,tm,om,l,m,n);
}

void rotate_matrix(im,om,axis,angle)
double im[][4];
double om[][4];
struct vector3 *axis;
double angle;
{
  double r1[4][4],r2[4][4],r3[4][4],rm[4][4];
  double a,b,c,v;
  double rad;
  unsigned short l,m,n;

  normalize(axis);
  a=axis->x;
  b=axis->y;
  c=axis->z;
  v=sqrt(b*b+c*c);

  r1[0][0]=1;
  r1[0][1]=0;
  r1[0][2]=0;
  r1[0][3]=0;
  r1[1][0]=0;
  r1[1][1]=1;
  r1[1][2]=0;
  r1[1][3]=0;
  r1[2][0]=0;
  r1[2][1]=0;
  r1[2][2]=1;
  r1[2][3]=0;
  r1[3][0]=0;
  r1[3][1]=0;
  r1[3][2]=0;
  r1[3][3]=1;

  if (v!=0.0) {
    r1[1][1]=c/v;
    r1[1][2]=b/v;
    r1[2][1]=-b/v;
    r1[2][2]=c/v;
  }

  r2[0][0]=v;
  r2[0][1]=0;
  r2[0][2]=a;
  r2[0][3]=0;
  r2[1][0]=0;
  r2[1][1]=1;
  r2[1][2]=0;
  r2[1][3]=0;
  r2[2][0]=-a;
  r2[2][1]=0;
  r2[2][2]=v;
  r2[2][3]=0;
  r2[3][0]=0;
  r2[3][1]=0;
  r2[3][2]=0;
  r2[3][3]=1;

  rad=MY_PI/180.0;
  r3[0][0]=cos(angle*rad);
  r3[0][1]=sin(angle*rad);
  r3[0][2]=0;
  r3[0][3]=0;
  r3[1][0]=-sin(angle*rad);
  r3[1][1]=cos(angle*rad);
  r3[1][2]=0;
  r3[1][3]=0;
  r3[2][0]=0;
  r3[2][1]=0;
  r3[2][2]=1;
  r3[2][3]=0;
  r3[3][0]=0;
  r3[3][1]=0;
  r3[3][2]=0;
  r3[3][3]=1;

  l=4;
  m=4;
  n=4;
  mult_matrix(r1,r2,rm,l,m,n);
  mult_matrix(rm,r3,rm,l,m,n);

  r2[0][2]=-a;
  r2[2][0]=a;

  if (v!=0.0) {
    r1[1][2]=-b/v;
    r1[2][1]=b/v;
  }

  mult_matrix(rm,r2,rm,l,m,n);
  mult_matrix(rm,r1,rm,l,m,n);
  mult_matrix(im,rm,om,l,m,n);
}

void tform_matrix(scale,translate,axis,angle,om)
struct vector3 *scale;
struct vector3 *translate;
struct vector3 *axis;
double angle;
double om[4][4];
{
  double im[4][4];
  double sc[4][4];
  double tm[4][4];
  double r1[4][4],r2[4][4],r3[4][4];
  double a,b,c,v;
  double rad;
  unsigned short l,m,n;

  init_matrix(om);

  sc[0][0]=scale->x;
  sc[0][1]=0;
  sc[0][2]=0;
  sc[0][3]=0;
  sc[1][0]=0;
  sc[1][1]=scale->y;
  sc[1][2]=0;
  sc[1][3]=0;
  sc[2][0]=0;
  sc[2][1]=0;
  sc[2][2]=scale->z;
  sc[2][3]=0;
  sc[3][0]=0;
  sc[3][1]=0;
  sc[3][2]=0;
  sc[3][3]=1;

  tm[0][0]=1;
  tm[0][1]=0;
  tm[0][2]=0;
  tm[0][3]=0;
  tm[1][0]=0;
  tm[1][1]=1;
  tm[1][2]=0;
  tm[1][3]=0;
  tm[2][0]=0;
  tm[2][1]=0;
  tm[2][2]=1;
  tm[2][3]=0;
  tm[3][0]=translate->x;
  tm[3][1]=translate->y;
  tm[3][2]=translate->z;
  tm[3][3]=1;

  normalize(axis);
  a=axis->x;
  b=axis->y;
  c=axis->z;
  v=sqrt(b*b+c*c);

  r1[0][0]=1;
  r1[0][1]=0;
  r1[0][2]=0;
  r1[0][3]=0;
  r1[1][0]=0;
  r1[1][1]=1;
  r1[1][2]=0;
  r1[1][3]=0;
  r1[2][0]=0;
  r1[2][1]=0;
  r1[2][2]=1;
  r1[2][3]=0;
  r1[3][0]=0;
  r1[3][1]=0;
  r1[3][2]=0;
  r1[3][3]=1;

  if (v!=0.0) {
    r1[1][1]=c/v;
    r1[1][2]=b/v;
    r1[2][1]=-b/v;
    r1[2][2]=c/v;
  }

  r2[0][0]=v;
  r2[0][1]=0;
  r2[0][2]=a;
  r2[0][3]=0;
  r2[1][0]=0;
  r2[1][1]=1;
  r2[1][2]=0;
  r2[1][3]=0;
  r2[2][0]=-a;
  r2[2][1]=0;
  r2[2][2]=v;
  r2[2][3]=0;
  r2[3][0]=0;
  r2[3][1]=0;
  r2[3][2]=0;
  r2[3][3]=1;

  rad=MY_PI/180.0;
  r3[0][0]=cos(angle*rad);
  r3[0][1]=sin(angle*rad);
  r3[0][2]=0;
  r3[0][3]=0;
  r3[1][0]=-sin(angle*rad);
  r3[1][1]=cos(angle*rad);
  r3[1][2]=0;
  r3[1][3]=0;
  r3[2][0]=0;
  r3[2][1]=0;
  r3[2][2]=1;
  r3[2][3]=0;
  r3[3][0]=0;
  r3[3][1]=0;
  r3[3][2]=0;
  r3[3][3]=1;

  l=4;
  m=4;
  n=4;
  mult_matrix(r1,r2,om,l,m,n);
  mult_matrix(om,r3,om,l,m,n);

  r2[0][2]=-a;
  r2[2][0]=a;

  if (v!=0.0) {
    r1[1][2]=-b/v;
    r1[2][1]=b/v;
  }

  mult_matrix(om,r2,om,l,m,n);
  mult_matrix(om,r1,om,l,m,n);
  mult_matrix(om,sc,om,l,m,n);
  mult_matrix(om,tm,om,l,m,n);
}


void vectorize(p1,p2,v)
struct vector3 *p1,*p2,*v;
{

  v->x=p2->x-p1->x;
  v->y=p2->y-p1->y;
  v->z=p2->z-p1->z;
}


double vect_length(v)
struct vector3 *v;
{
  double length;

  length=sqrt((v->x)*(v->x)+(v->y)*(v->y)+(v->z)*(v->z));
  return(length);
}


double dot_prod(v1,v2)
struct vector3 *v1,*v2;
{
  double dot;

  dot=(v1->x)*(v2->x)+(v1->y)*(v2->y)+(v1->z)*(v2->z);
  return(dot);
}


void cross_prod(v1,v2,v3)
struct vector3 *v1,*v2,*v3;
{

  v3->x=(v1->y)*(v2->z)-(v1->z)*(v2->y);
  v3->y=(v1->z)*(v2->x)-(v1->x)*(v2->z);
  v3->z=(v1->x)*(v2->y)-(v1->y)*(v2->x);
}

#undef MY_PI
