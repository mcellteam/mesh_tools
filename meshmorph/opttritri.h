#ifndef OPTTRITRI_H
#define OPTTRITRI_H

/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 *
 * Updated June 1999: removed the divisions -- a little faster now!
 * Updated October 1999: added {} to CROSS and SUB macros 
 *
 * int NoDivTriTriIsect(double V0[3],double V1[3],double V2[3],
 *                      double U0[3],double U1[3],double U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 */

#include <math.h>

#include "meshmorph.h"

#define FABS(x) (fabs(x))        /* implement as is fastest on your machine */

//#define MYEPSILON 1E-5

/* if USE_EPSILON_TEST is true then we do a check:
   if |dv|<EPSILON then dv=0.0;
   else no check is done (which is less robust)
   */
#define USE_EPSILON_TEST TRUE

/* sort so that a<=b */
#define SORT(a,b)       \
      if(a>b)    \
{          \
  double c; \
  c=a;     \
  a=b;     \
  b=c;     \
}


/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 */
#define EDGE_EDGE_TEST(V0,U0,U1)                      \
      Bx=U0->p[i0]-U1->p[i0];                                   \
By=U0->p[i1]-U1->p[i1];                                   \
Cx=V0->p[i0]-U0->p[i0];                                   \
Cy=V0->p[i1]-U0->p[i1];                                   \
f=Ay*Bx-Ax*By;                                      \
d=By*Cx-Bx*Cy;                                      \
if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
{                                                   \
  e=Ax*Cy-Ay*Cx;                                    \
  if(f>0)                                           \
  {                                                 \
    if(e>=0 && e<=f) return 1;                      \
  }                                                 \
  else                                              \
  {                                                 \
    if(e<=0 && e>=f) return 1;                      \
  }                                                 \
}

#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
{                                              \
  double Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
  Ax=V1->p[i0]-V0->p[i0];                            \
  Ay=V1->p[i1]-V0->p[i1];                            \
  /* test edge U0,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U0,U1);                    \
  /* test edge U1,U2 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U1,U2);                    \
  /* test edge U2,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U2,U0);                    \
}

#define POINT_IN_TRI(V0,U0,U1,U2)           \
{                                           \
  double a,b,c,d0,d1,d2;                     \
  /* is T1 completly inside T2? */          \
  /* check if V0 is inside tri(U0,U1,U2) */ \
  a=U1->p[i1]-U0->p[i1];                          \
  b=-(U1->p[i0]-U0->p[i0]);                       \
  c=-a*U0->p[i0]-b*U0->p[i1];                     \
  d0=a*V0->p[i0]+b*V0->p[i1]+c;                   \
  \
  a=U2->p[i1]-U1->p[i1];                          \
  b=-(U2->p[i0]-U1->p[i0]);                       \
  c=-a*U1->p[i0]-b*U1->p[i1];                     \
  d1=a*V0->p[i0]+b*V0->p[i1]+c;                   \
  \
  a=U0->p[i1]-U2->p[i1];                          \
  b=-(U0->p[i0]-U2->p[i0]);                       \
  c=-a*U2->p[i0]-b*U2->p[i1];                     \
  d2=a*V0->p[i0]+b*V0->p[i1]+c;                   \
  if(d0*d1>0.0)                             \
  {                                         \
    if(d0*d2>0.0) return 1;                 \
  }                                         \
}

int coplanar_tri_tri(vector3 & N,vector3 const * V0,vector3 const * V1,vector3 const * V2,
                     vector3 const * U0,vector3 const * U1,vector3 const * U2);



#define NEWCOMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,A,B,C,X0,X1) \
{ \
  if(D0D1>0.0f) \
  { \
    /* here we know that D0D2<=0.0 */ \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
  } \
  else if(D0D2>0.0f)\
  { \
    /* here we know that d0d1<=0.0 */ \
    A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
  } \
  else if(D1*D2>0.0f || D0!=0.0f) \
  { \
    /* here we know that d0d1<=0.0 or that D0!=0.0 */ \
    A=VV0; B=(VV1-VV0)*D0; C=(VV2-VV0)*D0; X0=D0-D1; X1=D0-D2; \
  } \
  else if(D1!=0.0f) \
  { \
    A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
  } \
  else if(D2!=0.0f) \
  { \
    A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
  } \
  else \
  { \
    /* triangles are coplanar */ \
    return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2); \
  } \
}



int NoDivTriTriIsect(vector3 const * V0,vector3 const * V1,vector3 const * V2,
                     vector3 const * U0,vector3 const * U1,vector3 const * U2);

int intersect_triangle3(vector3 const * orig, vector3 const * end, vector3 const * normal,
                        vector3 const * vert0, vector3 const * vert1, vector3 const * vert2,
                        result & r);
#endif
