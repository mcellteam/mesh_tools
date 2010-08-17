#include "opttritri.h"

#include "controls.h"

/** Determine whether two coplanar triangles intersect.
 * \param[in] N Not sure what this is.
 * \param[in] V0 First vertex of first triangle.
 * \param[in] V1 second vertex of first triangle.
 * \param[in] V2 Third vertex of first triangle.
 * \param[in] U0 First vertex of second triangle.
 * \param[in] U1 second vertex of second triangle.
 * \param[in] U2 Third vertex of second triangle.
 * \return 1 if triangles intersect; 0 otherwise.
 */

int coplanar_tri_tri(vector3 & N,vector3 const * V0,vector3 const * V1,vector3 const * V2,
                                 vector3 const * U0,vector3 const * U1,vector3 const * U2)
{
   vector3 A;
   short i0,i1;
   /* first project onto an axis-aligned plane, that maximizes the area */
   /* of the triangles, compute indices: i0,i1. */
   A.p[0]=FABS(N.p[0]);
   A.p[1]=FABS(N.p[1]);
   A.p[2]=FABS(N.p[2]);
   if(A.p[0]>A.p[1])
   {
      if(A.p[0]>A.p[2])
      {
          i0=1;      /* A[0] is greatest */
          i1=2;
      }
      else
      {
          i0=0;      /* A[2] is greatest */
          i1=1;
      }
   }
   else   /* A[0]<=A[1] */
   {
      if(A.p[2]>A.p[1])
      {
          i0=0;      /* A[2] is greatest */
          i1=1;
      }
      else
      {
          i0=0;      /* A[1] is greatest */
          i1=2;
      }
    }

    /* test all edges of triangle 1 against the edges of triangle 2 */
    EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V1,V2,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V2,V0,U0,U1,U2);

    /* finally, test if tri1 is totally contained in tri2 or vice versa */
    POINT_IN_TRI(V0,U0,U1,U2);
    POINT_IN_TRI(U0,V0,V1,V2);

    return 0;
}

/** Determine whether two triangles intersect.
 * \param[in] V0 First vertex of first triangle.
 * \param[in] V1 second vertex of first triangle.
 * \param[in] V2 Third vertex of first triangle.
 * \param[in] U0 First vertex of second triangle.
 * \param[in] U1 second vertex of second triangle.
 * \param[in] U2 Third vertex of second triangle.
 * \return 1 if triangles intersect; 0 otherwise.
 */

int NoDivTriTriIsect(vector3 const * V0,vector3 const * V1,vector3 const * V2,
                     vector3 const * U0,vector3 const * U1,vector3 const * U2)
//                     bool myflag)
{
  Controls & cs(Controls::instance());
  //vector3 E2;
  double d1,d2;
  double du0,du1,du2,dv0,dv1,dv2;
  double isect1[2], isect2[2];
  double du0du1,du0du2,dv0dv1,dv0dv2;
  short index;
  double vp0,vp1,vp2;
  double up0,up1,up2;
  double bb,cc,max;

  /* compute plane equation of triangle(V0,V1,V2) */
  vector3 E1(*V1-*V0);
  vector3 E2(*V2-*V0);
  vector3 N1(E1.cross(E2));
  d1=-N1.dot(*V0);
  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0=N1.dot(*U0)+d1;
  du1=N1.dot(*U1)+d1;
  du2=N1.dot(*U2)+d1;

  /* coplanarity robustness check */
#if USE_EPSILON_TEST==TRUE
  if(FABS(du0)<cs.get_my_double_epsilon()) du0=0.0;
  if(FABS(du1)<cs.get_my_double_epsilon()) du1=0.0;
  if(FABS(du2)<cs.get_my_double_epsilon()) du2=0.0;
  //if(FABS(du0)<MYEPSILON) du0=0.0;
  //if(FABS(du1)<MYEPSILON) du1=0.0;
  //if(FABS(du2)<MYEPSILON) du2=0.0;
#endif
  du0du1=du0*du1;
  du0du2=du0*du2;

  if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  E1 = *U1-*U0;
  E2 = *U2-*U0;
  vector3 N2(E1.cross(E2));
  d2=-N2.dot(*U0);
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  dv0=N2.dot(*V0)+d2;
  dv1=N2.dot(*V1)+d2;
  dv2=N2.dot(*V2)+d2;

#if USE_EPSILON_TEST==TRUE
  if(FABS(dv0)<cs.get_my_double_epsilon()) dv0=0.0;
  if(FABS(dv1)<cs.get_my_double_epsilon()) dv1=0.0;
  if(FABS(dv2)<cs.get_my_double_epsilon()) dv2=0.0;
  //if(FABS(dv0)<MYEPSILON) dv0=0.0;
  //if(FABS(dv1)<MYEPSILON) dv1=0.0;
  //if(FABS(dv2)<MYEPSILON) dv2=0.0;
#endif

  dv0dv1=dv0*dv1;
  dv0dv2=dv0*dv2;

  if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute direction of intersection line */
  vector3 D(N1.cross(N2));

  /* compute and index to the largest component of D */
  max=(double)FABS(D.p[0]);
  index=0;
  bb=(double)FABS(D.p[1]);
  cc=(double)FABS(D.p[2]);
  if(bb>max) max=bb,index=1;
  if(cc>max) max=cc,index=2;

  /* this is the simplified projection onto L*/
  vp0=V0->p[index];
  vp1=V1->p[index];
  vp2=V2->p[index];

  up0=U0->p[index];
  up1=U1->p[index];
  up2=U2->p[index];

  /* compute interval for triangle 1 */
  double a,b,myc,x0,x1;
  NEWCOMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,a,b,myc,x0,x1);

  /* compute interval for triangle 2 */
  double d,e,f,y0,y1;
  NEWCOMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,d,e,f,y0,y1);

  double xx,yy,xxyy,tmp;
  xx=x0*x1;
  yy=y0*y1;
  xxyy=xx*yy;

  tmp=a*xxyy;
  isect1[0]=tmp+b*x1*yy;
  isect1[1]=tmp+myc*x0*yy;

  tmp=d*xxyy;
  isect2[0]=tmp+e*xx*y1;
  isect2[1]=tmp+f*xx*y0;

  SORT(isect1[0],isect1[1]);
  SORT(isect2[0],isect2[1]);

//  if (myflag)
//  {
//    // DEBUG
//    printf("du0 = %g, du1 = %g, du2 = %g, dv0 = %g, dv1 = %g, dv2 = %g, isect1[0] = %g, isect1[1] = %g, isect2[0] = %g, isect2[1] = %g\n",du0,du1,du2,dv0,dv1,dv2,isect1[0],isect1[1],isect2[0],isect2[1]);
//    // DEBUG
//  }

  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
  return 1;
}

/** Determine if line segment and triangle intersect.
 * \param[in] orig One end of line segment.
 * \param[in] end Other end of line segment.
 * \param[in] normal Normal vector of triangle.
 * \param[in] vert0 First vertex of triangle.
 * \param[in] vert1 second vertex of triangle.
 * \param[in] vert2 Third vertex of triangle.
 * \param[out] r Record whether line intersects triange along line segment;
 * whether line intersects triangle; and whether line intersects edge of triangle.
 * \return 1 if triangles intersect; 0 otherwise.
 */

int intersect_triangle3(vector3 const * orig, vector3 const * end, vector3 const * normal,
			vector3 const * vert0, vector3 const * vert1, vector3 const * vert2,
                        result & r)
{
  Controls & cs(Controls::instance());
  // does line segment cross plane of face?
  // compute d = -normal dot triangle_vertex
  double d = -normal->dot(*vert0);
  // compute d_origin = normal dot origin + d
  double d_origin = normal->dot(*orig)+d;
  // compute d_end = normal dot end + d
  double d_end = normal->dot(*end)+d;
  // if both d_origin and d_end are nonzero and have same sign,
  // then line segment does not cross plane of triangle
  // and line segment cannot possibly intersect triangle
  if( ((d_origin > cs.get_my_double_epsilon()) && (d_end > cs.get_my_double_epsilon())) ||
      ((d_origin < cs.get_my_double_epsilon()) && (d_end < cs.get_my_double_epsilon())) )
    return 0;
 
  r.line_flag=true;

  double det,inv_det,u,v;

  /* compute line segment direction */
  vector3 dir(*orig-*end);

  /* find vectors for two edges sharing vert0 */
  vector3 edge1(*vert1-*vert0);
  vector3 edge2(*vert2-*vert0);

  /* begin calculating determinant - also used to calculate U parameter */
  vector3 pvec(dir.cross(edge2));

  /* if determinant is near zero, ray lies in plane of triangle */
  det = edge1.dot(pvec);

  /* calculate distance from vert0 to ray origin */
  vector3 tvec(*orig-*vert0);
  inv_det = 1.0 / det;

  vector3 qvec(tvec.cross(edge1));

  if (det > cs.get_my_double_epsilon())
  {
    u = tvec.dot(pvec);
    if (u < 0.0 || u > det)
      return 0;

    /* calculate V parameter and test bounds */
    v = dir.dot(qvec);
    if (v < 0.0 || u + v > det)
      return 0;

  }
  else if(det < -cs.get_my_double_epsilon())
  {
    /* calculate U parameter and test bounds */
    u = tvec.dot(pvec);
    if (u > 0.0 || u < det)
      return 0;

    /* calculate V parameter and test bounds */
    v = dir.dot(qvec) ;
    if (v > 0.0 || u + v < det)
      return 0;
  }
  else return 0;  /* ray is parallell to the plane of the triangle */

  if (fabs(u)<cs.get_my_double_epsilon() || fabs(v)<cs.get_my_double_epsilon() ||
      distinguishable(u+v,det)==false) r.poly_edge_flag=true;
  r.poly_flag=true;
  return 1;
}
