#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "tetgen2FV.h"
#include "c_3d.h"

	int circumcent(struct tet_mesh *tet_mesh)
	{
        struct node **nodes;
        struct tet_list *tet_head;
	struct tet_list *itet;
	double xo,yo,zo;
	double x,y,z;
        double r;
	double bx[4][4],by[4][4],bz[4][4];
	double c[4][4];
	double st[4][4];
	double tetmat[3][4];
        int i;

        tet_head=tet_mesh->tet_head;
        nodes=tet_mesh->nodes;

	for (itet=tet_head; itet!=NULL; itet=itet->next)
	{

          for (i=0;i<4;i++)
          {
	    tetmat[0][i]=nodes[itet->tet->node_index[i]]->x;
	    tetmat[1][i]=nodes[itet->tet->node_index[i]]->y;
	    tetmat[2][i]=nodes[itet->tet->node_index[i]]->z;
          }

	  st[0][0]=tetmat[0][0];
          st[0][1]=tetmat[1][0];
          st[0][2]=tetmat[2][0];
          st[0][3]=1.0;
          st[1][0]=tetmat[0][1];
          st[1][1]=tetmat[1][1];
          st[1][2]=tetmat[2][1];
          st[1][3]=1.0;
          st[2][0]=tetmat[0][2];
          st[2][1]=tetmat[1][2];
          st[2][2]=tetmat[2][2];
          st[2][3]=1.0;
          st[3][0]=tetmat[0][3];
          st[3][1]=tetmat[1][3];
          st[3][2]=tetmat[2][3];
          st[3][3]=1.0;

	  bx[0][0]=tetmat[0][0]*tetmat[0][0]+tetmat[1][0]*tetmat[1][0]+tetmat[2][0]*tetmat[2][0];
	  bx[0][1]=tetmat[1][0];
	  bx[0][2]=tetmat[2][0];
	  bx[0][3]=1.0;
	  bx[1][0]=tetmat[0][1]*tetmat[0][1]+tetmat[1][1]*tetmat[1][1]+tetmat[2][1]*tetmat[2][1];
	  bx[1][1]=tetmat[1][1];
	  bx[1][2]=tetmat[2][1];
	  bx[1][3]=1.0;
	  bx[2][0]=tetmat[0][2]*tetmat[0][2]+tetmat[1][2]*tetmat[1][2]+tetmat[2][2]*tetmat[2][2];
	  bx[2][1]=tetmat[1][2];
	  bx[2][2]=tetmat[2][2];
	  bx[2][3]=1.0;
	  bx[3][0]=tetmat[0][3]*tetmat[0][3]+tetmat[1][3]*tetmat[1][3]+tetmat[2][3]*tetmat[2][3];
	  bx[3][1]=tetmat[1][3];
	  bx[3][2]=tetmat[2][3];
	  bx[3][3]=1.0;

	  by[0][0]=tetmat[0][0]*tetmat[0][0]+tetmat[1][0]*tetmat[1][0]+tetmat[2][0]*tetmat[2][0];
          by[0][1]=tetmat[0][0];
          by[0][2]=tetmat[2][0];
          by[0][3]=1.0;
          by[1][0]=tetmat[0][1]*tetmat[0][1]+tetmat[1][1]*tetmat[1][1]+tetmat[2][1]*tetmat[2][1];
          by[1][1]=tetmat[0][1];
          by[1][2]=tetmat[2][1];
          by[1][3]=1.0;
          by[2][0]=tetmat[0][2]*tetmat[0][2]+tetmat[1][2]*tetmat[1][2]+tetmat[2][2]*tetmat[2][2];
          by[2][1]=tetmat[0][2];
          by[2][2]=tetmat[2][2];
          by[2][3]=1.0;
          by[3][0]=tetmat[0][3]*tetmat[0][3]+tetmat[1][3]*tetmat[1][3]+tetmat[2][3]*tetmat[2][3];
          by[3][1]=tetmat[0][3];
          by[3][2]=tetmat[2][3];
          by[3][3]=1.0;

	  bz[0][0]=tetmat[0][0]*tetmat[0][0]+tetmat[1][0]*tetmat[1][0]+tetmat[2][0]*tetmat[2][0];
          bz[0][1]=tetmat[0][0];
          bz[0][2]=tetmat[1][0];
          bz[0][3]=1.0;
          bz[1][0]=tetmat[0][1]*tetmat[0][1]+tetmat[1][1]*tetmat[1][1]+tetmat[2][1]*tetmat[2][1];
          bz[1][1]=tetmat[0][1];
          bz[1][2]=tetmat[1][1];
          bz[1][3]=1.0;
          bz[2][0]=tetmat[0][2]*tetmat[0][2]+tetmat[1][2]*tetmat[1][2]+tetmat[2][2]*tetmat[2][2];
          bz[2][1]=tetmat[0][2];
          bz[2][2]=tetmat[1][2];
          bz[2][3]=1.0;
          bz[3][0]=tetmat[0][3]*tetmat[0][3]+tetmat[1][3]*tetmat[1][3]+tetmat[2][3]*tetmat[2][3];
          bz[3][1]=tetmat[0][3];
          bz[3][2]=tetmat[1][3];
          bz[3][3]=1.0;


	  c[0][0]=tetmat[0][0]*tetmat[0][0]+tetmat[1][0]*tetmat[1][0]+tetmat[2][0]*tetmat[2][0];
          c[0][1]=tetmat[0][0];
          c[0][2]=tetmat[1][0];
          c[0][3]=tetmat[2][0];
          c[1][0]=tetmat[0][1]*tetmat[0][1]+tetmat[1][1]*tetmat[1][1]+tetmat[2][1]*tetmat[2][1];
          c[1][1]=tetmat[0][1];
          c[1][2]=tetmat[1][1];
          c[1][3]=tetmat[2][1];
          c[2][0]=tetmat[0][2]*tetmat[0][2]+tetmat[1][2]*tetmat[1][2]+tetmat[2][2]*tetmat[2][2];
          c[2][1]=tetmat[0][2];
          c[2][2]=tetmat[1][2];
          c[2][3]=tetmat[2][2];
          c[3][0]=tetmat[0][3]*tetmat[0][3]+tetmat[1][3]*tetmat[1][3]+tetmat[2][3]*tetmat[2][3];
          c[3][1]=tetmat[0][3];
          c[3][2]=tetmat[1][3];
          c[3][3]=tetmat[2][3];

	  xo=dmat_det_4d(bx)/(2.*dmat_det_4d(st));
	  yo=-dmat_det_4d(by)/(2.*dmat_det_4d(st));
	  zo=dmat_det_4d(bz)/(2.*dmat_det_4d(st));

	  x=xo-tetmat[0][0];
	  y=yo-tetmat[1][0];
	  z=zo-tetmat[2][0];
	  r=sqrt(x*x + y*y + z*z);

/*
	  printf ("%g %g %g   %g\n",xo,yo,zo,r);
*/

          itet->tet->cent.x=xo;
          itet->tet->cent.y=yo;
          itet->tet->cent.z=zo;
          itet->tet->r=r;
        }

	return(0);

	}



      double dmat_det_4d(double a[4][4])
      {
	double dmat;

      dmat = 
      a[0][0]*( 
        a[1][1] * ( a[2][2] * a[3][3] - a[2][3] * a[3][2] )
      - a[1][2] * ( a[2][1] * a[3][3] - a[2][3] * a[3][1] ) 
      + a[1][3] * ( a[2][1] * a[3][2] - a[2][2] * a[3][1] ) ) 
    - a[0][1] * ( 
        a[1][0] * ( a[2][2] * a[3][3] - a[2][3] * a[3][2] ) 
      - a[1][2] * ( a[2][0] * a[3][3] - a[2][3] * a[3][0] ) 
      + a[1][3] * ( a[2][0] * a[3][2] - a[2][2] * a[3][0] ) )
    + a[0][2] * (
        a[1][0] * ( a[2][1] * a[3][3] - a[2][3] * a[3][1] ) 
      - a[1][1] * ( a[2][0] * a[3][3] - a[2][3] * a[3][0] )
      + a[1][3] * ( a[2][0] * a[3][1] - a[2][1] * a[3][0] ) ) 
    - a[0][3] * ( 
        a[1][0] * ( a[2][1] * a[3][2] - a[2][2] * a[3][1] ) 
      - a[1][1] * ( a[2][0] * a[3][2] - a[2][2] * a[3][0] ) 
      + a[1][2] * ( a[2][0] * a[3][1] - a[2][1] * a[3][0] ) );

       return (dmat);
      }

