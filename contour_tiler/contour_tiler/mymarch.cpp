#include <stdio.h>
#include <math.h>
#include <memory.h>

#include "cubes.h"
#include "common.h"
#include "myutil.h"
#include "mymarch.h"


/* put the following variables as globle to save calling time */
static double Thres; 
static int Dimx,Dimy,Dimz,Detx,Dety,Detz;
static float *Xpos,*Ypos,*Zpos;

extern int STATIC_CHQ;

void my_clear_mymarch()
{
	Thres =0.0;
	Dimx=Dimy=Dimz=Detx=Dety=Detz = 0;
	Xpos = NULL; Ypos = NULL; Zpos = NULL; 

	calc_triangles(NULL, NULL, STATIC_CHQ, STATIC_CHQ, STATIC_CHQ, NULL);
	calc_lines(NULL, NULL, NULL);
}

void init_marching_cubes(double thres, float *xpos, float *ypos, float *zpos,
						 int dimx, int dimy, int dimz, int detx, int dety, int detz)
{
    Thres=thres;
    Xpos=xpos;
    Ypos=ypos;
    Zpos=zpos;
    Detx=detx;
    Dety=dety;
    Detz=detz;
    Dimx=dimx;
    Dimy=dimy;
    Dimz=dimz;
}

int count_triangles(float *ary0, float *ary1, int x,int y )
/* given a volume, count its triangles based on the thresholding value 
return : the number of triangles
*/ 
{ 
#define DIMENSION_NUM 3
	double ptval[8];
	int i,temp;
	unsigned int code;
	
	temp=y*Dimx+x;
	ptval[0]=ary0[temp];
	ptval[1]=ary0[temp+Detx];
	ptval[2]=ary1[temp+Detx];
	ptval[3]=ary1[temp];
	temp=(y+Dety)*Dimx+x;
	ptval[4]=ary0[temp];
	ptval[5]=ary0[temp+Detx];
	ptval[6]=ary1[temp+Detx];
	ptval[7]=ary1[temp];
	
	code=0;

	for (i=7; i>=0; i--)
	{
		code <<=1;
		if (ptval[i] < Thres) code |= 1;
	}
	
	return(cubes[code][0]);
}


int calc_triangles(float *ary0, float *ary1, int x,int y,int z, 
				   double triangles[12][SPACE_VERT_DIM]) 
				   /* given a volume, find its triangles based on the thresholding value 
				   triangles[12][6]: store the vertices, [0-2] for first poly, [3-5] for 2nd ..
				   return : the number of triangles
				   */ 
{ 
#define DIMENSION_NUM 3
	/* the vertice of edge */
	static int Edge_start[12]={0,1,2,3,4,5,6,7,0,1,3,2}; 
	static int Edge_end[12]=  {1,2,3,0,5,6,7,4,4,5,7,6};
	/* Some varables use "double" for the compatibility of Dan's library */
	double ptval[8],pts[8][DIMENSION_NUM]; 
	double edgevert[12][DIMENSION_NUM]; /* include normals */
	int x1,y1,z1,i,j,temp;
	unsigned int code;
	
	temp=y*Dimx+x;
	ptval[0]=ary0[temp];
	ptval[1]=ary0[temp+Detx];
	ptval[2]=ary1[temp+Detx];
	ptval[3]=ary1[temp];
	temp=(y+Dety)*Dimx+x;
	ptval[4]=ary0[temp];
	ptval[5]=ary0[temp+Detx];
	ptval[6]=ary1[temp+Detx];
	ptval[7]=ary1[temp];
	
	code=0;
	
	for (i=7; i>=0; i--)
	{
		code <<=1;
		if (ptval[i] < Thres) code |= 1;
	}
	
	if(code==0 || code==255)
	{
		return 0; /* no triangles generated */
	}
	/*
	printf("val= %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",ptval[0],ptval[1],
	ptval[2],ptval[3],ptval[4],ptval[5],ptval[6],ptval[7]);
	*/
	x1=x+Detx; y1=y+Dety; z1=z+Detz;
	pts[0][0]=Xpos[x];  pts[0][1]=Ypos[y]; pts[0][2]=Zpos[z]; 
	pts[1][0]=Xpos[x1]; pts[1][1]=Ypos[y]; pts[1][2]=Zpos[z]; 
	pts[2][0]=Xpos[x1]; pts[2][1]=Ypos[y]; pts[2][2]=Zpos[z1]; 
	pts[3][0]=Xpos[x];  pts[3][1]=Ypos[y]; pts[3][2]=Zpos[z1]; 
	pts[4][0]=Xpos[x];  pts[4][1]=Ypos[y1]; pts[4][2]=Zpos[z]; 
	pts[5][0]=Xpos[x1]; pts[5][1]=Ypos[y1]; pts[5][2]=Zpos[z]; 
	pts[6][0]=Xpos[x1]; pts[6][1]=Ypos[y1]; pts[6][2]=Zpos[z1]; 
	pts[7][0]=Xpos[x];  pts[7][1]=Ypos[y1]; pts[7][2]=Zpos[z1]; 
	/* normals 
	for(i=0; i<8; i++){
    double dist=0.0,ratio;
    dist=pts[i][3]*pts[i][3] + pts[i][4]*pts[i][4] + pts[i][5]*pts[i][5];
    dist=sqrt(dist);
    ratio=1.0/dist;
    pts[i][3]*=ratio;
    pts[i][4]*=ratio;
    pts[i][5]*=ratio;
	}
	*/
	
	/* edgevert[i], 0<=i<8, is between the node i & i+1 */ 
	/* This picture is drawed by Dan Schikore. 
	____________________             ____________________
	7/|                 /|6           /|        6        /|
	/  |               /  |          /  |               /  |
	/    |             /    |        /10  |             /11  |
	3/______|___________/2     |      /______|_____2_____/      |
	|       |          |       |     |       |          |       5
	|       |          |       |     |       7          |       |
	|       |          |       |     |       |          |       |
	|       |          |       |     |       |          |       |
	|       |__________|_______|     3       |_____4____1_______|
	|     4/           |      /5     |      /           |      /
	|    /             |    /        |    /8            |   9/
	|  /               |  /          |  /               |  /
    0|/_________________|/1           |/_________0_______|/
	Vertices                           Edges
	*/
	
	for(i=0; i<12 ; i++)
	{
		double ratio, ratio1;
		int start_v,end_v;
		start_v=Edge_start[i];
		end_v=Edge_end[i];
		ratio= (Thres-ptval[start_v])/(ptval[end_v]-ptval[start_v]);
	
		if(ratio >0.0 && ratio <1.0)
		{
			ratio1=1.0 - ratio;
			for(j=0; j<SPACE_VERT_DIM; j++) 
			{
				if(pts[start_v][j] ==  pts[end_v][j]) /* avoid rounding error */
					edgevert[i][j] = pts[start_v][j];
				else 
					edgevert[i][j] = pts[start_v][j]*ratio1 + pts[end_v][j]*ratio;
			}
		} 
	}

	j=cubes[code][0]*3;
	
	for (i=0; i<j; i++) 
	{
		memcpy(triangles[i], edgevert[cubes[code][i+1]],
			sizeof(triangles[0][0])*SPACE_VERT_DIM);
	}

	return cubes[code][0];
}


int calc_lines(double *box, double *vals, double lines[4][SPACE_VERT_DIM])
{
	static int cubes[16][5] = 
	{
        {0}, {1,0,3}, {1,0,1}, {1,1,3},
        {1,1,2},{2,0,1,2,3},{1,0,2},{1,2,3},
        {1,2,3},{1,0,2},{2,0,1,2,3},{1,1,2},
        {1,1,3}, {1,0,1}, {1,0,3}, {0}
	};
		
		static int Edge_start[4]={0,1,2,3};
		static int Edge_end[4]={1,2,3,0};
		double pts[4][4];
		double edgevert[4][PLANE_VERT_DIM]; /* include normals */
		int i,j;
		int line_num;
		unsigned int code;
		
		code=0;
	
		for (i=3; i>=0; i--)
		{
			code <<=1;
			if (vals[i] < Thres) code |= 1;
		}
		
		if(code==0 || code==15)
		{
			return 0; /* no line segment */
		}
		
		/*
		if(code==5 || code==0xa)
        printf("x0=%d y0=%d, 2 lines\n",x0,y0);
		*/
		
		pts[0][0]=box[0];  pts[0][1]=box[2];
		pts[1][0]=box[1];  pts[1][1]=box[2];
		pts[2][0]=box[1];  pts[2][1]=box[3];
		pts[3][0]=box[0];  pts[3][1]=box[3];
		
		/* gradient information */
		/*
		pts[0][2]=(vals[1]-vals[0]);
		pts[0][3]=(vals[3]-vals[0]);
		pts[1][2]=pts[0][2];
		pts[1][3]=(vals[2]-vals[1]);
		pts[2][2]=(vals[2]-vals[3]);
		pts[2][3]=pts[1][3];
		pts[3][2]=pts[2][2];
		pts[3][3]=pts[0][3];
		*/
		
		for(i=0; i<4; i++)
		{
			double ratio, ratio1;
			double v1,v2;
			int start_v,end_v;
			start_v=Edge_start[i];
			end_v=Edge_end[i];
			
			v1=vals[start_v];
			v2=vals[end_v];

			if((v1>=Thres && v2<Thres) || (v1<Thres && v2>=Thres))
			{
				ratio= (Thres-v1)/(v2-v1);
				ratio1=1.0 - ratio;
				
				for(j=0; j<4; j++) 
				{
					edgevert[i][j] = pts[start_v][j]*ratio1 + pts[end_v][j]*ratio;
				}
			}
		}

		line_num=cubes[code][0];
		j=line_num*2;
		
		for (i=0; i<j; i++) 
		{
			int temp;
			temp=cubes[code][i+1];
			lines[i][0]=edgevert[temp][0];
			lines[i][1]=edgevert[temp][1];
			/*    normals 
			val=edgevert[temp][2]*edgevert[temp][2] +
            edgevert[temp][3]*edgevert[temp][3];
			val=sqrt(val);
			lines[i][2]=edgevert[temp][2]*val;
			lines[i][3]=edgevert[temp][3]*val;
			*/
		}
		return line_num;
}



