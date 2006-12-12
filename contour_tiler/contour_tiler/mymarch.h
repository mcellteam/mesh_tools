#ifndef MYMARCH_H
#define MYMARCH_H

void init_marching_cubes(double thres, float *xpos, float *ypos, float *zpos,
						 int dimx, int dimy, int dimz, int detx, int dety, int detz);

int count_triangles(float *ary0, float *ary1, int x,int y );

int calc_triangles(float *ary0, float *ary1, int x,int y,int z,
				   double triangles[12][SPACE_VERT_DIM]);

int calc_lines(double *box, double *vals, double lines[4][SPACE_VERT_DIM]);

#endif
