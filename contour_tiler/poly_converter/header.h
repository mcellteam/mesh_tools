#ifndef _HEADER_H_
#define _HEADER_H_

typedef struct _V_
{
	float x;
	float y;
	float z;

	int isUsed;	//reused later on as the new index.
	int count;	//the total # of tris that share this vert.
	
} tempVert;

typedef struct _T_
{
	int v1;
	int v2;
	int v3;

	int type;

} tempTri;


static bool isSame(float one, float two)
{
	float TOLERANCE = 0.0001f;
	one -=two;

	if (one<0) one*= -1;

	if ((-1*TOLERANCE <= one) && (one <= TOLERANCE))
		return 1;
	return 0;
}

#endif

