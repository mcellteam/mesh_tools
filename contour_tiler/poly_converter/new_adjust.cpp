#ifdef _MSC_VER
#pragma warning (disable : 4786)
#endif //_MSC_VER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <map>

#include "header.h"

extern tempVert* verts;
extern tempTri* tris;
extern int new_verts, new_tris;
extern int** trisShared;

int* neighbors;
int usedNeighs, prevUsed, total_done;

std::map<int,int> myMap;
std::map<int,int>::const_iterator iter;

int isAligned (int ver1, int ver2)
{
	if (ver1 == 1)
		if (ver2 == 2)	return 1;
		else return 0;

	if (ver1 == 2)
		if (ver2 == 3)	return 1;
		else return 0;

	if (ver1 == 3)
		if (ver2 == 1)	return 1;
		else return 0;

	return -1; //its an error, but, just to make the compiler happy. :-(
}

void exchangeVerts(int tri, int ver1, int ver2)
{
	if (tris[tri].v1 == ver1)
	{
		tris[tri].v1 = ver2;
		if (tris[tri].v2 == ver2)
			tris[tri].v2 = ver1;
		else	tris[tri].v3 = ver1;
	}
	else if (tris[tri].v2 == ver1)
	{
		tris[tri].v2 = ver2;
		if (tris[tri].v1 == ver2)
			tris[tri].v1 = ver1;
		else	tris[tri].v3 = ver1;
	}
	else if (tris[tri].v3 == ver1)
	{
		tris[tri].v3 = ver2;
		if (tris[tri].v1 == ver2)
			tris[tri].v1 = ver1;
		else	tris[tri].v2 = ver1;
	}
}

//also, change the dammed order of the vertices for the triangle.
int triangle_angles(int one, int two, int ver1, int ver2)
{
	int v1, v2, c1, c2;

	v1= v2 = c1 = c2= -1;

	if (tris[one].v1 == ver1)	v1 =1;
	if (tris[one].v1 == ver2)	v2 =1;
	if (tris[one].v2 == ver1)	v1 =2;
	if (tris[one].v2 == ver2)	v2 =2;
	if (tris[one].v3 == ver1)	v1 =3;
	if (tris[one].v3 == ver2)	v2 =3;

	if (tris[two].v1 == ver1)	c1 =1;
	if (tris[two].v1 == ver2)	c2 =1;
	if (tris[two].v2 == ver1)	c1 =2;
	if (tris[two].v2 == ver2)	c2 =2;
	if (tris[two].v3 == ver1)	c1 =3;
	if (tris[two].v3 == ver2)	c2 =3;

	if ( (v1 == -1) || (v2 == -1) || (c1 == -1) || (c2 == -1) )
	{
		printf("some err in <triangle_angles> : %d %d %d %d\n", one, two, ver1, ver2);
		return 1; //wot 2 return... :-(
	}
	
	if (isAligned(v1, v2))
	{
		if (isAligned(c1, c2))
		{
			//problemo.
			exchangeVerts(two, ver1, ver2);
			return 0;
		}
		else
		{
			//no problemo
			return 1;
		}
	}
	else
	{
		if (isAligned(c1, c2))
		{
			//no problemo.
			return 1;
		}
		else
		{
			//problemo
			exchangeVerts(two, ver1, ver2);
			return 0;
		}
	}
}


void insert_tri(int tri)
{
	if (tris[tri].type == -1) return;

	iter = myMap.find(tri);
	if(iter == myMap.end()) //ie not found
	{
		myMap[tri] = tri;
		neighbors[usedNeighs++] = tri;
		total_done++;
	}
}

void align_us(int with, int what, int vert)
{
	int i, j, flag=-1;
	int v1[3], v2[3];
	
	if (tris[what].type != -1) return;

	v1[0] = tris[with].v1;	v1[1] = tris[with].v2;	v1[2] = tris[with].v3;
	v2[0] = tris[what].v1;	v2[1] = tris[what].v2;	v2[2] = tris[what].v3;
	
	for (i=0; i<3; i++)
	{
		if (v1[i] == vert) continue;

		for (j=0; j<3; j++)
		{
			if (v2[j] == vert) continue;

			if (v1[i] == v2[j])
				flag = v1[i];
		}
	}

	if (flag == -1)
		return;

	//then compare the two triangles.
	if (triangle_angles(with, what, vert, flag))
		tris[what].type = tris[with].type;
	else
		tris[what].type = !(tris[with].type);

	//Then insert this triangle into the NEIGHBORS array.
	insert_tri(what);
}

void orient_vert(int tri, int vert)
{
	int i;

	for (i=1; i<trisShared[vert][0]; i++)
	{
		if (tri != trisShared[vert][i])
			align_us(tri, trisShared[vert][i], vert);
	}
}

void correct_tri(int tri)
{
	orient_vert(tri, tris[tri].v1);
	orient_vert(tri, tris[tri].v2);
	orient_vert(tri, tris[tri].v3);			
}

//This one is called for each unconnected component. Assume that the first triangle is pointing outwards...
void getNextComponent()
{
	int i;

	for (i=0; i<new_tris; i++)
	{
		if (tris[i].type == -1)
			break;
	}
	
	tris[i].type =1;
	insert_tri(i);
	prevUsed =usedNeighs;
}

void start_fireworks()
{
	int i, j, lastone;
	int* tarray;

	neighbors = (int*) (malloc(sizeof(int) * new_tris));
	tarray = (int*) (malloc(sizeof(int) * new_tris));

	printf("\n<start_fireworks> started...\n");

	myMap.clear();
	lastone = usedNeighs= total_done =0;

	while (1)
	{
		prevUsed = usedNeighs;
		printf("still processing with %d Triangles\n", prevUsed);
		
		if (lastone == prevUsed)
			getNextComponent();
		else
			lastone = prevUsed;

		for (i=0; i<prevUsed; i++)
			correct_tri(neighbors[i]);

		if (total_done == new_tris)
		{
			printf("The reqd normal flipping is done.\n");
			break;
		}

		j=0;
		for (iter=myMap.begin(); iter!=myMap.end(); ++iter)
		{
			neighbors[j++] = (*iter).first;
		}
		usedNeighs =j;
	}

	free(neighbors);
	free(tarray);
	myMap.clear();

	printf("<start_fireworks> over...\n");
}

