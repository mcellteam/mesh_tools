#ifdef _MSC_VER
#pragma warning (disable : 4786)
#endif //_MSC_VER

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <map>

#include "header.h"

extern void start_fireworks();
extern char *ifname, *ofname;

#define TRIS_PER_VERT 50

int** trisShared;

std::multimap<int,int>myVertMap;
std::multimap<int,int>::const_iterator myIter;
typedef std::multimap<int,int>::value_type value_type;

extern tempVert* verts;
extern tempTri* tris;
extern int orig_verts, no_of_tris, new_verts, new_tris;

int isClose(float one)
{
	double TOLERNACE =0.0001;	//all verts within this tolerance are merged...
	
	if ( (-1*TOLERNACE<=one) && (one<=TOLERNACE) )
		return 1;
	return 0;
}

void init_the_maps()
{
	myVertMap.clear();
}

int insert_vertex(float vert[3])
{
	verts[orig_verts].x = vert[0];	verts[orig_verts].y = vert[1];	verts[orig_verts].z = vert[2];
	verts[orig_verts].count = verts[orig_verts].isUsed =0;
	orig_verts++;
	return (orig_verts-1);
}

int chq_vertex_in_map(float vert[3])
{
	int leng, potential;

	leng = (int)(1000*(sqrt(vert[0]*vert[0] + vert[1]*vert[1] + vert[2]*vert[2])));

	myIter = myVertMap.find(leng);
	if(myIter == myVertMap.end()) //ie not found
	{
		myVertMap.insert(value_type(leng, orig_verts));
		return (insert_vertex(vert));
	}
	else
	{
		for (myIter = myVertMap.lower_bound(leng); myIter != myVertMap.upper_bound(leng); ++myIter) 
		{
            potential = myIter->second;
			float dist = (float)( (verts[potential].x-vert[0])*(verts[potential].x-vert[0]) + 
						(verts[potential].y-vert[1])*(verts[potential].y-vert[1]) + 
						(verts[potential].z-vert[2])*(verts[potential].z-vert[2]));
		
			if (isClose(dist))
			{
				return potential;
			}
			else
			{
				myVertMap.insert(value_type(leng, orig_verts));
				return (insert_vertex(vert));
			}
		}
	}

	myVertMap.insert(value_type(leng, orig_verts));
	orig_verts++;
	return (orig_verts-1);
}

void reverse_ptrs()
{
	int i, v[3], c[3];

	trisShared = (int**)(malloc(sizeof(int) * new_verts));
	for (i=0; i<new_verts; i++)
	{
		trisShared[i] = (int*)(malloc(sizeof(int) * TRIS_PER_VERT));
		trisShared[i][0] = 0;
	}

	for (i=0; i<new_tris; i++)
	{
		v[0] = tris[i].v1;			v[1] = tris[i].v2;			v[2] = tris[i].v3;
		c[0] = trisShared[ v[0] ][0];	c[1] = trisShared[ v[1] ][0];	c[2] = trisShared[ v[2] ][0];

		trisShared[ v[0] ][ c[0] ] = i;		trisShared[ v[1] ][ c[1] ] = i;		trisShared[ v[2] ][ c[2] ] = i;
		
		trisShared[ v[0] ][ c[0] ] = c[0]+1;	trisShared[ v[1] ][ c[1] ] = c[1]+1;	trisShared[ v[2] ][ c[2] ] = c[2]+1;

		if (c[0] >= TRIS_PER_VERT) 
			printf("more than %d triangles share this vertex... %d for vert=%d\n", TRIS_PER_VERT, v[0], c[0]);
		if (c[1] >= TRIS_PER_VERT) 
			printf("more than %d triangles share this vertex... %d for vert=%d\n", TRIS_PER_VERT, v[1], c[1]);
		if (c[2] >= TRIS_PER_VERT) 
			printf("more than %d triangles share this vertex... %d for vert=%d\n", TRIS_PER_VERT, v[2], c[2]);
	}
}

void repeat_vertex(int vert)
{
	int leng, potential;

	leng = (int)(1000*(sqrt(verts[vert].x*verts[vert].x + verts[vert].y*verts[vert].y +verts[vert].z*verts[vert].z)));

	myIter = myVertMap.find(leng);
	if(myIter == myVertMap.end()) //ie not found
	{
		myVertMap.insert(value_type(leng, vert));
	}
	else
	{
		for (myIter = myVertMap.lower_bound(leng); myIter != myVertMap.upper_bound(leng); ++myIter) 
		{
            potential = myIter->second;
			float dist = (float)( (verts[potential].x-verts[vert].x)*(verts[potential].x-verts[vert].x) + 
						(verts[potential].y-verts[vert].y)*(verts[potential].y-verts[vert].y) + 
						(verts[potential].z-verts[vert].z)*(verts[potential].z-verts[vert].z));
		
			if (isClose(dist))
			{
				verts[vert].isUsed = potential;
				return;
			}
		}
		myVertMap.insert(value_type(leng, vert));
	}
}

void read_raw()
{
	FILE* fp;
	int i;
	
	fp = fopen(ifname, "r");
	if (fp == NULL)
	{
		printf("could not open %s\n", ifname);
		exit(0);
	}

	fscanf(fp, "%d %d", &orig_verts, &no_of_tris);
	myVertMap.clear();

	verts = (tempVert*)(malloc(sizeof(tempVert) * orig_verts));
	tris = (tempTri*)(malloc(sizeof(tempTri) * no_of_tris));
	new_verts=new_tris=0;

	for (i=0; i<orig_verts; i++)
	{
		fscanf(fp, "%f %f %f", &verts[i].x, &verts[i].y, &verts[i].z);
		verts[i].count =0;			verts[i].isUsed=0;
		repeat_vertex(i);
		if (i%1000 ==0)
			printf("%d read \n", i);
	}

	for (i=0; i<no_of_tris; i++)
	{
		fscanf(fp, "%d %d %d", &tris[i].v1, &tris[i].v2, &tris[i].v3);
		tris[i].type =-1;
	}

	fclose(fp);

	printf("read_file done for %d tris and %d verts\n", no_of_tris, orig_verts);
}

void process_crap()
{
	int i, j;
	tempVert* tverts;
	tempTri* ttris;

	//1. First, remove the close (or repeated) vertices from the data.
	for (i=0; i<no_of_tris; i++)
	{
		if (verts[ tris[i].v1 ].isUsed != 0)
			tris[i].v1 = verts[ tris[i].v1 ].isUsed;

		if (verts[ tris[i].v2 ].isUsed != 0)
			tris[i].v2 = verts[ tris[i].v2 ].isUsed;

		if (verts[ tris[i].v3 ].isUsed != 0)
			tris[i].v3 = verts[ tris[i].v3 ].isUsed;
	}

	for (i=0; i<no_of_tris; i++)
	{
		verts[tris[i].v1].count++;	verts[tris[i].v2].count++;	verts[tris[i].v3].count++;
	}

	for (i=0; i<orig_verts; i++)
	{
		if (verts[i].count !=0)
			verts[i].isUsed = new_verts++;
	}
	
	for (i=0; i<no_of_tris; i++)
	{
		tris[i].v1 = verts[tris[i].v1].isUsed;
		tris[i].v2 = verts[tris[i].v2].isUsed;
		tris[i].v3 = verts[tris[i].v3].isUsed;
	}

	//Then, realloc the memory for the verts
	tverts = (tempVert*) (malloc(sizeof(tempVert) * new_verts));
	for (i=0,j=0; i<orig_verts; i++)
	{
		if (verts[i].count !=0)
		{
			tverts[j].x = verts[i].x;	tverts[j].y = verts[i].y;	tverts[j].z = verts[i].z;
			tverts[j].count = verts[i].count;		j++;
		}
	}
	
	if (j != new_verts)
		printf("some major problem in re-verting %d and %d\n", j, new_verts);

	free(verts);
	verts = tverts;

	printf("reduced verts from %d to %d\n", orig_verts, new_verts);

	//2. Then, remove the invalid triangles from the data.
	for (i=0; i<no_of_tris; i++)
	{
		if ( (tris[i].v1 != tris[i].v2) && (tris[i].v1 != tris[i].v3) && 
			 (tris[i].v2 != tris[i].v3) )
		{
			tris[i].type = 0;
			new_tris++;
		}
	}

	//Then, realloc the memory for the triangles
	ttris = (tempTri*) (malloc(sizeof(tempTri) * new_tris));
	for (i=0,j=0; i<no_of_tris; i++)
	{
		if (tris[i].type ==0)
		{
			ttris[j].v1 = tris[i].v1;	ttris[j].v2 = tris[i].v2;	ttris[j].v3 = tris[i].v3;
			ttris[j].type = -1;				j++;
		}
	}
	
	if (j != new_tris)
		printf("some major problem in re-trianging %d and %d\n", j, new_tris);

	free(tris);
	tris = ttris;
	
	printf("reduced tris from %d to %d\n", no_of_tris, new_tris);

	//3. Then, Orient the normals uniformly...
	//reverse_ptrs();
	//start_fireworks();

	printf("Finished orienting the normals\n");
}

void write_raw_file()
{
	FILE* fp;
	int i;

	fp = fopen(ofname, "w");
	if (fp == NULL)
	{
		printf("could not open %s\n", ofname);
		exit(0);
	}

	fprintf(fp, "%d %d\n", new_verts, new_tris);

	for (i=0; i<new_verts; i++)
		fprintf(fp, "%f %f %f\n", verts[i].x, verts[i].y, verts[i].z);

	for (i=0; i<new_tris; i++)
		fprintf(fp, "%d %d %d\n", tris[i].v1, tris[i].v2, tris[i].v3);

	fclose(fp);
	printf("write_file done for %d tris and %d verts\n", no_of_tris, new_verts);
}

void write_stl_file()
{
	FILE* fp;
	int i;
	float n=0.0f;

	fp = fopen(ofname, "w");
	if (fp == NULL)
	{
		printf("could not open %s\n", ofname);
		exit(0);
	}

	fprintf(fp, "solid crap\n");
	
	for (i=0; i<no_of_tris; i++)
	{
		fprintf(fp, "facet normal %f %f %f\n", n, n, n);
		fprintf(fp, "outer loop\n");
		fprintf(fp, "vertex %f %f %f\n", verts[tris[i].v1].x, verts[tris[i].v1].y, verts[tris[i].v1].z);
		fprintf(fp, "vertex %f %f %f\n", verts[tris[i].v2].x, verts[tris[i].v2].y, verts[tris[i].v2].z);
		fprintf(fp, "vertex %f %f %f\n", verts[tris[i].v3].x, verts[tris[i].v3].y, verts[tris[i].v3].z);
		fprintf(fp, "endloop\n");
		fprintf(fp, "endfacet\n");
	}

	fclose(fp);
	printf("write_file done.\n");
}


void parse_raw(int flag)
{
	if ( (flag ==0) || (flag ==1) )
		read_raw();	
	
	if (flag>0)
	{
		process_crap();
		write_raw_file();
	}
	else
	{
		write_stl_file();
	}

	printf("all done, chq up %s\n", ofname);
}
