//this file converts a string of RAW files to a single RAW file.
//useful for converting a IV file to a RAW file.
//requires: common.h and common.cpp files (for the myPoint and triangle data structures.)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string>
#include <afx.h>

using namespace std ;


#include "common.h"

#define MAX_PTS 60000
#define MAX_TRIS 30000
#define MAX_FILES 50

char dirname[2000];
int total_files=0;
char output[2000];

FILE* fptrs[MAX_FILES];
int used=0;

myPoint allpoints[MAX_PTS];
triangle surface[MAX_TRIS];
int ptsUsed, trisUsed;

void open_all()
{
	CFileFind dirPtr;
	string str(dirname);
	int count=0;
	BOOL isWorking = true;

	str.append("\\*.raw");

	isWorking = dirPtr.FindFile(str.c_str(), 0);

	while (isWorking)
	{
		isWorking = dirPtr.FindNextFile();
		printf("%s\n", dirPtr.GetFileName());
		
		fptrs[count] = fopen(dirPtr.GetFilePath(), "r");
		if (fptrs[count] == NULL)
		{
			printf("Error in opening the file %s Exiting\n", dirPtr.GetFilePath());
			exit(0);
		}

		count++;
		if (count >= MAX_FILES)
		{
			printf("Can handle only 50 files at a time. Exiting\n");
			exit(0);
		}
	}
	
	dirPtr.Close();
	total_files = count;

	printf("filecount = %d opened\n", count);
}

void process_all()
{
	int i, j, k, counter;
	int tVerts=0, tTris=0;
	int t1, t2, t3;

	for (counter=0; counter<total_files; counter++)
	{
		tVerts = ptsUsed;
		tTris = trisUsed;

		//First the points and the triangles
		fscanf(fptrs[counter], "%d %d", &i, &j);
		
		//Then read the points
		for (k=0; k<i; k++)
		{
			fscanf(fptrs[counter], "%f %f %f", &allpoints[ptsUsed].x, &allpoints[ptsUsed].y, &allpoints[ptsUsed].z);
			ptsUsed++;
			if (ptsUsed >= MAX_PTS)
			{
				printf("points exceed the max size. pl. change the limit %d\n", ptsUsed);
				exit(0);
			}
		}

		//Then read the triangles
		for (k=0; k<j; k++)
		{
			fscanf(fptrs[counter], "%d %d %d", &t1, &t2, &t3);
			surface[trisUsed].v1 = tVerts + t1;		surface[trisUsed].v2 = tVerts + t2;		surface[trisUsed].v3 = tVerts + t3;
			trisUsed++;
			if (trisUsed>= MAX_TRIS)
			{
				printf("triangles exceed the max size. pl. change the limit %d\n", trisUsed);
				exit(0);
			}
		}

		printf("completed %d files\n", counter+1);
		fclose(fptrs[counter]);
	}

	printf("Finished reading all the data in\n");
}

void chq_me()
{
	int i, j, counter=0;

	for (i=0; i<ptsUsed; i++)
	{
		for (j=0; j<ptsUsed; j++)
		{
			if (i == j) continue;

			if (isEqual(allpoints[i], allpoints[j]))
				counter++;
		}
	}

	printf("In all the files, %d points are same\n", counter);
}

void write_all()
{
	int i;
	FILE* fp;
	
	fp = fopen(output, "w");
	if (fp == NULL)
	{
		printf("could not open %s for writing output\n", output);
		exit(0);
	}

	fprintf(fp, "%d %d\n", ptsUsed, trisUsed);

	for (i=0; i<ptsUsed; i++)
		fprintf(fp, "%f %f %f\n", allpoints[i].x, allpoints[i].y, allpoints[i].z);

	for (i=0; i<trisUsed; i++)
		fprintf(fp, "%d %d %d\n", surface[i].v1, surface[i].v2, surface[i].v3);

	fclose(fp);

	printf("Finished writing the data\n");
}

int main (int argc, char **argv)
{
	if (argc != 4) 
	{
		printf("Enter the directory name, the total number of files and the putput file name\n");
		printf("max %d files at a time !!!\n", MAX_FILES);
		return 1;
	}

	strcpy(dirname, argv[1]);
	total_files= atoi(argv[2]);
	strcpy(output, argv[3]);

	printf("dirname = %s\t and files = %d\t and output = %s\n", dirname, total_files, output);

	ptsUsed = trisUsed =0;

	open_all();
	process_all();
	chq_me();
	write_all();

	return 0;
}