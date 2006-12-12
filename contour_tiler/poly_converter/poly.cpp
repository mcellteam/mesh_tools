#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "header.h"

char *ifname, *ofname;
int type;
int maxInd;
tempTri* tris;
tempVert* verts;
int orig_verts, no_of_tris, new_verts, new_tris;
FILE* fp;

extern void ReadMeshFromWrl(FILE *fp);
extern void parse_raw(int flag);

extern void iv2raw();
extern void wrl2raw();
extern void vrml2raw();
extern int chq_vertex_in_map(float vert[3]);
extern void init_the_maps();
extern int no_of_tris, orig_verts;
/*
* From RAW to STL
*/

void raw_2_stl()
{
	int i;

	fp = fopen(ifname, "r");
	if (fp == NULL)
	{
		printf("error in opening the file. Exiting\n");
		exit(0);
	}

	if (fscanf(fp,"%d %d", &orig_verts, &no_of_tris) == EOF)
	{
		printf("Input file is not valid....Exiting...\n");
		exit(0);
	}

	tris = (tempTri*)(malloc(sizeof(tempTri) * no_of_tris));
	verts = (tempVert*)(malloc(sizeof(tempVert) *orig_verts));
	
	/*Now, read the dammed file and populate the DS's*/
	for (i=0; i<orig_verts; i++)
	{
		if (fscanf(fp,"%f %f %f", &verts[i].x, &verts[i].y, &verts[i].z) == EOF)		
		{
			printf("Input file has to have %d Vertices....Exiting...\n",orig_verts);
			exit(0);
		}
		
		if (!(i% 5000))
			printf("still working on vertices !!!! %d \n",i);
	}
	
	printf("Finished reading the Vertices.. Now reading the tempTris\n");
	
	for (i=0; i<no_of_tris; i++)
	{
		if (fscanf(fp,"%d %d %d", &tris[i].v1, &tris[i].v2, &tris[i].v3 ) == EOF)
		{
			printf("Input file has to have %d tempTris....Exiting...\n",no_of_tris);
			exit(0);
		}

		if (!(i% 5000))
			printf("still working on tempTris !!!! %d \n",i);
	}
	
	fclose(fp);
	printf("Finished reading the %d verts and %d tempTris\n", orig_verts, no_of_tris);
	
	fp = fopen(ofname, "w");
	if (fp == NULL)
	{
		printf("error in opening the output file. Exiting\n");
		exit(0);
	}

	fprintf(fp, "solid something\n");
	
	for (i=0; i<no_of_tris; i++)
	{
		fprintf(fp, "facet normal 0.00 0.00 0.00\n");
			fprintf(fp, "outer loop\n");
				fprintf(fp, "vertex %f %f %f\n", verts[ tris[i].v1 ].x, verts[ tris[i].v1 ].y, verts[ tris[i].v1 ].z);
				fprintf(fp, "vertex %f %f %f\n", verts[ tris[i].v2 ].x, verts[ tris[i].v2 ].y, verts[ tris[i].v2 ].z);
				fprintf(fp, "vertex %f %f %f\n", verts[ tris[i].v3 ].x, verts[ tris[i].v3 ].y, verts[ tris[i].v3 ].z);
			fprintf(fp, "endloop\n");
		fprintf(fp, "endfacet\n");
	}
	
	fprintf(fp, "endsolid\n");
	fclose(fp);

	printf("Finished writing the STL file. \n");
}


/*
* From VRML to RAW
*/
void parse_vrml()
{
	fp = fopen(ifname, "r");
	if (fp == NULL)
	{
		printf("error in opening the file. Exiting\n");
		exit(0);
	}

	ReadMeshFromWrl(fp);
}


/*
* From POLY to RAW
*/

float parse_Co_ord(char buff[1000], int str, int end)
{
	char *tbuff = (char*)(malloc(sizeof(char) * (end-str)));
	int i;
	float val;
	
	for (i=str; i<end; i++)
		tbuff[i-str] = buff[i];
	
	val = (float)(atof(tbuff));
	free(tbuff);
	
	return val;
}

//line is "1.463319e+001 -3.108064e+000 7.000000e+001"
int parse_line1(char buff[1000])
{
	int i, str, ord, len;
	float ver[3];
	
	len = strlen(buff);
	ord= str= 0;
	
	for (i= str; i<len; i++)
	{
		if (buff[i] == ' ')	
		{
			ver[ord] = parse_Co_ord(buff, str, i);
			ord++;
			str = i+1;
		}
	}
	
	ver[2] = parse_Co_ord(buff, str, len);

	return (chq_vertex_in_map(ver));
	
	//////////////////////////////////
	/*j = orig_verts;

	verts[j].x = ver[0];	verts[j].y = ver[1];	verts[j].z = ver[2];

	for (i=0; i<j; i++)
	{
		if ( (isSame(verts[i].x, ver[0])) && (isSame(verts[i].y, ver[1])) 
										  && (isSame(verts[i].z, ver[2])))
			return i;
	}
	orig_verts++;
	return (j);*/
	//////////////////////////////////

	//insert the hash table here
	/*e1.v1 = ver[0];		e1.v2 = ver[1];		e1.v3 = ver[2];		e1.key = (int)(e1.v1 + e1.v2 +e1.v3);
	val = myhash->Lookup(e1);

	if (val == -1)
	{
		j = orig_verts;
		verts[j].x = ver[0];	verts[j].y = ver[1];	verts[j].z = ver[2];
		orig_verts++;
		return j;
	}
	else
	{
		return val;
	}*/
	
}

void parse_poly()
{
	int i, ver[3];	
	char buff[5000];
	
	fp = fopen(ifname, "r");
	if (fp == NULL)
	{
		printf("error in opening the file. Exiting\n");
		exit(0);
	}
	
	/*First read the dammed file once to get the #tris and #points.*/
	i=0;
	
	while (1)
	{
		if (fgets(buff, 5000, fp) ==0)
			break;
		
		if (fgets(buff, 51000, fp) ==0)
			break;
		
		fgets(buff, 5000, fp);
		fgets(buff, 5000, fp);
		i++;
	}
	
	fclose(fp);
	
	no_of_tris =i;
	orig_verts=i*3;
	tris = (tempTri*)(malloc(sizeof(tempTri) * no_of_tris*3));
	verts = (tempVert*)(malloc(sizeof(tempVert) *orig_verts));
	orig_verts=0;
	init_the_maps();
	
	printf("File contains %d tempTris\n", no_of_tris);
	printf("Finished reading once and allocating. now going for the kill !! \n");
	
	/*Now, read the dammed file and populate the DS's*/
	fp = fopen(ifname, "r");
	if (fp == NULL)
	{
		printf("error in opening the file. Exiting\n");
		exit(0);
	}
	
	i=0;
	while (1)
	{
		if (fgets(buff, 5000, fp) ==0)
			break;
		
		if (fgets(buff, 5000, fp) ==0)
			break;
		ver[0] = parse_line1(buff);
		
		fgets(buff, 5000, fp);	ver[1] = parse_line1(buff);
		fgets(buff, 5000, fp);	ver[2] = parse_line1(buff);
		
		tris[i].v1 = ver[0];		tris[i].v2 = ver[1];		tris[i].v3 = ver[2];
		i++;
		
		if (i%500 == 0)
			printf("finished %d tempTris\n", i);
	}
	
	printf("phew!!. all done %d tempTris and %d verts\n", no_of_tris, orig_verts);
	parse_raw(2);
	exit(0);
}

/*
* From IPOLY to RAW
*/

void parse_ipoly()
{
	int i, crap[6];
	float temp[6];
	
	fp = fopen(ifname, "r");
	if (fp == NULL)
	{
		printf("error in opening the file. Exiting\n");
		exit(0);
	}
	
	if (fscanf(fp,"%d %d %d %d %d %d %d", &orig_verts, &crap[0], &no_of_tris,  &crap[1], 
										  &crap[2],  &crap[3], &crap[4]) == EOF)
	{
		printf("Input file is not valid....Exiting...\n");
		exit(0);
	}

	fscanf(fp,"%d %d %d", &crap[0],  &crap[1],  &crap[2]);
	
	tris = (tempTri*)(malloc(sizeof(tempTri) * no_of_tris));
	verts = (tempVert*)(malloc(sizeof(tempVert) *orig_verts));

	/*Now, read the dammed file and populate the DS's*/
	for (i=0; i<orig_verts; i++)
	{
		if (fscanf(fp,"%f %f %f %f %f %f", &temp[0], &temp[1], &temp[2], &temp[3], &temp[4],
										   &temp[5] ) == EOF)		
		{
			printf("Input file has to have %d Vertices....Exiting...\n",orig_verts);
			exit(0);
		}
		
		verts[i].x = temp[0];	verts[i].y = temp[1];	verts[i].z = temp[2];
		if (!(i% 5000))
			printf("still working on vertices !!!! %d \n",i);
	}
	
	printf("Finished reading the Vertices.. Now reading the tempTris\n");
	fscanf(fp,"%d %d", &crap[0],  &crap[1]);
	
	for (i=0; i<no_of_tris; i++)
	{
		fscanf(fp,"%d", &crap[0]);

		if (fscanf(fp,"%d %d %d", &tris[i].v1, &tris[i].v2, &tris[i].v3 ) == EOF)
		{
			printf("Input file has to have %d tempTris....Exiting...\n",no_of_tris);
			exit(0);
		}

		if (!(i% 5000))
			printf("still working on tempTris !!!! %d \n",i);
	}
	
	fclose(fp);
	printf("Finished reading the %d verts and %d tempTris\n", orig_verts, no_of_tris);
}

/*
* From GLX to RAW
*/

void fix_file()
{
	int i, temp=0;

	if (maxInd != orig_verts)
	{
		printf("The file is not even a 1-based raw file. Its something very very wierd\n");
		
		temp = maxInd - orig_verts + 1;
		for (i=0; i<no_of_tris; i++)
		{
			tris[i].v1 = tris[i].v1 - temp;
			tris[i].v2 = tris[i].v2 - temp;
			tris[i].v3 = tris[i].v3 - temp;
		}
	}
	else
	{
		printf("The file is a 1-based raw file.\n");

		for (i=0; i<no_of_tris; i++)
		{
			tris[i].v1--;		tris[i].v2--;		tris[i].v3--;
		}
	}

	printf("Hopefully, the file is now fixed.\n");
}

void parse_glx()
{
	int i;	
	float temp[6];	

	fp = fopen(ifname, "r");
	if (fp == NULL)
	{
		printf("error in opening the file. Exiting\n");
		exit(0);
	}

	if (fscanf(fp,"%d %d", &orig_verts, &no_of_tris) == EOF)
	{
		printf("Input file is not valid....Exiting...\n");
		exit(0);
	}

	tris = (tempTri*)(malloc(sizeof(tempTri) * no_of_tris*3));
	verts = (tempVert*)(malloc(sizeof(tempVert) *orig_verts));
	

	for (i=0; i<orig_verts; i++)
	{
		if (fscanf(fp,"%f %f %f %f %f %f", &temp[0], &temp[1], &temp[2], &temp[3], &temp[4],
										   &temp[5] ) == EOF)		
		{
			printf("Input file has to have %d Vertices....Exiting...\n",orig_verts);
			exit(0);
		}
		
		verts[i].x = temp[0];	verts[i].y = temp[1];	verts[i].z = temp[2];
		if (!(i% 5000))
			printf("still working on vertices !!!! %d \n",i);
	}
	
	printf("Finished reading the Vertices.. Now reading the tempTris\n");

	for (i=0; i<no_of_tris; i++)
	{
		if (fscanf(fp,"%d %d %d", &tris[i].v1, &tris[i].v2, &tris[i].v3 ) == EOF)
		{
			printf("Input file has to have %d tempTris....Exiting...\n",no_of_tris);
			exit(0);
		}

		if (maxInd < tris[i].v1) maxInd = tris[i].v1;
		if (maxInd < tris[i].v2) maxInd = tris[i].v2;
		if (maxInd < tris[i].v3) maxInd = tris[i].v3;

		if (!(i% 5000))
			printf("still working on tempTris !!!! %d \n",i);
	}
	
	fclose(fp);

	if (maxInd >= orig_verts)
	{
		printf("\n\n\nThe input file is invalid. The maximum vert index is %d, while the total pts is %d\n\n\n", maxInd, orig_verts);
		fix_file();
	}

	printf("Finished reading the file with %d tempTris and %d vertices \n", no_of_tris, orig_verts);
}

/*
* From FAT to RAW
*/
int  ReadHowmanyComponents()
{
   char  str[256];
   float  f0, f1, f2, f3, f4, f5, f6,f7,f8,f9,f10,f11;
   int    nvert, nface,component;

   fp = fopen(ifname, "r");
   if (fp==NULL) {
      return(0);
   }

   fgets(str, 256, fp);
   sscanf(str, "%d %d", &nvert, &nface);

   fgets(str, 256, fp);
   component = sscanf(str, "%f %f %f %f %f %f %f %f %f %f %f %f\n", 
                            &f0, &f1, &f2, &f3, &f4, &f5,
                            &f6, &f7, &f8, &f9, &f10, &f11);

   fclose(fp);
   return(component);
}

void parse_fat()
{
	int i, comps;	
	float temp[12];	

	fp = fopen(ifname, "r");
	if (fp == NULL)
	{
		printf("error in opening the file. Exiting\n");
		exit(0);
	}

	if (fscanf(fp,"%d %d", &orig_verts, &no_of_tris) == EOF)
	{
		printf("Input file is not valid....Exiting...\n");
		exit(0);
	}
	fclose(fp);

	tris = (tempTri*)(malloc(sizeof(tempTri) * no_of_tris*3));
	verts = (tempVert*)(malloc(sizeof(tempVert) *orig_verts*2));

	comps = ReadHowmanyComponents();

	fp = fopen(ifname, "r");
	if (fp==NULL) 
	{
		printf("error in opening the file. Exiting\n");
		exit(0);
	}
	fscanf(fp,"%d %d", &orig_verts, &no_of_tris);

	switch (comps)
	{
	case 6:
		for (i=0; i<orig_verts; i++)
		{
			if (fscanf(fp,"%f %f %f %f %f %f", &temp[0], &temp[1], &temp[2], &temp[3], &temp[4],
											   &temp[5] ) == EOF)		
			{
				printf("Input file has to have %d Vertices....Exiting...\n",orig_verts);
				exit(0);
			}
			
			verts[i].x = temp[0];	verts[i].y = temp[1];	verts[i].z = temp[2];
			verts[orig_verts+i].x = temp[3];	verts[orig_verts+i].y = temp[4];	verts[orig_verts+i].z = temp[5];
			if (!(i% 5000))
				printf("still working on vertices !!!! %d \n",i);
		}
	break;
	case 12:
		for (i=0; i<orig_verts; i++)
		{
			if (fscanf(fp,"%f %f %f %f %f %f %f %f %f %f %f %f", &temp[0], &temp[1], &temp[2], 
					      &temp[3], &temp[4], &temp[5], &temp[6], &temp[7], &temp[8], 
						  &temp[9], &temp[10], &temp[11] ) == EOF)		
			{
				printf("Input file has to have %d Vertices....Exiting...\n",orig_verts);
				exit(0);
			}
			
			verts[i].x = temp[0];	verts[i].y = temp[1];	verts[i].z = temp[2];
			verts[orig_verts+i].x = temp[3];	verts[orig_verts+i].y = temp[4];	verts[orig_verts+i].z = temp[5];

			if (!(i% 5000))
				printf("still working on vertices !!!! %d \n",i);
		}
	break;
	default:
		printf("The vertices should contain only \" X Y Z X Y Z NX NY NZ NX NY NZ \" or	\"X Y Z X Y Z\" i.e. only 12 or 6 components per line.\n");
		break;
	}

	printf("Finished reading the Vertices.. Now reading the tempTris\n");

	for (i=0; i<no_of_tris; i++)
	{
		if (fscanf(fp,"%d %d %d", &tris[i].v1, &tris[i].v2, &tris[i].v3 ) == EOF)
		{
			printf("Input file has to have %d tempTris....Exiting after %d...\n",no_of_tris, i);
			exit(0);
		}

		if (maxInd < tris[i].v1) maxInd = tris[i].v1;
		if (maxInd < tris[i].v2) maxInd = tris[i].v2;
		if (maxInd < tris[i].v3) maxInd = tris[i].v3;

		if (!(i% 5000))
			printf("still working on tempTris !!!! %d \n",i);
	}
	
	fclose(fp);

	if (maxInd >= orig_verts)
	{
		printf("\n\n\nThe input file is invalid. The maximum vert index is %d, while the total pts is %d\n\n\n", maxInd, orig_verts);
		fix_file();
	}

	printf("Finished reading the file with %d tempTris and %d vertices \n", no_of_tris, orig_verts);
}

void write_fat()
{
	int var, i;
	char buff[5000];

	for (var=0; var<2; var++)
	{
		sprintf(buff, "%d%s", var, ofname);
	
		fp = fopen(buff, "w");
		if (fp == NULL)
		{
			printf("error in opening the file. Exiting\n");
			exit(0);
		}
		
		fprintf(fp, "%d %d\n", orig_verts, no_of_tris);
		
		for (i=0; i<orig_verts; i++)
			fprintf(fp, "%f %f %f\n", verts[var*orig_verts +i].x, verts[var*orig_verts +i].y, verts[var*orig_verts +i].z);
		
		for (i=0; i<no_of_tris; i++)
			fprintf(fp, "%d %d %d\n", tris[i].v1, tris[i].v2, tris[i].v3);
		
		fclose(fp);
	}

	if(tris != NULL)
		free(tris);

	if (verts != NULL)
		free(verts);
}

/*
* From STL to RAW
*/

//line is "      vertex 1.463319e+001 -3.108064e+000 7.000000e+001"
int parse_line(char buff[1000])
{
	int i, j, str, ord, len;
	float ver[3];
	
	len = strlen(buff);
	j=0;
	for (i=0; i<len; i++)
	{
		if (buff[i] != ' ')
			break;
	}

	str=i+7; //i.e. the strlen of ("vertex")
	ord=0;
	
	for (i= str; i<len; i++)
	{
		if (buff[i] == ' ')	
		{
			ver[ord] = parse_Co_ord(buff, str, i);
			ord++;
			str = i+1;
		}
	}
	
	ver[2] = parse_Co_ord(buff, str, len);

	return (chq_vertex_in_map(ver));
	
	/*j = orig_verts;
	
	for (i=0; i<j; i++)
	{
		if ( (isSame(verts[i].x, ver[0])) && (isSame(verts[i].y, ver[1])) 
			&& (isSame(verts[i].z, ver[2])))
			return i;
	}
	
	verts[j].x = ver[0];	verts[j].y = ver[1];	verts[j].z = ver[2];
	orig_verts++;
	
	return (j);*/
	
	//insert the hash table here
	/*e1.v1 = ver[0];		e1.v2 = ver[1];		e1.v3 = ver[2];		e1.key = (int)(e1.v1 + e1.v2 +e1.v3);
	val = myhash->Lookup(e1);

	if (val == -1)
	{
		j = orig_verts;
		verts[j].x = ver[0];	verts[j].y = ver[1];	verts[j].z = ver[2];
		orig_verts++;
		return j;
	}
	else
	{
		return val;
	}*/
}

void parse_stl()
{
	int i, ver[3];	
	char buff[1000];
	
	fp = fopen(ifname, "r");
	if (fp == NULL)
	{
		printf("error in opening the file. Exiting\n");
		exit(0);
	}
	
	/*First read the dammed file once to get the #tris and #points.*/
	i=0;
	fgets(buff, 1000, fp);
	
	while (1)
	{
		if (fgets(buff, 1000, fp) ==0)
			break;
		
		if (fgets(buff, 1000, fp) ==0)
			break;
		
		fgets(buff, 1000, fp);
		fgets(buff, 1000, fp);
		fgets(buff, 1000, fp);
		
		fgets(buff, 1000, fp);
		fgets(buff, 1000, fp);
		i++;
	}
	
	fclose(fp);
	
	no_of_tris =i;
	orig_verts=i*3;
	tris = (tempTri*)(malloc(sizeof(tempTri) * no_of_tris*3));
	verts = (tempVert*)(malloc(sizeof(tempVert) *orig_verts));
	orig_verts=0;
	init_the_maps();
	
	printf("File contains %d tempTris\n", no_of_tris);
	printf("Finished reading once and allocating. now going for the kill !! \n");
	
	/*Now, read the dammed file and populate the DS's*/
	fp = fopen(ifname, "r");
	if (fp == NULL)
	{
		printf("error in opening the file. Exiting\n");
		exit(0);
	}
	
	i=0;
	fgets(buff, 1000, fp);
	while (1)
	{
		if (fgets(buff, 1000, fp) ==0)
			break;
		
		if (fgets(buff, 1000, fp) ==0)
			break;
		
		fgets(buff, 1000, fp);	ver[0] = parse_line(buff);
		fgets(buff, 1000, fp);	ver[1] = parse_line(buff);
		fgets(buff, 1000, fp);	ver[2] = parse_line(buff);
		
		fgets(buff, 1000, fp);
		fgets(buff, 1000, fp);
		
		tris[i].v1 = ver[0];		tris[i].v2 = ver[1];		tris[i].v3 = ver[2];
		i++;
		
		if (i%500 == 0)
			printf("finished %d tempTris\n", i);
	}
	
	printf("phew!!. all done %d tempTris and %d verts\n", no_of_tris, orig_verts);
}


/*
* Here, the common write_file routine.
*/
void write_raw()
{
	int i;
	
	fp = fopen(ofname, "w");
	if (fp == NULL)
	{
		printf("error in opening the file. Exiting\n");
		exit(0);
	}
	
	fprintf(fp, "%d %d\n", orig_verts, no_of_tris);
	
	for (i=0; i<orig_verts; i++)
		fprintf(fp, "%f %f %f\n", verts[i].x, verts[i].y, verts[i].z);
	
	for (i=0; i<no_of_tris; i++)
		fprintf(fp, "%d %d %d\n", tris[i].v1, tris[i].v2, tris[i].v3);
	
	fclose(fp);

	if(tris != NULL)
		free(tris);

	if (verts != NULL)
		free(verts);
}

void usage()
{
	printf("usage: poly -n <number> -in <input file> -out <output file> -h\n");
	printf("       <number> can be any of the following:\n");
	printf("      -1: for RAW to STL file\n");
	printf("       0: for fixing a RAW file\n");
	printf("       1: for POLY to RAW\n");
	printf("       2: for IPOLY to RAW\n");
	printf("       3: for GLX to RAW\n");
	printf("       4: for FAT to RAW\n");
	printf("       5: for STL to RAW\n");
	printf("       6: for VRML to RAW\n");
	printf("       7: for IV to RAW\n");
	printf("       8: for WRL to RAW\n");
	printf("       9: for RAW to STL\n");
	printf("\nInput and Output file names are the second and third arguments.\n");
	printf("\nOptional parameter <-h> prints out this usage message.\n");

}

void parse_args(int argc, char *argv[]) 
{
	int i;

	for (i=1; i<argc; i++)
	{
		if (!strcmp(argv[i],"-h") && !strcmp(argv[i],"-H"))
		{
			usage();
			exit(0);
		}
		
		if (!strcmp("-in",argv[i]) || !strcmp("-IN",argv[i]))
			ifname=argv[++i];

		if (!strcmp("-out",argv[i]) || !strcmp("-OUT",argv[i]))
			ofname=argv[++i];

		if (!strcmp("-n",argv[i]) || !strcmp("-N",argv[i]))
			type=atoi(argv[++i]);
	}
}

int main (int argc, char** argv)
{
	type = maxInd = -1;
	
	parse_args(argc, argv);

	if ((ifname == NULL) || (ofname == NULL))
	{
		usage();
		exit(0);
	}

	switch (type)
	{
		case -1:
			parse_raw(0);
			return(0);
			break;
		case 0:
			parse_raw(1); //to FIX a RAW file.
			return(0);
			break;

		case 1: //to convert the POLY file to RAW
			parse_poly();
		break;

		case 2: //to convert the IPOLY file to RAW
			parse_ipoly();
		break;

		case 3: //to convert the GLX file to RAW
			parse_glx();
		break;

		case 4: //to convert the FAT file to RAW
			parse_fat();
			write_fat();
			return(0);
		break;
		
		case 5: //to convert the STL file to RAW
			parse_stl();
		break;

		case 6: //to convert the VRML file to RAW
			vrml2raw();
			return(0);
		break;

		case 7: //to convert the IV file to RAW
			iv2raw();
			return(0);
		break;

		case 8: //to convert the WRL file to RAW
			wrl2raw();
			return(0);
		break;

		case 9: //to convert the RAW file to STL
			raw_2_stl();
			return(0);
			break;
		
		default:
			printf("unknown case in input argument.\n");
			usage();
		break;
	}
	
	write_raw();
}
