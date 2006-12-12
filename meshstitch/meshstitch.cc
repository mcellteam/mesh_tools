#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "classes.cc"
#include "functions.cc"

int main(int argc,char *argv[]){

	if (argc != 5)
	{
		printf("\nSyntax: meshstitch input_file1 input_file2 z_value orientation\n\n");
		printf("Detail: The z_value refers to the z axis coordinate of\n");
		printf("        the common contour between the two input meshes. \n");
		printf("        Orientation is '1' if the polygons to be removed\n");
		printf("        from input_file2 lie on the +z side of the z_value, and\n");
		printf("        and orientation is '-1' if the polygons to be removed\n");
		printf("        from input_file2 lie on the -z side of the z_value.\n");
		printf("        The polygons to be removed from input_file1 are assumed\n");
		printf("        to lie on the opposite of the z_value specified by orientation.\n");
		printf("Output: input_file1.clipped and input_file2.clipped.\n\n");
		return 1;
	}

	////////// declare variables /////////

	char *infile1,*infile2,outfile1[128],outfile2[128],line[2048];
	char *str,*eptr,name[32];
	FILE *F;
	void_list *q,*qq,*p,*pp,*prev,*ptr1,*ptr2;
	bool flag,new_flag;
	double epsilon = 1e-14;
	int max_vert1,max_vert2,max_face1,max_face2,i,j;

	// vertex linked lists for first and second file
	void_list *vlh1,*vl1,*vlh2,*vl2;
	Vertex *v1,*v2;
	vlh1 = vlh2 = NULL;

	// face linked lists for first and second file
	void_list *flh1,*fl1,*flh2,*fl2;
	Face *f1,*f2;
	flh1 = flh2 = NULL;

	// vertex with target z_value linked lists for first and second file
	void_list *v1_match_h,*v2_match_h;
	v1_match_h = v2_match_h = NULL;

	// shared vertex linked lists for first and second file
	void_list *v1_shared,*v1_shared_h,*v2_shared,*v2_shared_h;
	v1_shared_h = v2_shared_h = NULL;

	// candidate bad vertex linked lists for first and second file
	void_list *cbv1h,*cbv2h;
	cbv1h = cbv2h = NULL;

	// candidate bad face linked lists for first and second file
	void_list *cbf1h,*cbf2h;
	cbf1h = cbf2h = NULL;

	// bad vertex linked lists for first and second file
	void_list *bv1,*bv1h,*bv2,*bv2h;
	bv1h = bv2h = NULL;

	// bad face linked lists for first and second file
	void_list *bf1,*bf1h,*bf2,*bf2h;
	bf1h = bf2h = NULL;

	// extra vertex face search pool
	void_list *evfsp1,*evfsp2;
	evfsp1 = evfsp2 = NULL;

	// pointers to start of vertex and face linked lists for first and second file
	void_list *v1end,*v2end,*f1end,*f2end;
	v1end = v2end = f1end = f2end = NULL;

	// variables used to find extra vertices
	void_list *ev1h,*ev2h;
	ev1h = ev2h = NULL;

	// get z_value
	int z_value;
	z_value = (int) (strtod(argv[3],&eptr));

	// get target z_value
	int orientation;
	orientation = (int) strtod(argv[4],&eptr);
	int val;

	////////// get data /////////

	// open first file
	infile1 = argv[1];
	F = fopen(infile1,"r");
	if (!F) {
	printf("Couldn't open input file %s\n",infile1);
	return 1;
	}

	// for every line in first file
	for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) {

		// skip leading whitespace
		while (strchr(" \t,",*str)!=NULL) { str++;}

		// if first character is V for Vertex, add new linked list class instance
		if (strchr("V",*str)!=NULL){
			vl1 = new void_list();
			vl1->next = vlh1;
			v1 = new Vertex(str);
			vl1->data = (void*)v1;
			vlh1 = vl1;
		} 
		// if first character is F for Face, add new linked list class instance
		else if (strchr("F",*str)!=NULL){
			fl1 = new void_list();
			fl1->next = flh1;
			f1 = new Face(str);
			fl1->data = (void*)f1;
			flh1 = fl1;
		}
	}
	fclose(F);

	// opend second file
	infile2 = argv[2];
	F = fopen(infile2,"r");
	if (!F) {
	printf("Couldn't open input file %s\n",infile2);
	return 1;
	}

	// for every line in second file
	for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) {

		// skip leading whitespace
		while (strchr(" \t,",*str)!=NULL) { str++;}

		// if first character is V for Vertex, add new linked list class instance
		if (strchr("V",*str)!=NULL){
			vl2 = new void_list();
			vl2->next = vlh2;
			v2 = new Vertex(str);
			vl2->data = (void*)v2;
			vlh2 = vl2;
		} 
		// if first character is F for Face, add new linked list class instance
		else if (strchr("F",*str)!=NULL){
			fl2 = new void_list();
			fl2->next = flh2;
			f2 = new Face(str);
			fl2->data = (void*)f2;
			flh2 = fl2;
		}
	}
	fclose(F);

	// maked linked lists bidirectional
	vlh1 = addPrevious(vlh1);
	flh1 = addPrevious(flh1);
	vlh2 = addPrevious(vlh2);
	flh2 = addPrevious(flh2);

    // find maximum vertex index
	max_vert1=maxVert(vlh1);
	max_vert2=maxVert(vlh2);
	// find max face index
	max_face1=maxFace(flh1);
	max_face2=maxFace(flh2);

    // create array of pointers to vertex list
    Vertex* vert_array1[max_vert1+1];
    Vertex* vert_array2[max_vert2+1];

    // fill arrays
    for (p=vlh1;p!=NULL;p=p->next) {
		vert_array1[((Vertex*)p->data)->index]=(Vertex*)p->data;
    }
    for (p=vlh2;p!=NULL;p=p->next) {
		vert_array2[((Vertex*)p->data)->index]=(Vertex*)p->data;
    }

	////////// find shared vertices /////////

	// build linked list of first and second vertex elements with z_value matching target
	v1_match_h = matchZ(vlh1,z_value);
	v2_match_h = matchZ(vlh2,z_value);

	// maked linked lists bidirectional
	v1_match_h = addPrevious(v1_match_h);
	v2_match_h = addPrevious(v2_match_h);

	if(0){strcpy(name,"v1_match_h");printVertices(v1_match_h,name);}
	if(0){strcpy(name,"v2_match_h");printVertices(v2_match_h,name);}

	// find shared vertices
	findSharedVertices(v1_match_h,v2_match_h,v1_shared_h,v2_shared_h,epsilon);

	if(0){strcpy(name,"v1_shared_h");printVertices(v1_shared_h,name);}
	if(0){strcpy(name,"v2_shared_h");printVertices(v2_shared_h,name);}

	// get extra vertex face search pool
	evfsp1 = getFaceSearchPool(flh1,v1_match_h,z_value,vert_array1);
	evfsp2 = getFaceSearchPool(flh2,v2_match_h,z_value,vert_array2);

	if(0){strcpy(name,"evfsp1");printFaces(evfsp1,name);}
	if(0){strcpy(name,"evfsp2");printFaces(evfsp2,name);}

	// group vertices
	Group g1,g2;
	getContours(g1,evfsp1,vert_array1,max_vert1,z_value,v1_shared_h);
	getContours(g2,evfsp2,vert_array2,max_vert2,z_value,v2_shared_h);

	if(0){strcpy(name,"FILE 1");printContourFacesVertices(g1,name);}
	if(0){strcpy(name,"FILE 2");printContourFacesVertices(g2,name);}

	///// STRATEGY /////
	// (1) collect all vertices with z value half section
	//		thickness away from target z_value in orientation direction
	// (2) collect all polygons connected to candidate bad vertices from step 1
	// (3) for each collected polygon, if at least one vertex is a shared vertex from above
	// 		then any vertex on that polygon with z value half section
	//		thickness away from target z_value in orientation direction is bad
	// (4) for each collected polygon, if at least one vertex is a shared vertex
	// 		then polygon is bad

	// gather candidate bad vertices from first and second vertex linked list
	val = z_value+orientation*25*(-1);
	cbv1h = candidateBadVertices(vlh1,val);
	val = z_value+orientation*25;
	cbv2h = candidateBadVertices(vlh2,val);

	// contour_tiler occassionally creates an extra vertex, so need to find them
	// gather extra vertices from first and second match vertex linked lists
	extraVertices(g1,v1_shared_h);
	extraVertices(g2,v2_shared_h);

	///// find faces that have at least one vertex from cbv lists /////
	cbf1h = candidateBadFaces(flh1,cbv1h);
	cbf2h = candidateBadFaces(flh2,cbv2h);

	// maked linked lists bidirectional
	cbv1h = addPrevious(cbv1h);
	cbv2h = addPrevious(cbv2h);
	cbf1h = addPrevious(cbf1h);
	cbf2h = addPrevious(cbf2h);

	///// find bad vertices and faces from first vertex and first face linked lists
	identifyBadVerticesAndFaces(cbf1h,cbv1h,v1_shared_h,g1,bf1h,bv1h);

	// rescan candidate bad faces using known bad vertices
	// looking for bad faces that do not have face vertex that is shared
	// Apparently, all bad vertices are members of a bad faces with shared vertices,
	// so no further searching for bad vertices is performed in findBadFaces.
	bf1h = findBadFaces(cbf1h,bv1h,bf1h);

	///// find bad vertices and faces from second vertex and second face linked lists
	identifyBadVerticesAndFaces(cbf2h,cbv2h,v2_shared_h,g2,bf2h,bv2h);

	// rescan candidate bad faces using known bad vertices
	// looking for bad faces that do not have face vertex that is shared 
	// Apparently, all bad vertices are members of a bad faces with shared vertices,
	// so no further searching for bad vertices is performed in findBadFaces.
	bf2h = findBadFaces(cbf2h,bv2h,bf2h);

	if(0){strcpy(name,"bv1h");printVertices(bv1h,name);}
	if(0){strcpy(name,"bv2h");printVertices(bv2h,name);}

	if(0){strcpy(name,"vlh1");printVertices(vlh1,name);}
	if(0){strcpy(name,"vlh2");printVertices(vlh2,name);}

	if(0){strcpy(name,"bf1h");printFaces(bf1h,name);}
	if(0){strcpy(name,"bf2h");printFaces(bf2h,name);}

	if(0){strcpy(name,"flh1");printFaces(flh1,name);}
	if(0){strcpy(name,"flh2");printFaces(flh2,name);}

	// make linked lists bidirectional
	bv1h = addPrevious(bv1h);
	bv2h = addPrevious(bv2h);
	bf1h = addPrevious(bf1h);
	bf2h = addPrevious(bf2h);

	// need to identify the vertices on either side of the extra vertices
	lateralVertices(g1,bf1h,v1_shared_h,z_value,vert_array1);
	lateralVertices(g2,bf2h,v2_shared_h,z_value,vert_array2);

	convertLateral(vlh2,g1,epsilon);
	convertLateral(vlh1,g2,epsilon);

	val = z_value+orientation*25;
	gatherThirdDeleteFace(g1,flh2,val,vert_array2);
	val = z_value+orientation*25*(-1);
	gatherThirdDeleteFace(g2,flh1,val,vert_array1);

	if(0){strcpy(name,"FILE 1");printOtherContourData(g1,name);}
	if(0){strcpy(name,"FILE 2");printOtherContourData(g2,name);}

	if(0){printf("ADDED TO FILE 1\n");}
	addFacesVertices(g2,flh1,vlh1,max_vert1,max_face1);
	if(0){printf("ADDED TO FILE 2\n");}
	addFacesVertices(g1,flh2,vlh2,max_vert2,max_face2);

	// go through linked list backwards and add previous pointers
	vlh1 = addPrevious(vlh1);
	flh1 = addPrevious(flh1);
	vlh2 = addPrevious(vlh2);
	flh2 = addPrevious(flh2);

	///// remove bad vertices from first and second vertex list /////
	vlh1 = removeBadVertices(vlh1,bv1h);
	vlh2 = removeBadVertices(vlh2,bv2h);

	///// remove bad faces from first and second face list /////
	flh1 = removeBadFaces(flh1,bf1h);
	flh2 = removeBadFaces(flh2,bf2h);

	if(0){strcpy(name,"flh1");printFaces(flh1,name);}
	if(0){strcpy(name,"flh2");printFaces(flh2,name);}

	////////// write output to files //////////
	// set end pointers
    for (q=vlh1;q!=NULL;q=q->next) { v1end=q; if (q->next==NULL) break; }
    for (q=flh1;q!=NULL;q=q->next) { f1end=q; if (q->next==NULL) break; }
    for (q=vlh2;q!=NULL;q=q->next) { v2end=q; if (q->next==NULL) break; }
    for (q=flh2;q!=NULL;q=q->next) { f2end=q; if (q->next==NULL) break; }

	// open first file
	sprintf(outfile1,"%s.clipped",infile1);
	printf("\n%s written.\n",outfile1);
	F = fopen(outfile1,"w");
	if (!F) {
		printf("Couldn't open output file %s\n",outfile1);
		return 1;
	}

	// for every vertex in first file
	for (p=v1end;p!=NULL;p=p->previous) {
		sprintf(line,"Vertex %i  %.15g %.15g %.15g\n",
				((Vertex*)p->data)->index,
				((Vertex*)p->data)->x,
				((Vertex*)p->data)->y,
				((Vertex*)p->data)->z);
		fputs(line,F);
		if (p->previous==NULL) break;
	}

	// for every face in first file
	for (p=f1end;p!=NULL;p=p->previous) {
		sprintf(line,"Face %i  %i %i %i\n",
				((Face*)p->data)->index,
				((Face*)p->data)->v1,
				((Face*)p->data)->v2,
				((Face*)p->data)->v3);
		fputs(line,F);
		if (p->previous==NULL) break;
	}
	fclose(F);
	
	// open second file
	sprintf(outfile2,"%s.clipped",infile2);
	printf("%s written.\n\n",outfile2);
	F = fopen(outfile2,"w");
	if (!F) {
		printf("Couldn't open output file %s\n",outfile2);
		return 1;
	}

	// for every vertex in second file
	for (p=v2end;p!=NULL;p=p->previous) {
		sprintf(line,"Vertex %i  %.15g %.15g %.15g\n",
				((Vertex*)p->data)->index,
				((Vertex*)p->data)->x,
				((Vertex*)p->data)->y,
				((Vertex*)p->data)->z);
		fputs(line,F);
		if (p->previous==NULL) break;
	}

	// for every face in second file
	for (p=f2end;p!=NULL;p=p->previous) {
		sprintf(line,"Face %i  %i %i %i\n",
				((Face*)p->data)->index,
				((Face*)p->data)->v1,
				((Face*)p->data)->v2,
				((Face*)p->data)->v3);
		fputs(line,F);
		if (p->previous==NULL) break;
	}
	fclose(F);

	//cleanup
	cleanUp2V(vlh1);
	cleanUp2V(vlh2);
	cleanUp2F(flh1);
	cleanUp2F(flh2);
	cleanUp(bv1h);
	cleanUp(bv2h);
	cleanUp(bf1h);
	cleanUp(bf2h);
	cleanUp(v1_match_h);
	cleanUp(v2_match_h);
	cleanUp(v1_shared_h);
	cleanUp(v2_shared_h);
	cleanUp(cbf1h);
	cleanUp(cbf2h);
	cleanUp(cbv1h);
	cleanUp(cbv2h);
	cleanUp(evfsp1);
	cleanUp(evfsp2);
	return 0;
}
