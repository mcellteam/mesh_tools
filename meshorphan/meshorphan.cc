#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <dirent.h>
#include <errno.h>

using namespace std;

int distinguishable(double a,double b,double eps)
{
  double c;

  c=a-b;

  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;

  if (b<a) eps*=a;
  else eps*=b;
  return (c>eps);
}

class void_list
{
public:
  void_list *previous;
  void_list *next;
  void *data;
};

class Swap
{
public:
  int index_old;
  int index_new;
};

class Edge
{
public:
  int v1,v2; // vertex indices
  int c12,c21;  // count of times edge was traversed from v1 to v2 and v2 to v1
  void_list *f;
  bool crumb;
};

class Vertex
{
public:
  double x,y,z;
  int index;
  bool found;
  Vertex(char *triplet);
};

Vertex::Vertex(char *triplet)
{
	char val[80];
	char *eptr;
	int i;

	found = false;

	// get past 'Vertex'
  	while (strchr("Vertx",*triplet)!=NULL) {triplet++;}

	// grab vertex index
	while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
	i=0;
	while (strchr("0123456789+-eE.",*triplet)!=NULL)
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	index = (int) strtod(val,&eptr);
	if (val==eptr)
	{
		index=0;
		x=y=z=0;
		printf("Error in reading vertex index\n");
		return;
	}

	// grab x coord
	while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
	i=0;
	while (strchr("0123456789+-eE.",*triplet)!=NULL)
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	x = strtod(val,&eptr);
	if (val==eptr)
	{
		x=y=z=0;
		printf("Error in reading vertex\n");
		return;
	}

	// grab y coord
	while (strchr(" \t,",*triplet)!=NULL) triplet++;
	i=0;
	while (strchr("0123456789+-eE.",*triplet))
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	y = strtod(val,&eptr);
	if (val==eptr)
	{
		x=y=z=0;
		printf("Error in reading vertex\n");
		return;
	}

	// grab z coord
	while (strchr(" \t,",*triplet)!=NULL) triplet++;
	i=0;
	while (strchr("0123456789+-eE.",*triplet))
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	z = strtod(val,&eptr);
	if (val==eptr)
	{
		x=y=z=0;
		printf("Error in reading vertex\n");
		return;
	}
}

class Face
{
public:
	int index;	// Face index
	int v1,v2,v3;	// vertex indices
	Face(char *triplet);
};

Face::Face(char *triplet)
{
	char val[80];
	char *eptr;
	int i;

	// get past 'Face'
  	while (strchr("Face",*triplet)!=NULL) {triplet++;}

	// grab Face index
	while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
	i=0;
	while (strchr("0123456789+-eE.",*triplet)!=NULL)
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	index = (int) strtod(val,&eptr);
	if (val==eptr)
	{
		index=0;
		v1=v2=v3=0;
		printf("Error in reading face index\n");
		return;
	}

	// grab first vertex index
	while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
	i=0;
	while (strchr("0123456789+-eE.",*triplet)!=NULL)
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	v1 = (int) strtod(val,&eptr);
	if (val==eptr)
	{
		v1=v2=v3=0;
		printf("Error in reading vertex index\n");
		return;
	}

	// grab second vertex index
	while (strchr(" \t,",*triplet)!=NULL) triplet++;
	i=0;
	while (strchr("0123456789+-eE.",*triplet))
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	v2 = (int) strtod(val,&eptr);
	if (val==eptr)
	{
		v1=v2=v3=0;
		printf("Error in reading vertex index\n");
		return;
	}

	// grab third vertex index
	while (strchr(" \t,",*triplet)!=NULL) triplet++;
	i=0;
	while (strchr("0123456789+-eE.",*triplet))
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	v3 = (int) strtod(val,&eptr);
	if (val==eptr)
	{
		v1=v2=v3=0;
		printf("Error in reading vertex index\n");
		return;
	}
}

void_list * scanDir(void_list *files) {
	void_list *p;
	char *name;
    std::string str;
    std::string::size_type found;

	DIR *pdir;								// pointer to a directory data structure
	struct dirent *pent;					// pointer to dirent structure

    pdir = opendir("./");
    if (!pdir) {
        printf ("opendir() failure; could not open %s. terminating","./");
        exit(1);
    }
    errno = 0;
    while ((pent=readdir(pdir))){
		// copy char array to string
		str = pent->d_name;
		// if file of typ *.mesh
		found = str.find(".mesh",0);
        if (found != std::string::npos) {
			// save filename
			p = new void_list();
			p->next=files;
			name = new char[1024];
			strcpy(name,str.c_str());
			p->data=(void*)name;
			files=p;
			// print file found to screen
	        printf("file found: %s\n",name);
		}
    }
    if (errno) {
        printf ("readdir() failure; terminating");
        exit(1);
    }
    closedir(pdir);
	return files;
}

void getData(char *infile,void_list *&flh,void_list *&vlh){

	char line[2048],*str,*eptr;
	FILE *F;
	void_list *vl;
	Vertex *v;
	void_list *fl;
	Face *f;

	// open first file
	F = fopen(infile,"r");
	if (!F) {
	printf("Couldn't open input file %s\n",infile);
	return;
	}

	// for every line in file
	for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) {

		// skip leading whitespace
		while (strchr(" \t,",*str)!=NULL) { str++;}

		// if first character is V for Vertex, add new linked list class instance
		if (strchr("V",*str)!=NULL){
			vl = new void_list();
			vl->next = vlh;
			v = new Vertex(str);
			vl->data = (void*)v;
			vlh = vl;
		} 
		// if first character is F for Face, add new linked list class instance
		else if (strchr("F",*str)!=NULL){
			fl = new void_list();
			fl->next = flh;
			f = new Face(str);
			fl->data = (void*)f;
			flh = fl;
		}
	}
	fclose(F);
}

void_list* addPrevious(void_list* L) {
    // go through linked list backwards and add previous pointers
    void_list *p,*prev;
    prev = NULL;
    for (p=L;p!=NULL;p=p->next) {
        p->previous = prev;
        prev = p;
        if (p->next==NULL) break;
    }
    return L;
}

int maxVert(void_list *L){
    void_list *p;
    int max=0;
    for (p=L;p!=NULL;p=p->next) {
        if (max<((Vertex*)p->data)->index) {max=((Vertex*)p->data)->index;}
    }
    return max;
}

int maxFace(void_list *L){
    void_list *p;
    int max=0;
    for (p=L;p!=NULL;p=p->next) {
        if (max<((Face*)p->data)->index) {max=((Face*)p->data)->index;}
    }
    return max;
}

void writeData(char *name,void_list *vlh,void_list *flh,bool print_flag){
	void_list *q,*vend,*fend;
	Vertex *v;
	Face *f;
	FILE *F;
	char line[1024];
	int count=0;
	for (q=vlh;q!=NULL;q=q->next) {vend=q;}
	for (q=flh;q!=NULL;q=q->next) {fend=q;}

	// open file
	if (print_flag){
		F = fopen(name,"w");
		if (!F) {
			printf("Couldn't open output file %s\n",name);
			exit(0);
		}
	}
	
	// write out vertex linked list
	for (q=vend;q!=NULL;q=q->previous) {
		v = (Vertex*)q->data;
		if (v->found){
			if (print_flag){
		      	sprintf(line,"Vertex %i  %.15g %.15g %.15g\n",v->index,v->x,v->y,v->z);
				fputs(line,F);
			}
		} else {
	      	fprintf(stderr,"ORPHAN Vertex %i  %.15g %.15g %.15g\n",v->index,v->x,v->y,v->z);
			count++;
		}
	}
	// write out face linked list
	for (q=fend;q!=NULL;q=q->previous) {
		if (print_flag){
			f=(Face*)q->data;
	      	sprintf(line,"Face %i  %i %i %i\n",f->index,f->v1,f->v2,f->v3);
			fputs(line,F);
		}
	}
	if (print_flag){
		fclose(F);
		fprintf(stderr,"Removed %i orphaned vertices from %s\n",count,name);
	} else {
		fprintf(stderr,"Found %i orphaned vertices in %s\n",count,name);
	}

}

void printData(void_list *vlh,void_list *flh,bool print_flag){
	void_list *q,*vend,*fend;
	Vertex *v;
	Face *f;
	int count=0;
	for (q=vlh;q!=NULL;q=q->next) {vend=q;}
	for (q=flh;q!=NULL;q=q->next) {fend=q;}

	////////// write data to stdout /////////
	// write out vertex linked list
	for (q=vend;q!=NULL;q=q->previous) {
		v = (Vertex*)q->data;
		if (v->found){
      		if (print_flag){printf("Vertex %i  %.15g %.15g %.15g\n",v->index,v->x,v->y,v->z);}
		} else {
	      	fprintf(stderr,"ORPHAN Vertex %i  %.15g %.15g %.15g\n",v->index,v->x,v->y,v->z);
			count++;
		}
		
	}
	// write out face linked list
	for (q=fend;q!=NULL;q=q->previous) {
		if (print_flag){
			f=(Face*)q->data;
	      	printf("Face %i  %i %i %i\n",f->index,f->v1,f->v2,f->v3);
		}
	}
	if (print_flag){ fprintf(stderr,"Removed %i orphaned vertices\n",count);}
	else {fprintf(stderr,"Found %i orphaned vertices\n",count);}
}

void cleanup(void_list *f,void_list *v){
	void_list *p,*q;
	// faces
	p=f;
	while (p!=NULL) {
		q=p->next;
		delete (Face*)p->data;
		delete p;
		p=q;
	}
	// vertices
	p=v;
	while (p!=NULL) {
		q=p->next;
		delete (Vertex*)p->data;
		delete p;
		p=q;
	}
}

int main(int argc,char *argv[]){

	if (argc != 2)
	{
		printf("\nSyntax: meshorphan input_file|4all\n\n");
		printf("Description: Identifies and removes orphaned vertices.\n");
		printf("# removed vertices written to stderr.\n");
		printf("If 4all, then all mesh files in current directory are renumbered.\n");
		printf("Note 4all writes renumbered mesh back to original filename.\n\n");
		return 1;
	}

	////////// declare variables /////////
	void_list *q,*p,*files,*vend,*fend;
	Face *f;
	Vertex *v;
	files=NULL;
	bool flag,print_flag=true;
	int max_vertex;
	Vertex **verts;

	// vertex and face linked lists
	void_list *vlh,*flh;
	vlh=flh=NULL;

	// get files to analyze
	if(!strcmp(argv[1],"4all")){
		files=scanDir(files);
		flag=true;
	} else {
		files = new void_list();
		files->next=NULL;
		files->data=(void*)argv[1];
		flag=false;
	}

	// for all mesh files in current directory
	for(q=files;q!=NULL;q=q->next){
		fprintf(stderr,"\nremoving oprhaned vertices from %s...\n",(char*)q->data);
		////////// get data /////////
		getData((char*)q->data,flh,vlh);
		flh=addPrevious(flh);
		vlh=addPrevious(vlh);
        max_vertex=maxVert(vlh);
		verts = new Vertex*[max_vertex+1];
		// load arrays
		for (p=vlh;p!=NULL;p=p->next) {vend=p;}
		for (p=flh;p!=NULL;p=p->next) {fend=p;}
		for (p=vend;p!=NULL;p=p->previous) {verts[((Vertex*)p->data)->index]=(Vertex*)p->data;}
		////////// fix linked lists //////////
		// for each face
		for (p=flh;p!=NULL;p=p->next) {
			f=(Face*)p->data;
			verts[f->v1]->found=true;
			verts[f->v2]->found=true;
			verts[f->v3]->found=true;
		}
		if(flag){writeData((char*)q->data,vlh,flh,print_flag);}
		else {printData(vlh,flh,print_flag);}
		cleanup(flh,vlh);
		vlh=flh=NULL;
		delete[] verts;
	}
	return 0;
}
