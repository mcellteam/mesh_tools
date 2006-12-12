// #####################################################
// #####################################################

typedef  unsigned long int  u4;   /* unsigned 4-byte type */
typedef  unsigned     char  u1;   /* unsigned 1-byte type */

/* The mixing step */
#define mix(a,b,c) \
{ \
  a=a-b;  a=a-c;  a=a^(c>>13); \
  b=b-c;  b=b-a;  b=b^(a<<8);  \
  c=c-a;  c=c-b;  c=c^(b>>13); \
  a=a-b;  a=a-c;  a=a^(c>>12); \
  b=b-c;  b=b-a;  b=b^(a<<16); \
  c=c-a;  c=c-b;  c=c^(b>>5);  \
  a=a-b;  a=a-c;  a=a^(c>>3);  \
  b=b-c;  b=b-a;  b=b^(a<<10); \
  c=c-a;  c=c-b;  c=c^(b>>15); \
}

/* The whole new hash function */
u4 hash( register u1* k, u4 length, u4 initval)
//register u1 *k;        /* the key */
//u4           length;   /* the length of the key in bytes */
//u4           initval;  /* the previous hash, or an arbitrary value */
{
   register u4 a,b,c;  /* the internal state */
   u4          len;    /* how many key bytes still need mixing */

   /* Set up the internal state */
   len = length;
   a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
   c = initval;         /* variable initialization of internal state */

   /*---------------------------------------- handle most of the key */
   while (len >= 12)
   {
      a=a+(k[0]+((u4)k[1]<<8)+((u4)k[2]<<16) +((u4)k[3]<<24));
      b=b+(k[4]+((u4)k[5]<<8)+((u4)k[6]<<16) +((u4)k[7]<<24));
      c=c+(k[8]+((u4)k[9]<<8)+((u4)k[10]<<16)+((u4)k[11]<<24));
      mix(a,b,c);
      k = k+12; len = len-12;
   }

   /*------------------------------------- handle the last 11 bytes */
   c = c+length;
   switch(len)              /* all the case statements fall through */
   {
   case 11: c=c+((u4)k[10]<<24);
   case 10: c=c+((u4)k[9]<<16);
   case 9 : c=c+((u4)k[8]<<8);
      /* the first byte of c is reserved for the length */
   case 8 : b=b+((u4)k[7]<<24);
   case 7 : b=b+((u4)k[6]<<16);
   case 6 : b=b+((u4)k[5]<<8);
   case 5 : b=b+k[4];
   case 4 : a=a+((u4)k[3]<<24);
   case 3 : a=a+((u4)k[2]<<16);
   case 2 : a=a+((u4)k[1]<<8);
   case 1 : a=a+k[0];
     /* case 0: nothing left to add */
   }
   mix(a,b,c);
   /*-------------------------------------------- report the result */
   return c;
}

// #####################################################
// #####################################################

u4 computeHashValue(int key) {
	char cval[10];
	u4 result;
	sprintf(cval,"%i",key);
	result = hash((u1*)cval,(u4)strlen(cval),0xa205b064);
	return result;
}

int getEdgeHashVal(int va,int vb,int mask){
	u4 hashval1 = computeHashValue(va);
	u4 hashval2 = computeHashValue(vb);
	u4 result = hashval1+hashval2;
	return (int)(result&(mask-1));
}

// #####################################################
// #####################################################

Edge* findEdge(EdgeBlock *eb,int va,int vb){
	EdgeBlock *q;
	void_list *p;
	bool found=false;
	int hashval;

	// for each EdgeBlock
	q=eb;
	while(q!=NULL && !found){
		// compute hashval of va,vb
		hashval = getEdgeHashVal(va,vb,q->ht->s);
//		printf("EdgeBlock: va %i,vb %i,hashval %i\n",va,vb,hashval);
		// for each edge in hashtable element pointed to by hashval
		p = eb->ht->t[hashval];
		while (p!=NULL && !found){
			// if v1,v2
			if  ((((Edge*)p->data)->v1==va||((Edge*)p->data)->v2==va)&& 
				(((Edge*)p->data)->v1==vb||((Edge*)p->data)->v2==vb)){found=true;}
			else	{p=p->next;}
		}if(!found){q=q->next;}
	}
	if(found){return (Edge*)p->data;}
	else{return NULL;}
}

void updateEdge(Face *F,Edge *e,int va,int vb){
	void_list *qq;
	//add face to edge
	qq = new void_list();
	qq->next = e->f;
	qq->data = (void*)F;
	e->f = qq;
	if  (e->v1==va && e->v2==vb) {e->c12+=1;}
	else {e->c21+=1;}
	// add edge to face
	qq = new void_list();
	qq->next = F->e;
	qq->data = (void*)e;
	F->e = qq;
	
}

void createEdge(Face *F,EdgeBlock *eb,int va,int vb,Vertex **vert_array){
	EdgeBlock *q;
	Edge *e;
	bool found=false;
	void_list *p;
	Vertex *pa,*pb;
	double vax,vay,vaz,vbx,vby,vbz;

	// for each EdgeBlock
	q=eb;
	while(q!=NULL && !found){
		// is there an open Edge slot
		if (q->p<(q->n-1)){
			// if yes grab Edge pointer
			e=&q->e[q->p];
			q->p++;
			found = true;
		} else {q=q->next;}
	}

	// if no open Edge slot found in any EdgeBlock
	if (!found) {
		// create a new EdgeBlock
		q = new EdgeBlock(10000);
		// store pointer to last EdgeBlock
		q->next=eb;
		// create HashTable
		q->ht = new HashTable(q->n);
		// grab Edge pointer
		e=&q->e[q->p];
		q->p++;
	}

	// load Edge
	////////////////////
	e->v1 = va;
	e->v2 = vb;
	e->c12 = 1;
	e->c21 = 0;

	// get vertex info
	pa = vert_array[va];
	pb = vert_array[vb];
	vax = pa->x;
	vay = pa->y;
	vaz = pa->z;
	vbx = pb->x;
	vby = pb->y;
	vbz = pb->z;

	e->l = (vax-vbx)*(vax-vbx)+(vay-vby)*(vay-vby)+(vaz-vbz)*(vaz-vbz);

	// add face to edge
	// set link data to point to Face class instance
	e->f = new void_list();
	(e->f)->next = NULL;
	(e->f)->data = (void*)F;
	// add edge to face
	p = new void_list();
	p->next = F->e;
	p->data = (void*)e;
	F->e = p;
	// store edge pointer in hash table
	// compute hashval of va,vb
	int hashval = getEdgeHashVal(va,vb,q->ht->s);
	// load hash table
	p=new void_list();
	p->next = q->ht->t[hashval];
	p->data= (void*)e;
	q->ht->t[hashval] = p;
}

void checkEdge(Face *F,EdgeBlock *&eb,int va,int vb,Vertex **vert_array) {
	Edge *e;
	if(e=findEdge(eb,va,vb)){updateEdge(F,e,va,vb);}
	else {createEdge(F,eb,va,vb,vert_array);}
}

void getEdges(void_list *flh,EdgeBlock *eb,int num,Vertex **vert_array){
	void_list *q;
	Face *f;
	int i=0;
	// for each face
	for (q=flh;q!=NULL;q=q->next) {
		f=(Face*)q->data;
//		printf("face %i of %i\n",i++,num);
		checkEdge(f,eb,f->v1,f->v2,vert_array);
		checkEdge(f,eb,f->v2,f->v3,vert_array);
		checkEdge(f,eb,f->v3,f->v1,vert_array);
	}
}

void buildVertArray(void_list *vlh,Vertex **&vert_array,int num_verts){
	void_list *p;
	vert_array = new Vertex*[num_verts+1];
	int i;
	for (i=0;i<num_verts+1;i++){vert_array[i]=NULL;}
	for (p=vlh;p!=NULL;p=p->next){vert_array[((Vertex*)p->data)->index]=(Vertex*)p->data;}
	if(0) {
		for (p=vlh;p!=NULL;p=p->next){
			printf("Vertex %i %.15g %.15g %.15g\n",
			((Vertex*)p->data)->index,
			((Vertex*)p->data)->x,
			((Vertex*)p->data)->y,
			((Vertex*)p->data)->z
			);
		}
	}
	if (0) {
		Vertex *v;
		for(i=0;i<num_verts+1;i++){
			v=vert_array[i];
			if (v!=NULL){
				printf("i %i Vertex %i %.15g %.15g %.15g\n",i,
				v->index,v->x,v->y,v->z);
			}
		}
	}
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

void_list * addFace(void_list *as,Face *F){
	void_list *p;
	p=new void_list();
	p->next=as;
	p->data=(void*)F;
	as=p;
	return as;
}

void_list * getAssociatedFaces(Edge *e,void_list *flh){
	void_list *p,*as;
	as=NULL;
	Face* f;
	printf("Associated Faces and Vertices\n");
	for(p=flh;p!=NULL;p=p->next){
		f=(Face*)p->data;
		if( (f->v1==e->v1||f->v1==e->v2)||
			(f->v2==e->v1||f->v2==e->v2)||
			(f->v3==e->v1||f->v3==e->v2)
			){
			printf("Face %i %i %i %i\n",f->index,f->v1,f->v2,f->v3);
			as=addFace(as,f);
		}
	}
	return as;
}

int compare (const void* a, const void* b ) {
  return ( *(int*)a - *(int*)b );
}

void unique(int *array,int count){
	int i;
	for(i=0;i<count-1;i++){
		if(array[i]==array[i+1]){array[i]=0;}
	}
}

void getAssociatedVertices(void_list *af,void_list *vlh){

	// gather unique list of vertices in faces
	// count number of vertices in faces
	int count=0;
	void_list *p;
	Face *f;
	for(p=af;p!=NULL;p=p->next){count+=3;}
	// declare array
	int array[count];
	int i=0;
	for(p=af;p!=NULL;p=p->next){
		f=(Face*)p->data;
		array[i++]=f->v1;
		array[i++]=f->v2;
		array[i++]=f->v3;
	}
	// sort array
    qsort(array,count,sizeof(int),compare);
	unique(array,count);
	bool flag;
	Vertex *v;
	for(i=0;i<count;i++){
		if(array[i]){
			flag=false;
			p=vlh;
			while(p!=NULL && !flag){
				v=(Vertex*)p->data;
				if(v->index==array[i]){
					printf("Vertex %i %.15g %.15g %.15g\n",v->index,v->x,v->y,v->z);
					flag=true;
				}else{p=p->next;}
			}
		}
	}
}

void identifyBadFaces(EdgeBlock *eb,void_list *flh,void_list *vlh){
	EdgeBlock *q;
	Edge *e;
	void_list *p,*af;
	Face *f;
	int i,count=0;
	
	// for each EdgeBlock
	for(q=eb;q!=NULL;q=q->next){
		// for each edge in block
		for(i=0;i<q->n;i++){
			e=&q->e[i];
			// if edge traversal is not exactly once in each direction
			if ((e->c12+e->c21)==1) {
				printf("\nError: Edge is open.");
				printf(" va %i vb %i,",e->v1,e->v2);
				f = (Face*)(e->f)->data;
				printf(" Face %i %i %i %i\n",f->index,f->v1,f->v2,f->v3);
				count++;
			}
			else if (((e->c12==2)&&(e->c21==0))||((e->c12==0)&&(e->c21==2))){
				printf("\nError: Flip candidates.");
				printf(" va %i vb %i, c12 %i, c21 %i\n",e->v1,e->v2,e->c12,e->c21);
				for (p=e->f;p!=NULL;p=p->next) {
					f=(Face*)p->data;
					printf("     Face %i %i %i %i\n",f->index,f->v1,f->v2,f->v3);
					count++;
				}
			}
			else if ((e->c12&&(e->c21>1))||((e->c12>1)&&e->c21)) {
				printf("\nError: Topological Problem.");
				printf(" va %i vb %i, c12 %i, c21 %i\n",e->v1,e->v2,e->c12,e->c21);
				for (p=e->f;p!=NULL;p=p->next) {
					f=(Face*)p->data;
					printf("     Face %i %i %i %i\n",f->index,f->v1,f->v2,f->v3);
					count++;
				}
				af=getAssociatedFaces(e,flh);
				getAssociatedVertices(af,vlh);
			}
		}
	}
	printf("\n\n # BAD FACES = %i\n\n",count);
}

int threshold(Edge *e,double t){return(e->l>t);}

bool thresholdEdges(void_list *flh,double t){
	void_list *p,*q;
	Face *f;
	Edge *e1,*e2,*e3;
	bool flag = false;
	// for each face
	for(p=flh;p!=NULL;p=p->next){
		f=(Face*)p->data;
		q=f->e;
		// grab edge pointers
		e1=(Edge*)q->data;
		q=q->next;
		e2=(Edge*)q->data;
		q=q->next;
		e3=(Edge*)q->data;
		// threshold criteria
		if (threshold(e1,t)||threshold(e2,t)||threshold(e3,t)) {
			// mark face as subdivided
			f->subdivided=true;
			// mark edges as bisected
			e1->bisected=true;
			e2->bisected=true;
			e3->bisected=true;
			// set flag
			flag = true;
		}
	}
	return flag;
}

bool findMoreSubdivided(void_list *flh){
	void_list *p,*q;
	Face *f;
	Edge *e1,*e2,*e3;
	bool flag = false;
	// for each face
	for(p=flh;p!=NULL;p=p->next){
		f=(Face*)p->data;
		// if not subdivided
		if (!f->subdivided){
			q=f->e;
			// grab edge pointers
			e1=(Edge*)q->data;
			q=q->next;
			e2=(Edge*)q->data;
			q=q->next;
			e3=(Edge*)q->data;
			// if two or more edges are bisected
			if ( (e1->bisected && e2->bisected) ||
				(e2->bisected && e3->bisected) ||
				(e1->bisected && e3->bisected)) {
				// mark face as subdivided
				f->subdivided=true;
				// mark edges as bisected
				e1->bisected=true;
				e2->bisected=true;
				e3->bisected=true;
				// set flag
				flag = true;
			}
		}
	}
	return flag;
}

void findBisubdivided(void_list *flh){
	void_list *p,*q;
	Face *f;
	Edge *e1,*e2,*e3;
	// for each face
	for(p=flh;p!=NULL;p=p->next){
		f=(Face*)p->data;
		// if not subdivided
		if (!f->subdivided){
			q=f->e;
			// grab edge pointers
			e1=(Edge*)q->data;
			q=q->next;
			e2=(Edge*)q->data;
			q=q->next;
			e3=(Edge*)q->data;
			// if two or more edges are bisected
			if ( (e1->bisected && !e2->bisected && !e3->bisected) ||
				(!e1->bisected && e2->bisected && !e3->bisected) ||
				(!e1->bisected && !e2->bisected && e3->bisected)) {
				// mark face as subdivided
				f->bisubdivided=true;
			}
		}
	}
}

void createNewVertices(EdgeBlock *eb,void_list *&vlh,Vertex **vert_array,int num_verts){
	EdgeBlock *q;
	int i;
	Edge *e;
	void_list *p;
	double x,y,z;
	Vertex *v,*v1,*v2;

	// for each EdgeBlock
	for (q=eb;q!=NULL;q=q->next){
		// for each edge in EdgeBlock
		for (i=0;i<q->n;i++){
			e=&q->e[i];
			// if edge is bisected
			if(e->bisected){
				// get edge vertices
				v1=vert_array[e->v1];
				v2=vert_array[e->v2];
				// compute new vertex
				x=(v1->x+v2->x)/2;
				y=(v1->y+v2->y)/2;
				z=(v1->z+v2->z)/2;
				// create new vertex
				p = new void_list();
				p->next = vlh;
				v = new Vertex(++num_verts,x,y,z);
				p->data = (void*)v;
				vlh = p;
				// save new vertex pointer to edge
				e->v=v;
			}
		}
	}
}

int match(int i,Face *f,Edge *e){
	int va,vb,v1,v2;
	// designate intrinsic vertex indices
	if(i==0){
		va=f->v1;
		vb=f->v2;
	} else if(i==1){
		va=f->v2;
		vb=f->v3;
	} else if(i==2){
		va=f->v3;
		vb=f->v1;
	}
	// designate edge vertex indices
	v1=e->v1;
	v2=e->v2;
	// compare indices
	if((va==v1&&vb==v2)||(va==v2&&vb==v1)){return 1;}
	else {return 0;}
}

void matchEdges(int *&new_verts,Face *f){
	int i;
	void_list *q;
	Edge *e;
	// for each intrinsic edge
	for(i=0;i<3;i++){
		// for each extrinsic edge
		for(q=f->e;q!=NULL;q=q->next){
			e=(Edge*)q->data;
			// if edges match
			if (match(i,f,e)){
				if (e->bisected){
					// assign vert index
					new_verts[i]=((Vertex*)e->v)->index;
				} else {printf("edge is not bisected!\n");}
			}
		}
	}
}

void addFace(void_list *&flh,int &num_faces,int va,int vb,int vc){
	Face *f;
	void_list *p;
	p = new void_list();
	p->next = flh;
	p->previous = NULL;
	f = new Face(++num_faces,va,vb,vc);
	p->data = (void*)f;
	flh = p;
	// adjust previous
	p=p->next;
	p->previous=flh;
}

void_list* removeLink(void_list* L) {
    void_list* q;
    // and remove face from candidate face list
    if (L->previous!=NULL) {
        if (L->next!=NULL) {
            // if both previous and next exist
            (L->previous)->next = L->next;
            (L->next)->previous = L->previous;
        } else {
            // if previous exists and next does not
            (L->previous)->next = NULL;
        }
    } else {
        if (L->next!=NULL) {
            // if previous does not exist and next does
            (L->next)->previous = NULL;
        } // else { // if neither previous nor next exists }
    }
    // update pointer
    q=L->next;
    delete L;
    return q;
}


void createNewSubdividedFaces(void_list *&flh,int num_faces){
	void_list *p;
	Face *f;
	int *new_verts;
	new_verts = new int[3];

	// for each face
	p=flh;
	while(p!=NULL){
		f=(Face*)p->data;
		// if face is subdivided
		if(f->subdivided){
			// identify new vertices
			matchEdges(new_verts,f);
			// add new faces
			addFace(flh,num_faces,f->v1,new_verts[0],new_verts[2]);
			addFace(flh,num_faces,new_verts[0],f->v2,new_verts[1]);
			addFace(flh,num_faces,new_verts[0],new_verts[1],new_verts[2]);
			addFace(flh,num_faces,new_verts[2],new_verts[1],f->v3);
			// delete current face link
			delete f;
			// adjust list pointer
			if (p==flh) {flh=p->next;}
			// remove current face from face list
			p=removeLink(p);
		} else {p=p->next;}
	}
	delete[] new_verts;
}

bool refineMesh(EdgeBlock *eb,void_list *&flh,void_list *&vlh,
				double threshold,Vertex **vert_array,int num_faces,int num_verts){
	bool flag1,flag2 = true;
	flag1 = thresholdEdges(flh,threshold);
	while (flag2) {flag2=findMoreSubdivided(flh);}
	findBisubdivided(flh);
	createNewVertices(eb,vlh,vert_array,num_verts);
	createNewSubdividedFaces(flh,num_faces);
	return flag1;
}

void cleanup1(EdgeBlock *eb,Vertex **vert_array){
	EdgeBlock *pp,*qq;
	// vert_array
	delete[] vert_array;
	// EdgeBlocks
	// for each EdgeBlock
	pp=eb;
	while (pp!=NULL) {
		qq=pp->next;
		delete pp;
		pp=qq;
	}
}

void cleanup2(void_list *f,void_list *v){
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

void resetFaces(void_list *flh){
	void_list *p,*q,*pp;
	Face *f;
	// for each face
	for(pp=flh;pp!=NULL;pp=pp->next){
		f=(Face*)pp->data;
		// delete edge list
		p=f->e;
		while (p!=NULL) {
			q=p->next;
			delete p;
			p=q;
		}
		// reset values
		f->subdivided = false;
		f->bisubdivided = false;
		f->e=NULL;
	}
}

void printMesh(void_list *vlh,void_list *flh){
	void_list *vend,*fend,*q;
	Vertex *v;
	Face *f;

	flh=addPrevious(flh);
	vlh=addPrevious(vlh);

    ////////// add pointers to ends of linked lists /////////
    // read out vertex linked list
    for (q=vlh;q!=NULL;q=q->next) { vend = q; if (q->next==NULL) break; }

    // read out face linked list
    for (q=flh;q!=NULL;q=q->next) { fend = q; if (q->next==NULL) break; }

    ////////// write data to stdout /////////
	// write out vertex linked list
	for (q=vend;q!=NULL;q=q->previous) {
		v = (Vertex*)q->data;
		printf("Vertex %i  %.15g %.15g %.15g\n",v->index,v->x,v->y,v->z);
	}
	// write out face linked list
	for (q=fend;q!=NULL;q=q->previous) {
		f=(Face*)q->data;
		printf("Face %i  %i %i %i\n",f->index,f->v1,f->v2,f->v3);
	}
}
