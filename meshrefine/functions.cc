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

u4 computeHashValue(int key1, int key2) {
    char cval[10];
    u4 result;
    sprintf(cval,"%i",key1);
    result = hash((u1*)cval,(u4)strlen(cval),(u4)key2);
    return result;
}

// #####################################################
// #####################################################

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
            v = new Vertex(str,vl);
            vl->data = (void*)v;
            vlh = vl;
        }
        // if first character is F for Face, add new linked list class instance
        else if (strchr("F",*str)!=NULL){
            fl = new void_list();
            fl->next = flh;
            f = new Face(str,fl);
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

// #####################################################
// #####################################################

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

// #####################################################
// #####################################################

void addPointersToFaces(void_list *flh,hashtable_v &hmv){
    void_list *q;
    Face *f;
    int i=0;
    // for each face
    for (q=flh;q!=NULL;q=q->next) {
        f=(Face*)q->data;
        f->v1=hmv[(int)f->v1];
        f->v2=hmv[(int)f->v2];
        f->v3=hmv[(int)f->v3];
    }
}

void buildVertMap(void_list *vlh,hashtable_v &hmv){
    void_list *p;
    Vertex *v;
    // for each vertex
    for (p=vlh;p!=NULL;p=p->next){
        v=(Vertex*)p->data;
        hmv[v->index]=v;
    }
}

// #####################################################
// #####################################################

u4 keyPair(int a,int b){
    if (a<b){ return computeHashValue(a,b);}
    else {return computeHashValue(b,a);}
}

bool edgeMatch(Edge *e,int va,int vb) {
    if ( (e->v1->index==va && e->v2->index==vb) ||
        (e->v1->index==vb && e->v2->index==va) ){return true;}
    else {return false;}
}

void checkEdge(hashtable_t &hme,Face *f,Vertex *va,Vertex *vb,std::vector<Edge*> &e) {
    char s[32];
    Edge *ee=NULL,*et=NULL,*en=NULL;
    void_list *p;
    ///// find edge /////
    // if element exists given key
    if (hme.count(keyPair(va->index,vb->index))>0){
        ///// get Edge pointer /////
        // for each link in list
        for(p=hme[keyPair(va->index,vb->index)];p!=NULL;p=p->next){
            // does edge match
            et=(Edge*)p->data;
            if (edgeMatch(et,va->index,vb->index)){ee=et;}
        }
    }
    ///// if exists, then update /////
    if(ee!=NULL){ ee->update(f); }
    else { ///// edge pointer not found /////
        // new edge
        en = new Edge(f,va,vb);
        // add edge to vector
        e.push_back(en);
        // create new link in void_list
        p = new void_list();
        p->data=(void*)en;
        p->previous=NULL;
        // hash_map empty
        if (et==NULL){
            p->next=NULL;
        } else { // hash_map not empty 
            p->next=hme[keyPair(va->index,vb->index)];
            hme[keyPair(va->index,vb->index)]->previous=p;
        }
        // store edge pointer in hash table
        hme[keyPair(va->index,vb->index)]=p;
    }
}

void getEdges(void_list *flh,hashtable_t &hme,std::vector<Edge*> &e){
    void_list *q;
    Face *f;
    int i=0;
    // for each face
    for (q=flh;q!=NULL;q=q->next) {
        f=(Face*)q->data;
        checkEdge(hme,f,f->v1,f->v2,e);
        checkEdge(hme,f,f->v2,f->v3,e);
        checkEdge(hme,f,f->v3,f->v1,e);
    }
}
// #####################################################
// #####################################################

void_list* removeLink(void_list* L) {
    void_list *q;
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
// #####################################################
// #####################################################

void clearVertexAdjacencies(void_list *vlh){
    Vertex *v;
    // for each vertex
    for (void_list *p=vlh;p!=NULL;p=p->next) {
        v=(Vertex*)p->data;
        v->f.clear();
        v->v.clear();
        v->e.clear();
    }
}

void findVertexAdjacencies(void_list *flh,std::vector<Edge*> &e){
    Face *ff;
    std::vector<Edge*>::iterator i;
    // for each face, add face* to each face vertex
    for (void_list *p=flh;p!=NULL;p=p->next) {
        ff=(Face*)p->data;
        ff->v1->f.push_back(ff);
        ff->v2->f.push_back(ff);
        ff->v3->f.push_back(ff);
    }
    // for each edge
    for(i=e.begin();i!=e.end();i++){
        // add each vertex to other as adjacent
        (*i)->v1->v.push_back((*i)->v2);
        (*i)->v2->v.push_back((*i)->v1);
        // add edge* to each vertex
        (*i)->v1->e.push_back(*i);
        (*i)->v2->e.push_back(*i);
    }
}

// #####################################################
// #####################################################

void clearEdges(std::vector<Edge*> &e,hashtable_t &hme){
    // for each edge, delete pointer
    std::vector<Edge*>::iterator i;
    for(i=e.begin();i!=e.end();i++){delete (*i);}
    e.clear();
    // for each void_list, delete links
    void_list *p,*q;
    __gnu_cxx::hash_map<u4, void_list*,pair_hash,equ4,std::allocator<void_list*> >::iterator j;
//  ht_iterator j;
    for(j=hme.begin();j!=hme.end();j++){
        p=(*j).second;
        while (p!=NULL) {
            q=p->next;
            delete p;
            p=q;
        }
    }
    hme.clear();
}


void clearFaceEdges(void_list *flh){
    void_list *q;
    // for each face
    for (q=flh;q!=NULL;q=q->next) {((Face*)q->data)->clearEdges();}
}

int threshold(Edge *e,double t){return(e->l>t);}

bool thresholdEdges(void_list *flh,double t){
	void_list *p,*q;
	Face *f;
	bool flag = false;
	// for each face
	for(p=flh;p!=NULL;p=p->next){
		f=(Face*)p->data;
		// threshold criteria
		if (threshold(f->e1,t)||threshold(f->e2,t)||threshold(f->e3,t)) {
			// mark face as subdivided
			f->subdivided=true;
			// mark edges as bisected
			f->e1->bisected=true;
			f->e2->bisected=true;
			f->e3->bisected=true;
			// set flag
			flag = true;
		}
	}
	return flag;
}

bool findMoreSubdivided(void_list *flh){
	void_list *p,*q;
	Face *f;
	bool flag = false;
	// for each face
	for(p=flh;p!=NULL;p=p->next){
		f=(Face*)p->data;
		// if not subdivided
		if (!f->subdivided){
			// if two or more edges are bisected
			if ( (f->e1->bisected && f->e2->bisected) ||
				(f->e2->bisected && f->e3->bisected) ||
				(f->e1->bisected && f->e3->bisected)) {
				// mark face as subdivided
				f->subdivided=true;
				// mark edges as bisected
				f->e1->bisected=true;
				f->e2->bisected=true;
				f->e3->bisected=true;
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
	// for each face
	for(p=flh;p!=NULL;p=p->next){
		f=(Face*)p->data;
		// if not subdivided
		if (!f->subdivided){
			// if exactly one edge is bisected
			if ( (f->e1->bisected && !f->e2->bisected && !f->e3->bisected) ||
				(!f->e1->bisected && f->e2->bisected && !f->e3->bisected) ||
				(!f->e1->bisected && !f->e2->bisected && f->e3->bisected)) {
				// mark face as subdivided
				f->bisubdivided=true;
			}
		}
	}
}

void createNewVertices(std::vector<Edge*> &e,void_list *&vlh,int &max_verts){
	std::vector<Edge*>::iterator i;
	void_list *p;
	double x,y,z;
	Vertex *v;
	// for each edge
	for (i=e.begin();i!=e.end();i++){
		// if edge is bisected
		if((*i)->bisected){
			// compute new vertex
			x=((*i)->v1->x+(*i)->v2->x)/2.0;
			y=((*i)->v1->y+(*i)->v2->y)/2.0;
			z=((*i)->v1->z+(*i)->v2->z)/2.0;
			// create new vertex
			p = new void_list();
    		vlh->previous=p;
			p->next = vlh;
    		p->previous = NULL;
			v = new Vertex(++max_verts,x,y,z,p);
			p->data = (void*)v;
			vlh = p;
			// save new vertex pointer to edge
			(*i)->v=v;
		}
	}
}

void matchEdges(Vertex* new_verts[3],Face *f){
	void_list *q;
	// for each face edge,assign vert index
	new_verts[0]=f->e1->v;
	new_verts[1]=f->e2->v;
	new_verts[2]=f->e3->v;
}

void addFace(void_list *&flh,int num_faces,Vertex *va,Vertex *vb,Vertex *vc){
	void_list *p = new void_list();
    flh->previous=p;
	p->next = flh;
	p->previous = NULL;
	Face *f = new Face(num_faces,va,vb,vc,p);
	p->data = (void*)f;
	flh = p;
}

void createNewSubdividedFaces(void_list *&flh,int &max_faces){
	void_list *p;
	Face *f;
	Vertex* new_verts[3];
	// for each face
	p=flh;
	while(p!=NULL){
		f=(Face*)p->data;
		// if face is subdivided
		if(f->subdivided){
			// identify new vertices
			matchEdges(new_verts,f);
			// add new faces
			addFace(flh,++max_faces,f->v1,new_verts[0],new_verts[2]);
			addFace(flh,++max_faces,new_verts[0],f->v2,new_verts[1]);
			addFace(flh,++max_faces,new_verts[0],new_verts[1],new_verts[2]);
			addFace(flh,++max_faces,new_verts[2],new_verts[1],f->v3);
			// delete current face link
			delete f;
			// adjust list pointer
			if (p==flh) {flh=p->next;}
			// remove current face from face list
			p=removeLink(p);
		} else if (f->bisubdivided){
			// identify new vertices
			matchEdges(new_verts,f);
			// add new faces
			if(new_verts[0]!=NULL){
				addFace(flh,++max_faces,f->v1,new_verts[0],f->v3);
				addFace(flh,++max_faces,f->v3,new_verts[0],f->v2);
			} else if (new_verts[1]!=NULL) {
				addFace(flh,++max_faces,f->v2,new_verts[1],f->v1);
				addFace(flh,++max_faces,f->v1,new_verts[1],f->v3);
			} else if (new_verts[2]!=NULL) {
				addFace(flh,++max_faces,f->v3,new_verts[2],f->v2);
				addFace(flh,++max_faces,f->v2,new_verts[2],f->v1);
			}
			// delete current face link
			delete f;
			// adjust list pointer
			if (p==flh) {flh=p->next;}
			// remove current face from face list
			p=removeLink(p);
		} else {p=p->next;}
	}
}

bool refineMesh(std::vector<Edge*> &e,void_list *&flh,void_list *&vlh,double t,int &mf,int &mv){
	bool flag1,flag2 = true;
	flag1 = thresholdEdges(flh,t);
	while (flag2) {flag2=findMoreSubdivided(flh);}
	findBisubdivided(flh);
	createNewVertices(e,vlh,mv);
	createNewSubdividedFaces(flh,mf);
	return flag1;
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

void printMesh(void_list *vlh,void_list *flh){
	void_list *vend,*fend,*q;
	Vertex *v;
	Face *f;
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
		printf("Face %i  %i %i %i\n",f->index,f->v1->index,f->v2->index,f->v3->index);
	}
}

