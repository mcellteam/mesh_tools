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

int compare (const void* a, const void* b ) {
  return ( *(int*)a - *(int*)b );
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

int vertexInFace(Face* F,Vertex *V){
	return ((F->v1==V->index)||(F->v2==V->index)||(F->v3==V->index));
}

void_list * removeBadFaces(void_list *F,void_list *&B) {
	void_list *p,*q,*ptr1,*ptr2;
	bool flag;
	ptr1 = F;
	ptr2 = B;
	// for each face
	p=F;
	while (p!=NULL) {
		// search bad face list
		flag = true;
		q=ptr2;
		while(q!=NULL && flag) {
			// if face is bad
    	  	if ( ((Face*)p->data)->index == ((Face*)q->data)->index ) {
				// delete face link
				delete (Face*)p->data;
				// adjust list pointer
				if (p==ptr1) { ptr1 = p->next; }
				// remove face from face list
				p=removeLink(p);
				// adjust list pointer
				if (q==ptr2) { ptr2 = q->next; }
				// remove face from face list
				q=removeLink(q);
				flag = false;
			}else{ q=q->next;}
		}if(flag){ p=p->next;}
	}
	F=ptr1;
	B=ptr2;
	return F;
}

void_list * removeBadVertices(void_list *V,void_list *&B) {
	void_list *p,*q, *ptr1,*ptr2;
	bool flag;
	ptr1 = V;
	ptr2 = B;
	// for each vertex
	p=V;
	while (p!=NULL) {
		// search bad vertex list
		flag = true;
		q=ptr2;
		while(q!=NULL && flag) {
			// if vertex is bad
    	  	if ( ((Vertex*)p->data)->index == ((Vertex*)q->data)->index ) {
				// delete vertex link
				delete (Vertex*)p->data;
				// adjust list pointer
				if (p==ptr1) { ptr1 = p->next; }
				// remove vertex from vertex list
				p = removeLink(p);
				// adjust list pointer
				if (q==ptr2) { ptr2 = q->next; }
				// remove bad vertex from bad vertex list
				q = removeLink(q);
				flag = false;
			}else { q=q->next; }
		} if(flag){ p=p->next;}
	}
	V=ptr1;
	B=ptr2;
	return V;
}

void_list* matchZ(void_list *L,int z_value){
	// build linked list of first vertex elements with z_value matching target
	void_list *q,*t,*th;
	th = NULL;
	for (q=L;q!=NULL;q=q->next) {
      	if (((Vertex*)q->data)->z == z_value) {
			t = new void_list();
			t->next = th;
			t->data = (Vertex*)q->data;
			th = t;
		}
		if (q->next==NULL) break;
	}
	return th;
}

void_list* candidateBadVertices(void_list *V, int val){
	// gather candidate bad vertices from first vertex linked list
	void_list *q,*t,*th;
	th = NULL;
	// for each vertex
	for (q=V;q!=NULL;q=q->next) {
		// if z value matches target
		if (((Vertex*)q->data)->z == val) {
			// store vertex 
			t = new void_list();
			t->next = th;
			t->data = (Vertex*)q->data;
			th = t;
		}
		if (q->next==NULL) break;
	}
	return th;
}

void_list* candidateBadFaces(void_list *F,void_list *V){
	///// find faces that have at least one vertex from cbv lists /////
	void_list *q,*p,*t,*th;
	th = NULL;
	bool flag;
	// for each face from first face linked list
	for (q=F;q!=NULL;q=q->next) {
		// search cbv1 list 
		flag = true;
		p=V;
		while(p!=NULL && flag) {
			// if candidate bad vertex is at least one of three face vertices
    	  	if (vertexInFace((Face*)q->data,(Vertex*)p->data)) {
				// add face to candidate bad face list
				t = new void_list();
				t->next = th;
				t->data = (Face*)q->data;
				th = t;
				flag = false;
			}
			if (p->next==NULL) break;
			p=p->next;
		}
		if (q->next==NULL) break;
	}
	return th;
}

void_list* findBadFaces(void_list *&CBF,void_list* V,void_list* BF) {
	// scan candidate bad faces searching for known bad vertices
	void_list *q,*p,*t,*th,*ptr;
	ptr=CBF;
	th = BF;
	bool flag;
	// for each candidate bad face
	p=CBF;
	while (p!=NULL) {
		// search bad vertex list
		flag = true;
		q=V;
		while(q!=NULL && flag) {
			// if bad vertex is at least one of three face vertices
    	  	if (vertexInFace((Face*)p->data,(Vertex*)q->data)) {
				// add face to bad face list
				t = new void_list();
				t->next = th;
				t->data = (Face*)p->data;
				th = t;
				flag = false;
				// adjust list pointer
				if (p==ptr) { ptr = p->next; }
				// and remove face from candidate face list
				p=removeLink(p);
			} q=q->next;
		}
		if(flag){ p=p->next;}
	}
	CBF=ptr;
	return th;
}

void extraVertices(Group &g,void_list *S) {
	// find extra vertices
	void_list *p,*q,*pp;
	ExtraVertex *EV,*ev;
	bool found;
	int i;

	// for each group
	for(i=0;i<g.count;i++){
		// for each vertex in group
		for(q=g.c[i].verts;q!=NULL;q=q->next){
			found=false;
			// for each shared vertex
			p = S;
			while(p!=NULL && !found) {
				// if group vertex is shared
				if (((Vertex*)q->data)->index==((Vertex*)p->data)->index){
					found=true;
				}
				p=p->next;
			}
			// if group vertex is not shared
			if(!found){
				// then it is extra vertex
				EV = new ExtraVertex();
				EV->v = (Vertex*)q->data;
				pp = new void_list();
				pp->next = g.c[i].extra;
				pp->data = EV;
				g.c[i].extra = pp;
			}
		}
	}
}

void groupExtraVertices(Group &g){
	// group extra vertices in each contour into sites
	bool flag;
	int i,j,num_ev,next_group,g1,g2,g3,ga,gb,smallest,largest,*array,size,num;
	ExtraVertex *ev1,*ev2,*ev3,*eva,*evb;
	void_list *q,*p,*qq,*pp;
	///// handle sites with two or more extra vertices
	//for each contour group
	for(i=0;i<g.count;i++){
		flag=true;
		next_group=0;
		// while changes are made
		while(flag){
			flag=false;
			// for each face in contour group
			for(q=g.c[i].faces;q!=NULL;q=q->next){
				g1=g2=g3=ga=gb=0;
				// for each extra vertex in contour group
				for(p=g.c[i].extra;p!=NULL;p=p->next){
					// are any vertices extra?
					if( ((Face*)q->data)->v1 ==
						((Vertex*)((ExtraVertex*)p->data)->v)->index)
						{ev1=(ExtraVertex*)p->data;g1=ev1->g;}
					else if( ((Face*)q->data)->v2 ==
						((Vertex*)((ExtraVertex*)p->data)->v)->index)
						{ev2=(ExtraVertex*)p->data;g2=ev2->g;}
					else if( ((Face*)q->data)->v3 ==
						((Vertex*)((ExtraVertex*)p->data)->v)->index)
						{ev3=(ExtraVertex*)p->data;g3=ev3->g;}
				}
 				// were at least two vertices extra?
				if ((g1&&g2)||(g1&&g3)||(g2&&g3)){
					if(g1){ga=g1;eva=ev1;}
					if(g2){ if(!ga){ga=g2;eva=ev2;} else{gb=g2;evb=ev2;} }
					if(g3){gb=g3;evb=ev3;}
					// check to see either extra vertex is in a group
					if(ga<0 && gb<0){
						// if both ungrouped, then add both to same group
						eva->g=next_group;
						evb->g=next_group;
						next_group++;
						flag=true;
					}else if (ga==gb){
						// if same and greater than zero, then do nothing
					}else if (ga<0 && gb>0){
						eva->g=evb->g;
						flag=true;
					}else if (ga>0 && gb<0){
						evb->g=eva->g;
						flag=true;
					}else {
						flag=true;
						// if different and greater than zero
						// then identify smallest group #
						if (ga<gb) {smallest=ga;largest=gb;}
						else {smallest=gb;largest=ga;}
						// then set any extra vertex with largest group # to smallest group #
						for(p=g.c[i].extra;p!=NULL;p=p->next){
							if(((ExtraVertex*)p->data)->g==largest){((ExtraVertex*)p->data)->g=smallest;}
						}
					}
				}
			}
		}
		///// handle sites with only a single extra vertex
		flag=true;
		// while changes are made
		while(flag){
			flag=false;
			// for each face in contour group
			for(q=g.c[i].faces;q!=NULL;q=q->next){
				g1=g2=g3=ga=gb=0;
				// for each extra vertex in contour group
				for(p=g.c[i].extra;p!=NULL;p=p->next){
					// are any vertices extra?
					if( ((Face*)q->data)->v1 ==
						((Vertex*)((ExtraVertex*)p->data)->v)->index)
						{ev1=(ExtraVertex*)p->data;g1=ev1->g;}
					else if( ((Face*)q->data)->v2 ==
						((Vertex*)((ExtraVertex*)p->data)->v)->index)
						{ev2=(ExtraVertex*)p->data;g2=ev2->g;}
					else if( ((Face*)q->data)->v3 ==
						((Vertex*)((ExtraVertex*)p->data)->v)->index)
						{ev3=(ExtraVertex*)p->data;g3=ev3->g;}
				}
 				// were any vertices extra?
				if (g1||g2||g3){
					// check to see if extra vertex is in a group
					if(g1<0){
						// if both ungrouped, then add both to same group
						ev1->g=next_group;
						next_group++;
						flag=true;
					} else if(g2<0){
						// if both ungrouped, then add both to same group
						ev2->g=next_group;
						next_group++;
						flag=true;
					} else if(g3<0){
						// if both ungrouped, then add both to same group
						ev3->g=next_group;
						next_group++;
						flag=true;
					}
				}
			}
		}
	}

	///// count number of sites per contour
	// for each contour group
	for(i=0;i<g.count;i++){
		size=0;	
		num=0;
		// for each vertex in contour
		for(q=g.c[i].verts;q!=NULL;q=q->next){size++;}
		// create int array of size = number of vertices per contour+1
		array = new int[size+1];
		for(j=0;j<size+1;j++){
			array[j]=0;
		}
		// for each extra vertex in contour group
		for(p=g.c[i].extra;p!=NULL;p=p->next){
			// record group number in array
			array[((ExtraVertex*)p->data)->g]=1;
		}
		for(j=0;j<size+1;j++){
			if(array[j]){num++;}
		}
		// record number of sites in contour
		g.c[i].num=num;
		// cleanup
		delete[] array;
	}
}

void initializeSites(Group &g){
	int i,j,k;
	bool found;
	void_list *p,*q,*pp;
	ExtraVertex *ev;
	// for each contour
	for(i=0;i<g.count;i++){
		// initialize sites
		g.c[i].s = new Site[g.c[i].num];
	}
	// for each contour
	for(i=0;i<g.count;i++){
		// for each site in contour
		for(j=0;j<g.c[i].num;j++){
			// for each extra vertex
			for(p=g.c[i].extra;p!=NULL;p=p->next){
				ev=(ExtraVertex*)p->data;
				// if extra vertex group number matches site group
				if (ev->g==j) {g.c[i].s[j].n++;}
			}
			// allocate space for ExtraVertex*
			g.c[i].s[j].ev = new ExtraVertex*[g.c[i].s[j].n];
			// initialize ExtraVertex*
			for(k=0;k<g.c[i].s[j].n;k++){
				g.c[i].s[j].ev[k]=NULL;
			}
		}
	}
}

void orderSites(Group &g){
	int i,j,k;
	bool found;
	void_list *p,*q,*pp;
	ExtraVertex *ev,*sev;
	Vertex *target;
	// for each contour
	for(i=0;i<g.count;i++){
		// for each site in contour
		for(j=0;j<g.c[i].num;j++){
			///// identify site ev_0 and ev_n-1
			// for each extra vertex
			for(p=g.c[i].extra;p!=NULL;p=p->next){
				ev=(ExtraVertex*)p->data;
				found=false;
				// if extra vertex group number matches site group
				if (ev->g==j) {
					// if extra vertex laterals, la or lb, match site laterals l1 or l2
					// then assign extra vertex to site ev in appropriate place if ev==NULL
					if(((Vertex*)ev->la)->index==((Vertex*)g.c[i].s[j].l1)->index){
						if(g.c[i].s[j].ev[0]==NULL){
							g.c[i].s[j].ev[0]=ev;
							found=true;
							target=(Vertex*)ev->lb;
						}
					}
					if(!found && ((Vertex*)ev->la)->index==((Vertex*)g.c[i].s[j].l2)->index){
						if(g.c[i].s[j].ev[g.c[i].s[j].n-1]==NULL){
							g.c[i].s[j].ev[g.c[i].s[j].n-1]=ev;	
							found=true;
						}
					}
					if(!found && ((Vertex*)ev->lb)->index==((Vertex*)g.c[i].s[j].l1)->index){
						if(g.c[i].s[j].ev[0]==NULL){
							g.c[i].s[j].ev[0]=ev;
							found=true;
							target=(Vertex*)ev->la;
						}
					}
					if(!found && ((Vertex*)ev->lb)->index==((Vertex*)g.c[i].s[j].l2)->index){
						if(g.c[i].s[j].ev[g.c[i].s[j].n-1]==NULL){
							g.c[i].s[j].ev[g.c[i].s[j].n-1]=ev;	
							found=true;
						}
					}
				}
			}

			///// identify site ev_1 to ev_n-2
			// for site extra vertex slots ev_1 to ev_n-2
			for(k=1;k<g.c[i].s[j].n-1;k++){
				sev=g.c[i].s[j].ev[k];
				// for each extra vertex
				for(p=g.c[i].extra;p!=NULL;p=p->next){
					ev=(ExtraVertex*)p->data;
					// if ev->v matches target 
					if(((Vertex*)ev->v)->index==target->index){
						// record extra vertex in site
						g.c[i].s[j].ev[k]=ev;
						// set new target
						// if ev->la matches k-1
						if(((Vertex*)ev->la)->index==((Vertex*)((ExtraVertex*)g.c[i].s[j].ev[k-1])->v)->index){
							// then ev->lb is new target
							target=(Vertex*)ev->lb;
						}
						// else ev->lb is new target
						else {target=(Vertex*)ev->la;}
					}
				}
			}
		}
	}
}

void lateralVertices(Group &g,void_list *BF,void_list *S,int z_value,Vertex **vert_array) {
	// find lateral vertices
	void_list *p,*q,*pp;
	bool found;
	int i,j;
	Vertex *v1,*v2,*v3;
	ExtraVertex *ev;

	// group extra vertices in each contour into sites
	groupExtraVertices(g);

	// initialize sites
	initializeSites(g);

	///// load sites
	// for each contour group
	for(i=0;i<g.count;i++){
		// for each extra vertex in contour group
		for(p=g.c[i].extra;p!=NULL;p=p->next){
			ev=(ExtraVertex*)p->data;
			///// identify lateral vertices
			// for each group face
			for(q=g.c[i].faces;q!=NULL;q=q->next){
				v1=vert_array[((Face*)q->data)->v1];
				v2=vert_array[((Face*)q->data)->v2];
				v3=vert_array[((Face*)q->data)->v3];
				// if face contains extra vertex
				if(vertexInFace((Face*)q->data,(Vertex*)((ExtraVertex*)p->data)->v)){
					// if z_value matches and index is not same as extra vertex
					// remember need to account for 2 pairs of duplicate lateral vertices
					if((v1->z==z_value)&&(((Face*)q->data)->v1!=((Vertex*)ev->v)->index)){
						// this vertex is lateral
						if(ev->la==NULL){ev->la=v1;}
						else if (ev->lb==NULL&&(ev->la)->index!=((Face*)q->data)->v1){ev->lb=v1;}
					}
					if((v2->z==z_value)&&(((Face*)q->data)->v2!=(ev->v)->index)){
						// this vertex is lateral
						if(ev->la==NULL){ev->la=v2;}
						else if (ev->lb==NULL&&(ev->la)->index!=((Face*)q->data)->v2){ev->lb=v2;}
					}
					if((v3->z==z_value)&&(((Face*)q->data)->v3!=(ev->v)->index)){ 
						// this vertex is lateral
						if(ev->la==NULL){ev->la=v3;}
						else if (ev->lb==NULL&&(ev->la)->index!=((Face*)q->data)->v3){ev->lb=v3;}
					}
				}
			}

			// for each shared vertex
			found = true;
			q=S;
			while(q!=NULL && found) {
				// if lateral is shared
				if(((Vertex*)(ev->la))->index==((Vertex*)q->data)->index){found=false;}
				q=q->next;
			}
			// if lateral is shared
			if(!found){
				// then save extra vertex lateral as site lateral 1 or 2
				if(g.c[i].s[ev->g].l1==NULL)
				{g.c[i].s[ev->g].l1=(Vertex*)ev->la;}
				else {g.c[i].s[ev->g].l2=(Vertex*)ev->la;}
			}
			// for each shared vertex
			found = true;
			q=S;
			while(q!=NULL && found) {
				// if lateral is shared
				if(((Vertex*)(ev->lb))->index==((Vertex*)q->data)->index){found=false;}
				q=q->next;
			}
			// if lateral is shared
			if(!found){
				// then save extra vertex lateral as site lateral 1 or 2
				if(g.c[i].s[ev->g].l1==NULL){g.c[i].s[ev->g].l1=(Vertex*)ev->lb;}
				else {g.c[i].s[ev->g].l2=(Vertex*)ev->lb;}
			}
		}
	}

	// determine the sequence of extra vertices in sites
	orderSites(g);
}

void convertLateral(void_list *V,Group &g,double epsilon) {
	// convert lateral vertices to other object
	void_list *p;
	Vertex *o,*v1,*v2;
	bool f1,f2;
	int i,j;
	// for each contour
	for(i=0;i<g.count;i++){
		// for each site in contour
		for(j=0;j<g.c[i].num;j++){
			v1=g.c[i].s[j].l1;
			v2=g.c[i].s[j].l2;
			// for each vertex of other object
			f1=false;f2=false;
			p=V;
			while (p!=NULL && (!f1 || !f2)) {
				o=(Vertex*)p->data;
				// for both site lateral vertices
				if ( !distinguishable(o->x,v1->x,epsilon) &&
					!distinguishable(o->y,v1->y,epsilon) &&
					!distinguishable(o->z,v1->z,epsilon)
					){g.c[i].s[j].l1 = o;f1=true;}
				if ( !distinguishable(o->x,v2->x,epsilon) &&
					!distinguishable(o->y,v2->y,epsilon) &&
					!distinguishable(o->z,v2->z,epsilon)
					){g.c[i].s[j].l2= o;f2=true;}
				p=p->next;
			}
		}
	}
}

int faceContainsLaterals(Vertex *v1,Vertex *v2,Face *q,int val,Vertex **vert_array){
	return ((v1->index==q->v1||v1->index==q->v2||v1->index==q->v3)&&
			(v2->index==q->v1||v2->index==q->v2||v2->index==q->v3)
			&&( vert_array[q->v1]->z!=val && vert_array[q->v2]->z!=val && vert_array[q->v3]->z!=val)
			);
}

void gatherThirdDeleteFace(Group &g,void_list *&F,int val,Vertex **vert_array){
	void_list *p,*ptr;
	int i,j,orient;
	Vertex *v1,*v2;
	Face *f;
	bool found;
	ptr=F;
	// for each contour
	for(i=0;i<g.count;i++){
		// for each site in contour
		for(j=0;j<g.c[i].num;j++){
			v1=g.c[i].s[j].l1;
			v2=g.c[i].s[j].l2;
			// for each face in other object
			found = false;
			p=F;
			while(p!=NULL && !found){
				f=(Face*)p->data;
				// if face contains both lateral vertices and z != val
				if(faceContainsLaterals(v1,v2,f,val,vert_array)){
					found = true;
					// identify third vertex
					if((v1->index!=f->v1)&&(v2->index!=f->v1)){g.c[i].s[j].th=vert_array[f->v1];}
					else if((v1->index!=f->v2)&&(v2->index!=f->v2)){g.c[i].s[j].th=vert_array[f->v2];}
					else {g.c[i].s[j].th=vert_array[f->v3];}
					// identify orientation
					if(((v1->index==f->v1)&&(v2->index==f->v2))||((v1->index==f->v3)&&
					(v2->index==f->v1))||((v1->index==f->v2)&&(v2->index==f->v3))){g.c[i].s[j].orient= 1;}
					else {g.c[i].s[j].orient=0;}
					// delete face link
					delete f;
					// adjust list pointer
					if (p==ptr) { ptr = p->next; }
					// remove face
					p=removeLink(p);
				}else {p=p->next;}
			}
			// adjust pointer
			F=ptr;
		}
	}
}

void_list * addFaces(void_list *F,int v1,int v2, int th,int orient,int &MF) {
	// add faces
	void_list *n;
	Face *f;
	MF++;
	char buffer[128];
	sprintf (buffer, "Face 0 0 0 0");
	f = new Face(buffer);
	f->index = MF;
	if(orient) {f->v1=v1;f->v2=v2;}
	else {f->v1=v2;f->v2=v1;}
	f->v3 = th;
	n = new void_list();
	n->next = F;
	n->data = (void*)f;
	if(0){printf("Face %i %i %i %i\n",f->index,f->v1,f->v2,f->v3);}
	return n;
}

void_list * addVertex(void_list *V,Vertex *EV,int MV) {
	// add extra vertex
	void_list *t;
	Vertex *v;
	char buffer[128];
	sprintf (buffer, "Vertex 0 0 0 0");
	v = new Vertex(buffer);
	v->index = MV;
	v->x = EV->x;
	v->y = EV->y;
	v->z = EV->z;
	t = new void_list();
	t->next = V;
	t->data = (void*)v;
	if(0){printf("Vertex %i %.15g %.15g %.15g\n",v->index,v->x,v->y,v->z);}
	return t;
}


void addFacesVertices(Group &g,void_list *&F,void_list *&V,int max_vertex,int max_face){
	void_list *p,*ptr;
	int i,j,k,orient;
	Vertex *v1,*v2,*th;
	Face *f;
	bool found;
	ptr=F;

	// for each contour
	for(i=0;i<g.count;i++){
		// for each site in contour
		for(j=0;j<g.c[i].num;j++){
			th=g.c[i].s[j].th;
			// add face with l1, site ev0, and third
			orient=g.c[i].s[j].orient;
			v1=g.c[i].s[j].l1;
			v2=g.c[i].s[j].ev[0]->v;
			max_vertex++;
			F = addFaces(F,v1->index,max_vertex,th->index,orient,max_face);
			V = addVertex(V,v2,max_vertex);
			g.c[i].s[j].ev[0]->v=(Vertex*)V->data;
			// add face with site ev_n-1, l2, and third
			if (g.c[i].s[j].n>1){
				v1=g.c[i].s[j].ev[g.c[i].s[j].n-1]->v;
				v2=g.c[i].s[j].l2;
				max_vertex++;
				F = addFaces(F,max_vertex,v2->index,th->index,orient,max_face);
				V = addVertex(V,v1,max_vertex);
				g.c[i].s[j].ev[g.c[i].s[j].n-1]->v=(Vertex*)V->data;
			}else{
				v1=g.c[i].s[j].ev[g.c[i].s[j].n-1]->v;
				v2=g.c[i].s[j].l2;
				F = addFaces(F,max_vertex,v2->index,th->index,orient,max_face);
			}
			// add more faces
			if (g.c[i].s[j].n>1){
				for(k=0;k<g.c[i].s[j].n-2;k++){
					v1=g.c[i].s[j].ev[k]->v;
					v2=g.c[i].s[j].ev[k+1]->v;
					max_vertex++;
					F = addFaces(F,v1->index,max_vertex,th->index,orient,max_face);
					V = addVertex(V,v2,max_vertex);
					g.c[i].s[j].ev[k+1]->v=(Vertex*)V->data;
				}
				k=g.c[i].s[j].n-2;
				v1=g.c[i].s[j].ev[k]->v;
				v2=g.c[i].s[j].ev[k+1]->v;
				F = addFaces(F,v1->index,max_vertex,th->index,orient,max_face);
			}
		}
	}
}

void identifyBadVerticesAndFaces(void_list *&cbf,void_list *&cbv,void_list *v_shared,
								Group &g,void_list *&bfh,void_list *&bvh) {
	///// find bad vertices and faces from first vertex and first face linked lists
	void_list *ptr1,*ptr2,*p,*q,*pp;
	Face *temp;
	void_list *bf,*bv;
	bool flag;
	ptr1=cbf;
	ptr2=cbv;
	int i;
	// for each candidate bad face
	p=cbf;
	while (p!=NULL) {
		// search shared vertex list v_shared
		flag = true;
		q=v_shared;
		while(q!=NULL && flag) {
			// if shared vertex is at least one of three face vertices
    	  	if (vertexInFace((Face*)p->data,(Vertex*)q->data)) {
				// add face to bad face list
				bf = new void_list();
				bf->next = bfh;
				bf->data = (Face*)p->data;
				bfh = bf;
				flag = false;
				// store pointer
				temp = (Face*)p->data;
				// adjust list pointer
				if (p==ptr1) { ptr1 = p->next; }
				// and remove face from candidate face list
				p=removeLink(p);
				// search candidate bad vertex list for each face vertex
				pp=ptr2;
				while(pp!=NULL) {
					// if bad vertex is one of three face vertices
		    	  	if (vertexInFace(temp,(Vertex*)pp->data)) {
						// add vertex to bad vertex list
						bv = new void_list();
						bv->next = bvh;
						bv->data = (Vertex*)pp->data;
						bvh = bv;
						// adjust list pointer
						if (pp==ptr2) { ptr2 = pp->next; }
						// and remove vertex from candidate vertex list
						pp=removeLink(pp);
					}else { pp=pp->next;}
				}
			}
			if (q->next==NULL) break;
			q=q->next;
		}
		if(flag){
		// search extra vertex list
			// for each group
			for(i=0;i<g.count;i++){
				// for each extra vertex in group
				q=g.c[i].extra;
				while(q!=NULL && flag){
					// if extra vertex is at least one of three face vertices
			   	  	if (vertexInFace((Face*)p->data,(Vertex*)((ExtraVertex*)q->data)->v)) {
						// add face to bad face list
						bf = new void_list();
						bf->next = bfh;
						bf->data = (Face*)p->data;
						bfh = bf;
						flag = false;
						// store pointer
						temp = (Face*)p->data;
						// adjust list pointer
						if (p==ptr1) { ptr1 = p->next; }
						// and remove face from candidate face list
						p=removeLink(p);
						// search candidate bad vertex list for each face vertex
						pp=ptr2;
						while(pp!=NULL) {
							// if bad vertex is one of three face vertices
				    	  	if (vertexInFace(temp,(Vertex*)pp->data)) {
								// add vertex to bad vertex list
								bv = new void_list();
								bv->next = bvh;
								bv->data = (Vertex*)pp->data;
								bvh = bv;
								// adjust list pointer
								if (pp==ptr2) { ptr2 = pp->next; }
								// and remove vertex from candidate vertex list
								pp=removeLink(pp);
							}else { pp=pp->next;}
						}
					}
					if (q->next==NULL) break;
					q=q->next;
				}
			}
		}
		if(flag){ p=p->next;}
	}
	// adjust pointers
	cbf=ptr1;
	cbv=ptr2;
}

void_list* getFaceSearchPool(void_list *F,void_list *M,int z_value,Vertex** vert_array) {
	///// gather faces with at least two vertices from match list
	void_list *p,*q,*pp,*qq,*qh,*ph;
	qh=NULL;
	ph=NULL;
	// for each face
	for (p=F;p!=NULL;p=p->next) {
		// if at least two vertices match z_value
		if((( vert_array[((Face*)p->data)->v1]->z==z_value)&&
			( vert_array[((Face*)p->data)->v2]->z==z_value))||
			((vert_array[((Face*)p->data)->v1]->z==z_value)&&
			(vert_array[((Face*)p->data)->v3]->z==z_value))||
			((vert_array[((Face*)p->data)->v2]->z==z_value)&&
			(vert_array[((Face*)p->data)->v3]->z==z_value))){
			// add face pointer to list
			q = new void_list();
			q->next = qh;
			q->data = (void*)p->data;
			qh = q;
		}
	}
	return qh;
}

void findSharedVertices(void_list *&v1_match_h,void_list *&v2_match_h
						,void_list *&v1_shared_h,void_list *&v2_shared_h,double epsilon) {
	// find shared vertices
	void_list *q,*p,*v1_shared,*v2_shared,*ptr1,*ptr2;
	bool flag;
	ptr1=v1_match_h;
	ptr2=v2_match_h;
	p=ptr1;
	while (p!=NULL) {
		flag = true;
		q=ptr2;
		while (q!=NULL && flag) {
			if ( !distinguishable(((Vertex*)p->data)->x,((Vertex*)q->data)->x,epsilon) &&
				!distinguishable(((Vertex*)p->data)->y,((Vertex*)q->data)->y,epsilon) &&
				!distinguishable(((Vertex*)p->data)->z,((Vertex*)q->data)->z,epsilon)){
				// store vertex from first file
				v1_shared = new void_list();
				v1_shared->next = v1_shared_h;
				v1_shared->data = (Vertex*)p->data;
				v1_shared_h = v1_shared;
				// adjust list pointer
				if (p==ptr1) { ptr1 = p->next; }
				// and remove vertex from list
				p=removeLink(p);
				// store vertex from second file
				v2_shared = new void_list();
				v2_shared->next = v2_shared_h;
				v2_shared->data = (Vertex*)q->data;
				v2_shared_h = v2_shared;
				// adjust list pointer
				if (q==ptr2) { ptr2 = q->next; }
				// and remove vertex from candidate vertex list
				q=removeLink(q);
				flag = false;
			}else{q=q->next;}
		} if(flag){p=p->next;}
	}
	// adjust pointers
	v1_match_h=ptr1;
	v2_match_h=ptr2;
}

void getContours(Group &g,void_list *evfsp,
				Vertex **vert_array,int max_vert,int z_value,void_list *S){
	// Collect vertices and faces for each contour that has shared vertices
	// between the two mesh files
	int i,va,vb,ga,gb,next_group,smallest,largest;
	void_list *p,*q,*pp;
	bool flag=true;
	//create array of ints
	int array[max_vert+1];
	for(i=0;i<max_vert+1;i++){
		array[i]=0;
	}
	next_group=1;

	while (flag){
		flag=false;
		//for each face
		for(p=evfsp;p!=NULL;p=p->next){
			va=vb=0;
			// are any vertices at z_value?
			if (vert_array[((Face*)p->data)->v1]->z==z_value){
				va=((Face*)p->data)->v1;
			}
			if (vert_array[((Face*)p->data)->v2]->z==z_value){
				if(!va){va=((Face*)p->data)->v2;}
				else{vb=((Face*)p->data)->v2;}
			}
			if (vert_array[((Face*)p->data)->v3]->z==z_value){
				vb=((Face*)p->data)->v3;
			}
			if (va && vb) {
				// check to see if any are in a group
				ga=array[va];
				gb=array[vb];
				if(!ga && !gb){
					// if both zero, then add both to same group
					array[va]=array[vb]=next_group;
					next_group++;
					flag=true;
				}else if (ga==gb){
					// if same and nonzero, then do nothing
				}else if (!ga && gb){
					array[va]=array[vb];
					flag=true;
				}else if (ga && !gb){
					array[vb]=array[va];
					flag=true;
				}else {
					flag=true;
					// if different and nonzero
					// then identify smallest group #
					if (ga<gb) {smallest=ga;largest=gb;}
					else {smallest=gb;largest=ga;}
					// then set any vertex with other group #s to smallest group #
					for(i=0;i<max_vert+1;i++){
						if(array[i]==largest) array[i]=smallest;
					}
				}
			}
		}
	}

	///// check each group for shared vertices
	// create array of size next_group
	int shared[next_group];
	for(i=0;i<next_group;i++){
		shared[i]=0;
	}

	int num=1;
	// for each array
	for(i=0;i<max_vert+1;i++){
		// if group nonzero
		if(array[i]){
			flag = true;
			q=S;
			while (q!=NULL && flag) {
				// if vertex is shared
				if(((Vertex*)q->data)->index==i){
					// then set shared[group] to 1
					if(!shared[array[i]]){
					shared[array[i]]=num++;
					}
					flag = false;
				}
				q=q->next;
			}
		}
	}
	num--;

	// create contours
	g.count = num;
	g.c = new Contour[num];
	if(g.c[0].num!=0){printf("THIS SHOULD NOT HAVE PRINTED!!!!");}
	if(g.c[0].extra!=NULL){printf("THIS SHOULD NOT HAVE PRINTED!!!!");}

	// load vertices into contours
	// for each vertex 
	for(i=0;i<max_vert+1;i++){
		if(shared[array[i]]){
			p = new void_list();
			p->next = g.c[shared[array[i]]-1].verts;
			p->data = vert_array[i];
			g.c[shared[array[i]]-1].verts = p;
		}
	}
	
	// load face into contours
	//for each face
	for(p=evfsp;p!=NULL;p=p->next){
		flag = false;
		// for each vertex group
		i=0;
		while(i<num && !flag){
			// for each vertex in group
			q=g.c[i].verts;
			while(q!=NULL && !flag){
				// if vertex is part of face
				if(vertexInFace((Face*)p->data,(Vertex*)q->data)){
					flag=true;
					// then add face to group
					pp = new void_list();
					pp->next = g.c[i].faces;
					pp->data = p->data;
					g.c[i].faces = pp;
				}
				q=q->next;
			}
			i++;
		}
	}
}

void cleanUp(void_list *L){
	void_list *p,*q;
    p=L;
    while (p!=NULL) {
        q=p->next;
        delete p;
        p=q;
    }
}

void cleanUp2V(void_list *L){
	void_list *p,*q;
    p=L;
    while (p!=NULL) {
        q=p->next;
        delete (Vertex*)p->data;
        delete p;
        p=q;
    }
}

void cleanUp2F(void_list *L){
	void_list *p,*q;
    p=L;
    while (p!=NULL) {
        q=p->next;
        delete (Face*)p->data;
        delete p;
        p=q;
    }
}

void printVertices(void_list *L,char *str){
	void_list *q;
	for(q=L;q!=NULL;q=q->next){
        		printf("%s %i %.15g %.15g %.15g\n",str,
		                ((Vertex*)q->data)->index,
		                ((Vertex*)q->data)->x,
		                ((Vertex*)q->data)->y,
		                ((Vertex*)q->data)->z);
	}
}

void printFaces(void_list *L,char *str){
	void_list *q;
	for(q=L;q!=NULL;q=q->next){
        		printf("%s %i %i %i %i\n",str,
		                ((Face*)q->data)->index,
		                ((Face*)q->data)->v1,
		                ((Face*)q->data)->v2,
		                ((Face*)q->data)->v3);
	}
}

void printContourFacesVertices(Group &g,char *str){
	void_list *q;
	int i;
	printf("%s\n",str);
	for(i=0;i<g.count;i++){
		printf("\nContour %i\n",i);
		for(q=g.c[i].verts;q!=NULL;q=q->next){
			printf("vertex %i %.15g %.15g %.15g\n",
					((Vertex*)q->data)->index,
					((Vertex*)q->data)->x,
					((Vertex*)q->data)->y,
					((Vertex*)q->data)->z);
		}
		for(q=g.c[i].faces;q!=NULL;q=q->next){
			printf("face %i %i %i %i\n",
					((Face*)q->data)->index,
					((Face*)q->data)->v1,
					((Face*)q->data)->v2,
					((Face*)q->data)->v3);
		}
	}
}

void printOtherContourData(Group &g,char *str){
	void_list *q;
	int i,j,k;
	ExtraVertex *ev;
	printf("%s\n",str);
	for(i=0;i<g.count;i++){
		printf("\nContour group %i\n",i);
		printf("\nnumber of sites %i\n",g.c[i].num);
		for(q=g.c[i].extra;q!=NULL;q=q->next){
			ev=(ExtraVertex*)q->data;
        	printf("extra vertex group %i\n",ev->g);
        	printf("extra vertex %i %.15g %.15g %.15g\n",
					((Vertex*)ev->v)->index,
					((Vertex*)ev->v)->x,
					((Vertex*)ev->v)->y,
					((Vertex*)ev->v)->z);
			if(((ExtraVertex*)q->data)->la!=NULL){
        		printf("extra vertex la %i %.15g %.15g %.15g\n",
		               ((Vertex*)ev->la)->index,
		               ((Vertex*)ev->la)->x,
		               ((Vertex*)ev->la)->y,
		               ((Vertex*)ev->la)->z);
			}
			if(((ExtraVertex*)q->data)->lb!=NULL){
        		printf("extra vertex lb %i %.15g %.15g %.15g\n",
		               ((Vertex*)ev->lb)->index,
		               ((Vertex*)ev->lb)->x,
		               ((Vertex*)ev->lb)->y,
		               ((Vertex*)ev->lb)->z);
			}
		}
		for(j=0;j<g.c[i].num;j++){
        	printf("Site %i lateral1 %i %.15g %.15g %.15g\n",j,
		              ((Vertex*)g.c[i].s[j].l1)->index,
		              ((Vertex*)g.c[i].s[j].l1)->x,
		              ((Vertex*)g.c[i].s[j].l1)->y,
		              ((Vertex*)g.c[i].s[j].l1)->z);
        	printf("Site %i lateral2 %i %.15g %.15g %.15g\n",j,
		              ((Vertex*)g.c[i].s[j].l2)->index,
		              ((Vertex*)g.c[i].s[j].l2)->x,
		              ((Vertex*)g.c[i].s[j].l2)->y,
		              ((Vertex*)g.c[i].s[j].l2)->z);
        	printf("Site %i third %i %.15g %.15g %.15g\n",j,
		              ((Vertex*)g.c[i].s[j].th)->index,
		              ((Vertex*)g.c[i].s[j].th)->x,
		              ((Vertex*)g.c[i].s[j].th)->y,
		              ((Vertex*)g.c[i].s[j].th)->z);
        	printf("Site %i orientation %i\n",j,g.c[i].s[j].orient);
			for(k=0;k<g.c[i].s[j].n;k++){
				ev=g.c[i].s[j].ev[k];
	        	printf("Site %i ev %i vertex %i %.15g %.15g %.15g\n",j,k,
		              ((Vertex*)ev->v)->index,
		              ((Vertex*)ev->v)->x,
		              ((Vertex*)ev->v)->y,
		              ((Vertex*)ev->v)->z);
			}
		}
	}
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
