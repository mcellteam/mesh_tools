///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

void Object::scanFile(char *filename)
{
    char line[2048],*str;
    FILE *F;
    Vertex *vv;
    Face *ff;
	int vertex_num=0,polygon_num=0;
	std::vector<Vertex*> vp;
    // open file
    F = fopen(filename,"r");
    if (!F) { printf("Couldn't open input file %s\n",filename);return;}
    // for every line in file
    for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) {
        // skip leading whitespace
        while (strchr(" \t,",*str)!=NULL) { str++;}
        // if first character is V for Vertex, add new linked list class instance
        if (strchr("V",*str)!=NULL){
			vv=new Vertex(str);
			v.push_back(vv);
			vp.push_back(vv);
		}
        // if first character is F for Face, add new linked list class instance
        else if (strchr("F",*str)!=NULL){
			ff=new Face(str,vp);
			f.push_back(ff);
		}
    }
    fclose(F);
}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

std::string keyPair(int a,int b,int num_digits){
	char str[128],format[32];
	sprintf(format,"%%0%dd%%0%dd",num_digits,num_digits);
	if (a<b){ sprintf(str,format,a,b);}
	else { sprintf(str,format,b,a); }
	return str;
}

bool edgeMatch(Edge *e,int va,int vb) {
//	Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
//	e->getVertices(v1,v2,o1,o2);
    if ( (e->v1->index==va && e->v2->index==vb) ||
        (e->v1->index==vb && e->v2->index==va) ){return true;}
    else {return false;}
}

Edge* findEdge(Vertex* va,Vertex* vb,hashtable_t &hm,int num_digits){
	Edge *ee=NULL;
	std::string s = keyPair(va->index,vb->index,num_digits);
	// if element exists given key, then get Edge pointer
	if (hm.count(s)>0){ ee=hm[s]; }
	return ee;
}

void Edge::update(Face *f){
    //add face to edge
	if(f1==NULL) {f1=f;}
	else if (f2==NULL) {f2=f;}
	else { cout << "Error. Tried to add third face to edge.\n"
				<< "Face " << f->index 
				<< " " << f->v[0]->index
				<< " " << f->v[1]->index
				<< " " << f->v[2]->index
				<< endl;
			exit(1); 
	}
	// add edge pointer to face
	f->addEdge(this);
}

void Face::addEdge(Edge* ptr){
	if(e[0]==NULL){e[0]=ptr;}
	else if(e[1]==NULL){e[1]=ptr;}
	else if(e[2]==NULL){e[2]=ptr;}
	else { cout << "Error. Tried to add fourth edge to face.\n"
				<< "Face " << index 
				<< " " << (int)v[0]->index
				<< " " << (int)v[1]->index
				<< " " << (int)v[2]->index
				<< endl;
			exit(1); 
	}
}

void Object::createEdge(Face *ff,Vertex* va,Vertex* vb,hashtable_t &hm,int num_digits){
	// new edge
	Edge *en = new Edge(ff,va,vb);
//	Edge *en = new Edge(ff);
	// store edge pointer in hash table
	hm[keyPair(va->index,vb->index,num_digits)]=en;
	// add edge pointer to face
	ff->addEdge(en);
	// add edge pointer to object
	e.push_back(en);
}

void Object::checkEdge(Face *ff,Vertex *va,Vertex *vb,hashtable_t &hm,int num_digits) {
	Edge *ee=NULL;
    ee=findEdge(va,vb,hm,num_digits);
    if(ee!=NULL){ee->update(ff);}
    else {createEdge(ff,va,vb,hm,num_digits);}
}

int Object::setNumDigits(void){
	int max=0;
	// for each vertex in object
	for(std::vector<Vertex*>::iterator i=v.begin();i!=v.end();i++){
		if((*i)->index>max){max=(*i)->index;}
	}
	char str[64];
	sprintf(str,"%d",max);
	std::string s = str;
	return s.length();
}

void Object::clearEdges(void) {
	for(std::vector<Edge*>::iterator i=e.begin();i!=e.end();i++)
	{
		delete *i;
	}
	e.clear();
	for(std::vector<Face*>::iterator i=f.begin();i!=f.end();i++)
	{
		(*i)->e[0]=(*i)->e[1]=(*i)->e[2]=NULL;
	}
}

void Object::createEdges(void) {
	// clear existing edges if necessary
	if(e.empty()==false){clearEdges();}
	// determine number of digits in largest vertex index
	int num_digits = setNumDigits();
	// create map for finding edges
	hashtable_t hm;
	std::vector<Face*>::iterator i;
	// for each face
	for (i=f.begin();i!=f.end();i++) {
		checkEdge(*i,(*i)->v[0],(*i)->v[1],hm,num_digits);
        checkEdge(*i,(*i)->v[1],(*i)->v[2],hm,num_digits);
        checkEdge(*i,(*i)->v[2],(*i)->v[0],hm,num_digits);
	}
}
// #####################################################
// #####################################################

bool refineMesh(Object &o,double t){
	bool flag = true;
	fprintf(stderr,"Thresholding edges..............");
	fflush(stderr);
	flag = o.thresholdEdges(t);
	fprintf(stderr,"complete.\n");
	fflush(stderr);
	fprintf(stderr,"Creating new vertices...........");
	fflush(stderr);
	int n = o.createNewVertices();
	fprintf(stderr,"complete. %d new verts.\n",n);
	fflush(stderr);
	fprintf(stderr,"Creating new faces..............");
	fflush(stderr);
	n = o.createNewSubdividedFaces();
	fprintf(stderr,"complete. %d new faces.\n",n);
	fflush(stderr);
	fprintf(stderr,"Process new edges...............");
	fflush(stderr);
	n = o.processNewEdges();
	fprintf(stderr,"complete. %d new edges.\n",n);
	fflush(stderr);
	fprintf(stderr,"Check edges.....................");
	fflush(stderr);
	n = o.checkEdges();
	fprintf(stderr,"complete. %d border edges.\n",n);
	fflush(stderr);
	return flag;
}

int Object::checkEdges(void)
{
	int n=0;
	// for each edge
	for(std::vector<Edge*>::iterator i=e.begin();i!=e.end();i++)
	{
		// if edge is open
			if( ((*i)->f1==NULL && (*i)->f2!=NULL) ||
				((*i)->f1!=NULL && (*i)->f2==NULL) )
			{
				n++;
			}
	}
	return n;
}

Edge* Object::findMatchingBorderEdge2(Edge *ee)
{
	Edge *t=NULL;
	// for each border edge
	for(std::vector<BorderEdge*>::iterator i=be.begin();i!=be.end();i++)
	{
		// if border edge matches two vertices
			if( (*i)->e==ee )
			{
				return (*i)->e;
			}
	}
	return NULL;
}

Edge* Object::findMatchingBorderEdge(Vertex *v1,Vertex *v2)
{
	Edge *t=NULL;
	// for each border edge
	for(std::vector<BorderEdge*>::iterator i=be.begin();i!=be.end();i++)
	{
		// if border edge matches two vertices
			if( ((*i)->v1==v1 && (*i)->v2==v2) ||
				((*i)->v1==v2 && (*i)->v2==v1) )
			{
				return (*i)->e;
			}
	}
	return NULL;
}

int Object::processNewEdges(void)
{
	int n=0;
	// for each new face
	for(std::vector<Face*>::iterator i=nf.begin();i!=nf.end();i++)
	{
		// for each edge of new face
		for(int j=0;j<3;j++)
		{
			Vertex *v1 = (*i)->v[j];
			Vertex *v2 = (*i)->v[(j+1)%3];
			Edge* ee = findMatchingBorderEdge(v1,v2);
		    if(ee!=NULL)
			{
				// edge is border edge so add face to edge
				if(ee->f1==NULL) { ee->f1=*i; }
				else if (ee->f2==NULL) { ee->f2=*i; }
				else
				{
					cout << "Error! third face added to edge!\n";
					exit(0);
				}
			}
			else
			{
				n++;
				// new edge
				ee = new Edge(*i);
				// add edge pointer to face
				(*i)->addEdge(ee);
				// add edge pointer to object
				e.push_back(ee);
				// add edge pointer to border edges
				BorderEdge *bb = new BorderEdge(ee,v1,v2);
				be.push_back(bb);
			}
		}
	}
	return n;
}

bool Edge::threshold(double t){
	double l = getSqLength();
	return l>t;
}

bool Object::thresholdEdges(double t)
{
	bool flag = false;
	for(std::vector<Face*>::iterator i=f.begin();i!=f.end();i++)
	{
		// threshold criteria
		if ((*i)->e[0]->threshold(t)) {
			// mark face as subdivided
			(*i)->subdivided=true;
			// mark edges as bisected
			(*i)->e[0]->bisected=true;
			// set flag
			flag = true;
		}
		if ((*i)->e[1]->threshold(t)) {
			// mark face as subdivided
			(*i)->subdivided=true;
			// mark edges as bisected
			(*i)->e[1]->bisected=true;
			// set flag
			flag = true;
		}
		if ((*i)->e[2]->threshold(t)) {
			// mark face as subdivided
			(*i)->subdivided=true;
			// mark edges as bisected
			(*i)->e[2]->bisected=true;
			// set flag
			flag = true;
		}
	}
	return flag;
}

int Object::createNewVertices(void)
{
	int n=0;
	// for each edge
	for (std::vector<Edge*>::iterator i=e.begin();i!=e.end();i++){
		// if edge is bisected
		if((*i)->bisected==true){
			Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
			(*i)->getVertices(v1,v2,o1,o2);
			// compute new vertex
			double x=(v1->p[0]+v2->p[0])/2.0;
			double y=(v1->p[1]+v2->p[1])/2.0;
			double z=(v1->p[2]+v2->p[2])/2.0;
			// create new vertex
			Vertex *vv = new Vertex(v.size()+1,x,y,z);
			// save vertex in object
			v.push_back(vv);
			// save vertex in edge
			(*i)->v=vv;
			n++;
		}
	}
	return n;
}

int Face::matchEdges(Vertex* new_verts[3]){
	// count number of subdivided edges
	int j=0;
	for(int i=0;i<3;i++)
	{
		new_verts[i]=e[i]->v;
		if(new_verts[i]!=NULL){j++;}
	}
	return j;
}

void addFace(std::vector<Face*> &nf,int num_faces,Vertex *va,Vertex *vb,Vertex *vc){
	Face *f = new Face(num_faces,va,vb,vc);
	nf.push_back(f);
}

double dot(double a[3],double b[3])
{
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double getAspectRatio(Vertex *v0,Vertex *v1,Vertex *v2){

    // Make triangle edge vectors
    double va[3]={v1->p[0]-v0->p[0],v1->p[1]-v0->p[1],v1->p[2]-v0->p[2]};
    double vb[3]={v2->p[0]-v1->p[0],v2->p[1]-v1->p[1],v2->p[2]-v1->p[2]};
    double vc[3]={v0->p[0]-v2->p[0],v0->p[1]-v2->p[1],v0->p[2]-v2->p[2]};
    double vbase[3];
    double vopp[3];

  // Find length of longest edge
  double lmax=-DBL_MAX;
  double la=sqrt(dot(va,va));
  double lb=sqrt(dot(vb,vb));
  double lc=sqrt(dot(vc,vc));
  if (la>lmax)
  {
    lmax=la;
    vbase[0]=va[0];
    vbase[1]=va[1];
    vbase[2]=va[2];
    vc[0]=v2->p[0]-v0->p[0];
    vc[1]=v2->p[1]-v0->p[1];
    vc[2]=v2->p[2]-v0->p[2];
    vopp[0]=vc[0];
    vopp[1]=vc[1];
    vopp[2]=vc[2];
  }
  if (lb>lmax)
  {
    lmax=lb;
    vbase[0]=vb[0];
    vbase[1]=vb[1];
    vbase[2]=vb[2];
    va[0]=v0->p[0]-v1->p[0];
    va[1]=v0->p[1]-v1->p[1];
    va[2]=v0->p[2]-v1->p[2];
    vopp[0]=va[0];
    vopp[1]=va[1];
    vopp[2]=va[2];
  }
  if (lc>lmax)
  {
    lmax=lc;
    vbase[0]=vc[0];
    vbase[1]=vc[1];
    vbase[2]=vc[2];
    vb[0]=v1->p[0]-v2->p[0];
    vb[1]=v1->p[1]-v2->p[1];
    vb[2]=v1->p[2]-v2->p[2];
    vopp[0]=vb[0];
    vopp[1]=vb[1];
    vopp[2]=vb[2];
  }


  // Find shortest altitude
    double ll = sqrt(dot(vbase,vbase));
    vbase[0]=vbase[0]/ll;
    vbase[1]=vbase[1]/ll;
    vbase[2]=vbase[2]/ll;
    double dot_prod = dot(vbase,vopp);
    double alt[3]={vopp[0]-(dot_prod*vbase[0]),
                    vopp[1]-(dot_prod*vbase[1]),
                    vopp[2]-(dot_prod*vbase[2])};
    double amin=sqrt(dot(alt,alt));

    return lmax/amin;

}

double getMaxAR(Vertex *v0,Vertex *v1,Vertex *v2,Vertex *v3){
	// get max aspect ratio connecting v0 and v2
	// get aspect ratio of face v0,v1,v2
	double ar012 = getAspectRatio(v0,v1,v2);
	// get aspect ratio of face v0,v2,v3
	double ar023 = getAspectRatio(v0,v2,v3);
	// identify max aspect ratios
	if( ar012 > ar023 ){ return ar012;}
	else { return ar023;}
}

void Object::clearBorderEdges(void) {
	for(std::vector<BorderEdge*>::iterator i=be.begin();i!=be.end();i++)
	{
		delete *i;
	}
	be.clear();
}

int Object::createNewSubdividedFaces(void)
{
	int foo=1;
	nf.clear();
	clearBorderEdges();
	Vertex* new_verts[3];
//	std::vector<Face*> nf; // new faces
	// for each face
	for(std::vector<Face*>::iterator i=f.begin();i!=f.end();i++)
	{
		if(foo%1000==0){
			cout << foo << " of " << f.size() << endl;
		}
		foo++;
		// if face is subdivided
		if((*i)->subdivided==true){
			// identify new vertices
			int num_subdivided_edges = (*i)->matchEdges(new_verts);
			///// if one edge subdivided /////
			if(num_subdivided_edges==1){
				// add new faces
				if(new_verts[0]!=NULL){
					addFace(nf,++max_faces,(*i)->v[0],new_verts[0],(*i)->v[2]);
					addFace(nf,++max_faces,(*i)->v[2],new_verts[0],(*i)->v[1]);
				} else if (new_verts[1]!=NULL) {
					addFace(nf,++max_faces,(*i)->v[1],new_verts[1],(*i)->v[0]);
					addFace(nf,++max_faces,(*i)->v[0],new_verts[1],(*i)->v[2]);
				} else if (new_verts[2]!=NULL) {
					addFace(nf,++max_faces,(*i)->v[2],new_verts[2],(*i)->v[1]);
					addFace(nf,++max_faces,(*i)->v[1],new_verts[2],(*i)->v[0]);
				}
			} else if(num_subdivided_edges==2){
				if(new_verts[0]==NULL){
					// get max aspect ratio connecting new_verts[2] and f->v2
					double max1 = getMaxAR(new_verts[2],(*i)->v[0],(*i)->v[1],new_verts[1]);
					// get max aspect ratio connecting new_verts[1] and f->v1
					double max2 = getMaxAR(new_verts[1],new_verts[2],(*i)->v[0],(*i)->v[1]);
					// if connecting new_verts[2],f->v2 generates
					// faces with smallest aspect ratio
					if(max1<max2){
						// add new faces
						addFace(nf,++max_faces,(*i)->v[2],new_verts[2],new_verts[1]);
						addFace(nf,++max_faces,new_verts[2],(*i)->v[0],(*i)->v[1]);
						addFace(nf,++max_faces,new_verts[2],(*i)->v[1],new_verts[1]);
					}
					// else connecting new_verts[1],f->v1 generates
					// faces with smallest aspect ratio
					else {
						// add new faces
						addFace(nf,++max_faces,(*i)->v[2],new_verts[2],new_verts[1]);
						addFace(nf,++max_faces,new_verts[1],(*i)->v[0],(*i)->v[1]);
						addFace(nf,++max_faces,new_verts[2],(*i)->v[0],new_verts[1]);
					}
				} else if (new_verts[1]==NULL) {
					double max1 = getMaxAR(new_verts[0],(*i)->v[1],(*i)->v[2],new_verts[2]);
					double max2 = getMaxAR(new_verts[2],new_verts[0],(*i)->v[1],(*i)->v[2]);
					// if connecting new_verts[0],v3 generates
					// faces with smallest aspect ratio
					if(max1<max2){
						// add new faces
						addFace(nf,++max_faces,(*i)->v[0],new_verts[0],new_verts[2]);
						addFace(nf,++max_faces,new_verts[0],(*i)->v[1],(*i)->v[2]);
						addFace(nf,++max_faces,new_verts[0],(*i)->v[2],new_verts[2]);
					}
					// else connecting new_verts[2],v2 generates
					// faces with smallest aspect ratio
					else {
						// add new faces
						addFace(nf,++max_faces,(*i)->v[0],new_verts[0],new_verts[2]);
						addFace(nf,++max_faces,new_verts[0],(*i)->v[1],new_verts[2]);
						addFace(nf,++max_faces,(*i)->v[1],(*i)->v[2],new_verts[2]);
					}
				} else if (new_verts[2]==NULL) {
					double max1 = getMaxAR(new_verts[1],(*i)->v[2],(*i)->v[0],new_verts[0]);
					double max2 = getMaxAR(new_verts[0],new_verts[1],(*i)->v[2],(*i)->v[0]);
					// if connecting new_verts[1],v1 generates
					// faces with smallest aspect ratio
					if(max1<max2){
						// add new faces
						addFace(nf,++max_faces,(*i)->v[1],new_verts[1],new_verts[0]);
						addFace(nf,++max_faces,new_verts[1],(*i)->v[2],(*i)->v[0]);
						addFace(nf,++max_faces,new_verts[1],(*i)->v[0],new_verts[0]);
					}
					// else connecting new_verts[0],v3 generates
					// faces with smallest aspect ratio
					else {
						// add new faces
						addFace(nf,++max_faces,(*i)->v[1],new_verts[1],new_verts[0]);
						addFace(nf,++max_faces,new_verts[1],(*i)->v[2],new_verts[0]);
						addFace(nf,++max_faces,(*i)->v[2],(*i)->v[0],new_verts[0]);
					}
				}
			} else if(num_subdivided_edges==3){
				// add new faces
				addFace(nf,++max_faces,(*i)->v[0],new_verts[0],new_verts[2]);
				addFace(nf,++max_faces,new_verts[0],(*i)->v[1],new_verts[1]);
				addFace(nf,++max_faces,new_verts[0],new_verts[1],new_verts[2]);
				addFace(nf,++max_faces,new_verts[2],new_verts[1],(*i)->v[2]);
			} else {
				printf("\n\ncreateNewSubdividedFaces: Error: Face is subdivided, ");
				printf("but num_subdivided_edges = %d\n\n",num_subdivided_edges);
				exit(0);
			}
			// remove face from edge
			buildBorderEdges(*i);
			// free memory
			delete *i;
			// delete face* from vector
//			f.erase(i);
			*i=NULL;
//			i++;
//			df.push_back(*i);
		}// else {i++;}
	}
	// add new faces to old faces
//	f.insert(f.end(),nf.begin(),nf.end());
	return nf.size();
}

void Object::buildBorderEdges(Face *ff)
{
	// for each face edge
	for(int i=0;i<3;i++)
	{
//		// build Border edge
//		Edge* ee = findMatchingBorderEdge2(ff->e[i]);
//		if(ee==NULL)
//		{
//			BorderEdge *bb = new BorderEdge(ff->e[i]);
//			be.push_back(bb);
//		}
		// remove Face* from edge
		if(ff==ff->e[i]->f1){ ff->e[i]->f1=NULL; }
		else				{ ff->e[i]->f2=NULL; }
	}
}

/*
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
*/
/*
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
}*/

