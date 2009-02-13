double dot(double a[3],double b[3]){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

bool getPointEdgeDistance(double p[3],Face *f){
	// for each face edge
	for (int i=0;i<3;i++) {
		double Ax = f->v[i]->pN[0];
		double Ay = f->v[i]->pN[1];
		double Az = f->v[i]->pN[2];
		double Bx = f->v[(i+1)%3]->pN[0];
		double By = f->v[(i+1)%3]->pN[1];
		double Bz = f->v[(i+1)%3]->pN[2];
		double uDen = (Bx-Ax)*(Bx-Ax)+(By-Ay)*(By-Ay)+(Bz-Az)*(Bz-Az);
		if(uDen) {
			// u = AdotA-AdotB-AdotC+BdotC)/uDen
			double u =((p[0]-Ax)*(Bx-Ax)+(p[1]-Ay)*(By-Ay)+(p[2]-Az)*(Bz-Az))/uDen;
			// no need to check for u ==0 and u ==1, since current 
			// vertex/face plane coincidence was checked previously.
			// Closest point on face edge line to current vertex
			// occurs on face edge between face vertices
			if (u>0 && u<1) {
				double a=Ax+u*(Bx-Ax);
				double b=Ay+u*(By-Ay);
				double c=Az+u*(Bz-Az);
				if(!distinguishable(p[0],a) && 
					!distinguishable(p[1],b) && 
					!distinguishable(p[2],c)){ return true; }
			}
		}
	}
	return false;
}

void biggest(double *x, int& big) {
    if ( x[0] < x[1] ) {
        // x[0] < x[1]
        if ( x[0] < x[2] ) {
            // x[0] < x[1]
            // x[0] < x[2]
            if (x[1] < x[2]) {
                // x[0] < x[1]
                // x[1] < x[2]
                // x[0] is smallest
                // x[2] is biggest
                big = 2;
            } else {
                // x[0] < x[1]
                // x[2] <= x[1]
                // x[0] is smallest
                // x[1] is biggest or tied with x[2]
                big = 1;
            }
        } else {
            // x[0] < x[1]
            // x[2] <= x[0]
            // x[2] is smallest or tied with x[0]
            // x[1] is biggest
            big = 1;
        }
    } else {
        // x[1] <= x[0]
        if ( x[0] < x[2] ) {
            // x[1] <= x[0]
            // x[0] < x[2]
            // x1 is smallest or tied with x[0]
            // x2 is biggest
            big = 2;
        } else {
            big = 0;
        }
    }
}

void threeValueSort(double p1,double p2,double p3, double &biggest, double &smallest) {
    if ( p1 < p2 ) {
        // x[0] < x[1]
        if ( p2 < p3 ) {
            // x[0] < x[1]
            // x[1] < x[2]
            smallest = p1;
            biggest = p3;
		} else {
            // x[0] < x[1]
            // x[2] <= x[1]
            if (p1 < p3) {
                // x[0] < x[1]
                // x[0] < x[2]
                // x[2] <= x[1]
                smallest = p1;
                biggest = p2;
            } else {
                // x[0] < x[1]
                // x[2] <= x[1]
                // x[2] <= x[0]
                smallest = p3;
                biggest = p2;
            }
        }
    } else {
        // x[1] <= x[0]
        if ( p1 < p3 ) {
            // x[1] <= x[0]
            // x[0] < x[2]
            smallest = p2;
            biggest = p3;
        } else {
            // x[1] <= x[0]
            // x[2] <= x[0]
            if (p2 < p3) {
                // x[1] <= x[0]
                // x[2] <= x[0]
                // x[1] < x[2]
                smallest = p2;
                biggest = p1;
            } else {
                // x[1] <= x[0]
                // x[2] <= x[0]
                // x[2] <= x[1]
                smallest = p3;
                biggest = p1;
            }
        }
	}
}

// #####################################################
// #####################################################


void checkLineFaceIntersection(Face *f,double lp[2][3],bool &line_flag,
								bool &poly_flag, bool &poly_edge_flag) {
	//lp  = line_points
	// initialize flags
	line_flag=poly_flag=poly_edge_flag=false;
	// get face normal
	double pn[3];
	f->getNormal(pn);
	// use line points, polygon normal and one polygon vertex
	// to compute where line intersects plane of polygon (if not parallel)
	// pI=p1+u*(p2-p1) where u=dot(N,(p3-p1))/dot(N,(p2-p1))
	// pI = intersection of line and plane
	// p3 = any point on plane
	// p2,p1 = end points of line
	// N = plane normal
	// denominator of u
	double den = pn[0]*(lp[1][0]-lp[0][0])+pn[1]*(lp[1][1]-lp[0][1])+pn[2]*(lp[1][2]-lp[0][2]);
	// if line and polygon plane are not parallel
	if (den) {
		//pvc = polygon_vertex_coordinates
		Vertex *v0=f->v[0]; 
		// point of intersection
		double u = (pn[0]*(v0->pN[0]-lp[0][0]) 
			+ pn[1]*(v0->pN[1]-lp[0][1]) 
			+ pn[2]*(v0->pN[2]-lp[0][2]))/den;
		// if polygon cuts through line
		if (u > DOUBLE_EPSILON && u < DOUBLE_EPSILON_COMP) {
			line_flag = true;
		}
		// compute polygon double-area on each of three principal planes
		int tv[2],big;
		double pI[2],I[3];//,d[3];
		double area[3] = { fabs(pn[2]),// xy
							fabs(pn[0]),// yz
							fabs(pn[1])// zx
							};
		biggest(area,big);
		// pI = p1+u*(p2-p1)
		// where lp[1][] is p2 and lp[0][] is p1
		tv[0] = big;
		tv[1] = (big+1)%3;
		pI[0] = (1-u)*lp[0][big] + u*lp[1][big];
		pI[1] = (1-u)*lp[0][(big+1)%3] + u*lp[1][(big+1)%3];
		I[big]=pI[0];
		I[(big+1)%3]=pI[1];
		I[(big+2)%3] = (1-u)*lp[0][(big+2)%3] + u*lp[1][(big+2)%3];
		////////// is point of intersection on other polygon? //////////
		//pvc = polygon_vertex_coordinates
		Vertex *v1=f->v[1],*v2=f->v[2]; 
		// does point of intersection lie on polygon edge?
		if (getPointEdgeDistance(I,f)){poly_edge_flag = true;}
		// if point of intersection is not on polygon edge
		if (poly_edge_flag==false) {
			// compute three determinants
			double det[3];
			det[0] = (v0->pN[tv[0]]-pI[0])*(v1->pN[tv[1]]-pI[1])
					-(v1->pN[tv[0]]-pI[0])*(v0->pN[tv[1]]-pI[1]);
			det[1] = (v1->pN[tv[0]]-pI[0])*(v2->pN[tv[1]]-pI[1])
					-(v2->pN[tv[0]]-pI[0])*(v1->pN[tv[1]]-pI[1]);
			// proceed if determinants are DOUBLE_EPSILON away from zero
			if ((det[0]*det[1])>0){
				det[2]=(v2->pN[tv[0]]-pI[0])*(v0->pN[tv[1]]-pI[1])
						-(v0->pN[tv[0]]-pI[0])*(v2->pN[tv[1]]-pI[1]);
				if ((det[0]*det[2])>0){
					// line intersects polygon plane inside polygon
				    poly_flag = true;
				}
			}
		}
	}
}

// #####################################################
// #####################################################

int Container::checkEdgeEdgeIntersection(Face *cf,Face *of,bool share_edge_flag) {
	// cpvc = current_polygon_vertex_coordinates
	// opvc = other_polygon_vertex_coordinates
	// cv   = current_vertex
	// ov   = other_vertex
	Vertex *cv[2],*ov[2];
	// for each current face edge
	for (int i=0;i<3;i++) {
		// for each other face edge
		for (int j=0;j<3;j++) {
			cv[0] = cf->v[pairs[i][0]];
			cv[1] = cf->v[pairs[i][1]];
			ov[0] = of->v[pairs[j][0]];
			ov[1] = of->v[pairs[j][0]];
			// if the edges do not share a vertex
			if (cv[0]!=ov[0]&&cv[0]!=ov[1]&&cv[1]!=ov[0]&&cv[1]!=ov[1]) {
				// and the edges are not parallel
				bool parallel_flag = false;
				double x[3],y[3];
				for (int k=0;k<3;k++) {
					x[k] = cv[1]->pN[k]-cv[0]->pN[k];
					y[k] = ov[1]->pN[k]-ov[0]->pN[k];
				}
				double term1 = x[0]*y[0]+ x[1]*y[1]+ x[2]*y[2];
				if ( !distinguishable(term1*term1,
				(x[0]*x[0]+ x[1]*x[1]+ x[2]*x[2])*(y[0]*y[0]+ y[1]*y[1]+ y[2]*y[2])) ) {parallel_flag = true;}
				if (!parallel_flag) {
					// compute scalars
					double qDen = (cv[1]->pN[0]-cv[0]->pN[0])*(ov[1]->pN[1]-ov[0]->pN[1])
							-(cv[1]->pN[1]-cv[0]->pN[1])*(ov[1]->pN[0]-ov[0]->pN[0]);
					double qNum = (cv[1]->pN[0]-cv[0]->pN[0])*(cv[0]->pN[1]-ov[0]->pN[1])
							-(cv[1]->pN[1]-cv[0]->pN[1])*(cv[0]->pN[0]-ov[0]->pN[0]);
					double q = qNum/qDen;
					double uNum = ov[0]->pN[0]-cv[0]->pN[0]+q*(ov[1]->pN[0]-ov[0]->pN[0]);
					double uDen = cv[1]->pN[0]-cv[0]->pN[0];
					if(fabs(qDen)>DOUBLE_EPSILON && fabs(uDen)>DOUBLE_EPSILON) {
						double u = uNum/uDen;
						if ( 
							((u > DOUBLE_EPSILON && u < (1-DOUBLE_EPSILON)) && 
							(q > DOUBLE_EPSILON && q < (1-DOUBLE_EPSILON)) ) ||
							( share_edge_flag && ((u > DOUBLE_EPSILON && u < (1-DOUBLE_EPSILON)) || 
												(q > DOUBLE_EPSILON && q < (1-DOUBLE_EPSILON))) )
							) {
						    return(1);
						}
					}
				}
			}
		}
	}
	return(0);
}

// #####################################################
// #####################################################

int Container::checkFaceEdgeIntersection(Face *cf,Face *of){
	// NOTE THAT INTERSECTION CHECK IS ASSYMETRIC
	// I.E. CHECK IS IF CURRENT FACE EDGE INTERSECTS OTHER FACE ONLY.
	// NO CHECKING IS DONE IF CURRENT FACE IS INTERSECTED BY OTHER FACE EDGE.
	// THIS IS OK SINCE THE ROLES OF CURRENT AND OTHER FACE GET REVERSED LATER
	// IN THE CODE, SO ALL NECESSARY CHECKS ARE EVENTUALLY EXECUTED.
	// 
	// Actually, I just added code to make the test symmetric.
	double lp[2][3],*cpvc[3];
	bool line_flag=false, poly_flag=false, poly_edge_flag;
	//cpvc = current_polygon_vertex_coordinates
	// get face vertex coordinates
	cf->getVertexCoordinates(cpvc);
	// for each current polygon edge
	for (int i=0;i<3;i++) {
		if(!line_flag || !poly_flag) {
			lp[0][0] = cpvc[pairs[i][0]][0];
			lp[0][1] = cpvc[pairs[i][0]][1];
			lp[0][2] = cpvc[pairs[i][0]][2];
			lp[1][0] = cpvc[pairs[i][1]][0];
			lp[1][1] = cpvc[pairs[i][1]][1];
			lp[1][2] = cpvc[pairs[i][1]][2];
			checkLineFaceIntersection(of,lp,line_flag,poly_flag,poly_edge_flag);
		}
	}
	if (line_flag && (poly_flag || poly_edge_flag)) {return(1);}
	// get face vertex coordinates
	of->getVertexCoordinates(cpvc);
	// for each current polygon edge
	for (int i=0;i<3;i++) {
		if(!line_flag || !poly_flag) {
			lp[0][0] = cpvc[pairs[i][0]][0];
			lp[0][1] = cpvc[pairs[i][0]][1];
			lp[0][2] = cpvc[pairs[i][0]][2];
			lp[1][0] = cpvc[pairs[i][1]][0];
			lp[1][1] = cpvc[pairs[i][1]][1];
			lp[1][2] = cpvc[pairs[i][1]][2];
			checkLineFaceIntersection(cf,lp,line_flag,poly_flag,poly_edge_flag);
		}
	}
	// do polygons intersect?
	if (line_flag && (poly_flag || poly_edge_flag)) {return(1);}
	else {return(0);}
}

// #####################################################
// #####################################################

double Face::getAngle(Vertex *vv){
	Vertex *vA=vv,*vB=NULL,*vC=NULL;
	double AB[3],AC[3],abL,acL,costheta;
	// identify face vertices
	if (v[0]!=vv){vB=v[0];}
	else if (v[1]!=vv){vB=v[1];}
	else if (v[2]!=vv){vB=v[2];}
	if (v[0]!=vv && v[0]!=vB){vC=v[0];}
	else if (v[1]!=vv && v[1]!=vB){vC=v[1];}
	else if (v[2]!=vv && v[2]!=vB){vC=v[2];}
	// AB,AC
	AB[0]=vB->pN[0]-vA->pN[0];
	AB[1]=vB->pN[1]-vA->pN[1];
	AB[2]=vB->pN[2]-vA->pN[2];
	AC[0]=vC->pN[0]-vA->pN[0];
	AC[1]=vC->pN[1]-vA->pN[1];
	AC[2]=vC->pN[2]-vA->pN[2];
	// lengths
	acL=sqrt( dot(AC,AC) );
	abL=sqrt( dot(AB,AB) );
	costheta=( dot(AB,AC) )/abL/acL;
	cout.precision(12);
	if (costheta > 1) costheta=1;
	if (costheta < -1) costheta=-1;
	return acos(costheta);
}

void Vertex::getNormal(double *n) {
	f_iterator i;
	double t[3],theta,thetaT=0,L;
	n[0]=n[1]=n[2]=0;
	// for each adjacent face
	for (i=f.begin();i!=f.end();i++) {
		// get coordinates of polygon normal
		(*i)->getNormal(t);
		L=sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
		theta=(*i)->getAngle(this);
		thetaT+=theta;
		// and add to sum
		n[0] += t[0]/L*theta;
		n[1] += t[1]/L*theta;
		n[2] += t[2]/L*theta;
	}
	n[0] = n[0]/f.size()/thetaT;
	n[1] = n[1]/f.size()/thetaT;
	n[2] = n[2]/f.size()/thetaT;
}

// #####################################################
// #####################################################

double Edge::getAngle(void) {
	Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
	getVertices(v1,v2,o1,o2);
	// get outward normals of edge faces
	double n1[3],n2[3];
	f1->getNormal(n1);
	f2->getNormal(n2);
	// compute the cosine of angle between normals
	double normal_angle_cosine=dot(n1,n2)/sqrt(dot(n1,n1))/sqrt(dot(n2,n2));
	// compute angle between normals 
    if 		(normal_angle_cosine >= 1)	{
		return PI;
	}  else if (normal_angle_cosine <= -1)	{
		// normal_angle = PI;
		// gamma == 0 or 2PI
		return 0;
	} else {
		// normal_angle = acos(normal_angle_cosine);
	    // use the edge itself as a reference vector
	    double refvec[3] = {v2->pN[0]-o2->pN[0],v2->pN[1]-o2->pN[1],v2->pN[2]-o2->pN[2]};
	    // dot product of refvec and n1
	    double d = dot(refvec,n1);
		if(!d){
			return PI;
		} else {
	    return PI+d/fabs(d)*acos(normal_angle_cosine);
	}
		// I EXPECT 0 <= angle < 2*PI
	}
}

// #####################################################
// #####################################################

std::string keyPair(int a,int b,int num_digits){
	char str[128],format[32];
	sprintf(format,"%%0%dd%%0%dd",num_digits,num_digits);
	if (a<b){ sprintf(str,format,a,b);}
	else { sprintf(str,format,b,a); }
	return str;
}

bool edgeMatch(Edge *e,int va,int vb) {
	Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
	e->getVertices(v1,v2,o1,o2);
    if ( (v1->index==va && v2->index==vb) ||
        (v1->index==vb && v2->index==va) ){return true;}
    else {return false;}
}

Edge* Object::findEdge(Vertex* va,Vertex* vb,map_se &hm,int num_digits){
	Edge *ee=NULL;
	std::string s = keyPair(va->index,vb->index,num_digits);
	// if element exists given key, then get Edge pointer
	if (hm.count(s)>0){ ee=hm[s]; }
	return ee;
}

void Edge::getVertices(Vertex *&v1,Vertex *&v2,Vertex *&o1,Vertex *&o2){
	if(f1!=NULL && f2!=NULL){
		v1=vv1;
		v2=vv2;
		// find o1 on f1
		if(f1->v[0]!=vv1 && f1->v[0]!=vv2){o1=f1->v[0];}
		else if(f1->v[1]!=vv1 && f1->v[1]!=vv2){o1=f1->v[0];}
		else if(f1->v[2]!=vv1 && f1->v[2]!=vv2){o1=f1->v[0];}
		else {
			cout << "\n\nEdge::getVertices: o1 not identified!\n";
			cout << "	v1:\n";
			v1->printVertex(v1->o->name);
			cout << endl << "	v2:\n";
			v2->printVertex(v2->o->name);
			cout << endl << "	o1:\n";
			o1->printVertex(o1->o->name);
			cout << endl << "	o2:\n";
			o2->printVertex(o2->o->name);
			cout << endl << "	vv1:\n";
			vv1->printVertex(vv1->o->name);
			cout << endl << "	vv2:\n";
			vv2->printVertex(vv2->o->name);
			cout << endl;
			cout << endl;
			exit(0);
		}
		// find o2 on f2
		if(f2->v[0]!=vv1 && f2->v[0]!=vv2){o2=f2->v[0];}
		else if(f2->v[1]!=vv1 && f2->v[1]!=vv2){o2=f2->v[0];}
		else if(f2->v[2]!=vv1 && f2->v[2]!=vv2){o2=f2->v[0];}
		else {
			cout << "\n\nEdge::getVertices: o2 not identified!\n";
			cout << "	v1:\n";
			v1->printVertex(v1->o->name);
			cout << endl << "	v2:\n";
			v2->printVertex(v2->o->name);
			cout << endl << "	o1:\n";
			o1->printVertex(o1->o->name);
			cout << endl << "	o2:\n";
			o2->printVertex(o2->o->name);
			cout << endl << "	vv1:\n";
			vv1->printVertex(vv1->o->name);
			cout << endl << "	vv2:\n";
			vv2->printVertex(vv2->o->name);
			cout << endl;
			cout << endl;
			exit(0);
		}
	} else {
		// use f1 only
		v1=vv1;
		v2=vv2;
		// find o1 on f1
		if(f1->v[0]!=vv1 && f1->v[0]!=vv2){o1=f1->v[0];}
		else if(f1->v[1]!=vv1 && f1->v[1]!=vv2){o1=f1->v[0];}
		else if(f1->v[2]!=vv1 && f1->v[2]!=vv2){o1=f1->v[0];}
		else {
			cout << "\n\nEdge::getVertices: o1 not identified!\n";
			cout << "	v1:\n";
			v1->printVertex(v1->o->name);
			cout << endl << "	v2:\n";
			v2->printVertex(v2->o->name);
			cout << endl << "	o1:\n";
			o1->printVertex(o1->o->name);
			cout << endl << "	o2:\n";
			o2->printVertex(o2->o->name);
			cout << endl << "	vv1:\n";
			vv1->printVertex(vv1->o->name);
			cout << endl << "	vv2:\n";
			vv2->printVertex(vv2->o->name);
			cout << endl;
			cout << endl;
			exit(0);
		}
	}
}

void Edge::update(Face *f){
    //add face to edge
	if(f1==NULL) {f1=f;}
	else if (f2==NULL) {f2=f;}
	else { fvec.push_back(f); }
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

void Object::createEdge(Face *ff,Vertex* va,Vertex* vb,map_se &hm,int num_digits){
	// new edge
	Edge *en = new Edge(ff,va,vb);
	// store edge pointer in hash table
	hm[keyPair(va->index,vb->index,num_digits)]=en;
	// add edge pointer to face
	ff->addEdge(en);
	// add edge pointer to object
	e.push_back(en);
}

void Object::buildEdge(Face *ff,Vertex *va,Vertex *vb,map_se &hm,int num_digits) {
	Edge *ee=NULL;
    ee=findEdge(va,vb,hm,num_digits);
	// NOTE: 
    if(ee!=NULL){
		ee->update(ff);
	}
    else {
		createEdge(ff,va,vb,hm,num_digits);
	}
}

int Object::setNumDigits(void){
	int max=0;
	// for each vertex in object
	for(v_iterator i=v.begin();i!=v.end();i++){
		if((*i)->index>max){max=(*i)->index;}
	}
	char str[64];
	sprintf(str,"%d",max);
	std::string s = str;
	return s.length();
}

void Object::createEdges(void) {
	cout << "Create edges for object [" << name << "]...............................";
	cout.flush();
	// determine number of digits in largest vertex index
	int num_digits = setNumDigits();
	// create map for finding edges
	map_se hm;
	f_iterator i;
	// for each face
	for (i=f.begin();i!=f.end();i++) {
		buildEdge(*i,(*i)->v[0],(*i)->v[1],hm,num_digits);
        buildEdge(*i,(*i)->v[1],(*i)->v[2],hm,num_digits);
        buildEdge(*i,(*i)->v[2],(*i)->v[0],hm,num_digits);
	}
	cout << "complete.\n";cout.flush();
}

//void Object::verifyEdges(hashtable_t &hm,Face *ff){
void Object::verifyEdges(void){
	for(e_iterator i=e.begin();i!=e.end();i++){
		Edge *ee=*i;
		if(ee->valid()==false){
			cout << "Object::verifyEdges: Error! edge is invalid.\n";
			ee->printEdge(ee->vv1->o->name);
			cout << endl << endl;
			exit(0);
		}
		
	}

}

void Vertex::getAdjacentVertices(vector<Vertex*> &a){
	a.clear();
	e_iterator i;
	// for each adjacent edge
	for (i=e.begin();i!=e.end();i++) {
		Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
		(*i)->getVertices(v1,v2,o1,o2);
		// find vertex different from self and add different vertex to vector
		if (v1->index!=index){a.push_back(v1);}
		else if (v2->index!=index) {a.push_back(v2);}
		else { printf("Error. both vertices of edge are equal to current vertex.\n"); exit(1); }
	}
}

double Object::getMeanEdgeLength(void){
	e_iterator i;
	double L = 0;
	// for each edge
	for (i=e.begin();i!=e.end();i++) {
		L+=sqrt((*i)->getSqLength());
	}
	// compute mean edge length
	return (L/(double)e.size());
}

void Vertex:: getAdjacentFaces(hset_f &fset){
	e_iterator i;
	// for each adjacent edge
	for (i=e.begin();i!=e.end();i++) {
		// add edge faces to set
		fset.insert((*i)->f1);
		fset.insert((*i)->f2);
	}
}

bool Object::processEdge(Edge *ee,hmap_vd &hood,vec_e &bucket,Vertex *vv){
    bool empty = true;
    bool f1 =false,f2=false;
	Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
	ee->getVertices(v1,v2,o1,o2);
    Vertex *vp1=v1,*vp2=v2;
    // if v1 found in hashtable, then f1 = true;
    if(hood.find(vp1)!=hood.end()){f1=true;}
    // if v2 found in hashtable, then f1 = true;
    if(hood.find(vp2)!=hood.end()){f2=true;}
    // if neither v1 nor v2 are in struct list
    if(!f1&&!f2){
        // then add edge* to bucket
        bucket.push_back(ee);
        empty = false;
    } else if (f1 && f2){
    } else {
        if(f1){
            // if v1 found and v2 not found in hashtable
            // then add v2 to hashtable
            hood[vp2]=sqrt( (vp2->pN[0]-vv->pN[0])*(vp2->pN[0]-vv->pN[0])+
                            (vp2->pN[1]-vv->pN[1])*(vp2->pN[1]-vv->pN[1])+
                            (vp2->pN[2]-vv->pN[2])*(vp2->pN[2]-vv->pN[2]));
        }
        if(f2){
            // if v2 found and v1 not found in hashtable
            // then add v1 to hashtable
            hood[vp1]=sqrt( (vp1->pN[0]-vv->pN[0])*(vp1->pN[0]-vv->pN[0])+
                            (vp1->pN[1]-vv->pN[1])*(vp1->pN[1]-vv->pN[1])+
                            (vp1->pN[2]-vv->pN[2])*(vp1->pN[2]-vv->pN[2]));
        }
    }
    return empty;
}

void Object::collectFaces(hmap_vd &hood,set_v &disabled,vec_f &new_faces){
    new_faces.clear();
    // for each vertex in hood
    for(vd_iterator i=hood.begin();i!=hood.end();i++){
        // if vertex is thawed and not disabled
        if(!ifFrozen(hood,(*i).first) && disabled.find((*i).first)==disabled.end()){
            // for each adjacent face of thawed vertex
            for(f_iterator j=(*i).first->f.begin();j!=(*i).first->f.end();j++){
                // if any face vertex is thawed, then add face to collection
                if( !ifFrozen(hood,(*j)->v[0]) ||
                    !ifFrozen(hood,(*j)->v[1]) ||
                    !ifFrozen(hood,(*j)->v[2])){new_faces.push_back(*j);}
            }
            // add vertex to disabled list
            disabled.insert((*i).first);
        }
    }
}

bool Object::ifFrozen(hmap_vd &neighborhood,Vertex *vv){
    // if vertex is in hashtable
    if(neighborhood.find(vv)!=neighborhood.end()){
        if (neighborhood[vv]<NEIGHBORHOOD_RADIUS){return false;}
        else {return true;}
    }
    return true;
}

bool Object::thawedAndAble(hmap_vd &hood,set_v &disabled){
    // for each vertex in hood
    for(vd_iterator i=hood.begin();i!=hood.end();i++){
        // if vertex is thawed and not disabled
        if(!ifFrozen(hood,(*i).first) && disabled.find((*i).first)==disabled.end()){
            // then affirm that at least one vertex is thawed and able
            return true;
        }
    }
    return false;
}

void Object::newFindNeighborhoods(void){
    // for each vertex in object
    for (v_iterator i=v.begin();i!=v.end();i++) {
        // initialize vertex*->double hash table
        // represents a neighbor vertex and
        // the shortest cumulative edge length
        // to reach it from current vertex
        hmap_vd hood;
        hood.clear();
        // add current vertex to neighbor list
        // naturally, assign it zero length
        hood[*i]=0.0;
        // initialize set of vertices to constitute disabled list
        set_v disabled;
        disabled.clear();
        // init collection of neighborhood faces
        vec_f c;
        ///// initial round /////
        // init collection of new faces
        vec_f new_faces;
        collectFaces(hood,disabled,new_faces);
        // for each face in new collection
        for(f_iterator k=new_faces.begin();k!=new_faces.end();k++){
            // initialize container for edges
            // with neither vertex in hood
            vec_e bucket;
            bucket.clear();
            bool bucket_empty = true;
            // for each edge in face
            for(int j=0;j<3;j++){
                if(!processEdge((*k)->e[j],hood,bucket,*i)){
                    bucket_empty = false;
                }
            }
            if(!bucket_empty){
                // initialize another container for edges
                // with neither vertex in hood
                vec_e pail;
                pail.clear();
                bool pail_empty = true;
                // for each edge in bucket
                for(e_iterator j=bucket.begin();j!=bucket.end();j++){
                    if(!processEdge(*j,hood,pail,*i)){
                        pail_empty = false;
                    }
                }
                if(!pail_empty){
                    cout << "Error. Multiple rounds of bucket use required.\n";
                    exit(0);
                }
            }
        }

        ///// all subsequent rounds /////

        // while there are thawed vertices in neighbor list, hood,
        // that are also not disabled
        while(thawedAndAble(hood,disabled)){
            // init collection of new faces
            vec_f new_faces_too;
            collectFaces(hood,disabled,new_faces_too);
            // for each face in new collection
            for(f_iterator k=new_faces_too.begin();k!=new_faces_too.end();k++){
                // initialize container for edges
                // with neither vertex in hood
                vec_e bucket;
                bucket.clear();
                bool bucket_empty = true;
                // for each edge in face
                for(int j=0;j<3;j++){
                    if(!processEdge((*k)->e[j],hood,bucket,*i)){
                        bucket_empty = false;
                    }
                }
                if(!bucket_empty){
                    // initialize another container for edges
                    // with neither vertex in hood
                    vec_e pail;
                    pail.clear();
                    bool pail_empty = true;
                    // for each edge in bucket
                    for(e_iterator j=bucket.begin();j!=bucket.end();j++){
                        if(!processEdge(*j,hood,pail,*i)){
                            pail_empty = false;
                        }
                    }
                    if(!pail_empty){
                        cout << "Error. Multiple rounds of bucket use required.\n";
                        exit(0);
                    }
                }
                // add face to neighborhood, c
                c.push_back(*k);
            }
        }
        // copy local neighborhood to Object class neighborhood faces vector
        (*i)->nf.assign(c.begin(),c.end());
        // sort vectors
        sort((*i)->nf.begin(),(*i)->nf.end());
    }
}

void Object::findVertexAdjacencies(void){
	cout << "Finding vertex adjacencies for object [" << name << "].................";
	cout.flush();
	e_iterator i;
	f_iterator j;
	// for each edge, add edge* to both edge vertices
	for (i=e.begin();i!=e.end();i++){
		Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
		(*i)->getVertices(v1,v2,o1,o2);
		v1->e.push_back(*i);
		v2->e.push_back(*i);
	}
	// for each face, add face* to each face vertex
	for (j=f.begin();j!=f.end();j++){
		((*j)->v[0])->f.push_back(*j);
		((*j)->v[1])->f.push_back(*j);
		((*j)->v[2])->f.push_back(*j);
	}
	cout << "complete.\n";cout.flush();
}

// ######################################
// ######################################

void Object::boundObject(double* r) {
	v_iterator i;
	double xmin,xmax,ymin,ymax,zmin,zmax,x,y,z;
	//initialize mins and maxes
	xmin = v[0]->pN[0];
	xmax = v[0]->pN[0];
	ymin = v[0]->pN[1];
	ymax = v[0]->pN[1];
	zmin = v[0]->pN[2];
	zmax = v[0]->pN[2];
	// for each vertex in object
	for (i=v.begin();i!=v.end();i++) {
		///////// extract coordinates //////////
		x = (*i)->pN[0];
		y = (*i)->pN[1];
		z = (*i)->pN[2];
		if (x>xmax) {xmax = x;}
		else if (x<xmin) {xmin = x;}
		if (y>ymax) {ymax = y;}
		else if (y<ymin) {ymin = y;}
		if (z>zmax) {zmax = z;}
		else if (z<zmin) {zmin = z;}
	}
	r[0]=xmin;r[1]=xmax;
	r[2]=ymin;r[3]=ymax;
	r[4]=zmin;r[5]=zmax;
}

// #####################################################
// #####################################################

void Container::createEdges(void) {
	o_iterator i;
	// for each object, create edges
	for (i=o.begin();i!=o.end();i++) {
		(*i)->createEdges();
	}
}

void Container::findVertexAdjacencies(void) {
	o_iterator i;
	// for each object, find vertex adjacencies
	for (i=o.begin();i!=o.end();i++) {
		(*i)->findVertexAdjacencies();
	}
}

// #####################################################
// #####################################################

void Container::scanDir(const char *filename) {
	struct dirent *pent;			// pointer to dirent structure
    DIR *pdir = opendir(filename);	// pointer to a directory data structure
    if (!pdir) {printf("Error. Could not open %s.\n",filename);exit(1);}
	else { cout << "\nFolder found " << filename << endl << endl;}
    while ((pent=readdir(pdir))){
		// copy char array to string
		std::string str = pent->d_name;
		// if file of typ *.mesh
		std::string::size_type found = str.find(".mesh",0);
		// if found
        if (found != std::string::npos) {
			// save filename
			files.push_back(str);
			// update index
			num_files++;
			// print file found to screen
//			cout << "file found: " << str << "\n"; cout.flush();
		}
    }
    closedir(pdir);

}
// #####################################################
// #####################################################

Object* Container::processFile(std::string filename) {
	// create new Object
	Object *obj = new Object(filename);
	// scan file
	scanFile(obj,filename);
	// check object contents
	if(obj->v.empty()==true || obj->f.empty()==true){
		delete obj;
		cout << "\n Container::processFile: "
		<< "no valid mesh object found in "
		<< filename << ". Skipping file.\n";
		return NULL;
	} else {
		// save Object* in container
		o.push_back(obj);
		// save first vertex* in object
			if(obj->v.empty()==false){
			eg=obj->v.front();
		} else {
			cout << "\nContainer::processFile: Error. "
			<< "Object " << obj->name
			<< " contains no vertices.\n";
			exit(0);
		}
		// return
		return obj;
	}
}

void Container::update(Object *oo){
	// update Container
	num_obj++;
	num_vert+=oo->v.size();
	// DEBUG
//	cout << "\nContainer::update: "
//	<< "oo->f.size()=" << oo->f.size() << endl;
	// DEBUG
	num_face+=oo->f.size();
	num_edge+=oo->e.size();
	num_sep+=oo->num_sep;
//	if(oo->manifold==true){num_man_cc++;}
//	else {num_man_nn++;}
	if(oo->manifold==true){num_man[0]++;}
	else if(oo->manifold==false && oo->closed==true) {num_man[1]++;}
	else {num_man[2]++;}
	num_bou+=oo->num_bou;
	num_indistin+=oo->indistin_v.size();
	if     (oo->manifold==true && oo->consistent==true) {num_cons[0]++;}
	else if(oo->manifold==true && oo->consistent==false){num_cons[1]++;}
	else {num_cons[2]++;}
	num_vol+=oo->vol;
	if(oo->closed==true){num_clo_cc++;}
	else {num_clo_nn++;}
	num_bor+=oo->border.size();
	num_nonman_v+=oo->nonman_v.size();
	num_nonman_e+=oo->nonman_e.size();
	num_flip+=oo->flipped.size();
	if     (oo->consistent==true && oo->outward==true) {num_out[0]++;}
	else if(oo->consistent==true && oo->outward==false){num_out[1]++;}
	else {num_out[2]++;}
	num_orph+=oo->orphan.size();
	num_mis+=oo->missing_v.size();
	num_deg+=oo->degen.size();
	num_dupl_v+=oo->dupl_v_index.size();
	num_dupl_f+=oo->dupl_f_index.size();
}

// #####################################################
// #####################################################

void Container::scanFile (Object *obj,std::string filename) {
    char line[2048];
    // open file
    FILE *F = fopen(filename.c_str(),"r");
    if (!F) {
		cout <<"Couldn't open input file "
		<< filename << endl;
		exit(0);
	} else {
//		cout << "\n\n" << "/* ********************** "
//		<< "OBJECT ********************** */\n";
		//	print object name 
//		cout << "name: " << obj->name << endl;
//		cout.flush();
//		cout << "file found: " << filename << "\n"; cout.flush();
	}
    // for every line in file
    for (char *str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) {
        // skip leading whitespace
        while (strchr(" \t,",*str)!=NULL) { str++;}
        // if first character is V for Vertex, add new linked list class instance
        if (strchr("V",*str)!=NULL){
			Vertex *v=new Vertex(str,obj);
			obj->v.push_back(v);
			obj->vp.insert(std::make_pair(v->index,v));
			obj->found.insert(std::make_pair(v->index,false));
		}
        // if first character is F for Face, add new linked list class instance
        else if (strchr("F",*str)!=NULL){
			Face *f=new Face(str,obj);
			obj->f.push_back(f);
		}
    }
    fclose(F);
}

void Container::writeDistances(void) {
	char file[FILENAME_SIZE];
	// create output filename
	sprintf(file,"closest_point_distances.dat");
	// open output file
	std::ofstream newfile (file,std::ios::out);
	if(newfile.is_open()){
		newfile.precision(4);
		// for each object
		for(o_iterator i=o.begin();i!=o.end();i++) {
			// for each vertex in object
			for(v_iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
				// if vertex has a closest face
				if((*j)->cl!=NULL){
					// compute separation vector
					double s[3];
					for(int k=0;k<3;k++){ s[k]=(*j)->pC[k]-(*j)->pN[k]; }
					// print separation distance
					if((*i)->vertexIsNice(*j)){ newfile << sqrt(dot(s,s)) << endl;}
					else {newfile << -sqrt(dot(s,s)) << endl;}
				}
			}
		}
		newfile.close();
	}
}

// #####################################################
// #####################################################
void Space::initBoxes(double num_faces) {
	// subdivide space
	double world_volume = (world[1]-world[0])*(world[3]-world[2])*(world[5]-world[4]);
	double box_volume = FACES_PER_BOX * world_volume / num_faces;
//	cout << "Space::initBoxes: "
//	<< "FACES_PER_BOX = " << FACES_PER_BOX
//	<< ", world_volume = " << world_volume
//	<< ", num_faces = " << num_faces << "\n";
	space_length = pow ( fabs ( box_volume ), 1.0 / 3.0 );
//	cout << "Space::initBoxes: "
//	<< "space_length = " << space_length << "\n";
/*	cout << "\nSpace::initBoxes: "
	<< "world ["
	<< world[0] << ","
	<< world[1] << ","
	<< world[2] << ","
	<< world[3] << ","
	<< world[4] << ","
	<< world[5] << "]\n";
*/
	num_space[0] = (int) ceil( (world[1]-world[0])/space_length );
	num_space[1] = (int) ceil( (world[3]-world[2])/space_length );
	num_space[2] = (int) ceil( (world[5]-world[4])/space_length );
/*	cout << "Space::initBoxes: "
	<< "num_space ["
	<< num_space[0] << ","
	<< num_space[1] << ","
	<< num_space[2] << "]\n";
*/
	num_boxes = num_space[0]*num_space[1]*num_space[2];
	// allocate memory for boxes
	b.reserve(num_boxes);
	// store box limits in boxes class
	// for each box
	for (int z =0;z<num_space[2];z++) {
		for (int y =0;y<num_space[1];y++) {
			for (int x =0;x<num_space[0];x++) {
				b.push_back(new Box(x,y,z));
			}
		}
	}
}
// ######################################
// ######################################

void Container::boundWorld(Space &s) {
	o_iterator i;
	double xmin,xmax,ymin,ymax,zmin,zmax,range[6];
	//initialize mins and maxes
//	xmin = o[0]->v[0]->pN[0];
//	xmax = o[0]->v[0]->pN[0];
//	ymin = o[0]->v[0]->pN[1];
//	ymax = o[0]->v[0]->pN[1];
//	zmin = o[0]->v[0]->pN[2];
//	zmax = o[0]->v[0]->pN[2];
	xmin = eg->pN[0];
	xmax = eg->pN[0];
	ymin = eg->pN[1];
	ymax = eg->pN[1];
	zmin = eg->pN[2];
	zmax = eg->pN[2];
    ////////// loop through all objects //////////
    // for each object
	for (i=o.begin();i!=o.end();i++) {
		// get range of object vertices
		(*i)->boundObject(&range[0]);
		if (range[1]>xmax) {xmax = range[1];}
		if (range[0]<xmin) {xmin = range[0];}
		if (range[3]>ymax) {ymax = range[3];}
		if (range[2]<ymin) {ymin = range[2];}
		if (range[5]>zmax) {zmax = range[5];}
		if (range[4]<zmin) {zmin = range[4];}
    }
	if (xmin<0) {s.world[0]=xmin*1.01;} else {s.world[0]=xmin*0.99;}
	if (xmax<0) {s.world[1]=xmax*0.99;} else {s.world[1]=xmax*1.01;}
	if (ymin<0) {s.world[2]=ymin*1.01;} else {s.world[2]=ymin*0.99;}
	if (ymax<0) {s.world[3]=ymax*0.99;} else {s.world[3]=ymax*1.01;}
	if (zmin<0) {s.world[4]=zmin*1.01;} else {s.world[4]=zmin*0.99;}
	if (zmax<0) {s.world[5]=zmax*0.99;} else {s.world[5]=zmax*1.01;}
}

// #####################################################
// #####################################################

void Space::computeBoxesToCheck(Face *f,vec_b &bp) {
	vec_d xv,yv,zv;
	// identify face bounding box limits
	xv.push_back(f->v[0]->pN[0]);
	xv.push_back(f->v[1]->pN[0]);
	xv.push_back(f->v[2]->pN[0]);
	yv.push_back(f->v[0]->pN[1]);
	yv.push_back(f->v[1]->pN[1]);
	yv.push_back(f->v[2]->pN[1]);
	zv.push_back(f->v[0]->pN[2]);
	zv.push_back(f->v[1]->pN[2]);
	zv.push_back(f->v[2]->pN[2]);
	sort(xv.begin(),xv.end());
	sort(yv.begin(),yv.end());
	sort(zv.begin(),zv.end());
	// grab face 3D location range
	double dr[6];
	dr[0] = xv[0];  // -x
	dr[1] = xv[2];	//  x
	dr[2] = yv[0];	// -y
	dr[3] = yv[2];	//  y
	dr[4] = zv[0];  // -z
	dr[5] = zv[2];	//  z
/*	// DEBUG
	cout << "\nSpace::computeBoxesToCheck: "
	<< "data range ["
	<< dr[0] << ","
	<< dr[1] << ","
	<< dr[2] << ","
	<< dr[3] << ","
	<< dr[4] << ","
	<< dr[5] << "]\n";
*/	// DEBUG
	// collect boxes to check
	getBoxesFor3DLocations(dr,bp);
}

// #####################################################
// #####################################################

void Face::getNormal(double n[3]) {
	double uX, uY, uZ, vX, vY, vZ;
	// compute vectors 01 and 12
	if(v[1]==NULL){cout << "v[1]==NULL\n";cout.flush();}
	if(v[2]==NULL){cout << "v[2]==NULL\n";cout.flush();}
	uX = v[1]->pN[0]-v[0]->pN[0];
	uY = v[1]->pN[1]-v[0]->pN[1];
	uZ = v[1]->pN[2]-v[0]->pN[2];
	vX = v[2]->pN[0]-v[0]->pN[0];
	vY = v[2]->pN[1]-v[0]->pN[1];
	vZ = v[2]->pN[2]-v[0]->pN[2];
	// compute cross product (u x v)
	n[0] = uY*vZ-uZ*vY;
	n[1] = uZ*vX-uX*vZ;
	n[2] = uX*vY-uY*vX;
}

// #####################################################
// #####################################################
void Container::getExtraRay(Vertex *v,double lp[2][3],int index) {
	// get normal info
	double n[3];
	v->f[index]->getNormal(n);
	// compute centroid of first adjacent face
	double cx = (v->f[index]->v[0]->pN[0]+
				v->f[index]->v[1]->pN[0]+
				v->f[index]->v[2]->pN[0])/3.0;
	double cy = (v->f[index]->v[0]->pN[1]+
				v->f[index]->v[1]->pN[1]+
				v->f[index]->v[2]->pN[1])/3.0;
	double cz = (v->f[index]->v[0]->pN[2]+
				v->f[index]->v[1]->pN[2]+
				v->f[index]->v[2]->pN[2])/3.0;
	double L=sqrt( dot(n,n) );
	lp[0][0] = cx;
	lp[1][0] = lp[0][0]+n[0]/L*RAY_EPSILON;
	lp[0][1] = cy;
	lp[1][1] = lp[0][1]+n[1]/L*RAY_EPSILON;
	lp[0][2] = cz;
	lp[1][2] = lp[0][2]+n[2]/L*RAY_EPSILON;
}

int Container::findExtraPoint(Space &s,Vertex *v,double p[3],int index){
	double lp[2][3];
	int num_odd_objects;
	vec_i face_flags;
	vec_o tmp;
	f_iterator jj;
	vec_f crossed_faces,unique_faces;
	std::pair<f_iterator,vec_f::iterator> pp;
	std::pair<o_iterator,vec_o::iterator> ppp;
	getExtraRay(v,lp,index); // returns ray
	// look for intersected faces along ray, i.e. between face centroid and RAY ORIGIN
	collectNiceFaces(s,lp,unique_faces);
	findIntersectedFaces(lp,unique_faces,crossed_faces,face_flags);
	sort(crossed_faces.begin(),crossed_faces.end());
	// look for adjacent face in crossed_faces
	pp=equal_range(crossed_faces.begin(),crossed_faces.end(),v->f[index]);
	// if found, then remove
	if(pp.first!=pp.second){crossed_faces.erase(pp.first);}
	findOddMeshes(crossed_faces,face_flags,num_odd_objects,tmp);
	p[0]=lp[1][0];
	p[1]=lp[1][1];
	p[2]=lp[1][2];
	if (!tmp.empty()){
		ppp=equal_range(tmp.begin(),tmp.end(),v->o);
		if(ppp.first!=ppp.second){
			return 0;
		}
 	}
	return 1;
}

void Container::findCrossed1(Space &s,Vertex *v,double lp[2][3],vec_o &c){
	// find and return crossed objects between pN and extracellular point
	int num_odd_objects;
	vec_i face_flags;
	f_iterator i;
	vec_f crossed_faces,unique_faces;
	std::pair<f_iterator,vec_f::iterator> pp;	
	collectNiceFaces(s,lp,unique_faces);
	findIntersectedFaces(lp,unique_faces,crossed_faces,face_flags);
	// remove current vertex adjacent faces from crossed_faces
	// for each adjacent face
	sort(crossed_faces.begin(),crossed_faces.end());
	for(i=v->f.begin();i!=v->f.end();i++){
		pp=equal_range(crossed_faces.begin(),crossed_faces.end(),*i);
		// if adjacent face is found in crossed_faces, then remove from crossed_faces
		if(pp.first!=pp.second){crossed_faces.erase(pp.first);}
	}
	findOddMeshes(crossed_faces,face_flags,num_odd_objects,c);		
}

void Container::findCrossed2(Space &s,double lp[2][3],vec_o &c){
	// find and return crossed objects between pN and extracellular point
	int num_odd_objects;
	vec_i face_flags;
	vec_f crossed_faces,unique_faces;
	collectNiceFaces(s,lp,unique_faces);
	findIntersectedFaces(lp,unique_faces,crossed_faces,face_flags);
	findOddMeshes(crossed_faces,face_flags,num_odd_objects,c);		
}

void Container::collectCrossed(Space &s,Vertex *v,vec_o &cb){
	vec_o ca;
	std::pair<o_iterator,vec_o::iterator> pp;
	double lp[2][3],p[3];
	// find point, p, outside of current object
	unsigned int index = 0;
	while (!findExtraPoint(s,v,p,index)){
		index++;
		if (index>(v->f.size()-1)){
			cout << "\n\nEvery adjacent face failed!\n";
			cout << v->o->name << "->" << v->index 
			<< " current vertex [" << v->pN[0] << " "
			<< v->pN[1] << " "
			<< v->pN[2] << "]"
			<< endl;
			exit(0);
		}
	}
	// grab intersected objects between pN and RAY ORIGIN and return as ca
	lp[0][0]=v->pN[0];
	lp[0][1]=v->pN[1];
	lp[0][2]=v->pN[2];
	lp[1][0]=p[0];
	lp[1][1]=p[1];
	lp[1][2]=p[2];
	findCrossed1(s,v,lp,ca);
	// grab intersected objects along RAY as cb
	lp[0][0]=p[0];
	lp[0][1]=p[1];
	lp[0][2]=p[2];
	lp[1][0]=p[0];
	lp[1][1]=p[1];
	lp[1][2]=p[2];
	findClosestAxis(s,v,lp); // returns ray
	lp[0][0]=p[0];
	lp[0][1]=p[1];
	lp[0][2]=p[2];
	lp[1][1]=p[1];
	lp[1][2]=p[2];
	if(v->index==531 && !strcmp(v->o->name.c_str(),"a001_FILTERED_SMOOTH_SMOOTH")){
	}
	findCrossed2(s,lp,cb);
	// remove all meshes in ca from cb
	// since the object is not odd relative to pN
	sort(cb.begin(),cb.end());
	for (o_iterator ii=ca.begin();ii!=ca.end();ii++){
		pp=equal_range(cb.begin(),cb.end(),*ii);
		// if object in ca is found in cb, then remove from cb
		if(pp.first!=pp.second){cb.erase(pp.first);}
		// else add it to cb
		else {cb.push_back(*ii);}
	}
}

bool Container::updateNiceness(Vertex *v,vec_o &cb){
	std::pair<o_iterator,vec_o::iterator> pp;
	int old_nice = v->getVertexNiceness();
	// if vertex niceness changes then set flag=true
	bool flag = false;
	// if cb is not empty, then vertex is not nice
	if (!cb.empty()) {
		v->setVertexNiceness(1);
		// if vertex was nice
		if (!old_nice){
			flag = true;
			pp=equal_range(cb.begin(),cb.end(),v->o);
			// if vertex is inside self object
			if(pp.first!=pp.second){
				v->setVertexNiceness(2);
				cout << endl << endl 
				<< v->o->name << "->" << v->index 
				<< " vertex inside self [" << v->pN[0] << " "
				<< v->pN[1] << " "
				<< v->pN[2]
				<< "], cb.size() " << cb.size()
				<< endl;
				for (f_iterator jj=v->nf.begin();jj!=v->nf.end();jj++){
					cout << v->o->name << "->" << (*jj)->index 
					<< " adjacent face\n";
					(*jj)->printFace(v->o);
					cout << endl;
				}
			}
		}
	}else{ // else cb is empty, then vertex is nice
		// if vertex was nonnice, but not to self
		if (old_nice==1){	
			flag=true;
		}
		// if vertex was at least nonnice to self
		else if (old_nice==2){
			flag=true;
		}
		// update niceness
		v->setVertexNiceness(0);			
	}
	return flag;
}

bool Container::checkNiceness(Space &s,Vertex *v) {
	vec_o cb;
	// collect objects inside which vertex lies
	collectCrossed(s,v,cb);
	// update niceness of vertex based on cb
	return updateNiceness(v,cb);
}

void Container::findNice(Space &s) {
	cout << "Iteration 0: ";
	cout << "find nice vertices................";
	cout.flush();
	o_iterator i;
	v_iterator j;
    // for each object in container
	for (i=o.begin();i!=o.end();i++) {
        // for each vertex in object
		for (j=(*i)->v.begin();j!=(*i)->v.end();j++) {
			checkNiceness(s,*j);
		}
	}
	cout << "complete.\n";
	cout.flush();
}

// #####################################################
// #####################################################

void Container::findClosestAxis(Space &s,Vertex *v,double lp[2][3]) {
	// n=normal,r=ray
	// get normal info
	// identify nearest boundary
	double dis[6] = {
					fabs(v->pN[0]-s.world[0]),
					fabs(v->pN[0]-s.world[1]),
					fabs(v->pN[1]-s.world[2]),
					fabs(v->pN[1]-s.world[3]),
					fabs(v->pN[2]-s.world[4]),
					fabs(v->pN[2]-s.world[5])};
	int i=0;
	double min = dis[i];
	for(int j=1;j<6;j++){
		if(dis[j]<min){i=j;min=dis[j];}
	}
	// configure ray
	if		(i==0){
		lp[1][0] = lp[0][0]-2*(s.world[1]-s.world[0]);	// end x
		lp[1][1] = lp[0][1];							// end y
		lp[1][2] = lp[0][2];							// end z
	} else if (i==1){
		lp[1][0] = lp[0][0]+2*(s.world[1]-s.world[0]);	// end x
		lp[1][1] = lp[0][1];							// end y
		lp[1][2] = lp[0][2];							// end z
	} else if (i==2){
		lp[1][0] = lp[0][0];							// end x
		lp[1][1] = lp[0][1]-2*(s.world[3]-s.world[2]);	// end y
		lp[1][2] = lp[0][2];							// end z
	} else if (i==3){
		lp[1][0] = lp[0][0];							// end x
		lp[1][1] = lp[0][1]+2*(s.world[3]-s.world[2]);	// end y
		lp[1][2] = lp[0][2];							// end z
	} else if (i==4){
		lp[1][0] = lp[0][0];							// end x
		lp[1][1] = lp[0][1];							// end y
		lp[1][2] = lp[0][2]-2*(s.world[5]-s.world[4]);	// end z
	} else if (i==5){
		lp[1][0] = lp[0][0];							// end x
		lp[1][1] = lp[0][1];							// end y
		lp[1][2] = lp[0][2]+2*(s.world[5]-s.world[4]);	// end z
	}

}

// #####################################################
// #####################################################

inline bool comp(Face *a,Face* b) { return a->index == b->index; }

void Container::getBoxIndexRange(Space &s,double lp[2][3],int br[6]){
	// compute box index range that contains ray
	// note this range is zero lower-bounded (lowest range is zeroth box)
	// total range is 0..num_space[i]-1

	br[0] = s.location2Index(lp[0][0],"x");
	br[1] = s.location2Index(lp[1][0],"x");
	br[2] = s.location2Index(lp[0][1],"y");
	br[3] = s.location2Index(lp[1][1],"y");
	br[4] = s.location2Index(lp[0][2],"z");
	br[5] = s.location2Index(lp[1][2],"z");
}

void Container::collectNiceFaces(Space &s,double lp[2][3],vec_f &uf) {
	int br[6];
	vec_b b;
	b_iterator i;
	f_iterator j,new_end;
	
	bool flag = false;
	if(
		!distinguishable(lp[0][0],2719.54642515,1E-8) &&
		!distinguishable(lp[0][1],4249.96388355,1E-8) &&
		!distinguishable(lp[0][2],6242.92495541,1E-8)
		){flag=true;}

	// compute box index range that contains ray
	getBoxIndexRange(s,lp,br);

	if(flag){
		cout << "\nContainer::collectNiceFaces "
		<< " box index range ["
		<< br[0] << " "
		<< br[1] << " "
		<< br[2] << " "
		<< br[3] << " "
		<< br[4] << " "
		<< br[5] << "]\n";
	}

	////////// collect boxes to check //////////
	s.getBoxesFor3DIndices(br,b);
	////////// gather faces in boxes //////////
	// for each box
	for (i=b.begin();i!=b.end();i++) {
		// for each face in box
		for (j=(*i)->f.begin();j!=(*i)->f.end();j++) {
			uf.push_back(*j);
		}
	}
	// keep unique faces
	sort(uf.begin(),uf.end());
	new_end = unique(uf.begin(),uf.end());
	uf.assign(uf.begin(),new_end);
}

// #####################################################
// #####################################################

void Container::findIntersectedFaces(double lp[2][3],vector<Face*> &uf,
								vector<Face*> &cf,vector<int> &ff) {
	// uf = unique_faces
	// cf = crossed_faces
	// ff = face_flags
	f_iterator j;
	bool line_flag, poly_flag, poly_edge_flag;
	// for each unique polygon
	for (j=uf.begin();j!=uf.end();j++) {
		checkLineFaceIntersection(*j,lp,line_flag,poly_flag,poly_edge_flag);
		// does point intersect polygon
		if (line_flag && (poly_flag||poly_edge_flag)) {
			// add polygon_index to crossed array
			cf.push_back(*j);
			if (detect_polygon_edge_intersection) {
				// if intersection point falls on edge, signal with flag
				if ( poly_edge_flag ) {ff.push_back(1);}
				else {ff.push_back(0);}
			}
		}
	}
}

// #####################################################
// #####################################################

void Container::findOddMeshes(vec_f &cf,vec_i &ff,
							int& num_odd_objects,vec_o &tmp) {
	// find mesh objects crossed an odd number of times by ray
	vec_o ol; // object index list
	f_iterator i,j;
	int sum,k=0,L,parity,count;
    // for each crossed face
	for (i=cf.begin();i!=cf.end();i++) {
		if (detect_polygon_edge_intersection) {
			// if ray intersected face edge
			if (ff[k++]){
				//look for another face in same object with edge crossed
				sum=0;
				L=k;
				for (j=i+1;j!=cf.end();j++) {
					// if ray intersected face edge, and faces have same object index
					if (ff[L++] && ((*i)->index==(*j)->index)){sum++;}
				}
				// compute parity of sum
				parity=sum%2;
				// if even, add instance to object list
				if (!parity) {ol.push_back((*i)->v[0]->o);}
			}
		} else {ol.push_back((*i)->v[0]->o);}
	}
	///// sort object_index_list /////
	sort(ol.begin(),ol.end());
	///// count odd objects /////
	num_odd_objects=0;
	count=ol.size();
	k=0;
	tmp.clear();
	while (k<count) {
		if (k+1!=count) {
	    	// skip identical pairs of object indices
			if (ol[k]==ol[k+1]) {k++;k++;}
			// odd object
			else {
				tmp.push_back( ol[k]);
				num_odd_objects++;
				k++;
			}
		} else { // add remaining object to odd object list
			num_odd_objects++;
			tmp.push_back( ol[k]);
			k++;
		}
	}
}

// #####################################################
// #####################################################

bool Container::facesParallel(Face *cf,Face *of){
	// are current face and other face parallel
	// i.e. is angle between normals equal to zero?
	// i.e. is the square of the cosine of the angle equal to 1?
	double cn[3],on[3];
	// get face normals
	cf->getNormal(cn);
	of->getNormal(on);
	double term1 = cn[0]*on[0]+ cn[1]*on[1]+ cn[2]*on[2];
	if (!distinguishable(term1*term1,
	(cn[0]*cn[0]+cn[1]*cn[1]+cn[2]*cn[2])*(on[0]*on[0]+on[1]*on[1]+on[2]*on[2]))) {return true;}
	else {return false;}
}

bool Container::facesColinear(Face *cf,Face *of){
	// cpvi = current_polygon_vertex_indices
	// opvi = other_polygon_vertex_indices
	double on[3],*opvc[3],*cpvc[3];
	int x=2,y=2;
	// get face vertex coordinates
	cf->getVertexCoordinates(cpvc);
	of->getVertexCoordinates(opvc);
    // get face vertex indices
    int cpvi[3] = {cf->v[0]->index,cf->v[1]->index,cf->v[2]->index};
    int opvi[3] = {of->v[0]->index,of->v[1]->index,of->v[2]->index};
	// get face normal
	of->getNormal(on);
	// dot product of other face normal and 
	// line connecting point on other face to point on current face
	// try to choose vertices not shared by both face
	if 		((opvi[0]!=cpvi[0])&&(opvi[0]!=cpvi[1])&&(opvi[0]!=cpvi[2])) {x=0;}
	else if ((opvi[1]!=cpvi[0])&&(opvi[1]!=cpvi[1])&&(opvi[2]!=cpvi[2])) {x=1;}
	if 		((cpvi[0]!=opvi[0])&&(cpvi[0]!=opvi[1])&&(cpvi[0]!=opvi[2])) {y=0;}
	else if ((cpvi[1]!=opvi[0])&&(cpvi[1]!=opvi[1])&&(cpvi[1]!=opvi[2])) {y=1;}
	// if polygons colinear, then normal and line are orthogonal
	// and dot product will be ~zero
	if (fabs(on[0]*(opvc[x][0]-cpvc[y][0])
			+on[1]*(opvc[x][1]-cpvc[y][1])
			+on[2]*(opvc[x][2]-cpvc[y][2]))<DOUBLE_EPSILON) {return true;}
	else {return false;}
}

int Container::numUniqueVertices(Face *cf,Face *of,int single_shared_vert[2]){
	// how many vertices are shared between current and other face?
	int num_unique = 0;
	// for each current face vertex
	for (int j=0;j<3;j++) {
		for (int k=0;k<3;k++) {
			if (cf->v[j]!=of->v[k]) {num_unique++;} else {single_shared_vert[0]=j;single_shared_vert[1]=k;}
		}
	}
	return num_unique;
}

bool Container::checkFaceFaceIntersections(Face *cf,Face *of) {
	// cpvc = current_polygon_vertex_coordinates
	// opvc = other_polygon_vertex_coordinates
	// cpi = current_polygon_index
	// opi = other_polygon_index
	// cn   = current_normal
	// on   = other_normal
	// get face normals
	// get face vertex coordinates
	double *opvc[3],*cpvc[3];
	cf->getVertexCoordinates(cpvc);
	of->getVertexCoordinates(opvc);
	// are faces parallel?
	bool parallel_flag = facesParallel(cf,of);
	// are faces colinear?
	// i.e. is there a vertex on other face 
	// that lies in the plane of current face?
	bool colinear_flag = facesColinear(cf,of);
	// are faces coplanar?
	bool coplanar_flag=false;
	if (parallel_flag && colinear_flag) {coplanar_flag=true;}
	// get number of unique vertices between current and other face
	int single_shared_vert[2]={-1,-1};
	int num_unique = numUniqueVertices(cf,of,single_shared_vert);
	bool share_edge_flag=false,identical_flag=false,share_vert_flag=false;
	if		(num_unique == 8) {share_vert_flag = true;}
	else if (num_unique == 7) {share_edge_flag = true;}
	else if (num_unique == 6) {identical_flag = true;}
	////////// begin decision tree //////////
	if (coplanar_flag) {
		if (identical_flag) {
			// polygons are identical
			return true;
		} else {
			//do polygon edges intersect?	
			// if yes, intersect
			// if no, do not intersect
			if (checkEdgeEdgeIntersection(cf,of,share_edge_flag)) {return true;}
			else {return false;}
		}
	} else {
		if (share_vert_flag) {
			int m=0,n=1,p=0,q=1;
			// single vertex shared
			if (single_shared_vert[0]==0){m=1;n=2;}
			else if (single_shared_vert[0]==1){n=2;}
			if (single_shared_vert[1]==0){p=1;q=2;}
			else if (single_shared_vert[1]==1){q=2;}
			double lp[2][3] = {{cpvc[m][0],cpvc[m][1],cpvc[m][2]},
								{cpvc[n][0],cpvc[n][1],cpvc[n][2]}};
			bool line_flag=false, poly_flag=false, poly_edge_flag;
			checkLineFaceIntersection(of,lp,line_flag,poly_flag,poly_edge_flag);
			// do faces intersect?
			if (line_flag && (poly_flag || poly_edge_flag)) {return(1);}
			lp[0][0] = opvc[p][0];
			lp[0][1] = opvc[p][1];
			lp[0][2] = opvc[p][2];
			lp[1][0] = opvc[q][0];
			lp[1][1] = opvc[q][1];
			lp[1][2] = opvc[q][2];
			checkLineFaceIntersection(cf,lp,line_flag,poly_flag,poly_edge_flag);
			// do faces intersect?
			if (line_flag && (poly_flag || poly_edge_flag)) {return(1);}
			else {return false;}
		}else if (!share_edge_flag) {
			// do faces intersect?
			if (checkFaceEdgeIntersection(cf,of)) {return true;}
			else {return false;}
		} else { return false;}
	}
}

// #####################################################
// #####################################################

void Face::recordBoxes(vector<Box*> &ptr){
	b.assign(ptr.begin(),ptr.end());
}

void Space::deleteBoxes(void){
	b_iterator i;
	// for each box in space, clear vector of face*
	for (i=b.begin();i!=b.end();i++) {
		delete (*i);
	}
	b.clear();
}

void Space::clearBoxes(void){
    b_iterator i;
    // for each box in space, clear vector of face*
    for (i=b.begin();i!=b.end();i++) {(*i)->f.clear();}
}

void Space::recordFace(vector<Box*> &ptr,Face* f) {
	b_iterator i;
	// for each box, add face
	for (i=ptr.begin();i!=ptr.end();i++) {
		(*i)->f.push_back(f);
	}
}

void Container::assignFacesToBoxes(Space &s) {
	o_iterator i;
	f_iterator j;
	vec_b bp;
	// clear boxes
	s.clearBoxes();
	////////// identify in which boxes each face exists ////////
	// for each object in container
	for (i=o.begin();i!=o.end();i++) {
		// for each face in object
		for (j=(*i)->f.begin();j!=(*i)->f.end();j++) {
			// identify boxes  that overlap to any degree the face bounding box.
			// This conservative approach, i.e. always including the actual 
			// face-intersecting boxes plus others is meant to save time,
			// assuming the time to exclude the other boxes is greater than the
			// time to check intersection with other boxes later.
			bp.clear();
			s.computeBoxesToCheck(*j,bp);
			// check	
			if (bp.empty()) { 
				cout << "ERROR: NO BOXES FOR\n" << "Face " << (*j)->index << " " 
					<< ((*j)->v[0])->index << " " << ((*j)->v[1])->index << " "
					<< ((*j)->v[2])->index << endl;
					exit(1);
			}
			// record face in boxes class
			s.recordFace(bp,*j);
		}
	}
}

// #####################################################
// #####################################################

void Container::getSeparationDistances(Space &s){
	cout << "Iteration 0: ";
	cout << "get separation distances..........";
	cout.flush();
	o_iterator i;
	v_iterator j;
	// for each object in container
	for (i=o.begin();i!=o.end();i++) {
        // for each vertex in object
		for (j=(*i)->v.begin();j!=(*i)->v.end();j++) {
			////////// find closest point to current vertex //////////
			// false => just add value to table, do not touch existing elements
			//			there should be none anyway
			findClosest(s,*j);
		}
	}
	cout << "complete.\n";
	cout.flush();
}

// #####################################################
// #####################################################
void Container::getBoxes(vector<Box*> &bp,Vertex *v,int offset,Space &s){
	int cbi[3],br[6];
	// compute box index that contains current_vertex
	cbi[0] = s.location2Index(v->pN[0],"x");
	cbi[1] = s.location2Index(v->pN[1],"y");
	cbi[2] = s.location2Index(v->pN[2],"z");
	// box_range
	br[0]=cbi[0]-offset;
	br[1]=cbi[0]+offset;
	br[2]=cbi[1]-offset;
	br[3]=cbi[1]+offset;
	br[4]=cbi[2]-offset;
	br[5]=cbi[2]+offset;
	// handle case where vertex lies on subspace boundary
	if (!(br[0]*s.space_length-v->pN[0])){br[0]--;}
	if (!(br[2]*s.space_length-v->pN[1])){br[2]--;}
	if (!(br[4]*s.space_length-v->pN[2])){br[4]--;}
	// screen range
	br[0]=s.screenIndex(br[0],"x");
	br[1]=s.screenIndex(br[1],"x");
	br[2]=s.screenIndex(br[2],"y");
	br[3]=s.screenIndex(br[3],"y");
	br[4]=s.screenIndex(br[4],"z");
	br[5]=s.screenIndex(br[5],"z");
	// add box pointers to vector
	s.getBoxesFor3DIndices(br,bp);
}

bool Container::faceInNeighborhood(Face *f,Vertex *v){
	// if face is in different object than vertex, then return false
	if(f->v[0]->o!=v->o){return false;}
	// else if face and vertex are in same object
	std::pair<f_iterator,vec_f::iterator> p
		=equal_range(v->nf.begin(),v->nf.end(),f);
	// if face is in current vertex neighborhood
	if(p.first!=p.second){return true;}
	return false;
}

void Container::getCandidateFaces(vector<Box*> &bp,Vertex *v,hset_f &cf){
	// for each box in search
	for (b_iterator j=bp.begin();j!=bp.end();j++) {
		// for each face in box
		for (f_iterator k=(*j)->f.begin();k!=(*j)->f.end();k++) {
			// if face vertices are not in current vertex neighborhood, face is candidate
			if (!faceInNeighborhood(*k,v)){
				cf.insert(*k);
			}
		}
	}
}

bool Container::findClosest(Space &s,Vertex *v) {
	bool gate = false;
	// declare pair
	std::pair<double,double> p;
	// get vertex normal
	double n[3];
	v->getNormal(n);
	// get Box pointers for Face* collection
	vec_b bp;
	getBoxes(bp,v,NUM_ADJACENT_BOXES,s);
	// collect Face pointers
	hset_f cf;
	getCandidateFaces(bp,v,cf);
	// if candidate faces were found
	if (!cf.empty()){
		double squareD=0.0;
		// for each candidate face
		for (hf_iterator j=cf.begin();j!=cf.end();j++) {
			// if the closest point to current vertex was found on this face
			if(computeClosest(*j,v,squareD,n)){ gate=true;}
		}
		// if no closest point was found
		if(!gate){
			// reset pointer to closest face
			v->cl=NULL;
			// set closest point to current vertex location
			for (int i=0;i<3;i++) {v->pC[i]=v->pN[i];}
		}
	}
	return gate;
}

// #####################################################
// #####################################################

bool Container::getPlaneIntersection(Face *f,Vertex *v,double *n,double num,double den,Point &p){
	//compute point on face plane that is closest to current vertex
	// i.e. intersection of plane with face normal through current vertex
    bool line_flag, poly_flag=false,poly_edge_flag;
	double u=num/den;
	double lp[2][3] = {{v->pN[0],v->pN[1],v->pN[2]},
						{v->pN[0]+u*n[0],v->pN[1]+u*n[1],v->pN[2]+u*n[2]}};
	checkLineFaceIntersection(f,lp,line_flag,poly_flag,poly_edge_flag);
	// if intersection point is on face,then save point
	if (poly_flag) {p.add(lp[1][0],lp[1][1],lp[1][2]);}
	return (poly_flag || poly_edge_flag);
}

void Container::getEdgeIntersection(Vertex *v,double *P[3],Point &p){
	double a=dot(P[0],P[0]),
			b=dot(P[1],P[1]),
			c=dot(P[2],P[2]),
			d=dot(P[0],P[1]),
			e=dot(P[1],P[2]),
			f=dot(P[2],P[0]),
			g=dot(P[0],v->pN),
			h=dot(P[1],v->pN),
			i=dot(P[2],v->pN);
    // first pair of face vertices, P[0],P[1]
    double uDen = a-2*d+b;
    if(uDen) {
        double u = (a-d-g+h)/uDen;
        if (u>0 && u<1) {
            p.add(P[0][0]+u*(P[1][0]-P[0][0]),
                    P[0][1]+u*(P[1][1]-P[0][1]),
                    P[0][2]+u*(P[1][2]-P[0][2]));
        }
    }
    // second pair of face vertices, P[1],P[2]
    uDen = b-2*e+c;
    if(uDen) {
        double u = (b-e-h+i)/uDen;
        if (u>0 && u<1) {
            p.add(P[1][0]+u*(P[2][0]-P[1][0]),
                    P[1][1]+u*(P[2][1]-P[1][1]),
                    P[1][2]+u*(P[2][2]-P[1][2]));
        }
    }
    // third pair of face vertices, P[2],P[0]
    uDen = c-2*f+a;
    if(uDen) {
        double u = (c-f-i+g)/uDen;
        if (u>0 && u<1) {
            p.add(P[2][0]+u*(P[0][0]-P[2][0]),
                    P[2][1]+u*(P[0][1]-P[2][1]),
                    P[2][2]+u*(P[0][2]-P[2][2]));
        }
    }
}

bool Container::computeClosest(Face *f,Vertex *v,double &squareD,double vn[3]) {
	bool signal = false;
	// initialize point class instance with current vertex
	Point p(v->pN[0],v->pN[1],v->pN[2]);
	// get face vertex coordinates
	double *P[3] = {&(f->v[0])->pN[0],&(f->v[1])->pN[0],&(f->v[2])->pN[0]};

	// compute vector connecting arbitrary face vertex and current vertex
	double diff[3] = {P[0][0]-v->pN[0],P[0][1]-v->pN[1],P[0][2]-v->pN[2]};
	// compute indicators
	double num=dot(vn,diff);
	// if current vertex does not lie on face plane
	bool inside = false;
	if (!(num<DOUBLE_EPSILON)) {
		//compute point on face plane that is closest to current vertex
		// i.e. intersection of plane with face normal through current vertex
		// add to p if intersection is on face
		// return true if point lies inside or on an edge of face, false otherwise
		inside = getPlaneIntersection(f,v,vn,num,dot(vn,vn),p);
	}
	
	// if a point has been found, i.e. inside==true
	// then it is necessarily the closet, we are done
	// but if a point has not been found, then keep looking
	// closest point must lie on an edge or vertex
	// gather and compare
	if(!inside){
		// add each face vertex
		for (int i=0;i<3;i++) { p.add(P[i][0],P[i][1],P[i][2]); }
		// add points of minimum distance between current vertex and each face edge
		getEdgeIntersection(v,P,p);
	}
	// save closest point
	double c[3] = { p.a , p.b , p.c };
	// invert vertex normal if vertex is not nice
	double vn_copy[3];
	if (v->o->vertexIsNice(v)){vn_copy[0]=vn[0];vn_copy[1]=vn[1];vn_copy[2]=vn[2]; } 
	else 					  {vn_copy[0]=-vn[0];vn_copy[1]=-vn[1];vn_copy[2]=-vn[2]; } 
	// compute separation vector
	double sep_vec[3] = {c[0]-v->pN[0],c[1]-v->pN[1],c[2]-v->pN[2]};
	// compute cosine of angle between outward normal and separation vector
	// which is equal to dot product of vectors divided by vector magnitudes
	double cos_angle = dot(sep_vec,vn_copy)/sqrt(dot(vn_copy,vn_copy))/sqrt(dot(sep_vec,sep_vec));
	// is closest point located within angle window as defined in controls.cc?
	if ( cos_angle > 0){
		// compute square of separation distance
		double temp=0;
		for (int i=0;i<3;i++) {temp+=(c[i]-v->pN[i])*(c[i]-v->pN[i]);}
		// closest point must be within a specified angle of vertex normal
		// to avoid grabbing points on same object
		// alternatively, the point is allowed to be within neighborhood_radius of vertex
		// since that implies the point is not near vertex on surface
		if ( (cos_angle > CLOSEST_POINT_COSINE) || temp<(NEIGHBORHOOD_RADIUS*NEIGHBORHOOD_RADIUS)){
			// if square of separation distance is closer than square of SEARCH_RADIUS
			if (temp<(SEARCH_RADIUS*SEARCH_RADIUS)) {
				// if square of separation distance is less than
				// previously saved square of separation distance
				if (temp<squareD||!squareD) {
					// save 
					for (int i=0;i<3;i++) {v->pC[i]=c[i];}
					v->cl=f;
					squareD = temp;
					signal=true;
				}
			}
		}
	}
	return signal;
}

// #####################################################
// #####################################################

bool overlap(double a1,double a2,double b1,double b2){
	// if no overlap of f1 bb (a1,a2) and f2 bb (b1,b2)
	if (a2<b1 || b2<a1){return false;}
	else {return true;}
}

bool faceBBsOverlap(double bb[6],Face *f){
	double bb2[6];

	double *p1=f->v[0]->pN;
	double *p2=f->v[1]->pN;
	double *p3=f->v[2]->pN;

	threeValueSort(*p1,*p2,*p3,bb2[1],bb2[0]);
	// if overlap in x
	if(overlap(bb[0],bb[1],bb2[0],bb2[1])){
		p1++;p2++;p3++;
		threeValueSort(*p1,*p2,*p3,bb2[3],bb2[2]);
		// if overlap in y
		if(overlap(bb[2],bb[3],bb2[2],bb2[3])){
			p1++;p2++;p3++;
			threeValueSort(*p1,*p2,*p3,bb2[5],bb2[4]);
			// if overlap in z
			if(overlap(bb[4],bb[5],bb2[4],bb2[5])){
				return true;
			}else {return false;}
		}else {return false;}
	}else {return false;}
}

bool Box::faceIntersectionAlreadyKnown(Face *a,Face *b){
	// Face a is the intersected face
	// Face b is the intersecting face
	// 
	// get iterator pointing to location of Face a in it's object list
	ff_iterator i= a->findFaceInTable_intf();
	// if Face a is in it's object list
	if(i!=a->v[0]->o->intf.end()){
		// if Face b is in Face a's vector of intersecting faces
		if(find((*(*i).second).begin(),(*(*i).second).end(),b)!=(*(*i).second).end()){
			return true;
		} else {
			return false;
		}
	} else {
		return false;
	}
}


void Box::getFaceIntersection(Container *c) {
    // for each face in box
	for(f_iterator i=f.begin();i!=f.end();i++){
    	// for each face in box
		for(f_iterator j=f.begin();j!=f.end();j++){
			// if faces are different
//			if(*i!=*j && faceIntersectionAlreadyKnown(*i,*j)==false){
			if(*i!=*j){
				// if faces intersect
				if (c->checkFaceFaceIntersections(*i,*j)) {
					// if (*i) face is not already in object hashtable of intersected faces
					if((*i)->faceInTable_intf()==false){(*i)->addFaceToTable_intf();}
					// add (*j) face to vector of faces intersecting (*i) face
					// NOTE THE RISK OF DUPLICATE FACES IN THESE VECTORS
					(*i)->addFaceToVector(*j);
				}
			}
		}
	}
}

void Controls::store(Object *o){
//	cout << "\nControls::store: area.x.size()=" << area.x.size() << endl;cout.flush();
	area.x.insert(                  area.x.end(),         o->area.x.begin(),         o->area.x.end());
//	cout << "\nControls::store: BEFORE cs.aspect_ratio.x.size()=" << aspect_ratio.x.size() << endl;cout.flush();
//	cout << "\nControls::store: BEFORE o->aspect_ratio.x.size()=" << o->aspect_ratio.x.size() << endl;cout.flush();
	aspect_ratio.x.insert(  aspect_ratio.x.end(), o->aspect_ratio.x.begin(), o->aspect_ratio.x.end());
//	cout << "\nControls::store: AFTER cs.aspect_ratio.x.size()=" << aspect_ratio.x.size() << endl;cout.flush();
//	cout << "\nControls::store: AFTER o->aspect_ratio.x.size()=" << o->aspect_ratio.x.size() << endl;cout.flush();
	edge_length.x.insert(    edge_length.x.end(),  o->edge_length.x.begin(),  o->edge_length.x.end());
	edge_angle.x.insert(      edge_angle.x.end(),   o->edge_angle.x.begin(),   o->edge_angle.x.end());
	adjacent_face.x.insert(adjacent_face.x.end(),o->adjacent_face.x.begin(),o->adjacent_face.x.end());
}

// #####################################################
// #####################################################

void Object::evalAttributes(Space &s){
	cout << "Checking if object [" << name << "] is closed..........................";
	cout.flush();
	closed=isClosed();
	cout << "complete.\n";
	cout.flush();
	cout << "Checking if object [" << name << "] is manifold........................";
	cout.flush();
	manifold=isManifold();
	cout << "complete.\n";cout.flush();
	if(manifold){
		cout << "Checking if object [" << name << "] faces are consistently oriented....";
		cout.flush();
		consistent=isConsistent();
		cout << "complete.\n";cout.flush();
		if(consistent && closed){
			cout << "Checking if object [" << name << "] faces are oriented outward.........";
			cout.flush();
			outward=isOutward(s);
			cout << "complete.\n";cout.flush();
		}
	}
}

bool Edge::isConsistent(void){
	if(f2!=NULL){ // not a border edge
		bool forward=false;
		Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
		getVertices(v1,v2,o1,o2);
		// if match v1->v2
		if( (f1->v[0]==v1 && f1->v[1]==v2 )||
			(f1->v[2]==v1 && f1->v[0]==v2 )||
			(f1->v[1]==v1 && f1->v[2]==v2 )
			){forward=true;}
		// if match v1->v2
		if( (f2->v[0]==v1 && f2->v[1]==v2 )||
			(f2->v[2]==v1 && f2->v[0]==v2 )||
			(f2->v[1]==v1 && f2->v[2]==v2 )
			){
			if(forward){ return false; }
		} else {
			if(!forward){ return false; }
		}
	} 
	return true;
}

bool Object::rayIntersectsSide(char *str,double lp[2][3],double bb[6],double n[3]){
	double x1,y1,z1,a;
	if(!strcmp(str,"-x")){
		x1=bb[0];
		a=(x1-lp[0][0])/n[0];
		y1=a*n[1]+lp[0][1];
		z1=a*n[2]+lp[0][2];
		if(a>0 && y1>bb[2] && y1<bb[3] && z1>bb[4] && z1<bb[5] ){return true;}
	} else if(!strcmp(str,"+x")){
		x1=bb[1];
		a=(x1-lp[0][0])/n[0];
		y1=a*n[1]+lp[0][1];
		z1=a*n[2]+lp[0][2];
		if(a>0 && y1>bb[2] && y1<bb[3] && z1>bb[4] && z1<bb[5] ){return true;}
	} else if(!strcmp(str,"-y")){
		y1=bb[2];
		a=(y1-lp[0][1])/n[1];
		x1=a*n[0]+lp[0][0];
		z1=a*n[2]+lp[0][2];
		if(a>0 && x1>bb[0] && x1<bb[1] && z1>bb[4] && z1<bb[5] ){return true;}
	} else if(!strcmp(str,"+y")){
		y1=bb[3];
		a=(y1-lp[0][1])/n[1];
		x1=a*n[0]+lp[0][0];
		z1=a*n[2]+lp[0][2];
		if(a>0 && x1>bb[0] && x1<bb[1] && z1>bb[4] && z1<bb[5] ){return true;}
	} else if(!strcmp(str,"-z")){
		z1=bb[4];
		a=(z1-lp[0][2])/n[2];
		x1=a*n[0]+lp[0][0];
		y1=a*n[1]+lp[0][1];
		if(a>0 && x1>bb[0] && x1<bb[1] && y1>bb[2] && y1<bb[3] ){return true;}
	} else if(!strcmp(str,"+z")){
		z1=bb[5];
		a=(z1-lp[0][2])/n[2];
		x1=a*n[0]+lp[0][0];
		y1=a*n[1]+lp[0][1];
		if(a>0 && x1>bb[0] && x1<bb[1] && y1>bb[2] && y1<bb[3] ){return true;}
	}
	return false;
}

bool Object::rayIntersectsBB(double lp[2][3],Face *ff,double n[3]){
	///// compute face bounding box /////
	vec_d xv,yv,zv;
	// identify face bounding box limits
	xv.push_back(ff->v[0]->pN[0]);
	xv.push_back(ff->v[1]->pN[0]);
	xv.push_back(ff->v[2]->pN[0]);
	yv.push_back(ff->v[0]->pN[1]);
	yv.push_back(ff->v[1]->pN[1]);
	yv.push_back(ff->v[2]->pN[1]);
	zv.push_back(ff->v[0]->pN[2]);
	zv.push_back(ff->v[1]->pN[2]);
	zv.push_back(ff->v[2]->pN[2]);
	sort(xv.begin(),xv.end());
	sort(yv.begin(),yv.end());
	sort(zv.begin(),zv.end());
	// grab face 3D location range
	double bb[6];
	bb[0] = xv[0];  // -x
	bb[1] = xv[2];	//  +x
	bb[2] = yv[0];	// -y
	bb[3] = yv[2];	//  +y
	bb[4] = zv[0];  // -z
	bb[5] = zv[2];	//  +z

	///// for each side of bounding box /////
	if( rayIntersectsSide("-x",lp,bb,n) || 
		rayIntersectsSide("+x",lp,bb,n) || 
		rayIntersectsSide("-y",lp,bb,n) || 
		rayIntersectsSide("+y",lp,bb,n) || 
		rayIntersectsSide("-z",lp,bb,n) || 
		rayIntersectsSide("+z",lp,bb,n)
	){return true;}
		// if ray intersects side, then return true
	return false;
}

bool Object::isOutward(Space &s){
	// assuming object mesh is closed, manifold, and consistent...
	//
	bool line_flag=false, poly_flag=false, poly_edge_flag=false;
	f_iterator fff=f.begin();
	double count;
	do{
		// set origin face
		Face *ff = *fff;
		count=0;
		// compute normal of origin face
		double n[3];
		ff->getNormal(n);
		double L=sqrt( dot(n,n) );
		n[0]=n[0]/L;
		n[1]=n[1]/L;
		n[2]=n[2]/L;
		// ray origin = centroid of origin face
		double lp[2][3];
		lp[0][0] = (ff->v[0]->pN[0]+
					ff->v[1]->pN[0]+
					ff->v[2]->pN[0])/3.0;
		lp[0][1] = (ff->v[0]->pN[1]+
					ff->v[1]->pN[1]+
					ff->v[2]->pN[1])/3.0;
		lp[0][2] = (ff->v[0]->pN[2]+
					ff->v[1]->pN[2]+
					ff->v[2]->pN[2])/3.0;
		// ray end = point on normal advanced from origin
		//			 a distance equal to 2*world width along each axis
		double del[3]={fabs(s.world[1]-s.world[0]),fabs(s.world[3]-s.world[2]),fabs(s.world[5]-s.world[4])};
		int big;
		biggest(del,big);
		lp[1][0] = lp[0][0]+n[0]*2*del[big];
		lp[1][1] = lp[0][1]+n[1]*2*del[big];
		lp[1][2] = lp[0][2]+n[2]*2*del[big];
/*		// DEBUG
		cout.precision(12);
		cout << "\n\nObject::isOutward: "
		<< "ray  ["
		<< lp[0][0] << " "
		<< lp[0][1] << " "
		<< lp[0][2] << "]\n"
		<< "Object::isOutward: "
		<< "ray  ["
		<< lp[1][0] << " "
		<< lp[1][1] << " "
		<< lp[1][2] << "]\n";
*/ 		// DEBUG
		// for each face in object
		for(f_iterator i=f.begin();i!=f.end();i++){
			// if face in Object not same as origin face 
			if(*i!=ff){
				// if ray intersects face bounding box
				if(rayIntersectsBB(lp,*i,n)){
/*					// DEBUG
					cout << "Object::isOutward: "
					<< "ray intersects face BB:\n";
					(*i)->printFace((*i)->v[0]->o);
					cout << endl;
*/					// DEBUG
					// if hit, determine face-line intersection
					checkLineFaceIntersection(*i,lp,line_flag,poly_flag,poly_edge_flag);
					// check
					if(poly_edge_flag==true) {
						fff++;
						if(fff==f.end()){
							cout << "Object::isOutward: "
							<< "Error. Ray-tracing failed from all object <"
							<< name << "> faces.\n";
							exit(0);
						}
						break;
					}
					// does point intersect polygon
					if (poly_flag && line_flag) {count++;}
/*					// DEBUG
					cout << "Object::isOutward: "
					<< "poly_flag=" << poly_flag
					<< ", poly_edge_flag=" << poly_edge_flag
					<< ", count=" << count
					<< endl;
*/					// DEBUG
				}
			}
		}
	}while(poly_edge_flag==true);

	count = ceil(count);
	// if odd hit count then inward normal
/*	// DEBUG
	cout << "Object::isOutward: "
	<< "count=" << count
	<< ", static_cast<int>(count)%2="
	<< static_cast<int>(count)%2 << endl << endl; 
*/	// DEBUG
	if(static_cast<int>(count)%2){return false;}
	// if even hit count then outward normal
	else {return true;}
}

bool Object::isConsistent(void){
	// assuming object mesh is manifold...
	bool flag=true;
	// for each edge in object
	for(e_iterator i=e.begin();i!=e.end();i++){
		// if the one or two adjacent edge faces traverse the edge
		// in the same direction, either both v1->v2 or both v2->v1,
		// then the entire set of mesh faces is declared inconsistent
		// i.e. the first edge with a flipped face will trigger signal return
		// and no more edge checking is performed, but -p option
		// needs ALL flipped edges.
		if(!(*i)->isConsistent()){
			flag=false;
			// record offending edge
			flipped.push_back(*i);
		}
	}
	// update cumulative # flipped edges
//	cs.num_flip=cs.flipped.size();
	return flag;
}

bool Object::isClosed(void){
	bool flag=true;
	// for each edge in object
	for(e_iterator i=e.begin();i!=e.end();i++){
		if((*i)->f2==NULL){
			flag=false;
			// record offending edge
			border.push_back(*i);
		}
	}
	// update cumulative # border edges
//	cs.num_bor+=cs.border.size();
	return flag;
}

bool Edge::getStartingFace(Face* &sf){
	// if edge is manifold
	if(isManifold()==true){
		// grab starting face
		if(f1!=NULL){
			sf = f1;
		} else if (f2!=NULL) {
			sf = f2;
		} else {
			cout << "Error. Both edge faces are NULL!\n"; exit(0);
		}
		return true;
	} else {
		return false;
	}
}

Face* Edge::getNewFace(Face *old){
	// return the adjacent face that is
	// different from old adjacent face
	if(old!=f1 && old!=f2){
		cout << "Error. Neither edge adjacent face matches old face.\n";
		exit(0);
	}
	if(old==f1){return f2;}
	else {return f1;}
}

Edge* Face::getNewEdge(Edge *old,Vertex *vv){
	// return face edge that contains vertex vv
	// and is different from old edge
	if((e[0]!=old && e[1]!=old && e[2]!=old) ||
		(v[0]!=vv && v[1]!=vv && v[2]!=vv)){
		cout << "\n\nOld edge does not match any on face.\n";
		cout << "Current face:\n";
		printFace(v[0]->o);
		cout << endl;
		cout << "Old edge:\n";
		old->printEdge(old->f1->v[0]->o->name);
		cout << endl;
		cout << "Current vertex:\n";
		vv->printVertex(vv->o->name);
		cout << endl << endl;
		
		exit(0);
	}
	Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
	e[0]->getVertices(v1,v2,o1,o2);
	if(e[0]!=old && (v1==vv || v2==vv)){return e[0];}
	e[1]->getVertices(v1,v2,o1,o2);
	if(e[1]!=old && (v1==vv || v2==vv)){return e[1];}
	e[2]->getVertices(v1,v2,o1,o2);
	if(e[2]!=old && (v1==vv || v2==vv)){return e[2];}
	cout << "\n\nFace::getNewEdge: Error. No edge on face contains current vertex \n"
	<< "and is different from old vertex.\n\n";
	cout << "current vertex:\n";
	vv->printVertex(vv->o->name);
	cout << endl << endl;
	cout << "old edge:\n";
	old->printEdge(old->vv1->o->name);
	cout << endl << endl;
	cout << "current face:\n";
	printFace(v[0]->o);
	cout << endl << endl;
	exit(0);
}

bool Vertex::scanAdjFaces(Edge *se,Face *sf,bool &nonman){
	// collect touched edges
	vec_e te;
	te.push_back(se);
	// collect touched faces
	vec_f tf;
	tf.push_back(sf);
	// get new edge
	se=sf->getNewEdge(se,this);
	if (se->isManifold()==false){
		// vertex manifoldness cannot be determined
		// because it's complicated, so just return flag
		nonman=true;
		return false;
	}
	// while new edge has not already been touched
	while(find(te.begin(),te.end(),se)==te.end()){
		// keep new edge
		te.push_back(se);
		// get new face
		sf=se->getNewFace(sf);
		if (sf==NULL) {break;}
		// keep new face
		tf.push_back(sf);
		// get new edge
		se=sf->getNewEdge(se,this);
		if (se->isManifold()==false){
			// vertex manifoldness cannot be determined
			// because it's complicated, so just return flag
			nonman=true;
			return false;
		}
	}
	// if number of touched faces != number of vertex adjacent faces
	return tf.size()==f.size();
}

bool Vertex::isManifold(bool flag){
//	cout << "Vertex::isManifold: @start: flag=" << flag << endl;
	// try to "walk" around vertex adjacent faces
	// using adjacent face edges
	// if an edge is nonmanifold then abort mission
	// stop when the walk returns to starting adjacent face
	if(e.empty()==true){
		cout << "\n\nVertex::isManifold: Error."
		<< " Vertex was not an 'orphan', but has no edges.\n";
		cout << "\n\nVertex::isManifold: Confused vertex:\n";
		printVertex(o->name);
		cout << endl;
		exit(0);
	}
	// grab starting edge
	Edge *se=e.front();
	// grab starting face
	Face *sf=NULL,*cw=NULL,*ccw=NULL;
	// if edge is manifold
	if(se->isManifold()==true){
		// grab starting face
		if(se->f1!=NULL){ cw = se->f1; }
		if(se->f2!=NULL) { 
			if(cw==NULL) { cw = se->f2; }
			else { ccw = se->f2; }
		}
		if(cw==NULL){
			cout << "Error. Both edge faces are NULL!\n"; exit(0);
		}
	} else {
		// vertex manifoldness cannot be determined
		// because the starting edge is not manifold
		// thus determining manifoldness of vertex
		// is complicated, so just return flag
		return flag;
		// TODO: report number of vertices for which manifoldness
		// was not determined and alert user
	}

	// try clockwise
	sf=cw;
	bool nonman=false;
	// if fail i.e. not all adjacent faces touched
	if(scanAdjFaces(se,sf,nonman)==false){
		// if nonmanifold edge found, then bail
		if(nonman==true){return flag;}
		// if ccw face is not NULL
		if(ccw!=NULL){
			// try counter-clockwise
			sf=ccw;
			// if fail i.e. all adjacent faces still not touched
			if(scanAdjFaces(se,sf,nonman)==false){
				// record offending edge
				o->nonman_v.push_back(this);
				return false;
			} else {
				return flag;
			}
		} else {
			// record offending edge
			o->nonman_v.push_back(this);
			return false;
		}
	} else {
		// all adjacent faces touched
		return flag;
	}

/*	if(se->getStartingFace(sf)==false){
		// vertex manifoldness cannot be determined
		// because the starting edge is not manifold
		// thus determining manifoldness of vertex
		// is complicated, so just return flag
		return flag;
		// TODO: report number of vertices for which manifoldness
		// was not determined and alert user
	}
*/
}

bool Object::verticesManifold(bool flag){
	///// if vertices are manifold /////
	////// confirm that all adjacent faces /////
	///// are consecutively reachable by edge hops /////

	// sort orphan vertices
	sort(orphan.begin(),orphan.end());

	// for each vertex in object
	for(v_iterator i=v.begin();i!=v.end();i++){
		// if vertex is not an orphan
		if(binary_search(orphan.begin(),orphan.end(),*i)==false){
			flag = (*i)->isManifold(flag);
		}
	}
	return flag;
}

bool Edge::isManifold(void){
	// if three or more faces share an edge then NOT manifold
	// fvec stores adjacent faces beyond first and second faces
	// so if fvec is empty, then edge is manifold
	return fvec.empty()==true;
}

bool Object::edgesManifold(bool flag){
	// for each edge in object
	for(e_iterator i=e.begin();i!=e.end();i++){
		// if edge is NOT manifold
		if((*i)->isManifold()==false){
			flag=false;
			// record offending edge
			nonman_e.push_back(*i);
		}
	}
	return flag;
}

bool Object::isManifold(void){
	bool flag=true;
	// check edge manifoldness
	flag=edgesManifold(flag);
	// check vertex manifoldness
	flag=verticesManifold(flag);

	// update cumulative # nonmanifold edges and vertices
//	cs.num_nonman_e+=cs.nonman_e.size();
//	cs.num_nonman_v+=cs.nonman_v.size();
	return flag;
}

//void Controls::parse(int argc,char **argv,std::string message,char filename[FILENAME_SIZE]){
void Controls::parse(int argc,char **argv,std::string message){

	// if no arguments passed
	if(argc==1){
		cout << message << endl;
		exit(0);
	}

	int c;
	opterr=0;
	char *eptr=NULL;
	while((c=getopt(argc,argv,"ab:c:d:e:f:hipqz:")) != -1)
		switch(c)
			{
			case 'a':
				// determine attributes only
				attr=true;
				break;
			case 'b':
				// specify max edge length threshold (units of input file)
				signal[4]=true;
				thresholds[4]= strtod(optarg,&eptr);
				break;
			case 'c':
				// specify min edge length threshold (units of input file)
				signal[3]=true;
				thresholds[3]= strtod(optarg,&eptr);
				break;
			case 'd':
				// specify max edge angle threshold (degrees)
				signal[2]=true;
				thresholds[2]= strtod(optarg,&eptr);
				break;
			case 'e':
				// specify min edge angle threshold (degrees)
				signal[1]=true;
				thresholds[1]= strtod(optarg,&eptr);
				break;
			case 'f':
				// specify max aspect ratio threshold
				signal[0]=true;
				thresholds[0]= strtod(optarg,&eptr);
				break;
			case 'h':
				cout << message << endl;
				abort ();
			case 'i':
				interf=true;
				break;
			case 'p':
				print=true;
				sprintf(style,"detail");
				break;
			case 'q':
				print=true;
				sprintf(style,"cp");
				break;
			case 'z':
				zEM=strtod(optarg,&eptr);
				break;
			case '?':
				if (optopt == 'c')
					fprintf (stderr,"Option -%c requires an argument.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Unknown option '-%c.\n", optopt);
				else
					fprintf (stderr,
							"Unknown option character '\\x%x'.\n",
							optopt);
				exit(0);
			default:
				cout << message << endl;
				abort ();
			}

	if(interf==true && folder==false){
		cout << "\n\nError. The '-i' option applies only to folders.\n\n";
		exit(0);
	}

	if(optind<argc){
		for (int index=optind;index<argc;index++){
			// determine if file or folder was passed
			struct stat stat_p;
			stat(argv[index],&stat_p);
			if(stat_p.st_mode & S_IFDIR){ folder=true;}
			else{ folder=false;}
	
			// adjust input directory
			char filename[FILENAME_SIZE];
			strcpy(filename,argv[index]);
			if(folder){
		    	char *temp=strrchr(filename,'/');
				if(!temp) {strcat(filename,"/");}
				else if(*++temp) {strcat(filename,"/");}
			}
			inpath=filename;
		}
	} else {
		fprintf (stderr,"No input file argument found on command line.\n");
		exit(0);
	}
}

// ********************************************************
// ********************************************************

void Object::checkIntegrity(void){
	// check vertex integrity
	cout << "\nChecking vertex integrity for object [" << name << "]..................";
	cout.flush();
	orphanMissingContig();
	cout << "complete.\n";cout.flush();
	// check face integrity
	cout << "Checking face integrity for object [" << name << "]....................";
	cout.flush();
	degenContig();
	cout << "complete.\n";cout.flush();
}

bool Object::goodIntegrity(void){
	// integrity is good if there are
	// (1) no missing vertices, face references nonexistent vertex index
	// (2) no orphan vertices, vertices referenced by no face
	// (3) no degenerate face, face references same vertex more than once
	// Note: missing data, such as a face with only two vertex indices given
	// 	would have been caught during input file parsing.
	return ((missing_v.empty()==true) &&
			(orphan.empty()==true) &&
			(degen.empty()==true) &&
			(dupl_v_index.empty()==true) &&
			(dupl_f_index.empty()==true));
}

// VERTICES
void Object::orphanMissingContig(void){
	// keep unique vertices
	sort(missing_v.begin(),missing_v.end());
	i_iterator new_end_v = unique(missing_v.begin(),missing_v.end());
	missing_v.assign(missing_v.begin(),new_end_v);
	// keep unique faces
	sort(missing_f.begin(),missing_f.end());
	f_iterator new_end_f = unique(missing_f.begin(),missing_f.end());
	missing_f.assign(missing_f.begin(),new_end_f);

	// update cumulative # missing vertices
//	num_missing+=missing_v.size();

	///// orphan ////
	std::pair<ib_iterator,ib_iterator> pp;
	std::pair<iv_iterator,iv_iterator> qq;
	// for each element of map #1
	ib_iterator i=found.begin();
	while(i!=found.end()){
		// if this vertex was not used in any face
		if((*i).second==false){
			// find all elements in map with this key
			pp=found.equal_range((*i).first);
			// if no matching elements were found
			if(pp.first==pp.second){ cout << "Error. Weird. how can this be?\n";exit(0);}
			else { // orphan vertex (or vertices)
				// find matching elements in map #2
				qq=vp.equal_range((*i).first);
				// if no matching elements were found
				if(qq.first==qq.second){ cout << "Error. Weird. how can this be too?\n";exit(0);}
				else {
					// for each matching element
					for(iv_iterator n=qq.first;n!=qq.second;n++){
						// add vertex* to orphan list
						orphan.push_back((*n).second);
					}
				}
				i=pp.second;
			}
		} else {i++;}
	}

	// update cumulative # orphan vertices
//	num_orphan+=orphan.size();

	///// contiguousness /////
	contig_v=true;
	// since map #1 is automatically sorted by integer value
	int pace = 1;
	// for each other element in map
	for(i=found.begin();i!=found.end();i++){
		// if vertex indexing deviates 
		if((*i).first!=pace++){
			pace--;
			contig_v=false;
			// save five vertex indices to vec_cv
			if (pace==1){
				vec_cv[0]=vec_cv[1]=-1;
				vec_cv[2]=(*i).first; i++;
				vec_cv[3]=(*i).first; i++;
				vec_cv[4]=(*i).first;
				i--;i--;
			} else if (pace==2) {
				vec_cv[0]=-1;
				i--;
				vec_cv[1]=(*i).first; i++;
				vec_cv[2]=(*i).first; i++;
				vec_cv[3]=(*i).first; i++;
				vec_cv[4]=(*i).first;
				i--;i--;
			} else if (pace==static_cast<int>(v.size())-1) {
				i--;i--;
				vec_cv[0]=(*i).first; i++;
				vec_cv[1]=(*i).first; i++;
				vec_cv[2]=(*i).first; i++;
				vec_cv[3]=(*i).first; i++;
				vec_cv[4]=-1;
				i--;i--;
			} else if (pace==static_cast<int>(v.size())) {
				i--;i--;
				vec_cv[0]=(*i).first; i++;
				vec_cv[1]=(*i).first; i++;
				vec_cv[2]=(*i).first; i++;
				vec_cv[3]=vec_cv[4]=-1;
				i--;i--;
			} else {
				i--;i--;
				vec_cv[0]=(*i).first; i++;
				vec_cv[1]=(*i).first; i++;
				vec_cv[2]=(*i).first; i++;
				vec_cv[3]=(*i).first; i++;
				vec_cv[4]=(*i).first;
				i--;i--;
			}
			break;
		}
	}

	///// contiguousness /////
	// for each element in multimap(int->Vertex*)
	for(iv_iterator k=vp.begin();k!=vp.end();k++){
		// find # elements in map with this key
		if(vp.count((*k).first)>1){
			dupl_v_index.push_back((*k).second);
		}
	}
}

void Object::vertexAdjacentFaces(void){
//	mmap_iv af;
	// for each vertex in object
	for(v_iterator i=v.begin();i!=v.end();i++){
		int c=(*i)->f.size();
//		af.insert(std::make_pair(c,*i));
//		adjacent_face.n++;
		adjacent_face.sum+=c;
		adjacent_face.sum2+=c*c;
		adjacent_face.total+=c;
		if(c<adjacent_face.min) {adjacent_face.min=c;}
		if(c>adjacent_face.max) {adjacent_face.max=c;}
		// add to vector
		adjacent_face.x.push_back(c);
	}
	// build adjacent face histogram
	adjacent_face.createAdjacentFaceHistogram();
}

void Controls::vertexAdjacentFaces(void){
//	mmap_iv af;
	// for each element in vector
	for(d_iterator i=adjacent_face.x.begin();i!=adjacent_face.x.end();i++){
//		int c=(*i)->f.size();
//		af.insert(std::make_pair(c,*i));
//		adjacent_face.n++;
		adjacent_face.sum+=*i;
		adjacent_face.sum2+=(*i)*(*i);
		adjacent_face.total+=*i;
		if(*i<adjacent_face.min) {adjacent_face.min=*i;}
		if(*i>adjacent_face.max) {adjacent_face.max=*i;}
		// add to vector
//		adjacent_face.x.push_back(c);
	}
	// build adjacent face histogram
	adjacent_face.createAdjacentFaceHistogram();
}

// FACES
void Object::areaAspectRatio(Controls &cs){
	// for each face in object
	for(f_iterator i=f.begin();i!=f.end();i++){
//		area.n++;
//		aspect_ratio.n++;
		///// face area /////
		// compute face normal vector
		double n[3];
		(*i)->getNormal(n);
		// compute face area = half normal vector length
		double aa=sqrt(dot(n,n))/2.0;
		//
		area.sum+=aa;
		area.sum2+=aa*aa;
		area.total+=aa;
		if(aa<area.min) {area.min=aa;}
		if(aa>area.max) {area.max=aa;}
		// add to vector
		area.x.push_back(aa);
		///// aspect ratio /////
		double ar = (*i)->getAspectRatio();
		aspect_ratio.sum+=ar;
		aspect_ratio.sum2+=ar*ar;
		aspect_ratio.total+=ar;
		if(ar<aspect_ratio.min) aspect_ratio.min=ar;
		if(ar>aspect_ratio.max) aspect_ratio.max=ar;
		if(cs.signal[0]==true){
			// compare face aspect ratio to user-defined threshold
			if(ar>cs.thresholds[0]) bad_aspect[*i]=ar;
		}
		// add to vector
		aspect_ratio.x.push_back(ar);
	}
	// update cumulative surface area
//	a+=area.sum;
	// build face area histogram
	area.createHistogram();
	// build aspect ratio histogram
	aspect_ratio.createAspectRatioHistogram();
}

void Controls::areaAspectRatio(void){
	// for each element in area vector
	for(d_iterator i=area.x.begin();i!=area.x.end();i++){
		area.sum+=*i;
		area.sum2+=(*i)*(*i);
		area.total+=*i;
		if(*i<area.min) {area.min=*i;}
		if(*i>area.max) {area.max=*i;}
	}

	// for each element in aspect ratio vector
	for(d_iterator i=aspect_ratio.x.begin();i!=aspect_ratio.x.end();i++){
		aspect_ratio.sum+=*i;
		aspect_ratio.sum2+=(*i)*(*i);
		aspect_ratio.total+=*i;
		if(*i<aspect_ratio.min) aspect_ratio.min=*i;
		if(*i>aspect_ratio.max) aspect_ratio.max=*i;
	}
	// build face area histogram
	area.createHistogram();
	// build aspect ratio histogram
	aspect_ratio.createAspectRatioHistogram();
}

double Face::getAspectRatio(void){

	/* Make triangle edge vectors */
	double va[3]={v[1]->pN[0]-v[0]->pN[0],v[1]->pN[1]-v[0]->pN[1],v[1]->pN[2]-v[0]->pN[2]};
	double vb[3]={v[2]->pN[0]-v[1]->pN[0],v[2]->pN[1]-v[1]->pN[1],v[2]->pN[2]-v[1]->pN[2]};
	double vc[3]={v[0]->pN[0]-v[2]->pN[0],v[0]->pN[1]-v[2]->pN[1],v[0]->pN[2]-v[2]->pN[2]};
	double vbase[3]={0,0,0};
	double vopp[3]={0,0,0};

  /* Find length of longest edge */
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
	vc[0]=v[2]->pN[0]-v[0]->pN[0];
	vc[1]=v[2]->pN[1]-v[0]->pN[1];
	vc[2]=v[2]->pN[2]-v[0]->pN[2];
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
	va[0]=v[0]->pN[0]-v[1]->pN[0];
	va[1]=v[0]->pN[1]-v[1]->pN[1];
	va[2]=v[0]->pN[2]-v[1]->pN[2];
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
	vb[0]=v[1]->pN[0]-v[2]->pN[0];
	vb[1]=v[1]->pN[1]-v[2]->pN[1];
	vb[2]=v[1]->pN[2]-v[2]->pN[2];
    vopp[0]=vb[0];
    vopp[1]=vb[1];
    vopp[2]=vb[2];
  }

  /* Find shortest altitude */
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

void Object::degenContig(void){
	// buide map: face index (integer)->face*
	mmap_if tt;
	// for each face in object
	for(f_iterator i=f.begin();i!=f.end();i++){
		tt.insert(std::make_pair((*i)->index,*i));
	}
	contig_f=true;
	int pace=1;
	// for each face in map
	for(if_iterator i=tt.begin();i!=tt.end();i++){
		Face *ff=(*i).second;
		///// duplicity /////
		// find # elements in map with this key
		if(tt.count((*i).first)>1){
			dupl_f_index.push_back(ff);
		}
		///// degeneracy /////
		// if any two of the vertex indices are the same
		if( ff->v[0]==ff->v[1] || ff->v[0]==ff->v[2] || ff->v[1]==ff->v[2]){
			// then face is degenerate
			degen.push_back(ff);
//			num_degen++;
		}
		// update cumulative # degenerate faces
		///// contiguousness /////
		if(contig_f){
			if((*i).first!=pace++){
				// face indexing NOT contiguous
				contig_f=false;
				pace--;
				// save five face indices to vec_cf
				if (pace==1){
					vec_cf[0]=vec_cf[1]=-1;
					vec_cf[2]=(*i).first; i++;
					vec_cf[3]=(*i).first; i++;
					vec_cf[4]=(*i).first;
					i--;i--;
				} else if (pace==2) {
					vec_cf[0]=-1;
					i--;
					vec_cf[1]=(*i).first; i++;
					vec_cf[2]=(*i).first; i++;
					vec_cf[3]=(*i).first; i++;
					vec_cf[4]=(*i).first;
					i--;i--;
				} else if (pace==static_cast<int>(f.size())-1) {
					i--;i--;
					vec_cf[0]=(*i).first; i++;
					vec_cf[1]=(*i).first; i++;
					vec_cf[2]=(*i).first; i++;
					vec_cf[3]=(*i).first; i++;
					vec_cf[4]=-1;
					i--;i--;
				} else if (pace==static_cast<int>(f.size())) {
					i--;i--;
					vec_cf[0]=(*i).first; i++;
					vec_cf[1]=(*i).first; i++;
					vec_cf[2]=(*i).first; i++;
					vec_cf[3]=vec_cf[4]=-1;
					i--;i--;
				} else {
					i--;i--;
					vec_cf[0]=(*i).first; i++;
					vec_cf[1]=(*i).first; i++;
					vec_cf[2]=(*i).first; i++;
					vec_cf[3]=(*i).first; i++;
					vec_cf[4]=(*i).first;
					i--;i--;
				}
			}
		}
	}
}

void Object::setAll(Vertex *vv,hmap_fi &group,int &g){
	// for each vertex adjacent face
	for(f_iterator i=vv->f.begin();i!=vv->f.end();i++){
		group[*i]=g;
	}
	g++;
}

void Object::getGroups(Vertex *vv,hmap_fi &group,set_i &s){
	s.clear();
	// for each vertex adjacent face
	for(f_iterator i=vv->f.begin();i!=vv->f.end();i++){
		s.insert(group[*i]);
	}
}

int Object::getLowest(set_i &s){
	int i = 100000;
	// for each element in set
	for(std::set<int,lti>::iterator j=s.begin();j!=s.end();j++){
		// if element is lower than i and not 0
		if(*j<i && *j){i=*j;}
	}
	return i;
}

void Object::replaceGroups(Vertex *vv,hmap_fi &group,int z){
	// for each vertex adjacent face
	for(f_iterator i=vv->f.begin();i!=vv->f.end();i++){
		group[*i]=z;
	}
}

void Object::setZero(Edge *ee,hmap_fi &group,int z){
	fi_iterator i=group.find(ee->f1);
	if((*i).second==0){group[ee->f1]=z;}
	// if second adjacent face
	if(ee->f2!=NULL){
		i=group.find(ee->f2);
		if((*i).second==0){group[ee->f2]=z;}
	}
	// if more adjacent faces
	if(!ee->fvec.empty()){
		// for each adjacent face
		for(f_iterator j=ee->fvec.begin();j!=ee->fvec.end();j++){
			// if adjacent face has no group
			i=group.find(*j);
			if((*i).second==0){group[*j]=z;}
		}
	}
}

bool Object::removeSelectedFaces(vec_f &sf,vec_f &fv){
	bool flag=false;
	// for each selected face
	f_iterator i=sf.begin();
	while(i!=sf.end()){
		f_iterator j= find(fv.begin(),fv.end(),*i);
		// if selected face found in face vector
		if(j!=fv.end()){
			// remove face from face vector
			fv.erase(j);
			flag=true;
			i++;
		} else {
			// remove face from selected face vector
			//  so that the face is not used next round
			sf.erase(i);
		}
	}
	return flag;
}

void Object::getSelectedFaces(vec_f &sf,vec_f fv){
	// if selected faces vector is empty
	if(sf.empty()==true){
		// then grab first face from vector
		sf.push_back(fv.front());
	} else {
		// adjacent faces
		vec_f af;
		// for each selected face
		for(f_iterator i=sf.begin();i!=sf.end();i++){
			// grab adjacent faces of selected face
			if((*i)->e[0]->f1!=*i){af.push_back((*i)->e[0]->f1);}
			else 				  {af.push_back((*i)->e[0]->f2);}
			if((*i)->e[1]->f1!=*i){af.push_back((*i)->e[1]->f1);}
			else 				  {af.push_back((*i)->e[1]->f2);}
			if((*i)->e[2]->f1!=*i){af.push_back((*i)->e[2]->f1);}
			else 				  {af.push_back((*i)->e[2]->f2);}
		}
		// keep unique faces
		sort(af.begin(),af.end());
		f_iterator i=unique(af.begin(),af.end());
		sf.assign(af.begin(),i);
	}
}

int Object::countComponents(void){
	// initialize group number
	int g=0;
	///// create hashed map face*->integer (group #) /////
	hmap_fi group;
	// for each face
	for(f_iterator i=f.begin();i!=f.end();i++){
		group[*i]=g;
	}
	// increment group number
	g++;
	// set of integers
	set_i s;
	bool changes = true;
	while(changes==true){
		changes=false;
		// for each vertex in object
		for(v_iterator i=v.begin();i!=v.end();i++){
			// get vertex adjacent face groups
			getGroups(*i,group,s);
			std::set<int,lti>::iterator j=s.begin();
			// if no adjacent face group has been set
			if(s.size()==1 && *j==0){
				// set all adjacent faces to next available group #
				setAll(*i,group,g);
				changes=true;
			} else if(s.size()>1){ // more than one group
				// identify lowest group # larger than 0 in set
				int z=getLowest(s);
				// replace the larger group # with lowest group # in all faces
				replaceGroups(*i,group,z);
				changes=true;
			} 
		}
	}
	// analyze map
	s.clear();
	// for each element in group
	for(fi_iterator i=group.begin();i!=group.end();i++){
		s.insert((*i).second);
	}
	return s.size();
}

void Object::findIntersectingFaces(Container *c,Space &s){
	vec_f dummy;
	// DEBUG
//	cout << "\nObject::findIntersectingFaces: "
//	<< "# boxes = " << s.b.size() << endl;
//	int booyah = s.b.size();
//	int mofo = 1;
	// DEBUG
	// for each box
	for(b_iterator i=s.b.begin();i!=s.b.end();i++){
		// if box is not empty
		if((*i)->f.empty()==false){
			// check for intersections
		   	(*i)->getFaceIntersection(c);
			// intersecting faces end up in o->intf (hashtable_f_face)
			// DEBUG
//			cout << "Object::findIntersectingFaces: "
//			<< "box " << mofo++ << " of " << booyah << endl;
			// DEBUG
		}
	}
}


// EDGES
void Object::countBoundaries(void){
	num_bou=0;
	Boundary bb;
	// keep unique border edges
	sort(border.begin(),border.end());
	e_iterator new_end = unique(border.begin(),border.end());
	border.assign(border.begin(),new_end);
	// copy the border vector to ws
	vec_e ws;
	ws.assign(border.begin(),border.end());
	// while there are border edges in ws
	while(ws.empty()==false){
		// if the boundary is closed
		if(bb.open==false){
			// init boundary
			bb.init(ws.back());
			// pop last edge from ws
			ws.pop_back();
		}
		bool myfound = false;
		// for each edge in ws
		e_iterator i=ws.begin();
		while(i!=ws.end()){
			// if edge extends boundary
			if(bb.edgeExtendsBoundary(*i)){
				myfound = true;
				// remove edge from ws
				ws.erase(i);
				// if boundary is closed
				if(bb.closed()){
					// increment boundary count
					num_bou++;
					// break
					break;
				}
			} else { i++;}
		}
		// if no edges added to boundary then error
		if(myfound==false){
			cout << "Error: no edges in ws were added to boundary.\n";
			bb.print();
			exit(0);
		}
	}
	// update cumulative number boundaries
//	cs.num_bou+=num_bou;
}

void Object::processEdgeLengths(Controls &cs){
	// for each edge in object
	for(e_iterator i=e.begin();i!=e.end();i++){
		double l = (*i)->l;
//		edge_length.n++;
		edge_length.sum+=l;
		edge_length.sum2+=l*l;
		edge_length.total+=l;
		if(l<edge_length.min) edge_length.min=l;
		if(l>edge_length.max) edge_length.max=l;
		if(cs.signal[3]==true){
			if(l<cs.thresholds[3]) bad_length[*i]=l;
		}
		if(cs.signal[4]==true){
			if(l>cs.thresholds[4]) bad_length[*i]=l;
		}
		// check distinguishability
		Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
		(*i)->getVertices(v1,v2,o1,o2);
		if( !distinguishable(v1->pN[0],v2->pN[0]) &&
			!distinguishable(v1->pN[1],v2->pN[1]) &&
			!distinguishable(v1->pN[2],v2->pN[2]) 
			){ indistin.push_back(*i); }
		// add to vector
		edge_length.x.push_back(l);
	}
	edge_length.createHistogram();
}

void Controls::processEdgeLengths(void){
	// for each element in area vector
	for(d_iterator i=edge_length.x.begin();i!=edge_length.x.end();i++){
//		double l = (*i)->l;
//		edge_length.n++;
		edge_length.sum+=*i;
		edge_length.sum2+=(*i)*(*i);
		edge_length.total+=*i;
		if(*i<edge_length.min) edge_length.min=*i;
		if(*i>edge_length.max) edge_length.max=*i;
	}
	edge_length.createHistogram();
}

void Object::computeGenus(void){
	int num = v.size()-e.size()+f.size();
	if (num%2) {
		cout << "Round off error in genus computation. "
		<< "(#vertices + #faces - #edges) is not evenly divisible by 2.\n";
		cout << "v-e+f = " << num << endl;
		cout << v.size() << "-" << e.size() << "+" << f.size() << endl;
		exit(0);
	} else {
		genus = num_sep-num/2;
	}
}

void Object::computeEdgeAngles(Controls &cs){
	// for each edge in object
	for(e_iterator i=e.begin();i!=e.end();i++){
		// if edge has exactly two adjacent faces
		if((*i)->f2!=NULL && (*i)->fvec.empty()){
			double angle = (*i)->getAngle()*180/PI; // degrees
//			edge_angle.n++;
			edge_angle.sum+=angle;
			edge_angle.sum2+=angle*angle;
			edge_angle.total+=angle;
			if(angle<edge_angle.min) edge_angle.min=angle;
			if(angle>edge_angle.max) edge_angle.max=angle;
			if(cs.signal[1]==true){
				if(angle<cs.thresholds[1]) bad_angle[*i]=angle;
			}
			if(cs.signal[2]==true){
				if(angle>cs.thresholds[2]) bad_angle[*i]=angle;
			}
			// add to vector
			edge_angle.x.push_back(angle);
		}
	}
	edge_angle.createHistogram();
}

void Controls::computeEdgeAngles(void){
	// for each element in area vector
	for(d_iterator i=edge_angle.x.begin();i!=edge_angle.x.end();i++){
		edge_angle.sum+=*i;
		edge_angle.sum2+=(*i)*(*i);
		edge_angle.total+=*i;
		if(*i<edge_angle.min) edge_angle.min=*i;
		if(*i>edge_angle.max) edge_angle.max=*i;
	}
	edge_angle.createHistogram();
}

void Object::computeVolume(void){

	vol = 0.0;

	// for each face in object
	for(f_iterator i=f.begin();i!=f.end();i++){
		double x1=(*i)->v[0]->pN[0];
		double y1=(*i)->v[0]->pN[1];
		double z1=(*i)->v[0]->pN[2];
		double x2=(*i)->v[1]->pN[0];
		double y2=(*i)->v[1]->pN[1];
		double z2=(*i)->v[1]->pN[2];
		double x3=(*i)->v[2]->pN[0];
		double y3=(*i)->v[2]->pN[1];
		double z3=(*i)->v[2]->pN[2];
		/* compute determinant of oriented triangle */
		double det=x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1);
		vol+=det;
	}
  	vol=vol/6.0;
	// update cumulative volume
//	vol+=o->vol;
}

void Object::vertexDistin(void){
	///// check vertex distinguishability /////
	// multimap: double -> Vertex*
	mmap_dv mm;
	// random vector
	double rand_vec[3]={0.236416584579274058342,
						0.927225593011826276779,
						0.389099507126957178116};
	// for each vertex
	for(v_iterator i=v.begin();i!=v.end();i++){
		// add dot product of vertex and random vector to multimap
		mm.insert(std::make_pair(dot((*i)->pN,rand_vec),*i));
	}
	// for each multimap element
	dv_iterator j;
	for(dv_iterator i=mm.begin();i!=mm.end();i++){
		j=i;j++;
		if(j!=mm.end()){
			// if sequential pair of multimap elements is not distinguishable
			if( !distinguishable((*i).second->pN[0],(*j).second->pN[0]) &&
				!distinguishable((*i).second->pN[1],(*j).second->pN[1]) &&
				!distinguishable((*i).second->pN[2],(*j).second->pN[2]) 
				){
				indistin_v.push_back((*i).second);
				indistin_v.push_back((*j).second);
			}
		}
	}
}

void Object::evalCharacteristics(Container* c,Controls &cs,Space &s){
	// vertices
	cout << "Bound object [" << name << "]..........................................";
	cout.flush();
	boundObject(cs.bb);
	cout << "complete.\n";cout.flush();
	cout << "Check if vertices are distinguishable for object [" << name << "]......";
	cout.flush();
	vertexDistin();
	cout << "complete.\n";cout.flush();
	cout << "Analyze vertex adjacent faces for object [" << name << "]..............";
	cout.flush();
	vertexAdjacentFaces();
	cout << "complete.\n";cout.flush();
	// faces
	cout << "Identify separate components of object [" << name << "]................";
	cout.flush();
	num_sep=countComponents();
//	cs.num_sep+=num_sep;
	cout << "complete.\n";cout.flush();
	cout << "Compute face area and aspect ratio for object [" << name << "].........";
	cout.flush();
	areaAspectRatio(cs);
	cout << "complete.\n";cout.flush();
	// if not batch mode
	if(cs.interf==false){
		cout << "Find intersecting faces for object [" << name << "]....................";
		cout.flush();
		findIntersectingFaces(c,s);
		cout << "complete.\n";cout.flush();
	}
	// edges
	cout << "Identify boundaries for object [" << name << "]........................";
	cout.flush();
	countBoundaries();
	cout << "complete.\n";cout.flush();
	cout << "Analyze edge lengths for object [" << name << "].......................";
	cout.flush();
	processEdgeLengths(cs);
	cout << "complete.\n";cout.flush();

	if(manifold && consistent){
		// manifold and consistently oriented face normals
		cout << "Analyze edge angles for object [" << name << "]........................";
		cout.flush();
		computeEdgeAngles(cs);
		cout << "complete.\n";cout.flush();
		if(closed){
			// closed, manifold, and consistently oriented face normals
			cout << "Compute volume of object [" << name << "]..............................";
			cout.flush();
			computeVolume();
			cout << "complete.\n";cout.flush();
			// if number of components == 1
			// FUTURE IMPROVEMENT: only require orientable, not oriented
			// FUTURE IMPROVEMENT: separate components and compute genus of each component
			if(num_sep==1 && orphan.empty()){
				cout << "Compute genus of object [" << name << "]...............................";
				cout.flush();
				computeGenus();
				cout << "complete.\n";cout.flush();
			}
		}
	}

}

//void Controls::analyzeCumulative(Container* c){
void Controls::analyzeCumulative(void){
	cout << "\n" << "/* ********************** "
	<< "SET OF OBJECTS ********************** */\n\n";
	// vertices
	cout << "Analyze vertex adjacent faces for set of all objects...................";
	cout.flush();
	vertexAdjacentFaces();
	cout << "complete.\n";cout.flush();
	// faces
	cout << "Compute face area and aspect ratio for set of all objects..............";
	cout.flush();
	areaAspectRatio();
	cout << "complete.\n";cout.flush();
	// edges
	cout << "Analyze edge lengths for set of all objects............................";
	cout.flush();
	processEdgeLengths();
	cout << "complete.\n";cout.flush();
//	if(manifold && consistent){
		// manifold and consistently oriented face normals
		cout << "Analyze edge angles for set of all objects.............................";
		cout.flush();
		computeEdgeAngles();
		cout << "complete.\n";cout.flush();
//	}

}

void Object::analyze(Container *c,Controls &cs,Space &s){
	// single object
	// evaluate mesh attributes
	evalAttributes(s);
	// eval mesh characteristics
	if(cs.attr==false){evalCharacteristics(c,cs,s);}
}

void Controls::analyzeBatch(Container *c,Space &s){
	// find intersecting faces
	if(interf){
		cout << "Find intersecting faces of all objects...";cout.flush();
		Object *oo=c->o.front();
		oo->findIntersectingFaces(c,s);
		cout << "complete.\n";cout.flush();
	}
}

void Object::printChars(Controls &cs){
	cout << "\nMESH CHARACTERISTICS\n\n";
	cout << "    # vertices: " << v.size() << endl
	<< "    # faces: " << f.size() << endl
	<< "    # edges: " << e.size() << endl
	<< "    # components: " << num_sep << endl;
	if (manifold==false){
		cout << "    # boundaries: Since object is nonmanifold,\n"
		<< "    # boundaries: the number of boundaries may be underestimated.\n";
	}
	//////////////////// borders
	if (border.empty()){
		cout << "    # boundaries: none\n";
	} else {
		cout << "    # boundaries: " << border.size() << endl;
		//	if -p option, print offending
		if(cs.print){
			int j=1;
			// for each border edge
			for(e_iterator i=border.begin();i!=border.end();i++){
				cout << "    # boundaries: boundary edge " << j++ << endl;
				if(strcmp(cs.style,"cp")==false){
					(*i)->printEdgeCP();
				} else {
					(*i)->printEdge((*i)->f1->v[0]->o->name);
					cout << endl;
				}
			}
		}
	}
	///////////////////// indistinguishable vertices
	if (indistin_v.empty()){
		cout << "    # indistinguishable vertices: none\n";
	} else {
		cout << "    # indistinguishable vertices: " << indistin_v.size() << endl;
		//	if -p option, print offending
		if(cs.print){
			int j=1;
			// for each indistinguishable vertec
			for(v_iterator i=indistin_v.begin();i!=indistin_v.end();i++){
				cout << "    #  indistinguishable vertices: vertex " << j++ << endl;
				if(strcmp(cs.style,"cp")==false){
					(*i)->printVertexCP();
				} else {
					(*i)->printVertex((*i)->o->name);
					cout << endl;
				}
			}
		}
	}
	//////////////// if volume computed
	if(closed==true && manifold==true && consistent==true){
		cout << "    object volume: [(data units)^3]" << endl;
		cout << "    object volume: " << vol << endl;
	} else {
		cout << "    object volume: not computed, since ";
		if(closed==false){cout << "not closed,";}
		if(consistent==false){cout << "not consistent,";}
		if(manifold==false){cout << "not manifold";}
		cout << endl;
	}
	/////////////// if genus computed
	if(closed==true && manifold==true && consistent==true && num_sep==1 && orphan.empty()){
		cout << "    object genus: " << genus << endl;
	} else {
		cout << "    object genus: not computed, since ";
		if(closed==false){cout << "not closed,";}
		if(consistent==false){cout << "not consistent,";}
		if(manifold==false){cout << "not manifold,";}
		if(num_sep>1){cout << "#components=" << num_sep << ",";}
		if(!orphan.empty()){cout << "orphan vertices were found";}
		cout << endl;
	}
	//////////////// bounding box
	cout << "    bounding box: [data units]\n";
	cout << "    bounding box: [xmin,ymin,zmin][xmax,ymax,zmax]\n";
	cout << "    bounding box: ["
	<< cs.bb[0] << ","
	<< cs.bb[2] << ","
	<< cs.bb[4] << "]["
	<< cs.bb[1] << ","
	<< cs.bb[3] << ","
	<< cs.bb[5] << "]" << endl;
	////////////////// edges with indistinguishable vertices
	if(indistin.empty()==true){
		cout << "    # edges with indistinguishable vertices: none\n";
	} else {
		cout << "    # edges with indistinguishable vertices: "
		<< indistin.size() << endl;
	}
	//	if -p option, print offending
	if(cs.print){
		if(indistin.empty()==false){
			// for each afflicted edge
			for(e_iterator i=indistin.begin();i!=indistin.end();i++){
				if(strcmp(cs.style,"cp")==false){
					(*i)->printEdgeCP();
				} else {
					(*i)->printEdge((*i)->f1->v[0]->o->name);
					cout << endl;
				}
			}
		}
	}
	/////////////// intersecting faces
	// if not batch mode
	if(cs.interf==false){
		// intersecting faces
		if (intf.empty()){
			cout << "    # intersecting faces: none\n\n";
		} else {
			cout << "    # intersecting faces: " << intf.size() << endl;
			//	if -p option, print offending
			if(cs.print){
				int j=1;
				// for each intersected face
				for(ff_iterator i=intf.begin();i!=intf.end();i++){
					cout << "    # intersecting faces: intersected face " << j++ << endl;
					// print intersected face
					if(strcmp(cs.style,"cp")==false){
						(*i).first->printFaceCP();
					} else {
						(*i).first->printFace((*i).first->v[0]->o);
					}
					// keep unique list of intersecting faces
					sort((*(*i).second).begin(),(*(*i).second).end());
					f_iterator new_end = unique((*(*i).second).begin(),(*(*i).second).end());
					(*(*i).second).assign((*(*i).second).begin(),new_end);
					// print intersecting faces
					for(f_iterator k=(*(*i).second).begin();k!=(*(*i).second).end();k++){
						if(strcmp(cs.style,"cp")==false){
							(*k)->printFaceCP();
						} else {
							(*k)->printFace((*k)->v[0]->o);
							cout << endl;
						}
					}
				}
				cout << endl << endl;
			}
		}
	}
	//////////////// vertex adjacent faces
	cout << "    Vertex adjacent face statistics [faces]:" << endl;
	adjacent_face.printStats();
	cout << "    Vertex adjacent face histogram [faces]:" << endl;
	adjacent_face.printAdjacentFaceHistogram();
	cout << endl;
	///////////////// face area
	cout << "    Face area statistics [(data units)^2]:" << endl;
	cout << "       total    " << area.sum << endl;
	area.printStats();
	cout << "    Face area histogram [(data units)^2]:" << endl;
	area.printHistogram();
	cout << endl;
	/////////////////// face aspect ratio
	cout << "    Face aspect ratio statistics [unitless]:" << endl;
	aspect_ratio.printStats();
	cout << "    Face aspect ratio histogram [unitless]:" << endl;
	aspect_ratio.printAspectRatioHistogram();
	cout << "      (Aspect ratio is longest edge "
	<< "divided by shortest altitude)\n";
	// if face aspect ratio threshold specified
	if(cs.signal[0]==true){
		if(bad_aspect.empty()==true){
			cout << "    # faces with bad aspect ratios: none" << endl;
		} else {
			cout << "    # faces with bad aspect ratios: "
			<< bad_aspect.size() << endl;
		}
	}
	//	if -p option, print offending
	if(cs.print){
		if(bad_aspect.empty()==false){
			// print faces with aspect ratios that violate threshold
			Face_Pair fp(this);
			for(fd_iterator i=bad_aspect.begin();i!=bad_aspect.end();i++){
				cout << "    face aspect ratio: " << (*i).second << endl;
				if(strcmp(cs.style,"cp")==false){
					(*i).first->printFaceCP();
				} else {
					(*i).first->printFace((*i).first->v[0]->o);
					cout << endl;
				}
				fp.processBadFace((*i).first);
			}
		}
	}
	cout << "\n";
	////////////// edge length
	cout << "    Edge length statistics [data units]:" << endl;
	edge_length.printStats();
	cout << "    Edge length histogram [data units]:" << endl;
	edge_length.printHistogram();
	cout << endl;
	// if edge length threshold specified
	if(cs.signal[3]==true || cs.signal[4]==true){
		if(bad_length.empty()==true){
			cout << "    # edges with bad lengths: none" << endl;
		} else {
			cout << "    # edges with bad lengths: "
			<< bad_length.size() << endl;
		}
	}
	//	if -p option, print offending
	if(cs.print){
		if(bad_length.empty()==false){
			// print edges with lengthss that violate threshold
			for(ed_iterator i=bad_length.begin();i!=bad_length.end();i++){
				cout << "    edge length: " << (*i).second << endl;
				if(strcmp(cs.style,"cp")==false){
					(*i).first->printEdgeCP();
				} else {
					(*i).first->printEdge((*i).first->vv1->o->name);
					cout << endl;
				}
			}
		}
	}
	//////////////////// if edge angles computed
	if(manifold==true && consistent==true){
		cout << "    Edge angle statistics [degress]:" << endl;
		edge_angle.printStats();
		cout << "    Edge angle histogram [degress]:" << endl;
		edge_angle.printHistogram();
		cout << endl;
		// if edge angle threshold specified
		if(cs.signal[1]==true || cs.signal[2]==true){
			if(bad_angle.empty()==true){
				cout << "    # edges with bad angles: none" << endl;
			} else {
				cout << "    # edges with bad angles: "
				<< bad_angle.size() << endl;
			}
		}
		//	if -p option, print offending
		if(cs.print){
			if(bad_angle.empty()==false){
				// print edges with angles that violate threshold
				for(ed_iterator i=bad_angle.begin();i!=bad_angle.end();i++){
					cout << "    edge angle: " << (*i).second << endl;
					if(strcmp(cs.style,"cp")==false){
						(*i).first->printEdgeCP();
					} else {
						(*i).first->printEdge((*i).first->vv1->o->name);
						cout << endl;
					}
				}
			}
		}
	} else {
		cout << "    edge angles: not computed, since ";
		if(consistent==false){cout << "not consistent,";}
		if(manifold==false){cout << "not manifold";}
		cout << endl;
	}

}

void Object::printIntegrity(Controls &cs){
	cout << "\nMESH FILE INTEGRITY\n\n";
	// orphan vertices
	if(orphan.empty()){
		cout << "    # orphan vertices: none\n";
	} else {
		cout << "    # orphan vertices: " << orphan.size() << endl;
		//	if -p option, print offending
		if(cs.print){
			int j=1;
			// for each orphan vertex
			for(v_iterator i=orphan.begin();i!=orphan.end();i++){
				cout << "    # orphan vertices: orphan vertex " << j++ << endl;
				if(strcmp(cs.style,"cp")==false){
					(*i)->printVertexCP();
				} else {
					(*i)->printVertex((*i)->o->name);
					cout << endl;
				}
			}
		}
	}
	// missing vertices
	if(missing_v.empty()){
		cout << "    # missing vertices: none\n";
	} else {
		cout << "    # missing vertices: " << missing_v.size() << endl;
		//	if -p option, print offending
		if(cs.print){
			int j=1;
			// for each missing vertex
			for(i_iterator i=missing_v.begin();i!=missing_v.end();i++){
				cout << "    # missing vertices: #" << j++
				 << "-> missing vertex index " << *i << endl;
			}
			j=1;
			// for each missing face
			for(f_iterator i=missing_f.begin();i!=missing_f.end();i++){
				cout << "    # missing vertices: affected face " << j++ << endl;
				if((*i)->v[0]!=NULL){
					if(strcmp(cs.style,"cp")==false){
						(*i)->printFaceCP();
					} else {
						(*i)->printFace((*i)->v[0]->o);
						cout << endl;
					}
				} else if((*i)->v[1]!=NULL){
					if(strcmp(cs.style,"cp")==false){
						(*i)->printFaceCP();
					} else {
						(*i)->printFace((*i)->v[1]->o);
						cout << endl;
					}
				} else if((*i)->v[2]!=NULL){
					if(strcmp(cs.style,"cp")==false){
						(*i)->printFaceCP();
					} else {
						(*i)->printFace((*i)->v[2]->o);
						cout << endl;
					}
				}
			}
		}
	}
	// degenerate faces
	if(degen.empty()){
		cout << "    # degenerate faces: none\n";
	} else {
		cout << "    # degenerate faces: " << degen.size() << endl;
		//	if -p option, print offending
		if(cs.print){
			int j=1;
			// for each degenerate face
			for(f_iterator i=degen.begin();i!=degen.end();i++){
				cout << "    # degenerate faces: affected face " << j++ << endl;
				if(strcmp(cs.style,"cp")==false){
					(*i)->printFaceCP();
				} else {
					(*i)->printFace((*i)->v[0]->o);
					cout << endl;
				}
			}
		}
	}

	// duplicate indices
	if(dupl_v_index.empty()){
		cout << "    # duplicate vertex indices: none\n";
	} else {
		cout << "    # duplicat vertex indices: " << dupl_v_index.size() << endl;
		//	if -p option, print offending
		if(cs.print){
			int j=1;
			// for each vertex with a duplicate vertex index
			for(v_iterator i=dupl_v_index.begin();i!=dupl_v_index.end();i++){
				cout << "    # duplicate vertex indices: affected vertex " << j++ << endl;
				if(strcmp(cs.style,"cp")==false){
					(*i)->printVertexCP();
				} else {
					(*i)->printVertex((*i)->o->name);
					cout << endl;
				}
			}
		}
	}
	if(dupl_f_index.empty()){
		cout << "    # duplicate face indices: none\n";
	} else {
		cout << "    # duplicat face indices: " << dupl_f_index.size() << endl;
		//	if -p option, print offending
		if(cs.print){
			int j=1;
			// for each face with a duplicate face index
			for(f_iterator i=dupl_f_index.begin();i!=dupl_f_index.end();i++){
				cout << "    # duplicate face indices: affected face " << j++ << endl;
				if(strcmp(cs.style,"cp")==false){
					(*i)->printFaceCP();
				} else {
					(*i)->printFace((*i)->v[0]->o);
					cout << endl;
				}
			}
		}
	}


	// contiguous numbering
	if(contig_v){
		cout << "    contiguous vertex indexing from 1: yes\n";
	} else {
		cout << "    contiguous vertex indexing from 1: no\n";
		cout << "    contiguous vertex indexing from 1: bad index and +/- 2\n";
		cout << "    contiguous vertex indexing from 1: ";
		if(vec_cv[0]!=-1){ cout << vec_cv[0] << " ";}
		else 			{cout << "NA ";}
		if(vec_cv[1]!=-1){ cout << vec_cv[1] << " ";}
		else 			{cout << "NA ";}
		if(vec_cv[2]!=-1){ cout << vec_cv[2] << " ";}
		else 			{cout << "NA ";}
		if(vec_cv[3]!=-1){ cout << vec_cv[3] << " ";}
		else 			{cout << "NA ";}
		if(vec_cv[4]!=-1){ cout << vec_cv[4] << endl;}
		else 			{cout << "NA\n";}
	}
	if(contig_f){
		cout << "    contiguous face indexing from 1: yes\n";
	} else {
		cout << "    contiguous face indexing from 1: no\n";
		cout << "    contiguous face indexing from 1: bad index and +/- 2\n";
		cout << "    contiguous face indexing from 1: ";
		if(vec_cf[0]!=-1){ cout << vec_cf[0] << " ";}
		else 			{cout << "NA ";}
		if(vec_cf[1]!=-1){ cout << vec_cf[1] << " ";}
		else 			{cout << "NA ";}
		if(vec_cf[2]!=-1){ cout << vec_cf[2] << " ";}
		else 			{cout << "NA ";}
		if(vec_cf[3]!=-1){ cout << vec_cf[3] << " ";}
		else 			{cout << "NA ";}
		if(vec_cf[4]!=-1){ cout << vec_cf[4] << endl;}
		else 			{cout << "NA\n";}
	}
}

void Object::printAttr(Controls &cs){
	cout << "\nMESH ATTRIBUTES\n\n";
	// closed
	if(closed==true){
		cout << "    mesh is closed: yes\n";
	} else {
		cout << "    mesh is closed: no\n";
		//	if -p option, print offending
		if(cs.print){
			int j=1;
			cout << "    mesh is closed: # border edges - " << border.size() << endl;
			// for each border edge
			for(e_iterator i=border.begin();i!=border.end();i++){
				cout << "    mesh is closed: border edge # " << j++ << endl;
				if(strcmp(cs.style,"cp")==false){
					(*i)->printEdgeCP();
				} else {
					(*i)->printEdge((*i)->f1->v[0]->o->name);
					cout << endl;
				}
			}
		}
	}
	// manifold
	if(manifold==true){
		cout << "    mesh is manifold: yes\n";
	} else if(manifold==false && closed==true){
		cout << "    mesh is manifold: no\n";
		//	if -p option, print offending
		if(cs.print){
			if(nonman_v.empty()){
				cout << "    mesh is manifold: # nonmanifold vertices - none\n";
			} else {
				int j=1;
				cout << "    mesh is manifold: # nonmanifold vertices - " << nonman_v.size() << endl;
				// for each nonmanifold vertex
				for(v_iterator i=nonman_v.begin();i!=nonman_v.end();i++){
					cout << "    mesh is manifold: nonmanifold vertex # " << j++ << endl;
					if(strcmp(cs.style,"cp")==false){
						(*i)->printVertexCP();
					} else {
						(*i)->printVertex((*i)->o->name);
						cout << endl;
					}
				}
			}
			if(nonman_e.empty()){
				cout << "    mesh is manifold: # nonmanifold edges - none\n";
			} else {
				int j=1;
				cout << "    mesh is manifold: # nonmanifold edges - " << nonman_e.size() << endl;
				// for each nonmanifold edge
				for(e_iterator i=nonman_e.begin();i!=nonman_e.end();i++){
					cout << "    mesh is manifold: nonmanifold edge # " << j++ << endl;
					if(strcmp(cs.style,"cp")==false){
						(*i)->printEdgeCP();
					} else {
						(*i)->printEdge((*i)->f1->v[0]->o->name);
						cout << endl;
					}
				}
			}
		}
	} else {
		cout << "    mesh is manifold: undefined since mesh is open\n";
	}
	// consistent
	if(manifold==false){
			cout << "    mesh has consistently oriented face normals: undefined since not manifold\n";
	} else {
		if(consistent==true){
			cout << "    mesh has consistently oriented face normals: yes\n";
		} else {
			cout << "    mesh has consistently oriented face normals: no\n";
			//	if -p option, print offending
			if(cs.print){
				int j=1;
				cout << "    mesh has consistently oriented face normals: # flipped edges - " << flipped.size() << endl;
				// for each flipped edge
				// FUTURE IMPROVEMENT: if edge is nonmanifold then exclude from flipped list
				for(e_iterator i=flipped.begin();i!=flipped.end();i++){
					cout << "    mesh has consistently oriented face normals: flipped edge # " << j++ << endl;
					if(strcmp(cs.style,"cp")==false){
						(*i)->printEdgeCP();
					} else {
						(*i)->printEdge((*i)->f1->v[0]->o->name);
						cout << endl;
					}
				}
			}
		}
	}
	// outward
	if(manifold==false || consistent==false || closed==false){
		cout << "    mesh has outward oriented face normals: uncomputable since ";
		if(closed==false){cout << "not closed,";}
		if(consistent==false){cout << "not consistent,";}
		if(manifold==false){cout << "not manifold";}
		cout << endl;
	} else {
		if(outward==true){
			cout << "    mesh has outward oriented face normals: yes\n";
		} else {
			cout << "    mesh has outward oriented face normals: no\n";
		}
	}
}

void Object::print(Controls &cs){
/*	cout << "\n\n" << "************************ "
	<< "OBJECT *************************\n";
	//	print object name 
	cout << "name: " << name << endl;
*/
	//	print Integrity

	for(std::vector<Contour*>::iterator i=contours.begin();i!=contours.end();i++){
		(*i)->print(cs,name);
	}
}

void Controls::printCumulative(Container &c){
	// ISSUE WARNING IF ANY OBJECT FAILED INTEGRITY 
	// TEST OR LACK ATTRIBUTES
	// NOTE CUMULATIVE VOLUME ASSUMES ALL 
	// MESHES HAVE SAME ORIENTATION
	//	print Integrity
	printIntegrity(c);
	if(good_integrity==false){
		cout << "\n\nWarning: Attributes and "
		<< "characteristics were not evaluated,\n"
		<< " since mesh file failed the integrity check.\n\n";
	} else {
		//	print attributes
		printAttr(c);
		//	print characteristics
		if(attr==false) {
			printChars(c);
		}
	}
	cout << "/* ********************** "
	<< "END ********************** */\n\n";
}

void Controls::printIntegrity(Container &c){
	cout << "\nMESH SET INTEGRITY:\n\n";
	int a=c.countOrphan();
	if(a==0){
		cout << "    # orphan vertices: none\n";
	} else {
		cout << "    # orphan vertices: " << a << endl;
	}
	a=c.countMissing();
	if(a==0){
		cout << "    # missing vertices: none\n";
	} else {
		cout << "    # missin vertices: " << a << endl;
	}
	a=c.countDegen();
	if(a==0){
		cout << "    # degenerate faces: none\n";
	} else {
		cout << "    # degenerate faces: " << a << endl;
	}
	a=c.countDuplV();
	if(a==0){
		cout << "    # duplicate vertex indices: none\n";
	} else {
		cout << "    # duplicate vertex indices: " << a << endl;
	}
	a=c.countDuplF();
	if(a==0){
		cout << "    # duplicate face indices: none\n";
	} else {
		cout << "    # duplicate face indices: " << a << endl;
	}
	cout << endl;
}

void Controls::printAttr(Container &c){
	cout << "MESH SET ATTRIBUTES:\n\n";
	if(good_integrity==false){
		cout << "    Warning: These attribute summaries may be inaccurate,\n"
		<< "    since some mesh files failed the integrity check.\n";
	}
	// closed
	std::pair<int,int> pp=c.countClosed();
	if(pp.first==0){
		cout << "    # closed mesh objects: none\n";
	} else {
		cout << "    # closed mesh objects: " << pp.first << endl;
	}
	if(pp.second==0){
		cout << "    # open mesh objects: none\n";
		cout << "    # border edges: none\n";
	} else {
		cout << "    # open mesh objects: " << pp.second << endl;
		cout << "    # border edges: " << c.countBorder() << endl;
	}
	// manifold
	int val[3]={0,0,0};
	c.countManifold(val);
	if(val[0]==0){
		cout << "    # manifold mesh objects: none\n";
	} else {
		cout << "    # manifold mesh objects: " << val[0] << endl;
	}
	if(val[2]==0){
		cout << "    # mesh objects with undefined manifoldness: none\n";
	} else {
		cout << "    # mesh objects with undefined manifoldness: " << val[2] << endl;
	}
	if(val[1]==0){
		cout << "    # nonmanifold mesh objects: none\n";
		cout << "    # nonmanifold vertices: none\n";
		cout << "    # nonmanifold edges: none\n";
	} else {
		cout << "    # nonmanifold mesh objects: " << val[1] << endl;
		cout << "    # nonmanifold vertices: " << c.countNonmanV() << endl;
		cout << "    # nonmanifold edges: " << c.countNonmanE() << endl;
	}
	// consistent
	c.countConsistent(val);
	if(val[0]==0){
		cout << "    # mesh objects with consistently "
		<< "oriented face normals: none\n";
	} else {
		cout << "    # mesh objects with consistently "
		<< "oriented face normals: " << val[0] << endl;
	}
	if(val[1]==0){
		cout << "    # mesh objects with inconsistently "
		<< "oriented face normals: none\n";
	} else {
		cout << "    # mesh objects with inconsistently "
		<< "oriented face normals: " << val[1] << endl;
		cout << "       # flipped edges: " << c.countFlipped() << endl;
	}
	if(val[2]==0){
		cout << "    # mesh objects whose face normal "
		<< "orientation is undefined: none\n";
	} else {
		cout << "    # mesh objects whose face normal "
		<< "orientation is undefined: " << val[2] << endl;
	}
	// outward
	c.countOutward(val);
	if(val[0]==0){
		cout << "    # mesh objects with outward "
		<< "oriented face normals: none\n";
	} else {
		cout << "    # mesh objects with outward "
		<< "oriented face normals: " << val[0] << endl;
	}
	if(val[1]==0){
		cout << "    # mesh objects with inward "
		<< "oriented face normals: none\n";
	} else {
		cout << "    # mesh objects with inward "
		<< "oriented face normals: " << val[1] << endl;
	}
	if(val[2]==0){
		cout << "    # mesh objects whose face normal "
		<< "orientation is undefined: none\n";
	} else {
		cout << "    # mesh objects whose face normal "
		<< "orientation is undefined: " << val[2] << endl;
	}
	cout << endl;
}

void Controls::printChars(Container &c){
	cout << "MESH SET CHARACTERISTICS:\n\n";
	//	print characteristics
	if(good_integrity==false){
		cout << "    Warning: These characteristics "
		<< "summaries may be inaccurate,\n"
		<< "    since some mesh files failed the integrity check.\n";
	}
	cout << "    # objects: " << c.countObject() << endl;
	cout << "    # vertices: " << c.countVertex() << endl;
	cout << "    # faces: " << c.countFace() << endl;
	cout << "    # edges: " << c.countEdge() << endl;
	cout << "    # components: " << c.countComponents() << endl;
	int val[3]={0,0,0};
	c.countManifold(val);
	if (val[1]>0){
		cout << "    # boundaries: Since object is nonmanifold,\n"
		<< "    # boundaries: the number of boundaries may be underestimated.\n";
	}
	int a=c.countBoundaries();
	if(a==0){
		cout << "    # boundaries: none\n";
	} else {
		cout << "    # boundaries: " << a << endl;
	}
	a=c.countIndistin();
	if(a==0){
		cout << "    # indistinguishable vertices: none\n";
	} else {
		cout << "    # indistinguishable vertices: " << a << endl;
	}
	// volume
	c.countConsistent(val);
	if(val[0]>0){
		cout << "    object volume: [(data units)^3]" << endl;
		cout << "    object volume: " << c.countVol() << endl;
	} else {
		cout << "    object volume: not computed, since no mesh"
		<< "    object with consistently oriented face normals was found.\n";
	}
	// TODO
	// ADD GENUS STATS (HISTOGRAM TOO)
	// TODO
	// bounding box
	cout << "    bounding box: [data units]\n";
	cout << "    bounding box: [xmin,ymin,zmin][xmax,ymax,zmax]\n";
	cout << "    bounding box: ["
	<< bb[0] << ","
	<< bb[2] << ","
	<< bb[4] << "]["
	<< bb[1] << ","
	<< bb[3] << ","
	<< bb[5] << "]" << endl << endl;
	// vertex adjacent faces
	cout << "    Vertex adjacent face statistics [faces]:" << endl;
	adjacent_face.printStats();
	cout << "    Vertex adjacent face histogram [faces]:" << endl;
	adjacent_face.printAdjacentFaceHistogram();
	cout << endl;
	// face area
	cout << "    Face area statistics [(data units)^2]:" << endl;
	cout << "       total    " << area.sum << endl;
	area.printStats();
	cout << "    Face area histogram [(data units)^2]:" << endl;
	area.printHistogram();
	cout << endl;
	// face aspect ratio
	cout << "    Face aspect ratio statistics [unitless]:" << endl;
	aspect_ratio.printStats();
	cout << "    Face aspect ratio histogram [unitless]:" << endl;
	aspect_ratio.printAspectRatioHistogram();
	printf("      (Aspect ratio is longest edge divided by shortest altitude)\n");
	cout << endl;
	// edge length
	cout << "    Edge length statistics [data units]:" << endl;
	edge_length.printStats();
	// DEBUG
/*	bool ffound=false;
	for(std::vector<double>::iterator ww=edge_length.x.begin();ww!=edge_length.x.end();ww++){
		double diff = *ww-edge_length.max;
		if(fabs(diff)<1E-6){
			cout << "\nControls::printChars: max value=" << edge_length.max
			<< " found in edge_length.x\n";
			ffound=true;
		}
	}
	if(ffound==false){
		cout << "\nControls::printChars: max value=" << edge_length.max
		<< " NOT found in edge_length.x\n";
	}
*/	// DEBUG
	cout << "    Edge length histogram [data units]:" << endl;
	edge_length.printHistogram();
	cout << endl;
	// if edge angles computed, i.e. if at least one mesh was conistent
	if(val[0]>0){
		cout << "    Edge angle statistics [degrees]:" << endl;
		edge_angle.printStats();
/*		// DEBUG
		ffound=false;
		for(std::vector<double>::iterator ww=edge_angle.x.begin();ww!=edge_angle.x.end();ww++){
			double diff = *ww-edge_angle.max;
			if(fabs(diff)<1E-6){
				cout << "\nControls::printChars: max value=" << edge_angle.max
				<< " found in edge_angle.x\n";
				ffound=true;
			}
		}
		if(ffound==false){
			cout << "\nControls::printChars: max value=" << edge_angle.max
			<< " NOT found in edge_angle.x\n";
		}
*/		// DEBUG
		cout << "    Edge angle histogram [degrees]:" << endl;
		edge_angle.printHistogram();
		cout << endl;
	}
}

void Container::printBatch(Controls &cs){
	int num_intf = countIntFace();
	cout << "\n\n" << "/* ********************** BATCH RESULTS ********************** */\n";
	// intersecting faces
	if (num_intf==0){
		cout << "# intersecting faces: none\n";
	} else {
		cout << "# intersecting faces: " << num_intf << endl;
		//	if -p option, print offending
		if(cs.print){
			int j=1;
			// for each object
			for(o_iterator m=o.begin();m!=o.end();m++){
				// for each intersected face
				for(ff_iterator i=(*m)->intf.begin();i!=(*m)->intf.end();i++){
					cout << "# intersecting faces: intersected face " << j++ << endl;
					// print intersected face
					(*i).first->printFace((*i).first->v[0]->o);
					// print intersecting faces
					for(f_iterator k=(*(*i).second).begin();k!=(*(*i).second).end();k++){
						(*k)->printFace((*k)->v[0]->o);
						cout << endl;
					}
				}
			}
		}
	}
}

void Face::getzEM(Vertex *&aa,Vertex *&bb,double zEM){
	aa=bb=NULL;
	if(distinguishable(v[0]->pN[2],zEM)==false){aa=v[0];}
	if(distinguishable(v[1]->pN[2],zEM)==false){
		if(aa==NULL){aa=v[1];}else{bb=v[1];}
	}
	if(distinguishable(v[2]->pN[2],zEM)==false){
		if(bb!=NULL){
			cout << "\nFace::getzEM: "
			<< "Error. bb!=NULL\n";
			exit(0);
		}
		bb=v[2];
	}
	if(aa==bb){
		cout << "\nFace::getzEM: "
		<< "Error. aa==bb\n";
		exit(0);
	}
}

//void Object::addContour(vec_f &fg,double zEM){
void Object::addContour(vec_f &fg,double zEM){
	// instantiate contour
	Contour *cc;
	cc=new Contour();
	contours.push_back(cc);
	f_iterator i=fg.begin();
	Vertex *v1=NULL,*v2=NULL;
	(*i)->getzEM(v1,v2,zEM);
	Vertex *seed=v1;
	Vertex *cv=v2;
	if(cv==NULL || seed==NULL){
		cout << "Object::addContour: "
		<< " eiher cv or seed or both are NULL:\n";
		exit(0);	
	}
	cc->p.push_back(seed);
	cc->p.push_back(cv);
	fg.erase(i);
	if(fg.empty()==false){
		bool growing=true;
		while(growing==true){
			growing=false;
			// for each face in vector
			for(f_iterator j=fg.begin();j!=fg.end();j++){
				// if a face vertex is same as current vertex
				if((*j)->v[0]==cv||(*j)->v[1]==cv||(*j)->v[2]==cv){
					(*j)->getzEM(v1,v2,zEM);
					if(v1!=cv){cv=v1;}else{cv=v2;}
					if(cv!=seed){
						cc->p.push_back(cv);
						growing=true;
					}
					fg.erase(j);
					break;
				}
			}
		}
	}
	// DEBUG
//	cc->print(cs,name);
	// DEBUG
}

void Object::findContours(double zEM){
	//////// find vertices with z==zEM
//	vec_v vg;
	// for each vertex in object
//	for(v_iterator i=v.begin();i!=v.end();i++){
//		if(distinguishable((*i)->pN[2],zEM)==false){
//			vg.push_back(*i);
//		}
//	}
//	sort(vg.begin(),vg.end());
	//////// find faces with two vertices at z==zEM
	vec_f fg;
	// for each face in object
	for(f_iterator i=f.begin();i!=f.end();i++){
		int k=0;
		for(int j=0;j<3;j++){
//			if(binary_search(vg.begin(),vg.end(),(*i)->v[j])==true){k++;}
			if(distinguishable((*i)->v[j]->pN[2],zEM)==false){k++;}
		}
		if(k==2){fg.push_back(*i);}
	}
	// DEBUG
//	cout << "\nObject::findContours: "
//	<< " fg.size = " << fg.size() << endl;
	// DEBUG
	//////// get contours at z==zEM
	// while neither vector is empty
//	while(vg.empty()==false && fg.empty()==false){
	while(fg.empty()==false){
		// add a contour to object
//		addContour(fg,zEM);
		addContour(fg,zEM);
	}
}

