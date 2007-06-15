double dot(double a[3],double b[3]){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

bool getPointEdgeDistance(double p[3],double P[3][3]){
    int pairs[3][2];
    pairs[0][0] = 0;
    pairs[0][1] = 1;
    pairs[1][0] = 1;
    pairs[1][1] = 2;
    pairs[2][0] = 2;
    pairs[2][1] = 0;
    // cp=current_pairs
    double Ax,Ay,Az,Bx,By,Bz,uDen,AdotA_minus_AdotB,AdotB,u;
    // for each face edge
    for (int i=0;i<3;i++) {
        Ax = P[pairs[i][0]][0];
        Ay = P[pairs[i][0]][1];
        Az = P[pairs[i][0]][2];
        Bx = P[pairs[i][1]][0];
        By = P[pairs[i][1]][1];
        Bz = P[pairs[i][1]][2];
        AdotA_minus_AdotB = Ax*Ax+Ay*Ay+Az*Az-(Ax*Bx+Ay*By+Az*Bz);
        AdotB = Ax*Bx+Ay*By+Az*Bz;
        //uDen = AdotA-AdotB+BdotB-AdotB;
        uDen = AdotA_minus_AdotB+Bx*Bx+By*By+Bz*Bz-AdotB;
        cout.precision(15);
        if(uDen) {
            // u = AdotA-AdotB-AdotC+BdotC)/uDen
            u = (AdotA_minus_AdotB-(Ax*p[0]+Ay*p[1]+Az*p[2])
                +(Bx*p[0]+By*p[1]+Bz*p[2]))/uDen;
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

void checkLineFaceIntersection(Face *f,double lp[2][3],double pn[3],bool &line_flag,
                                bool &poly_flag, bool &poly_edge_flag,
                                bool d_polygon_edge_intersection) {
    //lp  = line_points
    // initialize flags
    line_flag=poly_flag=poly_edge_flag=false;
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
            //pvc = polygon_vertex_coordinates
            Vertex *v1=f->v[1],*v2=f->v[2];
            // get face vertex coordinates
            double pvc[3][3] = {{v0->pN[0],v0->pN[1],v0->pN[2]},
                                {v1->pN[0],v1->pN[1],v1->pN[2]},
                                {v2->pN[0],v2->pN[1],v2->pN[2]}};
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
            if (big == 0) {
                tv[0] = 0;
                tv[1] = 1;
                pI[0] = (1-u)*lp[0][0] + u*lp[1][0];
                pI[1] = (1-u)*lp[0][1] + u*lp[1][1];
                I[0]=pI[0];
                I[1]=pI[1];
                I[2] = (1-u)*lp[0][2] + u*lp[1][2];
            } else if (big == 1) {
                tv[0] = 1;
                tv[1] = 2;
                pI[0] = (1-u)*lp[0][1] + u*lp[1][1];
                pI[1] = (1-u)*lp[0][2] + u*lp[1][2];
                I[1]=pI[0];
                I[2]=pI[1];
                I[0] = (1-u)*lp[0][0] + u*lp[1][0];
            } else {
                tv[0] = 0;
                tv[1] = 2;
                pI[0] = (1-u)*lp[0][0] + u*lp[1][0];
                pI[1] = (1-u)*lp[0][2] + u*lp[1][2];
                I[0]=pI[0];
                I[2]=pI[1];
                I[1] = (1-u)*lp[0][1] + u*lp[1][1];
            }
            ////////// is point of intersection on other polygon? //////////
            // does point of intersection lie on polygon edge?
            if (d_polygon_edge_intersection) {
                if (getPointEdgeDistance(I,pvc)){poly_edge_flag = true;}
            }
            // if point of intersection is not on polygon edge
            if (!poly_edge_flag) {
                // compute three determinants
                double det[3];
                det[0] = (pvc[0][tv[0]]-pI[0])*(pvc[1][tv[1]]-pI[1])
                        -(pvc[1][tv[0]]-pI[0])*(pvc[0][tv[1]]-pI[1]);
                det[1] = (pvc[1][tv[0]]-pI[0])*(pvc[2][tv[1]]-pI[1])
                        -(pvc[2][tv[0]]-pI[0])*(pvc[1][tv[1]]-pI[1]);
                // proceed if determinants are DOUBLE_EPSILON away from zero
                if(fabs(det[0])>DOUBLE_EPSILON && fabs(det[1])>DOUBLE_EPSILON) {
                    if (((det[0]<0)&&(det[1]<0))||((det[0]>0)&&(det[1]>0))){
                        det[2]=(pvc[2][tv[0]]-pI[0])*(pvc[0][tv[1]]-pI[1])
                                -(pvc[0][tv[0]]-pI[0])*(pvc[2][tv[1]]-pI[1]);
                        if(fabs(det[2])>DOUBLE_EPSILON) {
                            if ( ( (det[0] < 0) && (det[1] < 0) && (det[2] < 0) ) ||
                                ( (det[0] > 0) && (det[1] > 0) && (det[2] > 0) ) ){
                                // line intersects polygon plane inside polygon
                                poly_flag = true;
                            }
                        }
                    }
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
	//
	double lp[2][3];
	bool line_flag=false, poly_flag=false, poly_edge_flag;
	// get face normal
	double n[3];
	of->getNormal(n);
	// for each current polygon edge
	for (int i=0;i<3;i++) {
		lp[0][0] = cf->v[pairs[i][0]]->pN[0];
		lp[0][1] = cf->v[pairs[i][0]]->pN[1];
		lp[0][2] = cf->v[pairs[i][0]]->pN[2];
		lp[1][0] = cf->v[pairs[i][1]]->pN[0];
		lp[1][1] = cf->v[pairs[i][1]]->pN[1];
		lp[1][2] = cf->v[pairs[i][1]]->pN[2];

		checkLineFaceIntersection(of,lp,n,line_flag,poly_flag,poly_edge_flag,true);
		if (line_flag && (poly_flag || poly_edge_flag)) {return(1);}
	}
	// get face normal
	cf->getNormal(n);
	// for each current polygon edge
	for (int i=0;i<3;i++) {
		lp[0][0] = of->v[pairs[i][0]]->pN[0];
		lp[0][1] = of->v[pairs[i][0]]->pN[1];
		lp[0][2] = of->v[pairs[i][0]]->pN[2];
		lp[1][0] = of->v[pairs[i][1]]->pN[0];
		lp[1][1] = of->v[pairs[i][1]]->pN[1];
		lp[1][2] = of->v[pairs[i][1]]->pN[2];

		checkLineFaceIntersection(cf,lp,n,line_flag,poly_flag,poly_edge_flag,true);
		if (line_flag && (poly_flag || poly_edge_flag)) {return(1);}
	}
	return(0);
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
	std::vector<Face*>::iterator i;
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
    if ( (e->v1->index==va && e->v2->index==vb) ||
        (e->v1->index==vb && e->v2->index==va) ){return true;}
    else {return false;}
}

Edge* Object::findEdge(Vertex* va,Vertex* vb,hashtable_t &hm,int num_digits){
	Edge *ee=NULL;
	std::string s = keyPair(va->index,vb->index,num_digits);
	// if element exists given key, then get Edge pointer
	if (hm.count(s)>0){ ee=hm[s]; }
	return ee;
}

void Edge::update(Face *f,Vertex *vc){
    //add face to edge
	if(f1==NULL) {f1=f;}
	else if (f2==NULL) {f2=f;}
	else {
		fvec.push_back(f);
	}
	if(o1==NULL) {o1=vc;}
	else if (o2==NULL) {o2=vc;}
	else { 
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

void Object::createEdge(Face *ff,Vertex* va,Vertex* vb,Vertex* vc,hashtable_t &hm,int num_digits){
	// new edge
	Edge *en = new Edge(ff,va,vb,vc,this);
	// store edge pointer in hash table
	hm[keyPair(va->index,vb->index,num_digits)]=en;
	// add edge pointer to face
	ff->addEdge(en);
	// add edge pointer to object
	e.push_back(en);
}

void Object::checkEdge(Face *ff,Vertex *va,Vertex *vb,Vertex *vc,hashtable_t &hm,int num_digits) {
	Edge *ee=NULL;
    ee=findEdge(va,vb,hm,num_digits);
    if(ee!=NULL){ee->update(ff,vc);}
    else {createEdge(ff,va,vb,vc,hm,num_digits);}
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

void Object::createEdges(void) {
	// determine number of digits in largest vertex index
	int num_digits = setNumDigits();
	// create map for finding edges
	hashtable_t hm;
	std::vector<Face*>::iterator i;
	// for each face
	for (i=f.begin();i!=f.end();i++) {
		checkEdge(*i,(*i)->v[0],(*i)->v[1],(*i)->v[2],hm,num_digits);
        checkEdge(*i,(*i)->v[1],(*i)->v[2],(*i)->v[0],hm,num_digits);
        checkEdge(*i,(*i)->v[2],(*i)->v[0],(*i)->v[1],hm,num_digits);
	}
}
/*
int Object::getMaxVertex(void){
	std::vector<Vertex*>::iterator i;
	// initialize max vertex index storage
	int a=v[0]->index;
	// for each vertex
	for (i=v.begin();i!=v.end();i++) {
		if ((*i)->index>a){a=(*i)->index;}
	}
	return a;
}*/

void Vertex::getAdjacentVertices(vector<Vertex*> &a){
	a.clear();
	std::vector<Edge*>::iterator i;
	// for each adjacent edge
	for (i=e.begin();i!=e.end();i++) {
		// find vertex different from self and add different vertex to vector
		if (((*i)->v1)->index!=index){a.push_back((*i)->v1);}
		else if (((*i)->v2)->index!=index) {a.push_back((*i)->v2);}
		else { printf("Error. both vertices of edge are equal to current vertex.\n"); exit(1); }
	}
}

double Object::getMeanEdgeLength(void){
	std::vector<Edge*>::iterator i;
	double L = 0;
	// for each edge
	for (i=e.begin();i!=e.end();i++) {
		L+=sqrt((*i)->getSqLength());
	}
	// compute mean edge length
	return (L/(double)e.size());
}

void Vertex:: getAdjacentFaces(hashset_f &fset){
	std::vector<Edge*>::iterator i;
	// for each adjacent edge
	for (i=e.begin();i!=e.end();i++) {
		// add edge faces to set
		fset.insert((*i)->f1);
		fset.insert((*i)->f2);
	}
}

bool Object::processEdge(Edge *ee,hashtable_v_double &hood,std::vector<Edge*> &bucket,Vertex *vp){
    bool empty = true;
    bool f1 =false,f2=false;
    Vertex *vp1=ee->v1,*vp2=ee->v2;
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
/*      bool u1=false,u2=false;
        double L = sqrt(ee->getSqLength());
        // if v2 + length < v1, then update v1
        if((hood[vp2]+L)<hood[vp1]){u1=true;}
        // if v1 + length < v2, then update v2
        if((hood[vp1]+L)<hood[vp2]){u2=true;}
        if(u1&&u2){cout << "Object::newFindNeighborhoods: Error. Something is weird.\n";exit(0);}
        if(u1){hood[vp1]=hood[vp2]+L;}
        if(u2){hood[vp2]=hood[vp1]+L;}
*/
    } else {
        if(f1){
            // if v1 found and v2 not found in hashtable
            // then add v2 to hashtable
            hood[vp2]=sqrt( (vp2->pN[0]-vp->pN[0])*(vp2->pN[0]-vp->pN[0])+
                            (vp2->pN[1]-vp->pN[1])*(vp2->pN[1]-vp->pN[1])+
                            (vp2->pN[2]-vp->pN[2])*(vp2->pN[2]-vp->pN[2]));
        }
        if(f2){
            // if v2 found and v1 not found in hashtable
            // then add v1 to hashtable
//          hood[vp1]=hood[vp2]+sqrt(ee->getSqLength());
            hood[vp1]=sqrt( (vp1->pN[0]-vp->pN[0])*(vp1->pN[0]-vp->pN[0])+
                            (vp1->pN[1]-vp->pN[1])*(vp1->pN[1]-vp->pN[1])+
                            (vp1->pN[2]-vp->pN[2])*(vp1->pN[2]-vp->pN[2]));
        }
    }
    return empty;
}

void Object::collectFaces(hashtable_v_double &hood,v_set &disabled,std::vector<Face*> &new_faces){
    new_faces.clear();
    // for each vertex in hood
    for(vdhm_iterator i=hood.begin();i!=hood.end();i++){
        // if vertex is thawed and not disabled
        if(!ifFrozen(hood,(*i).first) && disabled.find((*i).first)==disabled.end()){
            // for each adjacent face of thawed vertex
            for(std::vector<Face*>::iterator j=(*i).first->f.begin();j!=(*i).first->f.end();j++){
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

bool Object::ifFrozen(hashtable_v_double &neighborhood,Vertex *vv){
    // if vertex is in hashtable
    if(neighborhood.find(vv)!=neighborhood.end()){
//      cout << "Object::ifFrozen "
//      << "vertex cumul adge length " << neighborhood[vv]
//      << ", NEIGHBORHOOD_RADIUS " << NEIGHBORHOOD_RADIUS << endl;
        if (neighborhood[vv]<NEIGHBORHOOD_RADIUS){return false;}
        else {return true;}
    }
    /* else {
        cout << "\nObject::ifFrozen: Error. vertex not found in hashtable.\n";
        vv->printVertex(vv->o->name);
        cout << endl;
        exit(0);
    }*/
    return true;
}

bool Object::thawedAndAble(hashtable_v_double &hood,v_set &disabled){
    // for each vertex in hood
    for(vdhm_iterator i=hood.begin();i!=hood.end();i++){
        // if vertex is thawed and not disabled
//      bool mo = !ifFrozen(hood,(*i).first);
//      bool fo = disabled.find((*i).first)==disabled.end();
//      cout << "\nObject::thawedAndAble: "
//      << "thawed " << mo
//      << ", able " << fo
//      << endl;
        if(!ifFrozen(hood,(*i).first) && disabled.find((*i).first)==disabled.end()){
            // then affirm that at least one vertex is thawed and able
            return true;
        }
    }
    return false;
}

void Object::newFindNeighborhoods(void){
    // for each vertex in object
    for (std::vector<Vertex*>::iterator i=v.begin();i!=v.end();i++) {
        // initialize vertex*->double hash table
        // represents a neighbor vertex and
        // the shortest cumulative edge length
        // to reach it from current vertex
        hashtable_v_double hood;
        hood.clear();
        // add current vertex to neighbor list
        // naturally, assign it zero length
        hood[*i]=0.0;
        // initialize set of vertices to constitute disabled list
        v_set disabled;
        disabled.clear();
        // init collection of neighborhood faces
        std::vector<Face*> c;
//      c.assign((*i)->f.begin(),(*i)->f.end());

        ///// initial round /////
        // init collection of new faces
        std::vector<Face*> new_faces;
        collectFaces(hood,disabled,new_faces);
        // for each face in new collection
        for(std::vector<Face*>::iterator k=new_faces.begin();k!=new_faces.end();k++){
            // initialize container for edges
            // with neither vertex in hood
            std::vector<Edge*> bucket;
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
                std::vector<Edge*> pail;
                pail.clear();
                bool pail_empty = true;
                // for each edge in bucket
                for(std::vector<Edge*>::iterator j=bucket.begin();j!=bucket.end();j++){
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
//          c.push_back(*k);
        }

        ///// all subsequent rounds /////

        // while there are thawed vertices in neighbor list, hood,
        // that are also not disabled
        while(thawedAndAble(hood,disabled)){
            // init collection of new faces
            std::vector<Face*> new_faces_too;
            collectFaces(hood,disabled,new_faces_too);

/*          // print new face list, new_faces
            cout << "\n\n<<<<<<<<< new face list >>>>>>>>>>>\n";
            for(std::vector<Face*>::iterator k=new_faces.begin();k!=new_faces.end();k++){
                (*k)->printFace((*k)->v[0]->o->name);
                cout << endl;
            }
            // print hood vertex list, hood
            cout << "\n\n<<<<<<<<< hood vertex list >>>>>>>>>>>\n";
            for(vdhm_iterator j=hood.begin();j!=hood.end();j++){
                cout << " neighbor vertex " <<  (*j).first->index
                << ": length " << (*j).second << endl;
                (*j).first->printVertex((*j).first->o->name);
                cout << endl;
            }
            // print disabled vertex set, disabled
            cout << "\n\n<<<<<<<<< disabled vertex set >>>>>>>>>>>\n";
            for(vs_iterator j=disabled.begin();j!=disabled.end();j++){
                (*j)->printVertex((*j)->o->name);
                cout << endl;
            }
*/
            // for each face in new collection
            for(std::vector<Face*>::iterator k=new_faces_too.begin();k!=new_faces_too.end();k++){
                // initialize container for edges
                // with neither vertex in hood
                std::vector<Edge*> bucket;
                bucket.clear();
                bool bucket_empty = true;
                // for each edge in face
                for(int j=0;j<3;j++){
//                  if(!processEdge((*k)->e[j],hood,bucket)){
                    if(!processEdge((*k)->e[j],hood,bucket,*i)){
                        bucket_empty = false;
                    }
                }
                if(!bucket_empty){
                    // initialize another container for edges
                    // with neither vertex in hood
                    std::vector<Edge*> pail;
                    pail.clear();
                    bool pail_empty = true;
                    // for each edge in bucket
                    for(std::vector<Edge*>::iterator j=bucket.begin();j!=bucket.end();j++){
//                      if(!processEdge(*j,hood,pail)){
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

/*      // print hood
        if(!strcmp(name.c_str(),"d000_FILTERED_SMOOTH") && (*i)->index==18){
            cout << "\n\n****** current vertex *****" << endl;
            (*i)->printVertex((*i)->o->name);
            cout << endl;
            cout << "hash_map size " << hood.size() << endl;
            for(vdhm_iterator j=hood.begin();j!=hood.end();j++){
//          cout << "neighbor vertex " <<  (*j).first->index
//          << ": length " << (*j).second << endl;
//          (*j).first->printVertex((*j).first->o->name);
                (*j).first->printVertexCP();
                cout << endl;
            }
        }
*/
    }
}
/**/


void Object::findVertexAdjacencies(Controls &cs){
	std::vector<Edge*>::iterator i;
	std::vector<Face*>::iterator j;
	// for each edge, add edge* to both edge vertices
	for (i=e.begin();i!=e.end();i++){
		((*i)->v1)->e.push_back(*i);
		((*i)->v2)->e.push_back(*i);
	}
	// for each face, add face* to each face vertex
	for (j=f.begin();j!=f.end();j++){
		((*j)->v[0])->f.push_back(*j);
		((*j)->v[1])->f.push_back(*j);
		((*j)->v[2])->f.push_back(*j);
	}
	// if batch mode
	if(cs.interf || cs.sepdist){
		newFindNeighborhoods();
	}
}

// ######################################
// ######################################

void Object::boundObject(double* r) {
	std::vector<Vertex*>::iterator i;
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
	std::vector<Object*>::iterator i;
	// for each object, create edges
	for (i=o.begin();i!=o.end();i++) {
		(*i)->createEdges();
	}
}

void Container::findVertexAdjacencies(Controls &cs) {
	std::vector<Object*>::iterator i;
	// for each object, find vertex adjacencies
	for (i=o.begin();i!=o.end();i++) {
		(*i)->findVertexAdjacencies(cs);
	}
}

// #####################################################
// #####################################################

void Container::scanDir(Controls &cs,char *filename) {
	num_files = 0;
	struct dirent *pent;			// pointer to dirent structure

	if (cs.folder && !(cs.report && !cs.print)){
		cout << "\nFolder found " << filename << endl << endl;
	}

    DIR *pdir = opendir(filename);	// pointer to a directory data structure
    if (!pdir) {printf("Error. Could not open %s.\n",filename);exit(1);}
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
			if(!(cs.report && !cs.print)){
		        cout << "file found: " << str << "\n"; cout.flush();
			}
		}
    }
    closedir(pdir);

}
// #####################################################
// #####################################################

void Container::scanFiles(void) {
	Object *obj;
	char file[1024];
	// for each input file
	for (int count=0;count<num_files;count++) {
		// copy char array to string
        std::string str = files[count];
        // record object name
		std::string::size_type pos1 = str.find(".",0);
		if (!(pos1==std::string::npos)) {
			// ALLOCATE MEMORY FOR NEW OBJECT
			obj = new Object(str.substr(0,pos1));
		} else { cout << "Error! Object name was not found in " << str << "\n";exit(1);}
		// save pointer to object
		o.push_back(obj);
		// scan file
		sprintf(file,"%s",files[count].c_str());
		scanFile(obj,file);
	}
}

// #####################################################
// #####################################################

void Container::scanFile (Object *obj,char *filename) {
    char line[2048],*str;
    FILE *F;
    Vertex *v;
    Face *f;
//	int vertex_num=0,polygon_num=0;
//	std::vector<Vertex*> vp;
	hashtable_iv vp;
    // open file
    F = fopen(filename,"r");
    if (!F) { printf("Couldn't open input file %s\n",filename);exit(0);}
    // for every line in file
    for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) {
        // skip leading whitespace
        while (strchr(" \t,",*str)!=NULL) { str++;}
        // if first character is V for Vertex, add new linked list class instance
        if (strchr("V",*str)!=NULL){
//			vertex_num++;
			v=new Vertex(str,obj);
			obj->v.push_back(v);
//			vp.push_back(v);
			vp[v->index]=v;
		}
        // if first character is F for Face, add new linked list class instance
        else if (strchr("F",*str)!=NULL){
			f=new Face(str,vp,obj);
//			polygon_num++;
			obj->f.push_back(f);
		}
    }
    fclose(F);
}

void Container::writeDistances(int iteration) {
	char file[128];
	// create output filename
	if (append_iteration_number_to_distances) {
		sprintf(file,"%sclosest_point_distances_%i.dat",OUTPUT_DATA_DIR,iteration);
	} else {
		sprintf(file,"%sclosest_point_distances.dat",OUTPUT_DATA_DIR);
	}
	// open output file
	std::ofstream newfile (file,std::ios::out);
	if(newfile.is_open()){
		newfile.precision(4);
		// for each object
		for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
			// for each vertex in object
			for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
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
void Space::initBoxes(Container &c) {
	// subdivide space
	num_space[0] = (int) ceil( (c.world[1]-c.world[0])/SPACE_LENGTH );
	num_space[1] = (int) ceil( (c.world[3]-c.world[2])/SPACE_LENGTH );
	num_space[2] = (int) ceil( (c.world[5]-c.world[4])/SPACE_LENGTH );
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

void Container::boundWorld(void) {
	std::vector<Object*>::iterator i;
	double xmin,xmax,ymin,ymax,zmin,zmax,range[6];
	//initialize mins and maxes
	xmin = o[0]->v[0]->pN[0];
	xmax = o[0]->v[0]->pN[0];
	ymin = o[0]->v[0]->pN[1];
	ymax = o[0]->v[0]->pN[1];
	zmin = o[0]->v[0]->pN[2];
	zmax = o[0]->v[0]->pN[2];
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
	if (xmin<0) {world[0]=xmin*1.01;} else {world[0]=xmin*0.99;}
	if (xmax<0) {world[1]=xmax*0.99;} else {world[1]=xmax*1.01;}
	if (ymin<0) {world[2]=ymin*1.01;} else {world[2]=ymin*0.99;}
	if (ymax<0) {world[3]=ymax*0.99;} else {world[3]=ymax*1.01;}
	if (zmin<0) {world[4]=zmin*1.01;} else {world[4]=zmin*0.99;}
	if (zmax<0) {world[5]=zmax*0.99;} else {world[5]=zmax*1.01;}

}

void Space::boundWorld(Container& c) {
	std::vector<Object*>::iterator i;
	double xmin,xmax,ymin,ymax,zmin,zmax,range[6];
	//initialize mins and maxes
	xmin = c.o[0]->v[0]->pN[0];
	xmax = c.o[0]->v[0]->pN[0];
	ymin = c.o[0]->v[0]->pN[1];
	ymax = c.o[0]->v[0]->pN[1];
	zmin = c.o[0]->v[0]->pN[2];
	zmax = c.o[0]->v[0]->pN[2];
    ////////// loop through all objects //////////
    // for each object
	for (i=c.o.begin();i!=c.o.end();i++) {
		// get range of object vertices
		(*i)->boundObject(&range[0]);
		if (range[1]>xmax) {xmax = range[1];}
		if (range[0]<xmin) {xmin = range[0];}
		if (range[3]>ymax) {ymax = range[3];}
		if (range[2]<ymin) {ymin = range[2];}
		if (range[5]>zmax) {zmax = range[5];}
		if (range[4]<zmin) {zmin = range[4];}
    }
	if (xmin<0) {world[0]=xmin*1.01;} else {world[0]=xmin*0.99;}
	if (xmax<0) {world[1]=xmax*0.99;} else {world[1]=xmax*1.01;}
	if (ymin<0) {world[2]=ymin*1.01;} else {world[2]=ymin*0.99;}
	if (ymax<0) {world[3]=ymax*0.99;} else {world[3]=ymax*1.01;}
	if (zmin<0) {world[4]=zmin*1.01;} else {world[4]=zmin*0.99;}
	if (zmax<0) {world[5]=zmax*0.99;} else {world[5]=zmax*1.01;}

	if (print_flag) {
		cout << "\nworld bounds = [" 
			<< world[0] << " "
			<< world[1] << " "
			<< world[2] << " "
			<< world[3] << " "
			<< world[4] << " "
			<< world[5] << "]\n";
	}
}

// #####################################################
// #####################################################

void Space::computeBoxesToCheck(Face *f,std::vector<Box*> &bp) {
	std::vector<double> xv,yv,zv;
	double br[6];
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
	br[0] = xv[0];  // -x
	br[1] = xv[2];	//  x
	br[2] = yv[0];	// -y
	br[3] = yv[2];	//  y
	br[4] = zv[0];  // -z
	br[5] = zv[2];	//  z
	// collect boxes to check
	getBoxesFor3DLocations(br,bp);
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
	std::vector<int> face_flags;
	std::vector<Object*> tmp;
	std::vector<Face*>::iterator jj;
	std::vector<Face*> crossed_faces,unique_faces;
	std::pair<std::vector<Face*>::iterator,std::vector<Face*>::iterator> pp;
	std::pair<std::vector<Object*>::iterator,std::vector<Object*>::iterator> ppp;
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

void Container::findCrossed1(Space &s,Vertex *v,double lp[2][3],std::vector<Object*> &c){
	// find and return crossed objects between pN and extracellular point
	int num_odd_objects;
	std::vector<int> face_flags;
	std::vector<Face*>::iterator i;
	std::vector<Face*> crossed_faces,unique_faces;
	std::pair<std::vector<Face*>::iterator,std::vector<Face*>::iterator> pp;	
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

void Container::findCrossed2(Space &s,double lp[2][3],std::vector<Object*> &c){
	// find and return crossed objects between pN and extracellular point
	int num_odd_objects;
	std::vector<int> face_flags;
	std::vector<Face*> crossed_faces,unique_faces;
	collectNiceFaces(s,lp,unique_faces);
	findIntersectedFaces(lp,unique_faces,crossed_faces,face_flags);
	findOddMeshes(crossed_faces,face_flags,num_odd_objects,c);		
}

void Container::collectCrossed(Space &s,Vertex *v,std::vector<Object*> &cb){
	std::vector<Object*> ca;
	std::pair<std::vector<Object*>::iterator,std::vector<Object*>::iterator> pp;
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
	for (std::vector<Object*>::iterator ii=ca.begin();ii!=ca.end();ii++){
		pp=equal_range(cb.begin(),cb.end(),*ii);
		// if object in ca is found in cb, then remove from cb
		if(pp.first!=pp.second){cb.erase(pp.first);}
		// else add it to cb
		else {cb.push_back(*ii);}
	}
}

bool Container::updateNiceness(Vertex *v,std::vector<Object*> &cb){
	std::pair<std::vector<Object*>::iterator,std::vector<Object*>::iterator> pp;
	int old_nice = v->getVertexNiceness();
	// if vertex niceness changes then set flag=true
	bool flag = false;
	// if cb is not empty, then vertex is not nice
	if (!cb.empty()) {
		v->setVertexNiceness(1);
		// if vertex was nice
		if (!old_nice){
			flag = true;
//			nonnice++;
			pp=equal_range(cb.begin(),cb.end(),v->o);
			// if vertex is inside self object
			if(pp.first!=pp.second){
//				s_nonnice++;
				v->setVertexNiceness(2);
				cout << endl << endl 
				<< v->o->name << "->" << v->index 
				<< " vertex inside self [" << v->pN[0] << " "
				<< v->pN[1] << " "
				<< v->pN[2]
				<< "], cb.size() " << cb.size()
				<< endl;
				for (std::vector<Face*>::iterator jj=v->nf.begin();jj!=v->nf.end();jj++){
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
//			if(nonnice>0){nonnice--;}
		}
		// if vertex was at least nonnice to self
		else if (old_nice==2){
			flag=true;
//			if(nonnice>0){nonnice--;}
//			if(s_nonnice>0){s_nonnice--;}
		}
		// update niceness
		v->setVertexNiceness(0);			
	}
	return flag;
}

bool Container::checkNiceness(Space &s,Vertex *v) {
	std::vector<Object*> cb;
	// collect objects inside which vertex lies
	collectCrossed(s,v,cb);
	// update niceness of vertex based on cb
	return updateNiceness(v,cb);
}

void Container::findNice(Space &s) {
	cout << "Iteration 0: ";
	cout << "find nice vertices................";
	cout.flush();
	std::vector<Object*>::iterator i;
	std::vector<Vertex*>::iterator j;
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

void Container::collectNiceFaces(Space &s,double lp[2][3],std::vector<Face*> &uf) {
	int br[6];
	std::vector<Box*> b;
	std::vector<Box*>::iterator i;
	std::vector<Face*>::iterator j,new_end;
	
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
	s.getBoxesFor3DIndices(br,b,false);
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
	std::vector<Face*>::iterator j;
	bool line_flag, poly_flag, poly_edge_flag;
	// for each unique polygon
	for (j=uf.begin();j!=uf.end();j++) {
		// get face normal
		double n[3];
		(*j)->getNormal(n);
		checkLineFaceIntersection(*j,lp,n,line_flag,poly_flag,poly_edge_flag,false);
		// does point intersect polygon
		if (poly_flag) {
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

void Container::findOddMeshes(std::vector<Face*> &cf,std::vector<int> &ff,
							int& num_odd_objects,std::vector<Object*> &tmp) {
	// find mesh objects crossed an odd number of times by ray
	std::vector<Object*> ol; // object index list
	std::vector<Face*>::iterator i,j;
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
	// get number of unique vertices between current and other face
	int single_shared_vert[2]={-1,-1};
	int num_unique = numUniqueVertices(cf,of,single_shared_vert);
	bool share_edge_flag=false,identical_flag=false,share_vert_flag=false;
	if		(num_unique == 8) {share_vert_flag = true;}
	else if (num_unique == 7) {share_edge_flag = true;}
	else if (num_unique == 6) {identical_flag = true;}
	////////// begin decision tree //////////
	if (facesParallel(cf,of) && facesColinear(cf,of)) {
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
			// get face vertex coordinates
			double *opvc[3],*cpvc[3];
			cf->getVertexCoordinates(cpvc);
			of->getVertexCoordinates(opvc);
			int m=0,n=1,p=0,q=1;
			// single vertex shared
			if (single_shared_vert[0]==0){m=1;n=2;}
			else if (single_shared_vert[0]==1){n=2;}
			if (single_shared_vert[1]==0){p=1;q=2;}
			else if (single_shared_vert[1]==1){q=2;}
			double lp[2][3] = {{cpvc[m][0],cpvc[m][1],cpvc[m][2]},
								{cpvc[n][0],cpvc[n][1],cpvc[n][2]}};
			bool line_flag=false, poly_flag=false, poly_edge_flag;
			// get face normal
			double nn[3];
			of->getNormal(nn);
			checkLineFaceIntersection(of,lp,nn,line_flag,poly_flag,poly_edge_flag,true);
			// do faces intersect?
			if (line_flag && (poly_flag || poly_edge_flag)) {return(1);}
			lp[0][0] = opvc[p][0];
			lp[0][1] = opvc[p][1];
			lp[0][2] = opvc[p][2];
			lp[1][0] = opvc[q][0];
			lp[1][1] = opvc[q][1];
			lp[1][2] = opvc[q][2];
			// get face normal
			cf->getNormal(nn);
			checkLineFaceIntersection(cf,lp,nn,line_flag,poly_flag,poly_edge_flag,true);
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
	std::vector<Box*>::iterator i;
	// for each box in space, clear vector of face*
	for (i=b.begin();i!=b.end();i++) {
		delete (*i);
	}
	b.clear();
}

void Space::clearBoxes(void){
    std::vector<Box*>::iterator i;
    // for each box in space, clear vector of face*
    for (i=b.begin();i!=b.end();i++) {(*i)->f.clear();}
}

void Space::recordFace(vector<Box*> &ptr,Face* f) {
	std::vector<Box*>::iterator i;
	// for each box, add face
	for (i=ptr.begin();i!=ptr.end();i++) {(*i)->f.push_back(f);}
}

void Container::assignFacesToBoxes(Space &s) {
	std::vector<Object*>::iterator i;
	std::vector<Face*>::iterator j;
	std::vector<Box*> bp;
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
			// record boxes in face class
			(*j)->recordBoxes(bp);
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
	std::vector<Object*>::iterator i;
	std::vector<Vertex*>::iterator j;
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

struct vSortComp
{
    bool operator()(const vd_pair lhs,const vd_pair rhs)
    {
        //return (lhs->index < rhs->index);
        return (lhs.second > rhs.second);
    }
};
/*
void Container::getNonnice(hashset_v &target_vset) {
    // for each object* in container
    for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    	// for each vertex* in object
	    for (std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
			// if vertex is nonnice, then add to set
			if(!(*j)->o->vertexIsNice(*j)){target_vset.insert(*j);}
		}
	}
}*/
/*
void Container::getIntersectedVertices(hashset_v &target_vset) {
    // for each object* in container
    for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    	// for each face* in object
	    for (std::vector<Face*>::iterator j=(*i)->f.begin();j!=(*i)->f.end();j++) {
			// if current face has intersecting faces
			if((*j)->faceInTable_intf()){
				// then add face vertices to set
				target_vset.insert((*j)->v[0]);
				target_vset.insert((*j)->v[1]);
				target_vset.insert((*j)->v[2]);
			}
		}
	}
}*/

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
	if (!(br[0]*SPACE_LENGTH-v->pN[0])){br[0]--;}
	if (!(br[2]*SPACE_LENGTH-v->pN[1])){br[2]--;}
	if (!(br[4]*SPACE_LENGTH-v->pN[2])){br[4]--;}
	// screen range
	br[0]=s.screenIndex(br[0],"x");
	br[1]=s.screenIndex(br[1],"x");
	br[2]=s.screenIndex(br[2],"y");
	br[3]=s.screenIndex(br[3],"y");
	br[4]=s.screenIndex(br[4],"z");
	br[5]=s.screenIndex(br[5],"z");
	// add box pointers to vector
	s.getBoxesFor3DIndices(br,bp,false);
}

bool Container::faceInNeighborhood(Face *f,Vertex *v){
	// if face is in different object than vertex, then return false
	if(f->v[0]->o!=v->o){return false;}
	// else if face and vertex are in same object
	std::pair<std::vector<Face*>::iterator,std::vector<Face*>::iterator> p
		=equal_range(v->nf.begin(),v->nf.end(),f);
	// if face is in current vertex neighborhood
	if(p.first!=p.second){return true;}
	return false;
}

void Container::getCandidateFaces(vector<Box*> &bp,Vertex *v,hashset_f &cf){
	// for each box in search
	for (std::vector<Box*>::iterator j=bp.begin();j!=bp.end();j++) {
		// for each face in box
		for (std::vector<Face*>::iterator k=(*j)->f.begin();k!=(*j)->f.end();k++) {
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
	dd_pair p;
	// get vertex normal
	double n[3];
	v->getNormal(n);
	// get Box pointers for Face* collection
	std::vector<Box*> bp;
	getBoxes(bp,v,NUM_ADJACENT_BOXES,s);
	// collect Face pointers
	hashset_f cf;
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
	checkLineFaceIntersection(f,lp,n,line_flag,poly_flag,poly_edge_flag,true);
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

bool Face::getFaceIntersection(Container *c,bool just_check,std::vector<Face*> &int_f) {
	bool flag = false;
	// reset face element in table
//	clearFaceFromTable_intf();

	// for each box in which face lies, add faces in box to vector
	std::vector<Face*> of;
	for(std::vector<Box*>::iterator i=b.begin();i!=b.end();i++){
		of.insert(of.end(),(*i)->f.begin(),(*i)->f.end());
	}

	// sort and keep unique faces
	sort(of.begin(),of.end());
	std::vector<Face*>::iterator j = unique(of.begin(),of.end());
	of.assign(of.begin(),j);

	double bb[6];
	double *p1=v[0]->pN;
	double *p2=v[1]->pN;
	double *p3=v[2]->pN;
	threeValueSort(*p1,*p2,*p3, bb[1], bb[0]);
	p1++;p2++;p3++;
	threeValueSort(*p1,*p2,*p3, bb[3], bb[2]);
	p1++;p2++;p3++;
	threeValueSort(*p1,*p2,*p3, bb[5], bb[4]);

    // for each unique face
	j=of.begin();
	while(j!=of.end()){
		// if unique face is not same as current face
		if(*j!=this){
			if(faceBBsOverlap(bb,*j)){
				// if faces intersect
				if (c->checkFaceFaceIntersections(this,*j)) {
					if(!just_check){
//						c->ti++;
//						if (v[0]->o==(*j)->v[0]->o){c->si++;}
						// save intersecting face* to this face's intersecting face vector
						if(!faceInTable_intf()){addFaceToTable_intf();}
						addFaceToVector(*j);
					} else { int_f.push_back(*j); }
					// return
					flag = true;
				}
			}
		}
		j++;
	}
	return flag;
}

// #####################################################
// #####################################################

void Object::evalAttributes(Container *c,Controls &cs){
	closed=isClosed(cs);
	manifold=isManifold(cs);
	if(manifold){
		consistent=isConsistent(cs);
		if(consistent && closed){
			outward=isOutward(c);
		}
	}
}

bool Edge::isConsistent(void){
	if(f2!=NULL){ // not a border edge
		bool forward=false;
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

bool Object::rayIntersectsSide(char *str,double lp[3][3],double bb[6],double n[3]){
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

bool Object::rayIntersectsBB(double lp[3][3],Face *ff,double n[3]){
	///// compute face bounding box /////
	std::vector<double> xv,yv,zv;
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

bool Object::isOutward(Container *c){
	// assuming object mesh is closed, manifold, and consistent...
	//
	double count=0;
	// compute normal of first face in object
	Face *ff = f.front();
	double n[3];
	ff->getNormal(n);
	double L=sqrt( dot(n,n) );
	n[0]=n[0]/L;
	n[1]=n[1]/L;
	n[2]=n[2]/L;
	// ray origin = centroid of first adjacent face
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
	lp[1][0] = lp[0][0]+n[0]*2*fabs(c->world[1]-c->world[0]);
	lp[1][1] = lp[0][1]+n[1]*2*fabs(c->world[3]-c->world[2]);
	lp[1][2] = lp[0][2]+n[2]*2*fabs(c->world[5]-c->world[4]);
	// for each face in object
	for(std::vector<Face*>::iterator i=f.begin();i!=f.end();i++){
		// if face not same as first face in object
		if(*i!=ff){
			// if ray intersects face bounding box
			if(rayIntersectsBB(lp,ff,n)){
				// if hit, determine face-line intersection
				bool line_flag, poly_flag, poly_edge_flag;
				checkLineFaceIntersection(*i,lp,n,line_flag,poly_flag,poly_edge_flag,true);
				// does point intersect polygon
				if (poly_flag) {count++;}
				else if(poly_edge_flag) {count+=0.5;}	// expecting another edge hit
			}
		}
	}
	count = ceil(count);
	// if odd hit count then inward normal
	if(static_cast<int>(count)%2){return false;}
	// if even hit count then outward normal
	else {return true;}
}

bool Object::isConsistent(Controls &cs){
	// assuming object mesh is manifold...
	bool flag=true;
	// for each edge in object
	for(std::vector<Edge*>::iterator i=e.begin();i!=e.end();i++){
		// if the one or two adjacent edge faces traverse the edge
		// in the same direction, either both v1->v2 or both v2->v1,
		// then the entire set of mesh faces is declared inconsistent
		// i.e. the first edge with a flipped face will trigger signal return
		// and no more edge checking is performed, but -p option
		// needs ALL flipped edges.
		if(!(*i)->isConsistent()){
//			if(!cs.print){ return false;}
//			else {
				flag=false;
				// record offending edge
				cs.flipped.push_back(*i);
//			}
		}
	}
	// update cumulative # flipped edges
	cs.num_flip=cs.flipped.size();
	return flag;
}

bool Object::isClosed(Controls &cs){
	bool flag=true;
	// for each edge in object
	for(std::vector<Edge*>::iterator i=e.begin();i!=e.end();i++){
		if((*i)->f2==NULL){
//			if(!cs.print){ return false;}
//			else {
				flag=false;
				// record offending edge
				cs.border.push_back(*i);
//			}
		}
	}
	// update cumulative # border edges
	cs.num_bor+=cs.border.size();
	return flag;
}

bool Object::isManifold(Controls &cs){
	bool flag=true;

	///// if edges are manifold /////////////
	///// if one or two faces share an edge then manifold /////
	///// if three or more faces share an edge then NOT manifold /////
	// for each edge in object
	for(std::vector<Edge*>::iterator i=e.begin();i!=e.end();i++){
		if(!(*i)->fvec.empty()){
//			if(!cs.print){ return false;}
//			else {
				flag=false;
				// record offending edge
				cs.nonman_e.push_back(*i);
//			}
		}
	}

	///// if vertices are manifold /////
	////// confirm that all adjacent faces /////
	///// are consecutively reachable by edge hops /////
	// for each vertex in object
	for(std::vector<Vertex*>::iterator i=v.begin();i!=v.end();i++){
		// if vertex is not an orphan
		if(find(cs.orphan.begin(),cs.orphan.end(),*i)==cs.orphan.end()){
			std::vector<Face*> fc;
			fc.clear();
			// grab first adjacent face
			Face *f_f = (*i)->f.front();
			fc.push_back(f_f);
			Edge *ee=NULL,*eee=NULL;
			// if first face edge contains current vertex
			if(f_f->e[0]->v1==*i || f_f->e[0]->v2==*i){
				// then set ee
				ee = f_f->e[0];
			}
			// if second face edge contains current vertex
			if(f_f->e[1]->v1==*i || f_f->e[1]->v2==*i){
				if(ee==NULL){
					// then set ee
					ee = f_f->e[1];
				} else {
					// then set eee
					eee = f_f->e[1];
				}
			}
			// if third face edge contains current vertex
			if(f_f->e[2]->v1==*i || f_f->e[2]->v2==*i){
				if(eee==NULL){
					// then set eee
					eee = f_f->e[2];
				}
			}
	
			// initialize adjacent faces
			Face *ff=f_f,*fff=f_f;
			// while adjacent face disc has not been traced
			// i.e. until edge collection converges 
			while(ee!=eee){
				// identify new edge adjacent face that is not old adjacent edge face
				Face *ff_new=NULL,*fff_new=NULL;
				if(ee!=NULL){
					if(ee->f1!=ff){	ff_new = ee->f1;}
						else		{ ff_new = ee->f2;}
				} else {
					ff_new=NULL;
				}
				if(eee!=NULL){
					if(eee->f1!=fff){	fff_new = eee->f1;}
					else		{ fff_new = eee->f2;}
				} else {
					fff_new=NULL;
				}
				// collect faces
				if(ff_new!=fff_new){
					if(ff_new!=NULL){fc.push_back(ff_new);}
					if(fff_new!=NULL){fc.push_back(fff_new);}
				} else {
					if(ff_new!=NULL){fc.push_back(ff_new);}
					break;
				}
				////// identify new face adjacent edge that is not old adjacent edge /////
				// ee
				if(ff_new!=NULL){
					int tar = -1;
					// if first face edge contains current vertex
					if(ff_new->e[0]!=ee && (ff_new->e[0]->v1==*i || ff_new->e[0]->v2==*i)){
						tar=0;
					}
					// if second face edge contains current vertex
					else if(ff_new->e[1]!=ee && (ff_new->e[1]->v1==*i || ff_new->e[1]->v2==*i)){
						tar=1;
					}
					// if third face edge contains current vertex
					else {
						tar=2;
					}
					ee=ff_new->e[tar];
				} else {
					ee=NULL;
				}
				if(fff_new!=NULL){
					// eee
					int tar = -1;
					// if first face edge contains current vertex
					if(fff_new->e[0]!=eee && (fff_new->e[0]->v1==*i || fff_new->e[0]->v2==*i)){
						tar=0;
					}
					// if second face edge contains current vertex
					else if(fff_new->e[1]!=eee && (fff_new->e[1]->v1==*i || fff_new->e[1]->v2==*i)){
						tar=1;
					}
					// if third face edge contains current vertex
					else {
						tar=2;
					}
					eee = fff_new->e[tar];
				} else {
					eee=NULL;
				}
				// copy face info
				ff=ff_new;
				fff=fff_new;
			}
			// if the size of face collection != size of adjacent faces, then not manifold
			if(fc.size()!=(*i)->f.size()){
				flag=false;
				// record offending edge
				cs.nonman_v.push_back(*i);
			}
		}
	}
	// update cumulative # nonmanifold edges and vertices
	cs.num_nonman_e+=cs.nonman_e.size();
	cs.num_nonman_v+=cs.nonman_v.size();
	return flag;
}

bool Controls::gParse(char *str){
    char val[80];
    char *eptr;
    int i;
	if (strchr("-",*str)!=NULL){
		str++;
		if (strchr("g",*str)!=NULL){
			str++;
			if (strchr("=",*str)!=NULL){
				str++;
				// grab first threshold
				i=0;
				while (strchr("0123456789+-eE.",*str)!=NULL)
				{
				   val[i++] = *str++;
				}
				val[i]=0;
				t[0] = (int) strtod(val,&eptr);
				if (val==eptr)
				{
					t[0]=0;
					printf("Error in reading first threshold\n");
					exit(0);
				}

				// skip whitespace
				while (strchr(",",*str)!=NULL) { str++;}
					
			    // grab second threshold
			    i=0;
			    while (strchr("0123456789+-eE.",*str)!=NULL)
			    {
			        val[i++] = *str++;
			    }
			    val[i]=0;
			    t[1] = (int) strtod(val,&eptr);
			    if (val==eptr)
			    {
			        t[1]=0;
			        printf("Error in reading second threshold\n");
			        exit(0);
			    }

				// skip whitespace
				while (strchr(",",*str)!=NULL) { str++;}
				
			    // grab third threshold
			    i=0;
			    while (strchr("0123456789+-eE.",*str)!=NULL)
			    {
		        val[i++] = *str++;
			    }
			    val[i]=0;
			    t[2] = (int) strtod(val,&eptr);
			    if (val==eptr)
			    {
			        t[2]=0;
			        printf("Error in reading third threshold\n");
					exit(0);
				}
				report=true;
				return true;
			}else {return false;}
		}else {return false;}
	}else {return false;}
}

void Controls::parse(int argc,char **argv,char message[2048]){

	// if no arguments passed
	if(argc==1){
		cout << message << endl;
		exit(0);
	}
	// determine if file or folder was passed
	struct stat stat_p;
	stat(argv[1],&stat_p);
	if(stat_p.st_mode & S_IFDIR){ folder=true;}
	else{ folder=false;}

	for(int i=2;i<argc;i++){
		if(!strcmp(argv[i],"-a")){attr=true;}
		else if(!strcmp(argv[i],"-p")){print=true;}
		else if(!strcmp(argv[i],"-i")){interf=true;}
		else if(!strcmp(argv[i],"-s")){sepdist=true;}
		else {
			// check if NOT -g[t1,t2,t3]
			if(!gParse(argv[i])){
				// ERROR
				cout << "\n\nparse error. command line option = "
				<< argv[i] << endl << message << endl;
				exit(0);
			}
		}			
	}
	// enforce -a and -g mutual exclusivity
	if(attr==true && report==true){
		cout << "\n\nError. The '-a' and '-g' options are mutually exclusive,\n"
		<< "since '-g' incorporates mesh characteristics.\n\n";
		exit(0);
	}
	// enforce -i and -s application to folders
	if((interf==true || sepdist==true) && folder==false){
		cout << "\n\nError. The '-i' and '-s' options apply only to folders.\n\n";
		exit(0);
	}

}

// ********************************************************
// ********************************************************

void Object::checkIntegrity(Controls &cs){
	// check vertex integrity
	cs.orphanMissingContig(this);
	// check face integrity
	cs.DegenContig(this);
}

// VERTICES
void Controls::orphanMissingContig(Object *o){
	// build map #1: integer (vertex index) -> bool (flag)
	table_ib found;
	// build map #1: integer (vertex index) -> Vertex* (pointer)
	table_iv getPtr;
	//for each vertex in object
	for (std::vector<Vertex*>::iterator i=o->v.begin();i!=o->v.end();i++) {
		found[(*i)->index]=false;
		getPtr[(*i)->index]=*i;
	}
	///// missing /////
	// NOTE THIS IS REDUNDANT SINCE THIS INFORMATION IS GATHERED IN scanFile
	// ALSO NOTE MULTIPLE IDENTICAL FACE INDECES WITH MISSING VERTICES
	// WILL NOT BE HANDLED CORRECTLY. 
	//for each face in object
	for (std::vector<Face*>::iterator i=o->f.begin();i!=o->f.end();i++) {
//		cout << "Face " << (*i)->index << endl;
		// for each vertex index
		for(int j=0;j<3;j++){
//			int z=(*i)->v[j]->index;
			// if find index in map #1
			if((*i)->v[j]!=NULL){
				// set flag to true;
				int z=(*i)->v[j]->index;
				found[z]=true;
			} else {
				// missing vertex
				o->score=1;
				// add index to cs.missing_v
				missing_v.push_back(o->ii[(*i)->index]);
				// add face to cs.missing_f
				missing_f.push_back(*i);
			}
		}
	}

	// keep unique vertices
	sort(missing_v.begin(),missing_v.end());
	std::vector<int>::iterator new_end_v = unique(missing_v.begin(),missing_v.end());
	missing_v.assign(missing_v.begin(),new_end_v);
	// keep unique faces
	sort(missing_f.begin(),missing_f.end());
	std::vector<Face*>::iterator new_end_f = unique(missing_f.begin(),missing_f.end());
	missing_f.assign(missing_f.begin(),new_end_f);

	// update cumulative # missing vertices
	num_missing+=missing_v.size();
	///// orphan ////
	// for each element of map #1
	for(ib_iterator i=found.begin();i!=found.end();i++){
		if((*i).second==false){
			// orphan vertex
			// use map #2 to add vertex* to orphan vector
			orphan.push_back(getPtr[(*i).first]);
		}
	}
	// update cumulative # orphan vertices
	num_orphan+=orphan.size();
	///// contiguousness /////
	contig_v=true;
	// since map #1 is automatically sorted by integer value
	int pace = 1;
	// for each other element in map
	for(ib_iterator i=found.begin();i!=found.end();i++){
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
			} else if (pace==static_cast<int>(o->v.size())-1) {
				i--;i--;
				vec_cv[0]=(*i).first; i++;
				vec_cv[1]=(*i).first; i++;
				vec_cv[2]=(*i).first; i++;
				vec_cv[3]=(*i).first; i++;
				vec_cv[4]=-1;
				i--;i--;
			} else if (pace==static_cast<int>(o->v.size())) {
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
}

// FACES
void Controls::areaAspectRatio(Object *o){
	// for each face in object
	for(std::vector<Face*>::iterator i=o->f.begin();i!=o->f.end();i++){
		area.n++;
		aspect_ratio.n++;
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
		// add to vector
		aspect_ratio.x.push_back(ar);
	}
	// update cumulative surface area
	a+=area.sum;
	// build face area histogram
//	cout << "area [min,max]=[" << area.min << "," << area.max << "]"<< endl;
	area.createHistogram();
//	cout << "area [min,max]=[" << area.min << "," << area.max << "]"<< endl;
	// build aspect ratio histogram
	aspect_ratio.createHistogram();
}



double Face::getAspectRatio(void){

	/* Make triangle edge vectors */
	double va[3]={v[1]->pN[0]-v[0]->pN[0],v[1]->pN[1]-v[0]->pN[1],v[1]->pN[2]-v[0]->pN[2]};
	double vb[3]={v[2]->pN[0]-v[1]->pN[0],v[2]->pN[1]-v[1]->pN[1],v[2]->pN[2]-v[1]->pN[2]};
	double vc[3]={v[0]->pN[0]-v[2]->pN[0],v[0]->pN[1]-v[2]->pN[1],v[0]->pN[2]-v[2]->pN[2]};
	double vbase[3];
	double vopp[3];

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

void Controls::DegenContig(Object *o){
	// buide map: face index (integer)->face*
	table_if tt;
	// for each face in object
	for(std::vector<Face*>::iterator i=o->f.begin();i!=o->f.end();i++){
//		tt[(*i)->index]=*i;
		tt.insert(std::make_pair((*i)->index,*i));
	}
	contig_f=true;
	int pace=1;
	// for each face in map
	for(if_iterator i=tt.begin();i!=tt.end();i++){
		Face *f=(*i).second;
		///// degeneracy /////
		// if any two of the vertex indices are the same
		if( f->v[0]==f->v[1] || f->v[0]==f->v[2] || f->v[1]==f->v[2]){
			// then face is degenerate
			degen.push_back(f);
			num_degen++;
			o->score=1;
		}
		// update cumulative # degenerate faces
//		num_degen+=degen.size();
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
				} else if (pace==static_cast<int>(o->f.size())-1) {
					i--;i--;
					vec_cf[0]=(*i).first; i++;
					vec_cf[1]=(*i).first; i++;
					vec_cf[2]=(*i).first; i++;
					vec_cf[3]=(*i).first; i++;
					vec_cf[4]=-1;
					i--;i--;
				} else if (pace==static_cast<int>(o->f.size())) {
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

void Object::setAll(Edge *ee,hashtable_fi &group,int &g){
	group[ee->f1]=g;
	// if second adjacent face
	if(ee->f2!=NULL){ group[ee->f2]=g; }
	// if more adjacent faces
	if(!ee->fvec.empty()){
		// for each adjacent face
		for(std::vector<Face*>::iterator i=ee->fvec.begin();i!=ee->fvec.end();i++){
			// if adjacent face has group
			group[*i]=g;
		}
	}
	g++;
}

void Object::getGroups(Edge *ee,hashtable_fi &group,i_set &s){
	s.clear();
	s.insert(group[ee->f1]);
	// if second adjacent face
	if(ee->f2!=NULL){ s.insert(group[ee->f2]); }
	// if more adjacent faces
	if(!ee->fvec.empty()){
		// for each adjacent face
		for(std::vector<Face*>::iterator i=ee->fvec.begin();i!=ee->fvec.end();i++){
			// if adjacent face has group
			s.insert(group[*i]);
		}
	}
}

int Object::getLowest(i_set &s){
	int i = 100000;
	// for each element in set
	for(is_iterator j=s.begin();j!=s.end();j++){
		// if element is lower than i and not 0
		if(*j<i && *j){i=*j;}
	}
	return i;
}

void Object::replaceGroups(hashtable_fi &group,i_set &s,int z){
	// for each element of set
	for(is_iterator i=s.begin();i!=s.end();i++){
		// if element is not z
		if(*i!=z){
			// for each element in group
			for(fi_iterator j=group.begin();j!=group.end();j++){
				// if group is element then replace with z
				if((*j).second==*i){ (*j).second=z;
				}
			}
		}
	}
}

int Object::countComponents(void){
	int g=0;
	// cp edge vector to ws
	std::vector<Edge*> ws;
	ws.assign(e.begin(),e.end());
	///// create hashed map face*->integer (group #) /////
	hashtable_fi group;
	// for each face
	for(std::vector<Face*>::iterator i=f.begin();i!=f.end();i++){
		group[*i]=g;
	}
	g++;
	i_set s;
	while (!ws.empty()){
		// pop edge from ws
		Edge *ee=ws.back();
		ws.pop_back();
		// get edge groups
		getGroups(ee,group,s);
		is_iterator i=s.begin();
		// if no edge adjacent face has been set
		if(s.size()==1 && *i==0){
			// set all to next available group #
			setAll(ee,group,g);
		} else if(s.size()>1){ // more than one group
			// identify lowest group # larger than 0 in set
			int z=getLowest(s);
			// replace every other group # with lowest group #
			replaceGroups(group,s,z);
		} // else if all set and same then do nothing
	}
	// analyze map
	s.clear();
	// for each element in group
	for(fi_iterator i=group.begin();i!=group.end();i++){
		s.insert((*i).second);
	}
	return s.size();
}

void Controls::findIntersectingFaces(Container *c,Object *o){
	std::vector<Face*> dummy;
	// for each face in object
	for(std::vector<Face*>::iterator i=o->f.begin();i!=o->f.end();i++){
		// check for intersections
    	(*i)->getFaceIntersection(c,false,dummy);
		// intersecting faces end up in o->intf (hashtable_f_face)
	}
	// update cumulative # intersecting faces
	num_intf+=o->intf.size();
}
// EDGES
void Object::countBoundaries(Controls &cs){
	Edge *ee=NULL;
	Vertex* end=NULL,*begin=NULL;
	bool open_boundary = false;
	num_bou=0;
	// copy the border vector to ws
	std::vector<Edge*> ws;
	ws.assign(cs.border.begin(),cs.border.end());
	// while there are border edges in ws
	while(!ws.empty()){
		// if !open_boundary, i.e.e=NULL;end=NULL;begin=NULL;
		if(!open_boundary){
			// start a boundary
			open_boundary=true;
			// pop first edge from ws
			ee=ws.back();
			ws.pop_back();
			end=ee->v1;
			begin=ee->v2;
		}
		// else use open boundary, i.e. end,begin,e are valid
		// for each edge in ws
		for(std::vector<Edge*>::iterator i=ws.begin();i!=ws.end();i++){
			// if v1 or v2 are begin
			if((*i)->v1==begin || (*i)->v2==begin){
				// if v1 or v2 are end
				if((*i)->v1==end || (*i)->v2==end){
					// erase edge from ws
					ws.erase(i);
					ee=NULL;end=NULL;begin=NULL;
					num_bou++;
					break;
				} else {
					// erase edge from ws
					ws.erase(i);
					// begin=other vertex;
					if((*i)->v1==begin){ begin=(*i)->v2; }
					else { begin=(*i)->v1; }
					break;
				}
			}
		}
	}
	// update cumulative number boundaries
	cs.num_bou+=num_bou;
}

void Controls::processEdgeLengths(Object *o){
	// for each edge in object
	for(std::vector<Edge*>::iterator i=o->e.begin();i!=o->e.end();i++){
		double l = (*i)->l;
		edge_length.n++;
		edge_length.sum+=l;
		edge_length.sum2+=l*l;
		edge_length.total+=l;
		if(l<edge_length.min) edge_length.min=l;
		if(l>edge_length.max) edge_length.max=l;
		// check distinguishability
		if( !distinguishable((*i)->v1->pN[0],(*i)->v2->pN[0]) &&
			!distinguishable((*i)->v1->pN[1],(*i)->v2->pN[1]) &&
			!distinguishable((*i)->v1->pN[2],(*i)->v2->pN[2]) 
			){ indistin.push_back(*i); }
		// add to vector
		edge_length.x.push_back(l);
	}
	edge_length.createHistogram();
}

// OTHER
int Controls::isGood(Object *o){
	return	// closed
			o->closed &&
			// manifold
			o->manifold &&
			// consistently oriented faces
			o->consistent &&
			// outwardly oriented face normals
			o->outward &&
			// no intersecting faces
			!num_intf &&
			// no orphaned vertices
			orphan.empty() &&
			// no missing vertices
			missing_v.empty() &&
			// no degenerate faces
			degen.empty() &&
			// contiguous vertex indexing
			contig_v &&
			// contiguous face indexing
			contig_f &&
			// max aspect ratio < threshold
			(aspect_ratio.max < t[0]) &&
			// min edge angle > threshold
			(edge_angle.min > t[1]) &&
			// min edge length > threshold
			(edge_length.min > t[2]) &&
			// all edge vertices are distinguishable
			indistin.empty() &&
			// all vertices are distinguishable
			indistin_v.empty();
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
void Controls::computeEdgeAngles(Object *o){
	// for each edge in object
	for(std::vector<Edge*>::iterator i=o->e.begin();i!=o->e.end();i++){
		// if edge has exactly two adjacent faces
		if((*i)->f2!=NULL && (*i)->fvec.empty()){
			double angle = (*i)->getAngle();
			edge_angle.n++;
			edge_angle.sum+=angle;
			edge_angle.sum2+=angle*angle;
			edge_angle.total+=angle;
			if(angle<edge_angle.min) edge_angle.min=angle;
			if(angle>edge_angle.max) edge_angle.max=angle;
			// add to vector
			edge_angle.x.push_back(angle);
		}
	}
	edge_angle.createHistogram();
}

void Controls::computeVolume(Object *o){
	/*
	Calculate volume of a manifold surface mesh by summing over
	all triangles of the surface mesh, the quantity:

              |  x1   x2   x3  |
    (1/6) det |  y1   y2   y3  |
              |  z1   z2   z3  |

	where (xn,yn,zn) is vertex n of a surface triangle.  The surface normals of
	all triangles of the surface mesh must be consistently oriented.
	Volume will be positive if all normals face outward, and negative if
	all normals face inward.
	*/

	o->vol = 0.0;

	// for each face in object
	for(std::vector<Face*>::iterator i=o->f.begin();i!=o->f.end();i++){
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
		o->vol+=det;
	}
  	o->vol=o->vol/6.0;
	// update cumulative volume
	vol+=o->vol;
}

int Container::getScore(void){
	int t=0;
    // for each object in container
	for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
		t+=(*i)->score;
	}
	return t;
}

void Object::vertexDistin(Controls &cs){
	///// check vertex distinguishability /////
	// multimap: double -> Vertex*
	table_v mm;
	// random vector
	double rand_vec[3]={0.236416584579274058342,
						0.927225593011826276779,
						0.389099507126957178116};
	// for each vertex
	for(std::vector<Vertex*>::iterator i=v.begin();i!=v.end();i++){
		// add dot product of vertex and random vector to multimap
		mm.insert(std::make_pair(dot((*i)->pN,rand_vec),*i));
	}
	// for each multimap element
	tv_iterator j;
	for(tv_iterator i=mm.begin();i!=mm.end();i++){
		j=i;j++;
		if(j!=mm.end()){
			// if sequential pair of multimap elements is not distinguishable
			if( !distinguishable((*i).second->pN[0],(*j).second->pN[0]) &&
				!distinguishable((*i).second->pN[1],(*j).second->pN[1]) &&
				!distinguishable((*i).second->pN[2],(*j).second->pN[2]) 
				){
				cs.indistin_v.push_back((*i).second);
				cs.indistin_v.push_back((*j).second);
				cs.num_indistin+=2;
			}
		}
	}
}

void Object::evalCharacteristics(Container* c,Controls &cs){
	// vertices
	boundObject(cs.bb);
	vertexDistin(cs);
	// faces
	num_sep=countComponents();
	cs.num_sep+=num_sep;
	cs.areaAspectRatio(this);
	// if not batch mode
	if(cs.interf==false & cs.sepdist==false){
//		cout << "find intersecting faces.......";cout.flush();
		cs.findIntersectingFaces(c,this);
//		cout << "complete.\n";
	}
	// edges
	countBoundaries(cs);
	cs.processEdgeLengths(this);

	if(cs.report){ score = cs.isGood(this); }
	if(manifold && consistent){
		// manifold and consistently oriented face normals
//		cout << "analyze edge angles\n";
		cs.computeEdgeAngles(this);
		if(closed){
			// closed, manifold, and consistently oriented face normals
//			cout << "analyze volume\n";
			cs.computeVolume(this);
			// if number of components == 1
			// FUTURE IMPROVEMENT: only require orientable, not oriented
			// FUTURE IMPROVEMENT: separate components and compute genus of each component
			if(num_sep==1 && cs.orphan.empty()){
				computeGenus();
			}
		}
	}
}

void Object::analyze(Container *c,Controls &cs){
	// single object
	// evaluate mesh attributes
	evalAttributes(c,cs);
	// eval mesh characteristics
	if(!cs.attr){evalCharacteristics(c,cs);}
}

void Object::analyzeBatch(Container *c,Controls &cs,Space &s){
	// find intersecting faces
	if(cs.interf){
		cs.findIntersectingFaces(c,this);
	}
	// compute separation distance 
	if(cs.sepdist){
    	c->findNice(s);
	    c->getSeparationDistances(s);
	}
}

void Object::printChars(Controls &cs){
	cout << "\nMESH CHARACTERISTICS\n";
	cout << "number of vertices: " << v.size() << endl
	<< "number of faces: " << f.size() << endl
	<< "number of edges: " << e.size() << endl
	<< "number of components: " << num_sep << endl;
	if (manifold==false){
		cout << "number of boundaries: Since object is nonmanifold,\n"
		<< "number of boundaries: the number of boundaries may be underestimated.\n";
	}
	if (cs.border.empty()){
		cout << "number of boundaries: none\n";
	} else {
		cout << "number of boundaries: " << cs.border.size() << endl;
		//	if -p option, print offending
		if(cs.print){
			int j=1;
			// for each border edge
			for(std::vector<Edge*>::iterator i=cs.border.begin();i!=cs.border.end();i++){
				cout << "number of boundaries: boundary edge " << j++ << endl;
				(*i)->printEdge((*i)->o->name);
				cout << endl;
			}
		}
	}

	if (cs.indistin_v.empty()){
		cout << "number of indistinguishable vertices: none\n";
	} else {
		cout << "number of indistinguishable vertices: " << cs.indistin_v.size() << endl;
		//	if -p option, print offending
		if(cs.print){
			int j=1;
			// for each indistinguishable vertec
			for(std::vector<Vertex*>::iterator i=cs.indistin_v.begin();i!=cs.indistin_v.end();i++){
				cout << "number of  indistinguishable vertices: vertex " << j++ << endl;
				(*i)->printVertex((*i)->o->name);
				cout << endl;
			}
		}
	}

	// if volume computed
	if(closed==true && manifold==true && consistent==true){
		cout << "object volume: " << vol << endl;
	} else {
		cout << "object volume: not computed, since ";
		if(closed==false){cout << "not closed,";}
		if(consistent==false){cout << "not consistent,";}
		if(manifold==false){cout << "not manifold";}
		cout << endl;
	}
	// if genus computed
	if(closed==true && manifold==true && consistent==true && num_sep==1 && cs.orphan.empty()){
		cout << "object genus: " << genus << endl;
	} else {
		cout << "object genus: not computed, since ";
		if(closed==false){cout << "not closed,";}
		if(consistent==false){cout << "not consistent,";}
		if(manifold==false){cout << "not manifold,";}
		if(num_sep>1){cout << "#components=" << num_sep << ",";}
		if(!cs.orphan.empty()){cout << "orphan vertices were found";}
		cout << endl;
	}

	// bounding box
	cout << "bounding box: [xmin,xmax,ymin,ymax,zmin,zmax]\n";
	cout << "bounding box: ["
	<< cs.bb[0] << " , "
	<< cs.bb[1] << " , "
	<< cs.bb[2] << " , "
	<< cs.bb[3] << " , "
	<< cs.bb[4] << " , "
	<< cs.bb[5] << "]" << endl << endl;
	// face area
	cout << "Face area statistics:" << endl;
	cout << "   total    " << cs.area.total << endl;
	cs.area.printStats();
	cout << "Face area histogram:" << endl;
	cs.area.printHistogram();
	cout << endl;
	
	// face aspect ratio
	cout << "Face aspect ratio statistics:" << endl;
	cs.aspect_ratio.printStats();
	cout << "Face aspect ratio histogram:" << endl;
	cs.aspect_ratio.printHistogram();
	printf("  (Aspect ratio is longest edge divided by shortest altitude)\n");
	printf("\n");

	// edge length
	cout << "Edge length statistics:" << endl;
	cs.edge_length.printStats();
	cout << "Edge length histogram:" << endl;
	cs.edge_length.printHistogram();
	cout << endl;
	if(!cs.indistin.empty()){
		cout << "# edges with indistinguishable vertices: " << cs.indistin.size() << endl;
		// for each afflicted edge
		for(std::vector<Edge*>::iterator i=cs.indistin.begin();i!=cs.indistin.end();i++){
			(*i)->printEdge((*i)->o->name);
			cout << endl;
		}
	}

	// if edge angles computed
	if(manifold==true && consistent==true){
		cout << "Edge angle statistics:         (radians)" << endl;
		cs.edge_angle.printStats();
		cout << "Edge angle histogram:" << endl;
		cs.edge_angle.printHistogram();
		cout << endl;
	} else {
		cout << "edge angles: not computed, since ";
		if(consistent==false){cout << "not consistent,";}
		if(manifold==false){cout << "not manifold";}
		cout << endl;
	}
	// if not batch mode
	if(cs.interf==false & cs.sepdist==false){
		// intersecting faces
		if (intf.empty()){
			cout << "# intersecting faces: none\n\n";
		} else {
			cout << "# intersecting faces: " << intf.size() << endl;
			//	if -p option, print offending
			if(cs.print){
				int j=1;
				// for each intersected face
				for(ff_iterator i=intf.begin();i!=intf.end();i++){
					cout << "# intersecting faces: intersected face " << j++ << endl;
					// print intersected face
					(*i).first->printFace((*i).first->v[0]->o);
					// print intersecting faces
					for(std::vector<Face*>::iterator k=(*(*i).second).begin();k!=(*(*i).second).end();k++){
						(*k)->printFace((*k)->v[0]->o);
						cout << endl;
					}
				}
				cout << endl;
			}
		}
	}
}

bool Controls::binaryOutput(void){
	if (
	// tffttf folder -g -i
	(folder==true && attr==false && print==false && report==true && interf==true && sepdist==false) ||
	// tfftff folder -g
	(folder==true && attr==false && print==false && report==true && interf==false && sepdist==false) ||
	// ffftff file -g
	(folder==false && attr==false && print==false && report==true && interf==false && sepdist==false) 
	) { return true;}
	else {return false;}
}

void Controls::printIntegrity(void){
	cout << "\nMESH FILE INTEGRITY\n";
	// orphan vertices
	if(orphan.empty()){
		cout << "# orphan vertices: none\n";
	} else {
		cout << "# orphan vertices: " << orphan.size() << endl;
		//	if -p option, print offending
		if(print){
			int j=1;
			// for each orphan vertex
			for(std::vector<Vertex*>::iterator i=orphan.begin();i!=orphan.end();i++){
				cout << "# orphan vertices: orphan vertex " << j++ << endl;
				(*i)->printVertex((*i)->o->name);
				cout << endl;
			}
		}
	}
	// missing vertices
	if(missing_v.empty()){
		cout << "# missing vertices: none\n";
	} else {
		cout << "# missing vertices: " << missing_v.size() << endl;
		//	if -p option, print offending
		if(print){
			int j=1;
			// for each missing vertex
			for(std::vector<int>::iterator i=missing_v.begin();i!=missing_v.end();i++){
				cout << "# missing vertices: #" << j++
				 << "-> missing vertex index " << *i << endl;
			}
			j=1;
			// for each missing face
			for(std::vector<Face*>::iterator i=missing_f.begin();i!=missing_f.end();i++){
				cout << "# missing vertices: affected face " << j++ << endl;
				if((*i)->v[0]!=NULL){
					(*i)->printFace((*i)->v[0]->o);
					cout << endl;
				} else if((*i)->v[1]!=NULL){
					(*i)->printFace((*i)->v[1]->o);
					cout << endl;
				} else if((*i)->v[2]!=NULL){
					(*i)->printFace((*i)->v[2]->o);
					cout << endl;
				}
			}
		}
	}
	// degenerate faces
	if(degen.empty()){
		cout << "# degenerate faces: none\n";
	} else {
		cout << "# degenerate faces: " << degen.size() << endl;
		//	if -p option, print offending
		if(print){
			int j=1;
			// for each degenerate face
			for(std::vector<Face*>::iterator i=degen.begin();i!=degen.end();i++){
				cout << "# degenerate faces: affected face " << j++ << endl;
				(*i)->printFace((*i)->v[0]->o);
				cout << endl;
			}
		}
	}

	// contiguous numbering
	if(contig_v){
		cout << "contiguous vertex indexing from 1: yes\n";
	} else {
		cout << "contiguous vertex indexing from 1: no\n";
		cout << "contiguous vertex indexing from 1: bad index and +/- 2\n";
		cout << "contiguous vertex indexing from 1: ";
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
		cout << "contiguous face indexing from 1: yes\n";
	} else {
		cout << "contiguous face indexing from 1: no\n";
		cout << "contiguous face indexing from 1: bad index and +/- 2\n";
		cout << "contiguous face indexing from 1: ";
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
	cout << "\nMESH ATTRIBUTES\n";
	// closed
	if(closed==true){
		cout << "mesh is closed: yes\n";
	} else {
		cout << "mesh is closed: no\n";
		//	if -p option, print offending
		if(cs.print){
			int j=1;
			cout << "mesh is closed: # border edges - " << cs.border.size() << endl;
			// for each border edge
			for(std::vector<Edge*>::iterator i=cs.border.begin();i!=cs.border.end();i++){
				cout << "mesh is closed: border edge # " << j++ << endl;
				(*i)->printEdge((*i)->o->name);
				cout << endl;
			}
		}
	}
	// manifold
	if(manifold==true){
		cout << "mesh is manifold: yes\n";
	} else {
		cout << "mesh is manifold: no\n";
		//	if -p option, print offending
		if(cs.print){
			if(cs.nonman_v.empty()){
				cout << "mesh is manifold: # nonmanifold vertices - none\n";
			} else {
				int j=1;
				cout << "mesh is manifold: # nonmanifold vertices - " << cs.nonman_v.size() << endl;
				// for each nonmanifold vertex
				for(std::vector<Vertex*>::iterator i=cs.nonman_v.begin();i!=cs.nonman_v.end();i++){
					cout << "mesh is manifold: nonmanifold vertex # " << j++ << endl;
					(*i)->printVertex((*i)->o->name);
					cout << endl;
				}
			}
			if(cs.nonman_e.empty()){
				cout << "mesh is manifold: # nonmanifold edges - none\n";
			} else {
				int j=1;
				cout << "mesh is manifold: # nonmanifold edges - " << cs.nonman_e.size() << endl;
				// for each nonmanifold edge
				for(std::vector<Edge*>::iterator i=cs.nonman_e.begin();i!=cs.nonman_e.end();i++){
					cout << "mesh is manifold: nonmanifold edge # " << j++ << endl;
					(*i)->printEdge((*i)->o->name);
					cout << endl;
				}
			}
		}
	}
	// consistent
	if(manifold==false){
			cout << "mesh has consistently oriented face normals: uncomputable since notmanifold\n";
	} else {
		if(consistent==true){
			cout << "mesh has consistently oriented face normals: yes\n";
		} else {
			cout << "mesh has consistently oriented face normals: no\n";
			//	if -p option, print offending
			if(cs.print){
				int j=1;
				cout << "mesh has consistently oriented face normals: # flipped edges - " << cs.flipped.size() << endl;
				// for each flipped edge
				// FUTURE IMPROVEMENT: if edge is nonmanifold then exclude from flipped list
				for(std::vector<Edge*>::iterator i=cs.flipped.begin();i!=cs.flipped.end();i++){
					cout << "mesh has consistently oriented face normals: flipped edge # " << j++ << endl;
					(*i)->printEdge((*i)->o->name);
					cout << endl;
				}
			}
		}
	}
	// outward
	if(manifold==false || consistent==false || closed==false){
		cout << "mesh has outward oriented face normals: uncomputable since ";
		if(closed==false){cout << "not closed,";}
		if(consistent==false){cout << "not consistent,";}
		if(manifold==false){cout << "not manifold";}
		cout << endl;
	} else {
		if(outward==true){
			cout << "mesh has outward oriented face normals: yes\n";
		} else {
			cout << "mesh has outward oriented face normals: no\n";
		}
	}
}

void Object::print(Controls &cs,bool badmesh){
	// if not only -g option
	if(!cs.binaryOutput()){
		cout << "\n\n" << "/* ********************** OBJECT ********************** */\n";
		//	print object name 
		cout << "name: " << name << endl;
		//	print Integrity
		cs.printIntegrity();
		if(badmesh==true){
			cout << "\n\nWarning: Attributes and characteristics were not evaluated,\n"
			<< " since mesh file failed the integrity check.\n\n";
		} else {
			//	print attributes
			printAttr(cs);
			//	print characteristics
			if(cs.attr==false) { 
				printChars(cs);
			}
			// print goodness
			if(cs.report==true){
				if(cs.isGood(this)){
					cout << "mesh is good: yes\n";
				} else {
					cout << "mesh is good: no\n";
				}
			}
		}
	} else {
		// print goodness
		if(cs.folder==false){
			cout << cs.isGood(this) << endl;
		}
	}
}

void Controls::printCumulative(bool badmesh,Object *o){
	// ISSUE WARNING IF ANY OBJECT FAILED INTEGRITY TEST OR LACK ATTRIBUTES
	// NOTE CUMULATIVE VOLUME ASSUMES ALL MESHES HAVE SAME ORIENTATION
	// if not only -g option
	if(!binaryOutput()){
		cout << "\n\n" << "/* ********************** SET OF OBJECTS ********************** */\n";
		cout << "INTEGRITY:\n";
		cout << "number of orphan vertices: " << num_orphan << endl;
		cout << "number of missing vertices: " << num_missing << endl;
		cout << "number of degenerate faces: " << num_degen << endl;
		cout << "ATTRIBUTES:\n";
		if(badmesh==true){
			cout << "Warning: These attribute summaries are inaccurate,\n"
			<< " since some mesh files failed the integrity check.\n";
		}
		cout << "number of border edges: " << num_bor << endl;
		cout << "number of flipped edges: " << num_flip << endl;
		cout << "number of nonmanifold edges: " << num_nonman_e << endl;
		cout << "number of nonmanifold vertices: " << num_nonman_v << endl;
		//	print characteristics
		if(attr==false) {
			if(badmesh==true){
				cout << "Warning: These characteristics summaries are inaccurate,\n"
				<< " since some mesh files failed the integrity check.\n";
			}
			cout << "CHARACTERISTICS:\n";
			cout << "number of objects: " << num_obj << endl;
			cout << "number of vertices: " << num_vert << endl;
			cout << "number of faces: " << num_faces << endl;
			cout << "number of edges: " << num_edges << endl;
			cout << "mesh volume: " << vol << endl;
			cout << "mesh surface area: " << a << endl;
			cout << "number of components: " << num_sep << endl;
			cout << "number of boundaries: " << num_bou << endl;
			cout << "number of indistinguishable vertices: " << num_indistin << endl;
		}
		if(report==true){
			cout << "number of meshes that are NOT good: " << num_score << endl;
		}
		cout << endl << endl;
	} else {
		if(num_score || !isGood(o)){
			cout << "0" << endl;
		} else {
			cout << "1" << endl;
		}
	}
}

void Container::printBatch(Controls &cs){
	// if not only -g option
	if(!cs.binaryOutput()){
		cout << "\n\n" << "/* ********************** BATCH RESULTS ********************** */\n";
		// intersecting faces
		if (cs.num_intf==0){
			cout << "# intersecting faces: none\n";
		} else {
			cout << "# intersecting faces: " << cs.num_intf << endl;
			//	if -p option, print offending
			if(cs.print){
				int j=1;
				// for each object
				for(std::vector<Object*>::iterator m=o.begin();m!=o.end();m++){
					// for each intersected face
					for(ff_iterator i=(*m)->intf.begin();i!=(*m)->intf.end();i++){
						cout << "# intersecting faces: intersected face " << j++ << endl;
						// print intersected face
						(*i).first->printFace((*i).first->v[0]->o);
						// print intersecting faces
						for(std::vector<Face*>::iterator k=(*(*i).second).begin();k!=(*(*i).second).end();k++){
							(*k)->printFace((*k)->v[0]->o);
							cout << endl;
						}
					}
				}
			}
		}
		//	print sep distance
	    writeDistances(0);
	} else {
		// print goodness
		cout << cs.num_score << endl;
	}
}
