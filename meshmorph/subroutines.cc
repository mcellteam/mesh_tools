
bool getPointEdgeDistance(double p[3],double P[3][3],bool p_flag){
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
		if (p_flag){
			cout << "getPointEdgeDistance: "
			<< "Edge ["
			<< Ax << ","
			<< Ay << ","
			<< Az << "]["
			<< Bx << ","
			<< By << ","
			<< Bz << "]"
			<< "point ["
			<< p[0] << ","
			<< p[1] << ","
			<< p[2] << "]"
			<< endl;
			cout << "getPointEdgeDistance: "
			<< "uDen " << uDen << endl;
		}
		if(uDen) {
			// u = AdotA-AdotB-AdotC+BdotC)/uDen
			u = (AdotA_minus_AdotB-(Ax*p[0]+Ay*p[1]+Az*p[2])
				+(Bx*p[0]+By*p[1]+Bz*p[2]))/uDen;
			if (print_flag){
				cout << "getPointEdgeDistance: "
				<< "u " << u << endl;
			}
			// no need to check for u ==0 and u ==1, since current 
			// vertex/face plane coincidence was checked previously.
			// Closest point on face edge line to current vertex
			// occurs on face edge between face vertices
			if (u>0 && u<1) {
				double a=Ax+u*(Bx-Ax);
				double b=Ay+u*(By-Ay);
				double c=Az+u*(Bz-Az);
				if (p_flag){
					cout << "getPointEdgeDistance: "
					<< "distance " << sqrt((p[0]-a)*(p[0]-a)+(p[1]-b)*(p[1]-b)+(p[2]-c)*(p[2]-c)) << endl;
					cout << "getPointEdgeDistance: " << "!distinguishable(p[0],a) " << !distinguishable(p[0],a) << endl;
					cout << "getPointEdgeDistance: " << "!distinguishable(p[1],b) " << !distinguishable(p[1],b) << endl;
					cout << "getPointEdgeDistance: " << "!distinguishable(p[2],c) " << !distinguishable(p[2],c) << endl;
				}
//				if(sqrt((p[0]-a)*(p[0]-a)+(p[1]-b)*(p[1]-b)+(p[2]-c)*(p[2]-c))<DOUBLE_EPSILON){ return true; }
				if(!distinguishable(p[0],a) && !distinguishable(p[1],b) && !distinguishable(p[2],c)){ return true; }
			}
		}
	}
	return false;
}

// #####################################################
// #####################################################

//typedef  unsigned long int  u4;   /* unsigned 4-byte type */
//typedef  unsigned     char  u1;   /* unsigned 1-byte type */

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

void copyControlFile (void) {
	cout << "\ncopy control settings...";
	cout.flush();
    char line[1024],filein[128],fileout[128];
    // open input data file
	sprintf(filein,"controls.cc");
    std::ifstream inFile(filein);
    if (inFile.fail()) // if stream cannot be opened
    { cout << "Can't open " << filein ; exit(1); }
	// open output data file
	sprintf(fileout,"%s%s",OUTPUT_DATA_DIR,CONTROL_FILE);
    std::ofstream outFile(fileout);
    if (outFile.fail()) // if stream cannot be opened
    { cout << "Can't open " << fileout ; exit(1); }
	// foreach line in input file
    while (inFile.getline(line,1024)) {
		outFile << line << endl;
	}
	// close files
    inFile.close();
    outFile.close();
	cout << "complete.\n";
	cout.flush();
}
// #####################################################
// #####################################################

time_t recordTime(std::ofstream& myfile, time_t before, char string[128])
{
	time_t after,diff;
	after = time (NULL);
	diff = after-before;
	myfile << string << " " << diff << " seconds\n";
	myfile.flush();
	return after;
}

// #####################################################
// #####################################################

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

// #####################################################
// #####################################################

void threeValueSort(double x[3], double &biggest, double &smallest) {
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
                smallest = x[0];
                biggest = x[2];
            } else {
                // x[0] < x[1]
                // x[2] <= x[1]
                // x[0] is smallest
                // x[1] is biggest or tied with x[2]
                smallest = x[0];
                biggest = x[1];
            }
        } else {
            // x[0] < x[1]
            // x[2] <= x[0]
            // x[2] is smallest or tied with x[0]
            // x[1] is biggest
            smallest = x[2];
            biggest = x[1];
        }
    } else {
        // x[1] <= x[0]
        if ( x[0] < x[2] ) {
            // x[1] <= x[0]
            // x[0] < x[2]
            // x1 is smallest or tied with x[0]
            // x2 is biggest
            smallest = x[1];
            biggest = x[2];
        } else {
            // x[1] <= x[0]
            // x[2] <= x[0]
            if (x[1] < x[2]) {
                // x[1] <= x[0]
                // x[2] <= x[0]
                // x[1] < x[2]
                // x[1] is smallest
                // x[0] is biggest
                smallest = x[1];
                biggest = x[0];
            } else {
                // x[1] <= x[0]
                // x[2] <= x[1]
                // x[2] is smallest or tied with x[1]
                // x[0] is biggest or tied with x[1]
                smallest = x[2];
                biggest = x[0];
            }
        }
	}
}

// #####################################################
// #####################################################

bool intersectionInBounds(double Ai,double Bi,double Aj,double Bj,double pI[2]){
	return (// Ai<=pI<=Bi or Bi<=pI<=Ai
			((Ai<=pI[0] && pI[0]<=Bi) || (Bi<=pI[0] && pI[0]<=Ai)) &&
			// Aj<=pI<=Bj or Bj<=pI<=Aj
			((Aj<=pI[1] && pI[1]<=Bj) || (Bj<=pI[1] && pI[1]<=Aj))
			);
}

void checkLineFaceIntersection(Face *f,double lp[2][3],bool &line_flag,
								bool &poly_flag, bool &poly_edge_flag,
								bool d_polygon_edge_intersection) {
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
				if (getPointEdgeDistance(I,pvc,print_flag)){poly_edge_flag = true;}
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
			checkLineFaceIntersection(of,lp,line_flag,poly_flag,poly_edge_flag,true);
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
			checkLineFaceIntersection(cf,lp,line_flag,poly_flag,poly_edge_flag,true);
		}
	}
	// do polygons intersect?
	if (line_flag && (poly_flag || poly_edge_flag)) {return(1);}
	else {return(0);}
}

// #####################################################
// #####################################################

double dot(double a[3],double b[3]){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

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

double Vertex::getSeparationForceEnergy(double force[3],bool flag){
	// compute separation vector
	double s[3];
	for(int i=0;i<3;i++){ s[i]=pC[i]-pN[i]; }
	// compute separation distance //////////
	double sd = sqrt(dot(s,s));
	// compute separation error (signed value)
	// NOTE THIS ASSUMES VERTEX IS INSIDE SAME OBJECT
	// ON WHICH CLOSEST POINT WAS FOUND, OTHERWISE
	// OTHERWISE SE!=sd+TARGET_SEPARATION
	// also note that se=sd+TARGET_SEPARATION is correct
	// for nonnice vertex, since pC-pN is oriented so as
	// to move vertex out of violated object
	double se;
	if(!o->vertexIsNice(this)){se=sd+TARGET_SEPARATION;}
	else{se=sd-TARGET_SEPARATION;}
	// if pC==pN, i.e. curretn vertex is on surface of neighbor object
	// then use current vertex outward normal as separation vector
	// recompute separation distance since it is used as normalization factor
	if (!sd) {
		getNormal(s);
		sd = sqrt(dot(s,s));
	}
	if(flag){
		// spring_force = spring constant * stretch
		// let scaled_spring force = spring_force/separation_distance
		double scaled_spring_force = (SEPARATION_WEIGHT/100.0)*se/sd;
		// force cartesian component = spring_force * unit vector
		// where unit vector = cartesian separation_distance component / separation_distance
		// force cartesian component = scaled_spring_force * cartesian separation_distance component
		force[0]+=scaled_spring_force*s[0];
		force[1]+=scaled_spring_force*s[1];
		force[2]+=scaled_spring_force*s[2];
	}
	// return energy
/*	double mofo = (SEPARATION_WEIGHT/100.0)/2.0*se*se;
//	if (mofo==140){
	if (!distinguishable(mofo,140,1E-8)){
		cout << "Vertex::getSeparationForceEnergy: "
		<< o->name << "->" << index
		<< " sd " << sd
		<< ", se " << se
		<< endl;
	}*/
	return (SEPARATION_WEIGHT/100.0)/2.0*se*se;
}

double Edge::getStretchForceEnergy(Vertex *v,double force[3],bool flag){
	// get pointer to adjacent vertex
	Vertex *av;
	if(v1==v){av=v2;}
	else {av=v1;}
	// compute separation vector
	double s[3];
	for(int i=0;i<3;i++){ s[i]=av->pN[i]-v->pN[i]; }
	// compute separation distance //////////
	double sd = sqrt(dot(s,s));
	// compute separation error (signed value)
	double se = sd-l;
	if(flag){
		// force contribution
		// spring_force = spring constant * stretch
		// scaled_spring force = spring_force/edge length
		double scaled_spring_force = (EDGE_STRETCH_WEIGHT/100.0)*se/sd;
		// force cartesian component = spring_force * unit vector
		// where unit vector = (adjacent vertex position - target vertex position) / edge length
		// force cartesian component = scaled_spring_force * cartesian component difference
		force[0]+=scaled_spring_force*s[0];
		force[1]+=scaled_spring_force*s[1];
		force[2]+=scaled_spring_force*s[2];
	}
	// energy contribution
	return (EDGE_STRETCH_WEIGHT/100.0)/2.0*se*se;
}

double Vertex::getEdgeStretchForceEnergy(double force[3],bool flag){
	double energy=0;
	// for each adjacent edge of current vertex
	for (std::vector<Edge*>::iterator j=e.begin();j!=e.end();j++) {
		// compute force and energy of edge
		energy+=(*j)->getStretchForceEnergy(this,force,flag);
	}
	return energy;
}

double Vertex::getEdgeAngleFaceIntersectionForceEnergy(double force[3],bool flag){
	double energy=0;
	// for each adjacent face of current vertex
	for (std::vector<Face*>::iterator k=f.begin();k!=f.end();k++) {
		// identify which face edge has current vertex as o1 or o2
		Edge *ae; // ae = pointer to associated edge
		if      ((*k)->e[0]->o1==this || (*k)->e[0]->o2==this){ae=(*k)->e[0];}
		else if ((*k)->e[1]->o1==this || (*k)->e[1]->o2==this){ae=(*k)->e[1];}
		else												  {ae=(*k)->e[2];}
		// force and contribution
		if (ae->o1==this){energy+=ae->getForceEnergy(0,force,flag);}
		else			 {energy+=ae->getForceEnergy(1,force,flag);}
		if(flag){
			// if adjacent face has intersection force
			if((*k)->faceInTable_iv()){
				// add intersection_force
				(*k)->getForceFromTable(force);
			}
		}
	}
	return energy;
}

void Vertex::getForceEnergy(double force[3]) {
	////////// get force and energy contributions from membrane separation //////////
	// true -> compute force
	if (cl!=NULL) { getSeparationForceEnergy(force,true); }
	////////// get force and energy contributions from adjacent edge stretching //////////
	getEdgeStretchForceEnergy(force,true);
	////////// get force contribution from edge angles and face intersections //////////
	getEdgeAngleFaceIntersectionForceEnergy(force,true);
}

// #####################################################
// #####################################################
void Vertex::computeNewCoords(Container *c,double pH[3],double gain) {
	// compute flip of associated edges
//	computeEdgeFlip();
	// get force
	double force[3]={0,0,0};
	getForceEnergy(force);
	// compute new vertex coordinates
//	for (int k=0;k<3;k++) { pH[k] = c->gain*force[k]+pN[k]; }
	for (int k=0;k<3;k++) { pH[k] = gain*force[k]+pN[k]; }
	// update container quantities
	c->force+=sqrt(force[0]*force[0]+force[1]*force[1]+force[2]*force[2]);
}
// #####################################################
// #####################################################
/*
double Edge::getAngle(void) {
	// get outward normals of edge faces
	double n1[3],n2[3];
	f1->getNormal(n1);
	f2->getNormal(n2);
	// compute the cosine of angle between normals
	double normal_angle_cosine=dot(n1,n2)/sqrt(dot(n1,n1))/sqrt(dot(n2,n2));
	// compute angle between normals 
    double normal_angle;
    if 		(normal_angle_cosine > 1)	{normal_angle = 0;}
    else if (normal_angle_cosine < -1)	{normal_angle = PI;}
    else 								{normal_angle = acos(normal_angle_cosine);}
    // use the edge itself as a reference vector
    double refvec[3] = {v1->pN[0]-v2->pN[0],v1->pN[1]-v2->pN[1],v1->pN[2]-v2->pN[2]};
    // compute the cross product of the normal vectors, i.e. n1 x n2
    double cross[3] = {
		n1[1]*n2[2]-n1[2]*n2[1],
		n1[2]*n2[0]-n1[0]*n2[2],
		n1[0]*n2[1]-n1[1]*n2[0]};
    // dot product of cross and refvec
    double d = dot(cross,refvec);
    // if cross and refvec point in same direction
    // i.e. if dotproduct of cross and refvec is positive
    // then edge angle = pi-normal_angle
    // else if cross and refvec point in opposite directions
    // i.e. if dotproduct of cross and refvec is negative
	// then edge angle = pi+normal_angle

    double angle;
	if (fabs(normal_angle-PI)<DOUBLE_EPSILON) {angle = 0;}
	else if (d>DOUBLE_EPSILON) {angle = PI-normal_angle;}
	else if (d<-DOUBLE_EPSILON) {angle = PI+normal_angle;}
	else {angle = PI;}
    return angle;
	// I EXPECT 0 <= angle < 2*PI
}*/

double Edge::getAngle(void) {
	// get outward normals of edge faces
	double n1[3],n2[3];
	f1->getNormal(n1);
	f2->getNormal(n2);
	// compute the cosine of angle between normals
	double normal_angle_cosine=dot(n1,n2)/sqrt(dot(n1,n1))/sqrt(dot(n2,n2));
	// compute angle between normals 
//    double normal_angle;
    if 		(normal_angle_cosine >= 1)	{
		//normal_angle == 0;
		// therefore gamma == 0
/*	if(v1->index==421 && v2->index==491 && 
		!strcmp(v1->o->name.c_str(),"a001_FILTERED_PRIMPED_CLIPPED_HEALED")){
		cout << "\nEdge::getAngle angle should be PI"<< endl;
	}*/
		return PI;
	}  else if (normal_angle_cosine <= -1)	{
		// normal_angle = PI;
		// gamma == 0 or 2PI
/*	if(v1->index==421 && v2->index==491 && 
		!strcmp(v1->o->name.c_str(),"a001_FILTERED_PRIMPED_CLIPPED_HEALED")){
		cout << "\nEdge::getAngle angle shoudl be 0" << endl;
	}*/
		return 0;
	} else {
		// normal_angle = acos(normal_angle_cosine);
	    // use the edge itself as a reference vector
	    double refvec[3] = {v2->pN[0]-o2->pN[0],v2->pN[1]-o2->pN[1],v2->pN[2]-o2->pN[2]};
	    // dot product of refvec and n1
	    double d = dot(refvec,n1);
//		double gamma = -d/abs(d)*acos(normal_angle_cosine);
//	    return PI-gamma;
//	if(v1->index==421 && v2->index==491 && 
		if(!d){
			return PI;
		} else {
//		!strcmp(v1->o->name.c_str(),"a001_FILTERED_PRIMPED_CLIPPED_HEALED")){
//		cout << "\nEdge::getAngle normal_angle_cosine " << normal_angle_cosine << endl;
//		cout << "\nEdge::getAngle acos " << acos(normal_angle_cosine) << endl;
//		cout << "\nEdge::getAngle d " << d << endl;
//		cout << "\nEdge::getAngle d/fabs(d) " << d/fabs(d) << endl;
	    return PI+d/fabs(d)*acos(normal_angle_cosine);
	}
		// I EXPECT 0 <= angle < 2*PI
	}
}

int Edge::computeFlip(void) {
	double n1[3],n2[3],cross[3],vec[3],force_dir[3];
	// indicates direction of force relative to normal,(same,+1 or opposite,-1)
	// get polygon outward normals
	f1->getNormal(n1);
	f2->getNormal(n2);
	////////// compute flip //////////
	// normals cross product = n1 x n2
	cross[0] = n1[1]*n2[2]-n1[2]*n2[1];
	cross[1] = n1[2]*n2[0]-n1[0]*n2[2];
	cross[2] = n1[0]*n2[1]-n1[1]*n2[0];
	// vector: edge to other_vertex_indices[0]
	vec[0] = o1->pN[0]-v1->pN[0];
	vec[1] = o1->pN[1]-v1->pN[1];
	vec[2] = o1->pN[2]-v1->pN[2];
	// force direction = normals cross product x vec
	force_dir[0] = cross[1]*vec[2]-cross[2]*vec[1];
	force_dir[1] = cross[2]*vec[0]-cross[0]*vec[2];
	force_dir[2] = cross[0]*vec[1]-cross[1]*vec[0];
	// if dot product of force direction and normal < 0 , then flip = -1
	if ( (force_dir[0]*n1[0]+force_dir[1]*n1[1]+force_dir[2]*n1[2]) < 0){return -1;}
	else {return 1;}
}

double Edge::getForceEnergy(int i,double force[3],bool flag) {
	// get polygon outward normals
	double n[3];
	if (i){f2->getNormal(n);}
	else  {f1->getNormal(n);}
	// compute normal length
	double L = sqrt(dot(n,n));
	// compute cosine of angle between normals
	double angle=getAngle();
	// record if angle is smallest so far
//	c->checkAngle(angle); 
	// determine force direction
	// when angle is less than PI, then angle_error is negative
	// so force will be opposite direction of face normal
	// when angle is greater than PI, then angle_error is positive
	// so force will be in same direction as face normal
	double angle_error = angle-PI;
	if(flag){
		// force 
//		double force_magn = (ANGLE_STRETCH_WEIGHT/100.0)*angle_error/L*o->getEdgeFlip(this);
		double force_magn = (ANGLE_STRETCH_WEIGHT/100.0)*angle_error/L;
		force[0]+=force_magn*n[0];
		force[1]+=force_magn*n[1];
		force[2]+=force_magn*n[2];
	}
	// energy
	return (ANGLE_STRETCH_WEIGHT/100.0)/2.0*angle_error*angle_error;
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
	else { cout << "Error. Tried to add third face to edge.\n"
				<< "Face " << f->index 
				<< " " << f->v[0]->index
				<< " " << f->v[1]->index
				<< " " << f->v[2]->index
				<< endl;
			exit(1); 
	}
	if(o1==NULL) {o1=vc;}
	else if (o2==NULL) {o2=vc;}
	else { cout << "Error. Tried to add third other vertex to edge.\n"
				<< "Vertex " << vc->index 
				<< " " << (int)vc->pN[0]
				<< " " << (int)vc->pN[1]
				<< " " << (int)vc->pN[2]
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
    // clean up hashmap
//    ht_iterator j;
//    for(j=hm.begin();j!=hm.end();j++){
//    delete (*j).second;
//    }
}

void checkEdges(Container &c){
	std::vector<Object*>::iterator i;
	std::vector<Edge*>::iterator j;
	// for each object
	for(i=c.o.begin();i!=c.o.end();i++){
		// for each edge in object
		for(j=(*i)->e.begin();j!=(*i)->e.end();j++){
			if ((*j)->o==NULL){
				cout << "checkEdges: Object pointer is NULL!\n";
				std::string s= "NULL OBJECT";
				(*j)->printEdge(s);
				exit(0);
			}
		}
	}
}

int Object::getMaxVertex(void){
	std::vector<Vertex*>::iterator i;
	// initialize max vertex index storage
	int a=v[0]->index;
	// for each vertex
	for (i=v.begin();i!=v.end();i++) {
		if ((*i)->index>a){a=(*i)->index;}
	}
	return a;
}

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
//		cout << sqrt((*i)->getSqLength()) << endl;
		L+=sqrt((*i)->getSqLength());
	}
	// compute mean edge length
	return (L/(double)e.size());
}

//void printNeighborhood(std::string str,int i,std::vector<Vertex*> &n,int iter){
void printNeighborhood(std::string str,int i,std::vector<Face*> &nf,int iter){
//	std::vector<Vertex*>::iterator j;
//	// for each vertex in neighborhood
//	for(j=n.begin();j!=n.end();j++){
	// for each face in neighborhood
	cout << "\n\nprintNeighborhood: iter " << iter
	<< ", current vertex " << str << "->" << i << endl;
	for(std::vector<Face*>::iterator j=nf.begin();j!=nf.end();j++){
		(*j)->printFace((*j)->v[0]->o->name);
		cout << endl;
	}
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

/**/

//neighbor::neighbor(Vertex *a, double b){
//	v=a;
//	l=b;
//}
/*
bool Object::processEdge(Edge *ee,hashtable_v_double &hood,std::vector<Edge*> &bucket){
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
		bool u1=false,u2=false;
		double L = sqrt(ee->getSqLength());
		// if v2 + length < v1, then update v1
		if((hood[vp2]+L)<hood[vp1]){u1=true;}
		// if v1 + length < v2, then update v2
		if((hood[vp1]+L)<hood[vp2]){u2=true;}
		if(u1&&u2){cout << "Object::newFindNeighborhoods: Error. Something is weird.\n";exit(0);}
		if(u1){hood[vp1]=hood[vp2]+L;}
		if(u2){hood[vp2]=hood[vp1]+L;}
	} else {
		if(f1){
			// if v1 found and v2 not found in hashtable
			// then add v2 to hashtable
			hood[vp2]=hood[vp1]+sqrt(ee->getSqLength());
		}
		if(f2){
			// if v2 found and v1 not found in hashtable
			// then add v1 to hashtable
			hood[vp1]=hood[vp2]+sqrt(ee->getSqLength());
		}
	}
	return empty;
}
*/
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
/*		bool u1=false,u2=false;
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
//			hood[vp1]=hood[vp2]+sqrt(ee->getSqLength());
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
//		cout << "Object::ifFrozen "
//		<< "vertex cumul adge length " << neighborhood[vv]
//		<< ", NEIGHBORHOOD_RADIUS " << NEIGHBORHOOD_RADIUS << endl;
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
//		bool mo = !ifFrozen(hood,(*i).first);
//		bool fo = disabled.find((*i).first)==disabled.end();
//		cout << "\nObject::thawedAndAble: "
//		<< "thawed " << mo
//		<< ", able " << fo
//		<< endl;
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
//		c.assign((*i)->f.begin(),(*i)->f.end());

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
//			c.push_back(*k);
		}

		///// all subsequent rounds /////

		// while there are thawed vertices in neighbor list, hood,
		// that are also not disabled
		while(thawedAndAble(hood,disabled)){
			// init collection of new faces
			std::vector<Face*> new_faces_too;
			collectFaces(hood,disabled,new_faces_too);

/*			// print new face list, new_faces
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
//					if(!processEdge((*k)->e[j],hood,bucket)){
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
//						if(!processEdge(*j,hood,pail)){
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

/*		// print hood
		if(!strcmp(name.c_str(),"d000_FILTERED_SMOOTH") && (*i)->index==18){
			cout << "\n\n****** current vertex *****" << endl;
			(*i)->printVertex((*i)->o->name);
			cout << endl;
			cout << "hash_map size " << hood.size() << endl;
			for(vdhm_iterator j=hood.begin();j!=hood.end();j++){
//			cout << "neighbor vertex " <<  (*j).first->index
//			<< ": length " << (*j).second << endl;
//			(*j).first->printVertex((*j).first->o->name);
				(*j).first->printVertexCP();
				cout << endl;
			}
		}
*/
	}
}
/**/

void Object::findNeighborhoods(void){
	vector<Vertex*> newverts,adj,bank;
	std::vector<Vertex*>::iterator i,k,m;
	Bit b;
	hashset_f fset;
	////////// compute neighborhood radius in edge lengths //////////
//	double nr = 2.0*(PI-asin(TARGET_SEPARATION/NEIGHBORHOOD_RADIUS))*NEIGHBORHOOD_RADIUS;
//	int radius = (int) ceil(nr/getMeanEdgeLength());
	int radius = (int) ceil(NEIGHBORHOOD_RADIUS/getMeanEdgeLength());
//	cout << radius << " " << name << endl; 
	////////// collect neighborhood vertices //////////
	// for each vertex in object
	for (i=v.begin();i!=v.end();i++) {
		// initialize
		b.init(v.size());
		newverts.clear();
		bank.clear();
		fset.clear();
		// add self to neighborhood
		bank.push_back(*i);
		b.addToAdjacent(bank);
		////// build neighborhood //////
		for (int j=0;j<radius;j++) {
			b.getNewVertices(newverts,v);
			// for each new vertex
			bank.clear();
			for (k=newverts.begin();k!=newverts.end();k++) {
				// get adjacent vertices
				(*k)->getAdjacentVertices(adj);
				// if not last round of vertex searching
				if (j<(radius-1)){
					// get adjacent faces of current vertex
					(*k)->getAdjacentFaces(fset);
				}
				// add all adjacent vertices of new vertex to bank
				bank.insert(bank.end(),adj.begin(),adj.end());
			}
			// add adjacent vertices of new vertex to adjacent bit map
			b.addToAdjacent(bank);
		}
		// copy fset to neighborhood faces vector
		(*i)->nf.assign(fset.begin(),fset.end());
		// trim excess capacity from vector
		std::vector<Face*>((*i)->nf.begin(),(*i)->nf.end()).swap((*i)->nf);
		// sort vectors
		sort((*i)->nf.begin(),(*i)->nf.end());
		b.clear();
	}
}

void Object::findVertexAdjacencies(void){
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
	// find neighborhoods
//	findNeighborhoods();
	newFindNeighborhoods();
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

void Container::writeObjectList(void) {
	char file[128];
	std::vector<Object*>::iterator i;
	// for each object, accumulate number of vertices and faces
	for (i=o.begin();i!=o.end();i++) {
		object_count++;
		vertex_count += (*i)->v.size();
		face_count += (*i)->f.size();
		edge_count += (*i)->e.size();
	}
	// open file
	sprintf(file,"%s%s",OUTPUT_DATA_DIR,OBJECT_LIST_FILE);
	Olist.open(file);
	// add stuff
	Olist << "Input data directory = " << INPUT_DATA_DIR << "\n"
			<< "Total number of input files = " << num_files << "\n"
			<< "Total number of (objects, vertices, faces, edges) = (" 
			<< object_count << ","
			<< vertex_count << ","
			<< face_count << ","
			<< edge_count << ")\n\n";
	//
    Olist.width(15);
	Olist << left << "Object name";	
    Olist.width(15);
	Olist << left << "#vertices";
    Olist.width(15);
	Olist << left << "#faces";
    Olist.width(15);
	Olist << left << "#edges" << endl;

	// for each object, write name and index
	for (i=o.begin();i!=o.end();i++) {
	    Olist.width(15);
    	Olist << left << (*i)->name;
	    Olist.width(15);
		Olist << left << (*i)->v.size();
	    Olist.width(15);
		Olist << left << (*i)->f.size();
	    Olist.width(15);
		Olist << left << (*i)->e.size() << endl;
	}
	Olist.close();
}

// #####################################################
// #####################################################

void Container::statusFileInit(void) {
	char file[128];
	sprintf(file,"%s%s",OUTPUT_DATA_DIR,CONT_LOG_FILE);
	Cfile.open(file);
	Cfile << "\nLEGEND_______________________\n"
	<< "'Iteration'\n"
	<< "0th iteration is the original data.\n"
	<< "For 1st iteration vertices have moved once or not at all.\n"
	<< "--\n"
	<< "'Nonnice'\n"
	<< "Nonnice vertices are inside of another object, possibly the parent object.\n"
	<< "--\n"
	<< "'Self-Nonnice'\n"
	<< "Self-nonnice vertices are inside of their parent object.\n"
	<< "--\n"
	<< "'Int. Faces'\n"
	<< "Total number of pairs of intersecting faces in all objects.\n"
	<< "--\n"
	<< "'Self-Int. Faces'\n"
	<< "Number of pairs of self-intersecting faces, i.e. faces share a common parent object.\n"
	<< "--\n"
	<< "'Force'\n"
	<< "Cumulative force on all vertices in all objects.\n"
	<< "--\n"
	<< "'Energy'\n"
	<< "Cumulative potential energy in all vertices in all objects.\n"
	<< "--\n"
	<< "'Mean Disp.'\n"
	<< "Mean displacement of all vertices in all objects chosen to move during the iteration.\n"
	<< "--\n"
	<< "'Min. Edge Angle (rad.)'\n"
	<< "The minimum edge angle of all edges in all objects.\n"
	<< "Minimum edge angle is the angle between the two faces that share the edge.\n"
	<< "___________________________________________\n\n\n\n";

/*    Cfile.width(12);
    Cfile << left << "Iteration";
    Cfile.width(10);
    Cfile << left << "Nonnice";
    Cfile.width(10);
    Cfile << left << "S.Nonnice";
    Cfile.width(10);
    Cfile << left << "I.Faces";
    Cfile.width(10);
    Cfile << left << "S.Faces";*/

	Cfile << "I            S\n";
	Cfile << "t            .\n";
	Cfile << "e            n            S\n";
	Cfile << "r            o            .\n";
	Cfile << "a            n            F\n";
	Cfile << "t            n            a\n";
	Cfile << "i            i            c\n";
	Cfile << "o            c            e\n";
	Cfile << "n   Nonnice  e   I.Faces  s   ";
    Cfile.width(10);
    Cfile << left << "N";
    Cfile.width(16);
    Cfile << left << "Energy";
    Cfile.width(11);
    Cfile << left << "Disp_Mean";
    Cfile.width(31);
    Cfile << "Disp_Min:Disp_Max ";
    Cfile << left << "Min. Edge Angle (rad.)";
    Cfile.width(20);
    Cfile << left << "duration (sec)" << endl;
	Cfile.flush();
}

// #####################################################
// #####################################################

void Container::fileOutit(void) {
	Cfile.close();
}

// #####################################################
// #####################################################

void Container::updateFile(int iteration,bool flag,double tim) {
	char name[1024];
	cout << "Iteration " << iteration << ": ";
	cout << "update log files..................";
	cout.flush();
	Cfile.precision(10);
    Cfile.width(4);
    Cfile << left << iteration;
    Cfile.width(9);
    Cfile << left << nonnice;
    Cfile.width(4);
    Cfile << left << s_nonnice;
    Cfile.width(9);
    Cfile << left << ti/2;
    Cfile.width(4);
    Cfile << left << si/2;
	Cfile.precision(8);
    Cfile.width(10);
    Cfile << left << N;
    Cfile.width(16);
    Cfile << left << energy;
    if (!flag){Cfile.width(42);Cfile << left << "NA";}
    else           {
        Cfile.width(11);
        Cfile << left << md[1];
//        Cfile << left << "(";
        Cfile.width(31);
//        Cfile << left << d_min;
//        Cfile << left << ":";
//        Cfile << left << d_max;
		sprintf(name,"(%.6g:%.6g)",d_min,d_max);
//		sprintf(name,"e %.6g, ",c.energy);
		Cfile << name;
//        Cfile << "(" << d_min << ":" << d_max << ")";
//        Cfile.width(3);
//        Cfile << left << ")";
    }
	Cfile.precision(10);
    Cfile << left << min_edge_angle;
    Cfile.width(20);
    Cfile << left << tim << endl;
	Cfile.flush();
	cout << "complete.\n";
	cout.flush();
}

// #####################################################
// #####################################################

void Container::createEdges(void) {
	if(!write_verbose_init){
		cout << "Creating edges.................................";
		cout.flush();
	}
	std::vector<Object*>::iterator i;
	// for each object, create edges
	for (i=o.begin();i!=o.end();i++) {
		if(write_verbose_init){
			cout << "Creating edges for " << (*i)->name << "...";
			cout.flush();
		}
		(*i)->createEdges();
		if(write_verbose_init){
			cout << "complete.\n";
			cout.flush();
		}
	}
	if(!write_verbose_init){
		cout << "complete.\n";
		cout.flush();
	}
}

void Container::findVertexAdjacencies(void) {
	if(!write_verbose_init){
		cout << "Finding vertex adjacencies.....................";
		cout.flush();
	}
	std::vector<Object*>::iterator i;
	// if 
	if (write_neighborhood_to_file){
		// open dat file
		char log_file[50];
		sprintf(log_file,"%s%s",OUTPUT_DATA_DIR,NEIGHBORHOOD_FILE);
		std::ofstream myfile (log_file);
		myfile.close();
	}
	// for each object, find vertex adjacencies
	if(write_verbose_init){ cout << endl;} 
	for (i=o.begin();i!=o.end();i++) {
		if(write_verbose_init){
			cout << "Finding adjacencies for " << (*i)->name << "...";
			cout.flush();
		}
		(*i)->findVertexAdjacencies();
		if(write_verbose_init){
			cout << "complete.\n";
			cout.flush();
		}
	}
	if(!write_verbose_init){
		cout << "complete.\n";
		cout.flush();
	}
}

// #####################################################
// #####################################################

void Vertex::computeEdgeFlip(void) {
	std::vector<Edge*> ee;
	// for each adjacent face of vertex*
	for (std::vector<Face*>::iterator i=f.begin();i!=f.end();i++) {
		// for each edge of face, if edge does not refer to this vertex, then add edge to vector
		if((*i)->e[0]->v1!=this && (*i)->e[0]->v2!=this){ee.push_back((*i)->e[0]);}
		if((*i)->e[1]->v1!=this && (*i)->e[1]->v2!=this){ee.push_back((*i)->e[1]);}
		if((*i)->e[2]->v1!=this && (*i)->e[2]->v2!=this){ee.push_back((*i)->e[2]);}
	}
	// clear edge hashtable in object
	o->clearFlipTable();
	// for each collected edge, compute flip and store in hashtable
	for (std::vector<Edge*>::iterator j=ee.begin();j!=ee.end();j++) {
		o->addEdgeFlip(*j,(*j)->computeFlip());
	}
}
// #####################################################
// #####################################################

void Container::scanDir(void) {
	num_files = 0;
    std::string str;
    std::string::size_type found;
	DIR *pdir;						// pointer to a directory data structure
	struct dirent *pent;			// pointer to dirent structure
//	cout << "\n\nInput Data Directory: " << INPUT_DATA_DIR << "\n"
//	<< "Reading Directory..............................";cout.flush();
	if(write_verbose_init){ cout << endl;cout.flush();}
    pdir = opendir(INPUT_DATA_DIR);
    if (!pdir) {printf("Error. Could not open %s.\n",INPUT_DATA_DIR);exit(1);}
    while ((pent=readdir(pdir))){
		// copy char array to string
		str = pent->d_name;
		// if file of typ *.mesh
		found = str.find(".mesh",0);
		// if found
        if (found != std::string::npos) {
			// save filename
			files.push_back(str);
			// update index
			num_files++;
			// print file found to screen
			if(write_verbose_init){
		        cout << "file found: " << str << "\n"; cout.flush();
			}
		}
    }
    closedir(pdir);
	if(write_verbose_init){ cout << endl;cout.flush();}
//	cout << "complete.\n";

}
// #####################################################
// #####################################################

void Container::scanFiles(void) {
    std::string str;
    std::string::size_type pos1;
	Object *obj;
	char file[1024];
	cout << "\n\nInput Data Directory: " << INPUT_DATA_DIR << "\n"
	<< "Reading Directory..............................";cout.flush();
	// for each input file
	for (int count=0;count<num_files;count++) {
		// copy char array to string
        str = files[count];
        // record object name
		pos1 = str.find(".",0);
		if (!(pos1==std::string::npos)) {
			// ALLOCATE MEMORY FOR NEW OBJECT
			obj = new Object(str.substr(0,pos1));
		} else { cout << "Error! Object name was not found in " << str << "\n";exit(1);}
		// save pointer to object
		o.push_back(obj);
		// scan file
		sprintf(file,"%s%s",INPUT_DATA_DIR,files[count].c_str());
		if(write_verbose_init){
			cout << "loading file " << files[count] << "...";
		}
		cout.flush();
		scanFile(obj,file);
		if(write_verbose_init){
			cout << "complete.\n";
			cout.flush();
		}
	}
//	cout << "\n";
	cout << "complete.\n";
}

// #####################################################
// #####################################################

void Container::scanFile (Object *obj,char *filename) {
    char line[2048],*str;
    FILE *F;
    Vertex *v;
    Face *f;
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
			vertex_num++;
			v=new Vertex(str,obj);
			obj->v.push_back(v);
			vp.push_back(v);
			if (print_flag) { cout.precision(15);cout<<filename<<": Vertex "<<v->index 
								<<" "<<v->pN[0]<<" "<<v->pN[1]<<" "<<v->pN[2]<<"\n";
			}
		}
        // if first character is F for Face, add new linked list class instance
        else if (strchr("F",*str)!=NULL){
			f=new Face(str,vp);
			polygon_num++;
			obj->f.push_back(f);
			if (print_flag) { cout.precision(15);cout<<filename<<": Face "<<f->index 
					<<" "<<f->v[0]->index<<" "<<f->v[1]->index<<" "<<f->v[2]->index<<"\n";
			}
		}
    }
    fclose(F);
}


// #####################################################
// #####################################################

void Container::buildMeshAfter(int iteration) {
	char file[128];
	std::vector<Object*>::iterator i;
	std::vector<Vertex*>::iterator j;
	std::vector<Face*>::iterator k;
	// for each object
	for(i=o.begin();i!=o.end();i++) {
		// create output filename
		if (append_iteration_number_to_mesh) {
			sprintf(file,"%s%s_%s_%i.mesh",OUTPUT_DATA_DIR,(*i)->name.c_str(),OUTPUT_SUFFIX,iteration);
		} else {
			sprintf(file,"%s%s_%s.mesh",OUTPUT_DATA_DIR,(*i)->name.c_str(),OUTPUT_SUFFIX);
		}
		// open output file
		std::ofstream newfile (file,std::ios::out);
		if(newfile.is_open()){
			newfile.precision(12);
			// for each vertex in object
			for(j=(*i)->v.begin();j!=(*i)->v.end();j++) {
				// print index and final coordinates
				newfile << "Vertex "
						<< (*j)->index << " "
						<< (*j)->pN[0] << " "
						<< (*j)->pN[1] << " "
						<< (*j)->pN[2] << "\n";
			}
			// for each face in object
			for(k=(*i)->f.begin();k!=(*i)->f.end();k++) {
				newfile << "Face "
						<< (*k)->index << " "
						<< (*k)->v[0]->index << " "
						<< (*k)->v[1]->index << " "
						<< (*k)->v[2]->index << "\n";
			}
		}
		newfile.close();
	}
}

void Container::writeDistancesNOCP(int iteration) {
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
				if((*j)->cl==NULL){
//				if((*j)->cl!=NULL &&
//					( (*j)->pN[0] > 3000 && (*j)->pN[0] < 7000 &&
//						(*j)->pN[1] > 3000 && (*j)->pN[1] < 6500 &&
//						(*j)->pN[2] > 3500 && (*j)->pN[2] < 7800)
//					){
					// compute separation vector
					// print separation distance
						newfile
						<< (*j)->pN[0] << " "
						<< (*j)->pN[1] << " "
						<< (*j)->pN[2] << " 1 0 0 1\n";
				}
			}
		}
		newfile.close();
	}
}

void Container::writeDistances50(int iteration) {
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
//				if((*j)->cl!=NULL &&
//					( (*j)->pN[0] > 3000 && (*j)->pN[0] < 7000 &&
//						(*j)->pN[1] > 3000 && (*j)->pN[1] < 6500 &&
//						(*j)->pN[2] > 3500 && (*j)->pN[2] < 7800)
//					){
					// compute separation vector
					double s[3];
					for(int k=0;k<3;k++){ s[k]=(*j)->pC[k]-(*j)->pN[k]; }
					// print separation distance
					double doo = sqrt(dot(s,s));
					if (doo>50.0){
						newfile
						<< (*j)->pN[0] << " "
						<< (*j)->pN[1] << " "
						<< (*j)->pN[2] << " 1 0 0 1\n";
					}
				}
			}
		}
		newfile.close();
	}
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
//				if((*j)->cl!=NULL &&
//					( (*j)->pN[0] > 3000 && (*j)->pN[0] < 7000 &&
//						(*j)->pN[1] > 3000 && (*j)->pN[1] < 6500 &&
//						(*j)->pN[2] > 3500 && (*j)->pN[2] < 7800)
//					){
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
void Space::initBoxes(void) {
	// subdivide space
	num_space[0] = (int) ceil( (world[1]-world[0])/SPACE_LENGTH );
	num_space[1] = (int) ceil( (world[3]-world[2])/SPACE_LENGTH );
	num_space[2] = (int) ceil( (world[5]-world[4])/SPACE_LENGTH );
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
void Space::fileInit(void) {
	char file[128];
	char str[15];
	sprintf(file,"%s%s",OUTPUT_DATA_DIR,SPACE_LOG_FILE);
	Sfile.open(file);
    Sfile.width(15);
	Sfile << left << "#boxes(x,y,z)";
    Sfile.width(15);
	Sfile << left << "total #boxes";	
    Sfile.width(30);
	Sfile << left << "world bounds (x,y,z)[min max]" << endl;
	sprintf(str,"%d,%d,%d",num_space[0],num_space[1],num_space[2]);
    Sfile.width(15);
	Sfile << left << str;
    Sfile.width(15);
	Sfile << left << num_boxes;
	Sfile << "[" << world[0] << " " << world[1] << " " << world[2] << " "
		<< world[3] << " " << world[4] << " " << world[5] << "]\n";
	Sfile.close();
}

// #####################################################
// #####################################################

void Face::getNormal(double n[3]) {
	double uX, uY, uZ, vX, vY, vZ;
	// compute vectors 01 and 12
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

double Vertex::getSqSepDist(void){
	return (pC[0]-pN[0])*(pC[0]-pN[0])
			+(pC[1]-pN[1])*(pC[1]-pN[1])
			+(pC[2]-pN[2])*(pC[2]-pN[2]);
}

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
//	if(
//		!distinguishable(lp[0][0],2719.54642515,1E-8) &&
//		!distinguishable(lp[0][1],4249.96388355,1E-8) &&
//		!distinguishable(lp[0][2],6242.92495541,1E-8)
//		){
/*		cout << "\nContainer::collectCrossed: "
		<< " lp [" << lp[0][0]
		<< " " << lp[0][1]
		<< " " << lp[0][2]
		<< " " << lp[1][0]
		<< " " << lp[1][1]
		<< " " << lp[1][2] << "]\n";*/
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
			nonnice++;
			pp=equal_range(cb.begin(),cb.end(),v->o);
			// if vertex is inside self object
			if(pp.first!=pp.second){
				s_nonnice++;
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
					(*jj)->printFace(v->o->name);
					cout << endl;
				}
			}
		}
	}else{ // else cb is empty, then vertex is nice
		// if vertex was nonnice, but not to self
		if (old_nice==1){	
			flag=true;
			if(nonnice>0){nonnice--;}
		}
		// if vertex was at least nonnice to self
		else if (old_nice==2){
			flag=true;
			if(nonnice>0){nonnice--;}
			if(s_nonnice>0){s_nonnice--;}
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
/*
void Container::findClosestAxis(Space &s,Vertex *v,double lp[2][3]) {
	// n=normal,r=ray
	// get normal info
	double n[3];
	v->getNormal(n);
	double L=sqrt( dot(n,n) );
	// identify nearest boundary
	// ray is always in x
	lp[0][0] = v->pN[0]+n[0]/L;			// start x
	lp[1][0] = lp[0][0]+2*s.world[1];	// end x
	lp[0][1] = v->pN[1]+n[1]/L;			// start y
	lp[1][1] = lp[0][1];				// end y
	lp[0][2] = v->pN[2]+n[2]/L;			// start z
	lp[1][2] = lp[0][2];				// end z
}*/

void Container::findClosestAxis(Space &s,Vertex *v,double lp[2][3]) {
	// n=normal,r=ray
	// get normal info
//	double n[3];
//	v->getNormal(n);
//	double L=sqrt( dot(n,n) );
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
//	lp[0][0] = v->pN[0]+n[0]/L;			// start x
//	lp[0][1] = v->pN[1]+n[1]/L;			// start y
//	lp[0][2] = v->pN[2]+n[2]/L;			// start z
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

/*	if(v->index==531 && !strcmp(v->o->name.c_str(),"a001_FILTERED_SMOOTH_SMOOTH")){
		cout << "\n\nContainer::findClosestAxis: "
		<< "i " << i
		<< ", min " << min
		<< "\nContainer::findClosestAxis: "
		<< "dis [" << dis[0]
		<< " " << dis[1]
		<< " " << dis[2]
		<< " " << dis[3]
		<< " " << dis[4]
		<< " " << dis[5] << "]"
		<< "\nContainer::findClosestAxis: "
		<< " lp [" << lp[0][0]
		<< " " << lp[0][1]
		<< " " << lp[0][2]
		<< " " << lp[1][0]
		<< " " << lp[1][1]
		<< " " << lp[1][2] << "]"
		<< "\nContainer::findClosestAxis: "
		<< "world [" << s.world[0]
		<< " " << s.world[1]
		<< " " << s.world[2]
		<< " " << s.world[3]
		<< " " << s.world[4]
		<< " " << s.world[5] << "]\n";
	}*/
}

// #####################################################
// #####################################################

inline bool comp(Face *a,Face* b) { return a->index == b->index; }

void Container::getBoxIndexRange(Space &s,double lp[2][3],int br[6]){
	// compute box index range that contains ray
	// note this range is zero lower-bounded (lowest range is zeroth box)
	// total range is 0..num_space[i]-1
/*	// +x
	br[0] = s.location2Index(lp[0][0],"x");
	br[1] = s.num_space[0]-1;
	br[2] = s.location2Index(lp[0][1],"y");
	br[3] = br[2];
	br[4] = s.location2Index(lp[0][2],"z");
	br[5] = br[4];*/

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
		checkLineFaceIntersection(*j,lp,line_flag,poly_flag,poly_edge_flag,false);
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
//							int& num_odd_objects,std::vector<std::string> &tmp) {
							int& num_odd_objects,std::vector<Object*> &tmp) {
	// find mesh objects crossed an odd number of times by ray
//	std::vector<std::string> ol; // object index list
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
//				if (!parity) {ol.push_back((*i)->o->name);}
				if (!parity) {ol.push_back((*i)->v[0]->o);}
			}
//		} else {ol.push_back((*i)->o->name);}
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
//			if (!strcmp(ol[k].c_str(),ol[k+1].c_str())) {k++;k++;}
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
//	double X[3],Y[3],on[3],*opvc[3],*cpvc[3];
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
	double cn[3],on[3];
	cf->getNormal(cn);
	of->getNormal(on);
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
			checkLineFaceIntersection(of,lp,line_flag,poly_flag,poly_edge_flag,true);
			// do faces intersect?
			if (line_flag && (poly_flag || poly_edge_flag)) {return(1);}
			lp[0][0] = opvc[p][0];
			lp[0][1] = opvc[p][1];
			lp[0][2] = opvc[p][2];
			lp[1][0] = opvc[q][0];
			lp[1][1] = opvc[q][1];
			lp[1][2] = opvc[q][2];
			checkLineFaceIntersection(cf,lp,line_flag,poly_flag,poly_edge_flag,true);
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
	cout << "assign faces to boxes..........................";
	cout.flush();
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
	cout << "complete.\n";
	cout.flush();
}

// #####################################################
// #####################################################

void checkFaces(Container &c) {
	std::vector<Object*>::iterator i;
	std::vector<Vertex*>::iterator j;
	// for each object in container
	for (i=c.o.begin();i!=c.o.end();i++) {
    	// for each face in object
		std::vector<Face*>::iterator w;
	    for (w=(*i)->f.begin();w!=(*i)->f.end();w++) {
			// if any two vertices have same location
			if (
				(
				!(distinguishable((*w)->v[0]->pN[0],(*w)->v[1]->pN[0])) &&
				!(distinguishable((*w)->v[0]->pN[1],(*w)->v[1]->pN[1])) &&
				!(distinguishable((*w)->v[0]->pN[2],(*w)->v[1]->pN[2])) 
				) ||
				(
				!(distinguishable((*w)->v[1]->pN[0],(*w)->v[2]->pN[0])) &&
				!(distinguishable((*w)->v[1]->pN[1],(*w)->v[2]->pN[1])) &&
				!(distinguishable((*w)->v[1]->pN[2],(*w)->v[2]->pN[2])) 
				) ||
				(
				!(distinguishable((*w)->v[2]->pN[0],(*w)->v[0]->pN[0])) &&
				!(distinguishable((*w)->v[2]->pN[1],(*w)->v[0]->pN[1])) &&
				!(distinguishable((*w)->v[2]->pN[2],(*w)->v[0]->pN[2])) 
				)
				){
			}
		}
	}
}


void Container::getSeparationDistances(Space &s,Monitor& stats){
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
			findClosest(s,*j,stats,false);
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

void Container::getNonnice(hashset_v &target_vset) {
    // for each object* in container
    for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    	// for each vertex* in object
	    for (std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
			// if vertex is nonnice, then add to set
			if(!(*j)->o->vertexIsNice(*j)){target_vset.insert(*j);}
		}
	}
}

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
	else {
		std::pair<std::vector<Face*>::iterator,std::vector<Face*>::iterator> q
			=equal_range(v->f.begin(),v->f.end(),f);
		// if face is in adjacent face vector
		// since total neighborhood = adjacent faces + non-adjacent faces
		if(q.first!=q.second){return true;}
		return false;
	}
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

void Monitor::validateVertex(int index,Container &c){
	///// check for discrepencies in separation errors /////
	// for each pair in topN
	for(tv_iterator i=topN.begin();i!=topN.end();i++){
		// if vertex index matches target
		if((*i).second->index==index){
			// if the se stored in topN does not match the computed se
			// then topN is stale, i.e. not up-to-date
			if( distinguishable((*i).first,fabs(c.computeSeparationError((*i).second,(*i).second->getSqSepDist())))) {
				cout << "\n\nMonitor::validateMultimap: "
				<< "Error! topN-stored se "
				<< (*i).first
				<< " does not match current value "
				<< fabs(c.computeSeparationError((*i).second,(*i).second->getSqSepDist()))
				<< ".\n";
				(*i).second->printVertex((*i).second->o->name);
				cout << endl << endl;
				exit(0);
			}
		}
	}
}

//void Container::validateMultimap(table_v &topN){
void Monitor::validateMultimap(void){
	///// check for duplicate vertex* /////
	// instantiate set of vertex*
	v_set uniq_set;
	// for each pair in topN
	for(tv_iterator i=topN.begin();i!=topN.end();i++){
		// load v_set
		uniq_set.insert((*i).second);
	}
	// if the two sets are not the same size
	// then topN likely contains duplicate vertex* entries
	if(uniq_set.size()!=topN.size()){
		cout << "\n\nMonitor::validateMultimap: "
		<< "Error! topN likely contains duplicate vertex* entries.\n";
		exit(0);
	}
	///// check for presence of vertices with no closest point /////
	// for each pair in topN
	for(tv_iterator i=topN.begin();i!=topN.end();i++){
		if((*i).second->cl==NULL){
			cout << "\n\nMonitor::validateMultimap: "
			<< "Error! topN contains vertex* with no closest point.\n";
			exit(0);
		}
	}
}

void Monitor::validateOld(void){
	tv_iterator j;
	// for each pair in old
	for(td_iterator i=old.begin();i!=old.end();i++){
		if((*i).second==0){cout << "Monitor::validateOld: Error. (*i).second==0.\n";exit(0);}
		if(!entryInTopN((*i).first,(*i).second,j)){
			cout << "\n\nMonitor::validateOld: "
			<< "Error! old contains entry not found in TopN.\n";
			exit(0);
		}
	}
}

bool Monitor::entryInTopN(Vertex *v,double vd_old,tv_iterator &j){
//	// topN has ascending order by key (double, squared virtual displacement)
	// topN has descending order by key (double, squared virtual displacement)
	std::pair<tv_iterator,tv_iterator> p;
	p=topN.equal_range(vd_old);
	// if not found, then return false
	if(p.first==p.second){return false;}
	// else was found
	else {
		// for each entry in range
		for(tv_iterator i=p.first;i!=p.second;i++){ 
			// if the data element on entry is equal to v
			if(v==(*i).second){
				// set j and return true
				j = i;
				return true;
			}
		}
		cout << "\n\nError. Matching multimap element not found.\n";
		cout << "old vd " << vd_old << endl;
		exit(0);
	}
}

void Monitor::validateTopN(char* str){
	for(tv_iterator i=topN.begin();i!=topN.end();i++){
		if((*i).second->cl==NULL){
			std::string s = str;
			cout << "ERROR at " << s
			<< ", vertex " << (*i).second->index << " in topN has no closest point.\n";
			exit(0);
		}
	}
}

/*
void Monitor::updateSets(Vertex *v,double new_sqD,bool flag){
	// if vertex* found in old remove first
	if(old.find(v)!=old.end()){ updateTopN(v,old[v],new_sqD,flag);}
	// else just add to topN (0.0 passed as dummy value,false flag prevents usage)
	else { updateTopN(v,0.0,new_sqD,false); }
	updateOld(v,new_sqD);
}*/

void Monitor::updateTopN(Vertex *v,double vd_old,double vd_new,bool flag){
	if (flag){
		tv_iterator t;
		// if find table entry with old se
		if (entryInTopN(v,vd_old,t)){
			// remove entry
			topN.erase(t);
//			cout << "\nfound\n";
		} else {
//			cout << "\nnot found\n";
		}
	}
	// add new se entry to table
	topN.insert(std::make_pair(vd_new,v));
}

double Container::computeSeparationError(Vertex* v,double old_se){
	// compute separation error (signed value)
	if(!v->o->vertexIsNice(v)){return sqrt(old_se)+TARGET_SEPARATION;}
	else{return sqrt(old_se)-TARGET_SEPARATION;}
}

bool Container::findClosest(Space &s,Vertex *v,Monitor& stats,bool flag) {
	bool gate = false;
	// declare pair
	dd_pair p;
	// get Box pointers for Face* collection
	std::vector<Box*> bp;
	getBoxes(bp,v,NUM_ADJACENT_BOXES,s);
	// collect Face pointers
	hashset_f cf;
	getCandidateFaces(bp,v,cf);
	// if candidate faces were found
	if (!cf.empty()){
		double squareD=0.0;
		// get vertex normal
		double n[3];
		v->getNormal(n);
		// for each candidate face
		for (hf_iterator j=cf.begin();j!=cf.end();j++) {
			// if the closest point to current vertex was found on this face
			if(computeClosest(*j,v,squareD,n)){ gate=true;}
		}
		// if closest point was found
		if(gate){
			if(flag){
				// update sets with vertex squared virtual displacement
				stats.updateSets(v,getVertexSqD(v),flag);
			} 
		} else { // vertex has no closest point
			// reset pointer to closest face
			v->cl=NULL;
			// remove vertex* from topN
			stats.topN.erase(stats.old[v]);
			// remove vertex* from old
			stats.old.erase(v);
			// set closest point to current vertex location
			for (int i=0;i<3;i++) {v->pC[i]=v->pN[i];}
		}
	}
	return gate;
}

// #####################################################
// #####################################################

/*
void Container::getPlaneIntersection(Face *f,Vertex *v,double *n,double num,double den,Point &p){
	//compute point on face plane that is closest to current vertex
	// i.e. intersection of plane with face normal through current vertex
    bool line_flag, poly_flag=false,poly_edge_flag;
	double lp[2][3],intersect[3],u=num/den;
	for (int i=0;i<3;i++) {intersect[i]=v->pN[i]+u*n[i];}
	lp[0][0] = v->pN[0];
	lp[0][1] = v->pN[1];
	lp[0][2] = v->pN[2];
	lp[1][0] = intersect[0];
	lp[1][1] = intersect[1];
	lp[1][2] = intersect[2];
	checkLineFaceIntersection(f,lp,line_flag,poly_flag,poly_edge_flag,false);
	// if intersection point is on face,then save point
//	if (poly_flag) {p.add(v->pN[0]+u*n[0],v->pN[1]+u*n[1],v->pN[2]+u*n[2]);}
	if (poly_flag) {p.add(intersect[0],intersect[1],intersect[2]);}
}*/

bool Container::getPlaneIntersection(Face *f,Vertex *v,double *n,double num,double den,Point &p){
	//compute point on face plane that is closest to current vertex
	// i.e. intersection of plane with face normal through current vertex
    bool line_flag, poly_flag=false,poly_edge_flag;
	double u=num/den;
	double lp[2][3] = {{v->pN[0],v->pN[1],v->pN[2]},
						{v->pN[0]+u*n[0],v->pN[1]+u*n[1],v->pN[2]+u*n[2]}};
	checkLineFaceIntersection(f,lp,line_flag,poly_flag,poly_edge_flag,true);
	// if intersection point is on face,then save point
	if (poly_flag) {p.add(lp[1][0],lp[1][1],lp[1][2]);}
	return (poly_flag || poly_edge_flag);
}

/*
void Container::getEdgeIntersection(Vertex *v,double *P[3],Point &p){
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
		if(uDen) {
			// u = AdotA-AdotB-AdotC+BdotC)/uDen
			u = (AdotA_minus_AdotB-(Ax*v->pN[0]+Ay*v->pN[1]+Az*v->pN[2])
				+(Bx*v->pN[0]+By*v->pN[1]+Bz*v->pN[2]))/uDen;
			// no need to check for u ==0 and u ==1, since current 
			// vertex/face plane coincidence was checked previously.
			// Closest point on face edge line to current vertex
			// occurs on face edge between face vertices
			if (u>0 && u<1) {
				p.add(Ax+u*(Bx-Ax),Ay+u*(By-Ay),Az+u*(Bz-Az));
			}
		}
	}
}*/

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

/*
bool Container::computeClosest(Face *f,Vertex *v,double &squareD,double vn[3]) {
	bool signal = false;
	// initialize point class instance with current vertex
	Point p(v->pN[0],v->pN[1],v->pN[2]);
	// get face vertex coordinates
	double *P[3] = {&(f->v[0])->pN[0],&(f->v[1])->pN[0],&(f->v[2])->pN[0]};
	// add each face vertex
	for (int i=0;i<3;i++) { p.add(P[i][0],P[i][1],P[i][2]); }
	// add points of minimum distance between current vertex and each face edge
	getEdgeIntersection(v,P,p);
	// compute vector connecting arbitrary face vertex and current vertex
	double diff[3] = {P[0][0]-v->pN[0],P[0][1]-v->pN[1],P[0][2]-v->pN[2]};
	// compute indicators
	double num=dot(vn,diff);
	// if current vertex does not lie on face plane
	if (!(num<DOUBLE_EPSILON)) {
		//compute point on face plane that is closest to current vertex
		// i.e. intersection of plane with face normal through current vertex
		// add to pp if intersection is on face
		getPlaneIntersection(f,v,vn,num,dot(vn,vn),p);
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
}*/

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

void Container::checkEdgeAngles(void) {
    // for each object* in container
    for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    	// for each Edge* in object
	    for (std::vector<Edge*>::iterator j=(*i)->e.begin();j!=(*i)->e.end();j++) {
			// record if angle is smallest so far
			checkAngle((*j)->getAngle()); 
		}
	}
}

void Vertex::assignHolding(double pH[3]){
	std::vector<Edge*>::iterator i;
	Vertex *vp;
	double d=0;
	///// if no adjacent vertex has same position as pH /////
	// for each adjacent edge
	for(i=e.begin();i!=e.end();i++){
		// vp=edge vertex different from this vertex
		if ((*i)->v1==this){vp=(*i)->v2;}
		else {vp=(*i)->v1;}
		// if non-self edge vertex is indistinguishable from pH, then displace pH
		if ( !distinguishable(vp->pN[0],pH[0])&&
			!distinguishable(vp->pN[1],pH[1])&&
			!distinguishable(vp->pN[2],pH[2])){ d=VERTEX_EPSILON;}
	}
	pN[0]=pH[0]+d;
	pN[1]=pH[1]+d;
	pN[2]=pH[2]+d;
}

bool Container::checkForIntersections(Vertex *v,Space &s,bool sw,hashtable_f &nb){
	// NOTE THIS FUNCTION IS NEVER CALLED WITH sw===false
	// INSTEAD getFaceIntersection IS CALLED FROM computeIntersectionForce
	std::vector<Face*> dummy;
	// for each adjacent face of current vertex
	for(std::vector<Face*>::iterator i=v->f.begin();i!=v->f.end();i++){
		//  if face intersects any other face, then return true
		if(sw)	{if((*i)->getFaceIntersectionCheck(this,s,nb)){return true;}}
		else	{if((*i)->getFaceIntersection(this,false,dummy)){return true;}}
	}
	// no adjacent faces of current vertex intersect any other face
	return false;
}

void Container::collectEdgeAngles(Vertex *v,Monitor& stats){
	stats.e_angle.clear();
	///// collect set of all unique edges of adjacent faces /////
	// for each adjacent face of current vertex
	for(std::vector<Face*>::iterator i=v->f.begin();i!=v->f.end();i++){
		stats.e_angle[(*i)->e[0]] = (*i)->e[0]->getAngle();
		stats.e_angle[(*i)->e[1]] = (*i)->e[1]->getAngle();
		stats.e_angle[(*i)->e[2]] = (*i)->e[2]->getAngle();
	}
}

bool Container::angleChangeIsWrong(double old_angle,double new_angle){
	if (old_angle < PI){
		// angle should increase towards PI
		// angle increase is correct
		// if angle increases, return false
		// if angle decreases, return true
		if ( new_angle>old_angle ) 	{return false;}
		else 						{return true;} 
	} else {
		// assume old_angle > PI, i.e. not exactly PI
		// angle should decrease towards PI
		// angle decrease is correct
		// if angle decreases, return false
		// if angle increases, return true
		if ( new_angle<old_angle ) 	{return false;}
		else 						{return true;} 
	}
}

bool Container::checkForSmallAngles(Monitor &stats){
	// for each element in hashtable (edge*->double)
	for(edhm_iterator j=stats.e_angle.begin();j!=stats.e_angle.end();j++){
		// if edge angle is less than threshold
		// and if angle change is in wrong direction, then return true
		double new_angle=(*j).first->getAngle();
		// small new_angles are acceptable if old_angle was also small and angle is improving
		// if new angle is small
		if (fabs(new_angle)<EDGE_ANGLE_THRESHOLD || fabs(2*PI-new_angle)<EDGE_ANGLE_THRESHOLD) {
			// if old angle is not small or angle change is wrong, then return true
			if(
				 !(fabs((*j).second)<EDGE_ANGLE_THRESHOLD || fabs(2*PI-(*j).second)<EDGE_ANGLE_THRESHOLD)
				|| angleChangeIsWrong((*j).second,new_angle) ){return true;}
		}
	}
	// no adjacent edges of current vertex violate edge angle threshold
	return false;
}

void Face::updateBoxes(hashtable_f &nb){
	// b = list of box* in which this face previously lay
	// nb = multimap of current vertex adjacent faces and the boxes in which they now lie
	//
	// sort old list
	sort(b.begin(),b.end());

	///// copy list of box*s in which this face lies to vector /////
	std::vector<Box*> nw;
	// grab list of box* in which this face now lies
	std::pair<tf_iterator,tf_iterator> np=nb.equal_range(this);
	// for each box* in new list (nb)
	for(tf_iterator i=np.first;i!=np.second;i++){
		nw.push_back((*i).second);
	}

	///// remove box* in common between old and new lists /////
	// for each box* in new list (nb)
	tf_iterator i=np.first;
	while(i!=np.second){
		// if new box* found in old list (b)
		std::pair<std::vector<Box*>::iterator,std::vector<Box*>::iterator> p;
		p=equal_range(b.begin(),b.end(),(*i).second);
		// if found, then remove box* from old and new list
		if(p.first!=p.second){
			b.erase(p.first);
			tf_iterator k = i;k++;
			nb.erase(i);
			i=k;
		}
		else {i++;}
	}

	///// remove face* from remaining box* in old list /////
	if (!b.empty()){
		// for each remaining Box* in b (the old list)
		for (std::vector<Box*>::iterator j=b.begin();j!=b.end();j++){
			// sort Face* vector
			sort((*j)->f.begin(),(*j)->f.end());
			// look for this face in Box* face list
			std::pair<std::vector<Face*>::iterator,std::vector<Face*>::iterator> q;
			q=equal_range((*j)->f.begin(),(*j)->f.end(),this);
			// if found, then remove
			if(q.first!=q.second){(*j)->f.erase(q.first);}
		}
	}

	///// add face* to remaining box* in new list /////
	np=nb.equal_range(this);
	// for each box* in new list (nb)
	for(i=np.first;i!=np.second;i++){
		// add this face to Box* face vector
		(*i).second->f.push_back(this);
	}

	// copy new list to old list
	b.assign(nw.begin(),nw.end());
}

void Vertex::getBoundingBox(double bb[6]){
	// assume bb=bounding box = [xmn xmax ymin ymax zmin zmax]
	// for each adjacent face
	for(std::vector<Face*>::iterator i=f.begin();i!=f.end();i++){
		// for each vertex of face
		for(int j=0;j<3;j++){
			Vertex *vv=(*i)->v[j];
			if      (vv->pN[0]<bb[0]){bb[0]=vv->pN[0];}
			else if (vv->pN[0]>bb[1]){bb[1]=vv->pN[0];}
			if      (vv->pN[1]<bb[2]){bb[2]=vv->pN[1];}
			else if (vv->pN[1]>bb[3]){bb[3]=vv->pN[1];}
			if      (vv->pN[2]<bb[4]){bb[4]=vv->pN[2];}
			else if (vv->pN[2]>bb[5]){bb[5]=vv->pN[2];}
		}
	}
}

void Container::getAffectedVerticesAndEdgesBefore(Space &s,Vertex *v,double temp[3],Monitor &stats){
	double bb[6] = {temp[0],temp[0],temp[1],temp[1],temp[2],temp[2]};
	// grab bounding box of set of all adjacent faces to vertex
	v->getBoundingBox(bb);
	bb[0] -= NUM_ADJACENT_BOXES*SPACE_LENGTH;
	bb[1] += NUM_ADJACENT_BOXES*SPACE_LENGTH;
	bb[2] -= NUM_ADJACENT_BOXES*SPACE_LENGTH;
	bb[3] += NUM_ADJACENT_BOXES*SPACE_LENGTH;
	bb[4] -= NUM_ADJACENT_BOXES*SPACE_LENGTH;
	bb[5] += NUM_ADJACENT_BOXES*SPACE_LENGTH;

	std::vector<Box*> bp;
	std::vector<Box*>::iterator i;
	std::vector<Face*>::iterator j;
	s.getBoxesFor3DLocations(bb,bp);
//	if (bp.empty()){cout << "bp is empty!\n";}
	stats.av.clear();
	// for each box
	for(i=bp.begin();i!=bp.end();i++){
		// for each face in box
		for(j=(*i)->f.begin();j!=(*i)->f.end();j++){
			// if face vertex has a closest point and 
			// closest point lies on an adjacent face of active vertex, v
			// then add vertex to affected vertex
			if( (*j)->v[0]->cl!=NULL && find(v->f.begin(),v->f.end(),(*j)->v[0]->cl)!=v->f.end() ){
				stats.av.insert((*j)->v[0]);
			}
			if( (*j)->v[1]->cl!=NULL && find(v->f.begin(),v->f.end(),(*j)->v[1]->cl)!=v->f.end() ){
				stats.av.insert((*j)->v[1]);
			}
			if( (*j)->v[2]->cl!=NULL && find(v->f.begin(),v->f.end(),(*j)->v[2]->cl)!=v->f.end() ){
				stats.av.insert((*j)->v[2]);
			}
		}
	}
}

void Container::getAffectedVerticesAndEdgesAfter(Space &s,Vertex *v,double temp[3],Monitor &stats){
	double bb[6] = {temp[0],temp[0],temp[1],temp[1],temp[2],temp[2]};
	v->getBoundingBox(bb);
	bb[0] -= NUM_ADJACENT_BOXES*SPACE_LENGTH;
	bb[1] += NUM_ADJACENT_BOXES*SPACE_LENGTH;
	bb[2] -= NUM_ADJACENT_BOXES*SPACE_LENGTH;
	bb[3] += NUM_ADJACENT_BOXES*SPACE_LENGTH;
	bb[4] -= NUM_ADJACENT_BOXES*SPACE_LENGTH;
	bb[5] += NUM_ADJACENT_BOXES*SPACE_LENGTH;

	std::vector<Box*> bp;
	s.getBoxesFor3DLocations(bb,bp);
	// for each box
	for(std::vector<Box*>::iterator i=bp.begin();i!=bp.end();i++){
		// for each face in box
		for(std::vector<Face*>::iterator j=(*i)->f.begin();j!=(*i)->f.end();j++){
			// store face vertices
			stats.av.insert((*j)->v[0]);
			stats.av.insert((*j)->v[1]);
			stats.av.insert((*j)->v[2]);
		}
	}

	// collect edges from current vertex adjacent faces
	stats.ae.clear();
	for(std::vector<Face*>::iterator i=v->f.begin();i!=v->f.end();i++){
		// store face edges
		stats.ae.insert((*i)->e[0]);
		stats.ae.insert((*i)->e[1]);
		stats.ae.insert((*i)->e[2]);
	}
}

void Monitor::initTable(void){
	nb.clear();
	// for each adjacent face
//	for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
		// create Box* vector
//        nb[*k]= new std::vector<Box*>();
//	}
}

void Monitor::clearTable(void){
	// for each adjacent face
//	for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
		// delete Box* vector
//        delete nb[*k];
//	}
	nb.clear();
}

void Monitor::saveOld(void){
	old.clear();
	// find vertex in topN
	for(tv_iterator i=topN.begin();i!=topN.end();i++){
		// add to old
		old[(*i).second]=(*i).first;
	}
}

void Monitor::updateOld(Vertex *v,double new_vd){
	old[v]=new_vd;
}

void Container::collectAdjacentFaceNormals(table_fd& adj_n,Vertex *v){
	for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
		// get vertex normal
		double *n = new double[3];
		v->getNormal(n);
		adj_n[*k]=n;
	}
}

void Container::freeAdjacentFaceNormals(table_fd& adj_n){
	for(fd_iterator j=adj_n.begin();j!=adj_n.end();j++){
		delete[] (*j).second;
	}
}

void Container::updateClosest(Space &s,Vertex *v,Monitor &stats){
	double dummy[3] = {0.0,0.0,0.0};
	// collect adjacent face normals
	table_fd adj_n;
	collectAdjacentFaceNormals(adj_n,v);
	hashset_v vset2; // all vertices within NUM_ADJACENT_BOXES of adjacent faces
					// that do not require full search for closest point
	///// sort affected vertices according to whether full closest point search is required /////
	///// if full search is required add to set_nice_changed which also need full search /////
	// for each collected vertex
	for(hv_iterator kk=stats.av.begin();kk!=stats.av.end();kk++){
		// if collected vertex is different than active vertex, v
		if (*kk!=v){
			// if collected vertex has a closest point and 
			// closest point lies on an adjacent face of active vertex, v
			if( (*kk)->cl!=NULL && find(v->f.begin(),v->f.end(),(*kk)->cl)!=v->f.end() ){
				// add vertex to set_nice_changed
				stats.set_nice_changed.insert(*kk);
			} else {
				vset2.insert(*kk);
			}
		}
	}
	// for each collected vertex requiring full closest point search
	for(hv_iterator kk=stats.set_nice_changed.begin();kk!=stats.set_nice_changed.end();kk++){
		// then update collected vertex's closest point
		// true => remove existing matching element, if any
		if(findClosest(s,*kk,stats,true)){
			// if closest point changed then update global energy
			// false -> do not compute force, hence dummy
			if (*kk!=v){
				energy=energy-stats.v_energy[*kk]+(*kk)->getSeparationForceEnergy(dummy,false);
			}
		} else {
			// vertex lost it's closest point
			if (*kk!=v){
				energy=energy-stats.v_energy[*kk]
						+(SEPARATION_WEIGHT/100.0)/2.0*SEARCH_RADIUS*SEARCH_RADIUS;
			}
		}
	}
	// for each collected vertex not requiring full search
	for(hv_iterator kk=vset2.begin();kk!=vset2.end();kk++){
		///// check all adjacent faces to see /////
		///// if affected vertex's closest point has changed /////
		bool gate=false;
		double squareD = 0.0;
		// compute current square of separation distance for collected vertex
		if((*kk)->cl!=NULL){ squareD = (*kk)->getSqSepDist(); }
		// for each adjacent face of vertex *v
		for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
			// if face vertices are not in current vertex neighborhood
			if (!faceInNeighborhood(*k,*kk)){
				// if the closest point to current vertex was found on this face
				if(computeClosest(*k,*kk,squareD,adj_n[*k])){ gate=true; }
				//if(computeClosest(*k,*kk,squareD,adj_n[*k],false)){ gate=true; }
			}
		}
		// if new closest point was found
		if(gate){
			// update global energy
			// false -> do not compute force, hence dummy
			energy=energy-stats.v_energy[*kk]+(*kk)->getSeparationForceEnergy(dummy,false);
			// update sets with vertex squared virtual displacement
			stats.updateSets(*kk,getVertexSqD(*kk),true);
		}
	}
	// free adjacent face normal memory
	freeAdjacentFaceNormals(adj_n);
}


void Container::removeOldIntersections(Vertex *v,hashset_v &vset2){
	// for every adjacent face
	for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
		// if adjacent face has intersecting faces
		if((*k)->faceInTable_intf()){
			// for each intersecting face of adjacent face
			std::vector<Face*> *fv=(*k)->getIntersectingFaces();
			for(std::vector<Face*>::iterator p=(*fv).begin();p!=(*fv).end();p++){
				// add intersecting face vertices to vset2
				vset2.insert((*p)->v[0]);
				vset2.insert((*p)->v[1]);
				vset2.insert((*p)->v[2]);
				// remove adjacent face from intersecting face's vector
				(*p)->removeFaceFromVector(*k);
				// update si for removal of both adjacent and intersecting faces
				if((*k)->v[0]->o==(*p)->v[0]->o){si-=2;}
				// update ti for removal of both adjacent and intersecting faces
				ti-=2;
			}
		}
		// clear adjacent face's intersecting face vector
		(*k)->clearFaceFromTable_intf();
		// update adjacent face intersect force, i.e. set force to zero
		(*k)->clearFaceFromTable_iv();
	}
}

void Container::updateEdgeAngles(Vertex *v){
	// for every adjacent face
	for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
		// search for new mininmum edge angle
		checkAngle( (*k)->e[0]->getAngle() ); 
		checkAngle( (*k)->e[1]->getAngle() ); 
		checkAngle( (*k)->e[2]->getAngle() ); 
	}
}

void Container::updateNewIntersections(Vertex *v,hashset_v &vset2){
	// for every adjacent face
	for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
		// if adjacent face is currently intersected
		if ((*k)->computeIntersectionForce(this)){
			// for each intersecting face of adjacent face
			std::vector<Face*> *fv=(*k)->getIntersectingFaces();
			for(std::vector<Face*>::iterator p=(*fv).begin();p!=(*fv).end();p++){
				// add intersecting face vertices to vset2
				vset2.insert((*p)->v[0]);
				vset2.insert((*p)->v[1]);
				vset2.insert((*p)->v[2]);
				// if intersected face was previously not intersected
				if(!(*p)->faceInTable_intf()){
					// create new face vector in table for intersected face
					(*p)->addFaceToTable_intf();
					// add adjacent face* to intersected face
					(*p)->addFaceToVector(*k);
					// update si
					if((*k)->v[0]->o==(*p)->v[0]->o){si++;}
					// update ti
					ti++;
					// update intersected face intersection force
					(*p)->calculateIntersectionForce(this);
				} else {
					// add adjacent face* to intersected face
					(*p)->addFaceToVector(*k);
					// update si
					if((*k)->v[0]->o==(*p)->v[0]->o){si++;}
					// update ti
					ti++;
				}
			}
		}
	}
}

void Container::updateAdjacentFaceBoxes(Vertex *v,Monitor &stats){
	// for each adjacent face, update Box*s
	for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
		(*k)->updateBoxes(stats.nb);
	}
}

void Container::getNiceCheckSet(Vertex *v,Monitor &stats){
	// store vertices for which niceness may have changed
	stats.set_nice_check.clear();
	stats.set_nice_check.insert(v);
	removeOldIntersections(v,stats.set_nice_check);
	updateNewIntersections(v,stats.set_nice_check);
}

void Container::getNiceSet(Space &s,Monitor &stats){
	stats.set_nice_changed.clear();
	// update niceness and save to set
	for(hv_iterator i=stats.set_nice_check.begin();i!=stats.set_nice_check.end();i++){
		if(checkNiceness(s,*i)){stats.set_nice_changed.insert(*i);}
	}
}

double Container::getVertexSqD(Vertex *v){
	// compute flip of associated edges
//	computeEdgeFlip(v);
	// compute new vertex coords
	double pH[3]; // new holding position coordinates (x,y,z)
	v->computeNewCoords(this,pH,gain);
	return (pH[0]-v->pN[0])*(pH[0]-v->pN[0])+
				(pH[1]-v->pN[1])*(pH[1]-v->pN[1])+
				(pH[2]-v->pN[2])*(pH[2]-v->pN[2]);
}

void Monitor::updateEnergyMap(Vertex *v,double new_energy){
	v_energy[v]=new_energy;
}

void Monitor::updateSets(Vertex *v,double new_sqD,bool flag){
	// if vertex* found in old remove first
	if(old.find(v)!=old.end()){ updateTopN(v,old[v],new_sqD,flag);}
	// else just add to topN
	else { updateTopN(v,0.0,new_sqD,false); }
	updateOld(v,new_sqD);
}

void Container::computeGlobalEnergy(void){
	energy=0;
	double dummy[3] = {0.0,0.0,0.0};
	// for each object in container
	for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++){
		// for each vertex in object
		for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++){
			// false -> do not compute force, hence dummy
			energy+=(*j)->getSeparationForceEnergy(dummy,false);
/*			if (energy != energy){
				cout << "\nthis vertex started the problem\n";
				(*j)->printVertex((*i)->name);
				cout << endl; exit(0);
			}*/
		}
		// for each edge in object
		for(std::vector<Edge*>::iterator j=(*i)->e.begin();j!=(*i)->e.end();j++){
			// choice of v1 is arbitrary, v2 could have been used
			// false -> do not compute force, hence dummy

/*			if(
				(*j)->v1->index==421 && 
				(*j)->v2->index==491 && 
				!strcmp((*i)->name.c_str(),"a001_FILTERED_PRIMPED_CLIPPED_HEALED")
				){
				cout << "\ngetStretchForceEnergy " 
				<< (*j)->getStretchForceEnergy((*j)->v1,dummy,false) << endl;
				cout << "\ngetForceEnergy " 
				<< (*j)->getForceEnergy(0,dummy,false) << endl;
			}*/

			energy+=(*j)->getStretchForceEnergy((*j)->v1,dummy,false)
								+ (*j)->getForceEnergy(0,dummy,false);
/*			if (energy != energy){
				cout << "\nthis edge started the problem\n";
				(*j)->printEdge((*i)->name);
				cout << endl; exit(0);
			}*/
		}
	}
}

void Container::getVertexAndEdgeEnergy(Monitor &stats){
	stats.v_energy.clear();
	stats.e_energy.clear();
	double dummy[3] = {0.0,0.0,0.0};
	// for each affected vertex
	for(hv_iterator i=stats.av.begin();i!=stats.av.end();i++){
		// load hashtable with vertex energy
		// false -> do not compute force, hence dummy
		if((*i)->cl!=NULL){
			stats.v_energy[*i]=(*i)->getSeparationForceEnergy(dummy,false);
		} else {
			stats.v_energy[*i]=(SEPARATION_WEIGHT/100.0)/2.0*SEARCH_RADIUS*SEARCH_RADIUS;
		}
	}
	// for each affected edge
	for(he_iterator i=stats.ae.begin();i!=stats.ae.end();i++){
		// load hashtable with edge energy
		// choice of v1 is arbitrary, v2 could have been used
		// false -> do not compute force, hence dummy
		stats.e_energy[*i]=(*i)->getStretchForceEnergy((*i)->v1,dummy,false)
							+ (*i)->getForceEnergy(0,dummy,false);
	}
}

void Container::updateMovedVertexEnergy(Vertex *v,Monitor &stats){
	double dummy[3] = {0.0,0.0,0.0};
	///// vertex contribution /////
	// false -> do not compute force, hence dummy
	energy=energy-stats.v_energy[v]+(v)->getSeparationForceEnergy(dummy,false);
	///// adjacent edge contribution /////
	// for each adjacent edge of current vertex
	for (std::vector<Edge*>::iterator i=v->e.begin();i!=v->e.end();i++) {
		// compute force and energy of edge
		energy=energy-stats.e_energy[*i]+(*i)->getStretchForceEnergy(v,dummy,false);
	}
	///// opposing edge contribution /////
	// for each adjacent face of current vertex
	for (std::vector<Face*>::iterator i=v->f.begin();i!=v->f.end();i++) {
		// identify which face edge has current vertex as o1 or o2
		Edge *ae; // ae = pointer to associated edge
		if      ((*i)->e[0]->o1==v || (*i)->e[0]->o2==v){ae=(*i)->e[0];}
		else if ((*i)->e[1]->o1==v || (*i)->e[1]->o2==v){ae=(*i)->e[1];}
		else											{ae=(*i)->e[2];}
		// energy contribution
		int tag=1;
		if (ae->o1==v){tag=0;}
		energy=energy-stats.e_energy[ae]+ae->getForceEnergy(tag,dummy,false);
	}
}

void Container::updateVertexVD(Vertex *v, Monitor &stats){
	// set for storing nearby vertices to update
	v_set nearby;
	// for each collected edge from current vertex adjacent faces
	for(he_iterator i=stats.ae.begin();i!=stats.ae.end();i++){
		// if edge is not adjacent to current vertex
//		if(v->e.find(*i)==v->e.end()){
		if(
			find(v->e.begin(),v->e.end(),*i)==v->e.end()
			){
			// insert all four edge vertices (v1,v2,o1,o2) into set
			 nearby.insert((*i)->v1);
			 nearby.insert((*i)->v2);
			 nearby.insert((*i)->o1);
			 nearby.insert((*i)->o2);
		}
	}
	// update collected vertices including current vertex
	// for each collected vertex
	for(vs_iterator i=nearby.begin();i!=nearby.end();i++){
		// if vertex has a closest point, then update sets
		if((*i)->cl!=NULL){
			stats.updateSets(*i,getVertexSqD(*i),true);
		}
	}
}

bool Container::assignNewVertexCoords(Space &s,Vertex *v,double pH[3],Monitor &stats) {
	// build face* hashtable
	stats.nb.clear();
	// store current vertex position
	double pO[3]={v->pN[0],v->pN[1],v->pN[2]};
	// collect vertices whose closest point is on adjacent face
	// within NUM_ADJACENT_BOXES of adjacent faces
	getAffectedVerticesAndEdgesBefore(s,v,pO,stats);
	// collect edge angles
	collectEdgeAngles(v,stats);
	// set current position to holding position
	v->assignHolding(pH);
	// collect all vertices and edges within NUM_ADJACENT_BOXES of adjacent faces
	getAffectedVerticesAndEdgesAfter(s,v,pO,stats);
	// load vertex and edge energy hashmaps
	getVertexAndEdgeEnergy(stats);
	// if no faces intersect and no edges have small angles
	bool int_flag = (!checkForIntersections(v,s,true,stats.nb) || INTERSECTION_WEIGHT==100.0);
	bool angle_flag = !checkForSmallAngles(stats);
	if (int_flag && angle_flag){
//	if ( (!checkForIntersections(v,s,true,stats.nb) || INTERSECTION_WEIGHT==100.0) &&
//		!checkForSmallAngles(stats)){
		// for each adjacent face, update Box*s
		updateAdjacentFaceBoxes(v,stats);
		// for each adjacent face, update edge angles
		updateEdgeAngles(v);
		// store vertices for which niceness may have changed
		getNiceCheckSet(v,stats);
		// detect niceness changes and collect changed vertices
		getNiceSet(s,stats);
		// update closest point for affected vertices around adjacent faces
		updateClosest(s,v,stats);
		// update sets with squared virtual displacement of nearby vertics
		updateVertexVD(v,stats);
		// update sets with this vertex squared virtual displacement
//		stats.updateSets(v,getVertexSqD(v),true);
		// update global energy due to this vertex
		updateMovedVertexEnergy(v,stats);
		// clear face* hashtable
		stats.clearTable();
		return true;
	} else { // move vertex back
//		if(!int_flag) {cout << "Vertex not moved: intersections detected. ";}
//		if(!angle_flag) {cout << "Vertex not moved: small angles detected. ";}
		v->pN[0]=pO[0];
		v->pN[1]=pO[1];
		v->pN[2]=pO[2];
		stats.clearTable();
		return false;
	}
}

double getMinDistance(Vertex *p,Vertex *e1,Vertex *e2){
	double Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz;
	double AdotA,AdotB,AdotC,BdotC,BdotB;
	double uDen,uNum,u,Ix,Iy,Iz,min;
	Ax = e1->pN[0];
	Ay = e1->pN[1];
	Az = e1->pN[2];
	Bx = e2->pN[0];
	By = e2->pN[1];
	Bz = e2->pN[2];
	Cx = p->pN[0];
	Cy = p->pN[1];
	Cz = p->pN[2];
	AdotA = Ax*Ax+Ay*Ay+Az*Az;
	AdotB = Ax*Bx+Ay*By+Az*Bz;
	AdotC = Ax*Cx+Ay*Cy+Az*Cz;
	BdotC = Bx*Cx+By*Cy+Bz*Cz;
	BdotB = Bx*Bx+By*By+Bz*Bz;
	uDen = AdotA+BdotB-2*AdotB;
	if(uDen) {
		uNum = AdotA-AdotB-AdotC+BdotC;
		u = uNum/uDen;
		// no need to check for u ==0 and u ==1, since current 
		// vertex/face plane coincidence was checked previously.
		// Closest point on face edge line to current vertex
		// occurs on face edge between face vertices
		Ix=Ax+u*(Bx-Ax);
		Iy=Ay+u*(By-Ay);
		Iz=Az+u*(Bz-Az);
	} else {cout << "What the...\n";exit(0);}
	min = sqrt((Ix-Cx)*(Ix-Cx)+(Iy-Cy)*(Iy-Cy)+(Iz-Cz)*(Iz-Cz));
	return min;
}

void checkAspectRatio(Container &c) {
	std::vector<Object*>::iterator i;
	std::vector<Face*>::iterator j;
	std::vector<double> L,S;
	double ar;
    // for each object* in container
    for (i=c.o.begin();i!=c.o.end();i++) {
    	// for each Face* in object
	    for (j=(*i)->f.begin();j!=(*i)->f.end();j++) {
			L.clear();
			L.push_back(sqrt((*j)->e[0]->getSqLength()));
			L.push_back(sqrt((*j)->e[1]->getSqLength()));
			L.push_back(sqrt((*j)->e[2]->getSqLength()));
			sort(L.begin(),L.end());
			// L[2] is longest edge
			S.clear();
			S.push_back(getMinDistance((*j)->v[0],(*j)->v[1],(*j)->v[2]));
			S.push_back(getMinDistance((*j)->v[1],(*j)->v[0],(*j)->v[2]));
			S.push_back(getMinDistance((*j)->v[2],(*j)->v[0],(*j)->v[1]));
			sort(S.begin(),S.end());
			// S[0] is shortest triangle altitude
			ar = L[2]/S[0];
			if (ar>10000) {cout << "aspect ratio " << ar << endl;(*j)->printFace((*j)->v[0]->o->name);cout << endl << endl;}
		}
	}
}

// #####################################################
// #####################################################

bool facesInSameNeighborhood(Container *c,Face *cf, Face *of){
	///// if any vertex of other face is in /////
	///// neighborhood of any current face vertex /////
	// for each vertex in current face
	for (int i=0;i<3;i++){
		// if face vertices are in current vertex neighborhood
		if (c->faceInNeighborhood(of,cf->v[i])){return true;}
	}
	return false;
}

void Face::calculateIntersectionForce(Container *c) {
	// clear force
	this->clearFaceFromTable_iv();
	// if this face has intersecting faces
	if(this->faceInTable_intf()){
		// add this face to force table
		this->addFaceToTable_iv();
		// for each intersecting face
		std::vector<Face*> *fv=this->getIntersectingFaces();
		for(std::vector<Face*>::iterator i=(*fv).begin();i!=(*fv).end();i++){
			double n1[3],n2[3],L,fvec[3],R[3],T[3],P[3],m1,m2;
			// get current face normal
			getNormal(n1);
			m1=sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);
			n1[0]=n1[0]/m1;
			n1[1]=n1[1]/m1;
			n1[2]=n1[2]/m1;
			// get intersecting face normal
			(*i)->getNormal(n2);
			m2=sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]);
			n2[0]=n2[0]/m2;
			n2[1]=n2[1]/m2;
			n2[2]=n2[2]/m2;
			// compute resultant vector
			R[0]=n1[0]+n2[0];
			R[1]=n1[1]+n2[1];
			R[2]=n1[2]+n2[2];
			// compute cross product, T, as N1 X N2
			T[0] = n1[1]*n2[2]-n1[2]*n2[1];
			T[1] = n1[2]*n2[0]-n1[0]*n2[2];
			T[2] = n1[0]*n2[1]-n1[1]*n2[0];
			// compute cross product, P, as T X R
			P[0] = T[1]*R[2]-T[2]*R[1];
			P[1] = T[2]*R[0]-T[0]*R[2];
			P[2] = T[0]*R[1]-T[1]*R[0];
			// P length
			L=sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
			if (facesInSameNeighborhood(c,this,*i)){
				// if P and N1 point in opposite directions
				// i.e. if dot product of P and N1 is less than 0
				if((P[0]*n1[0]+P[1]*n1[1]+P[2]*n1[2])<0){
					// then invert direction of P
					P[0]=-P[0];
					P[1]=-P[1];
					P[2]=-P[2];
				}
			} else {
				// if P and N1 point in the same direction
				// i.e. if dot product of P and N1 is greater than 0
				if((P[0]*n1[0]+P[1]*n1[1]+P[2]*n1[2])>0){
					// then invert direction of P
					P[0]=-P[0];
					P[1]=-P[1];
					P[2]=-P[2];
				}
			}
			// compute force
			fvec[0]=P[0]/L*INTERSECTION_WEIGHT/(SEPARATION_WEIGHT+EDGE_STRETCH_WEIGHT
					+ANGLE_STRETCH_WEIGHT+INTERSECTION_WEIGHT) ;
			fvec[1]=P[1]/L*INTERSECTION_WEIGHT/(SEPARATION_WEIGHT+EDGE_STRETCH_WEIGHT
					+ANGLE_STRETCH_WEIGHT+INTERSECTION_WEIGHT) ;
			fvec[2]=P[2]/L*INTERSECTION_WEIGHT/(SEPARATION_WEIGHT+EDGE_STRETCH_WEIGHT
					+ANGLE_STRETCH_WEIGHT+INTERSECTION_WEIGHT) ;
			// set intersection force
			this->addForceToFace(fvec);
		}
	}
}


bool Face::computeIntersectionForce(Container *c) {
	std::vector<Face*> dummy;
	if(getFaceIntersection(c,false,dummy)){
		calculateIntersectionForce(c);
		return true;
	}
	return false;
}

void Container::computeFaceIntersectionForce(void) {
	cout << "Iteration 0: ";
	cout << "compute face intersection force...";
	cout.flush();
	// for each object
	for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
		// for each face in object
		for (std::vector<Face*>::iterator j=(*i)->f.begin();j!=(*i)->f.end();j++) {
			(*j)->computeIntersectionForce(this);
		}
	}
	cout << "complete.\n";
	cout.flush();
}

// #####################################################
// #####################################################
/*
bool Face::getFaceIntersectionCheck(Container *c,Space &s,hashtable_f &nb) {
	// get boxes in which this face lies
	std::vector<Box*> bp;
	s.computeBoxesToCheck(this,bp);
	// add box*s to multimap nb
	for(std::vector<Box*>::iterator i=bp.begin();i!=bp.end();i++){
		nb.insert(std::make_pair(this,*i));
	}
	// grab all box*s for this face
	std::pair<tf_iterator,tf_iterator> p;
	p=nb.equal_range(this);
	// if none found, then error
	if(p.first==p.second){
		cout << "ERROR: NO BOXES FOR\n" << "Face " << this->index << " " 
			<< (this->v[0])->index << " " << (this->v[1])->index << " "
			<< (this->v[2])->index << endl;
			exit(1);
	}
	// collect all faces in chosen boxes
	std::vector<Face*> of;
	for(tf_iterator i=p.first;i!=p.second;i++){
		of.insert(of.end(),(*i).second->f.begin(),(*i).second->f.end());
	}
	// sort and keep unique faces
	sort(of.begin(),of.end());
	std::vector<Face*>::iterator j;
	j = unique(of.begin(),of.end());
	of.assign(of.begin(),j);
    // for each unique face
	j=of.begin();
	while(j!=of.end()){

		// if unique face is not same as current face
		// and unique face and current face are from the same object
		if(*j!=this && (*j)->v[0]->o==this->v[0]->o){
			// if faces intersect
			if (c->checkFaceFaceIntersections(this,*j)) {
//				cout << endl << endl;
//				cout << "adjacent face\n";
//				this->printFace(this->v[0]->o->name);
//				cout << endl;
//				cout << "intersecting face\n";
//				(*j)->printFace((*j)->v[0]->o->name);
//				cout << endl;
				return true;
			}
		}
		j++;
	}
	return false;
}*/

bool Face::getFaceIntersectionCheck(Container *c,Space &s,hashtable_f &nb) {
	// get boxes in which this face now lies
	std::vector<Box*> bp;
	s.computeBoxesToCheck(this,bp);
	// add box*s to multimap nb
	for(std::vector<Box*>::iterator i=bp.begin();i!=bp.end();i++){
		nb.insert(std::make_pair(this,*i));
	}
	// collect all unique faces in chosen boxes
	f_set of;
	for(std::vector<Box*>::iterator i=bp.begin();i!=bp.end();i++){
		of.insert((*i)->f.begin(),(*i)->f.end());
	}
    // for each unique face
	for(fs_iterator i=of.begin();i!=of.end();i++){
//		// and unique face and current face are from the same object
//		if(*i!=this && (*i)->v[0]->o==this->v[0]->o){
		// if unique face is not same as current face
		if(*i!=this){
			// if faces intersect
			if (c->checkFaceFaceIntersections(this,*i)) { 
				// if this face previously had no intersections, then return true
				if(!this->faceInTable_intf()){return true;}
				// else face was previously intersected
				else {
					// grab previous intersecting faces of this face
					std::vector<Face*> *fv=this->getIntersectingFaces();
					bool same=false;
					// for each previous intersecting face
					for(std::vector<Face*>::iterator p=(*fv).begin();p!=(*fv).end();p++){
						// if current intersecting face is of same object
						//  as previous intersecting face
						if((*p)->v[0]->o==(*i)->v[0]->o){same=true;}
					}
					// if intersecting objects have changed from before to now, then return true
					if(!same){return true;}
				}
			}
		}
	}
	return false;
}

bool Face::getFaceIntersection(Container *c,bool just_check,std::vector<Face*> &int_f) {
	std::vector<Face*>::iterator j;
	bool flag = false;
	// reset face element in table
	clearFaceFromTable_intf();
	// for each box in which face lies, add faces in box to vector
	std::vector<Face*> of;
	for(std::vector<Box*>::iterator i=b.begin();i!=b.end();i++){
		of.insert(of.end(),(*i)->f.begin(),(*i)->f.end());
	}
	// sort and keep unique faces
	sort(of.begin(),of.end());
	j = unique(of.begin(),of.end());
	of.assign(of.begin(),j);
    // for each unique face
	j=of.begin();
	while(j!=of.end()){
		// if unique face is not same as current face
		if(*j!=this){
			// if faces intersect
			if (c->checkFaceFaceIntersections(this,*j)) {
				if(!just_check){
					c->ti++;
					if (v[0]->o==(*j)->v[0]->o){c->si++;}
					// save intersecting face* to this face's intersecting face vector
					if(!faceInTable_intf()){addFaceToTable_intf();}
					addFaceToVector(*j);
				} else { int_f.push_back(*j); }
				// return
				flag = true;
			}
		}
		j++;
	}
	return flag;
}

void Monitor::loadTopN(Container *c){
	topN.clear();
	// for each object in model
	for(std::vector<Object*>::iterator j=c->o.begin();j!=c->o.end();j++){
		// for each vertex in object
		for(std::vector<Vertex*>::iterator i=(*j)->v.begin();i!=(*j)->v.end();i++){
			// if vertex has a closest point
			if((*i)->cl!=NULL){
				// compute flip of associated edges
//				c->computeEdgeFlip(*i);
				// compute new vertex coords
				double pH[3]; // new holding position coordinates (x,y,z)
				// load energy map
				(*i)->computeNewCoords(c,pH,c->gain);
				// compute virtual displacement
				double vd = (pH[0]-(*i)->pN[0])*(pH[0]-(*i)->pN[0])+
							(pH[1]-(*i)->pN[1])*(pH[1]-(*i)->pN[1])+
							(pH[2]-(*i)->pN[2])*(pH[2]-(*i)->pN[2]);
				// add virtual squared displacement to topN
				if((*i)->cl==NULL){
					cout << "WHAT THE! It just had a closest and now it doesn't. What changed?\n";
					exit(0);
				}
				topN.insert(std::make_pair(vd,*i));
				
			}
		}
	}
}

void Monitor::printVertexSelect(Container &c,int interval){
	int zero=0;
	char file[1024];
	// open log file
	sprintf(file,"%s%s.%d",OUTPUT_DATA_DIR,VERTEX_SELECTION_FILE,interval);
	std::ofstream myfile (file);
	for(std::vector<Object*>::iterator i=c.o.begin();i!=c.o.end();i++){
		for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++){
			vhm_iterator t = touch_map.find(*j);
			// if vertex was touched
			if(t!=touch_map.end()){
				myfile << touch_map[*j] << endl;
			} else {
				myfile << zero << endl;
			}
		}
	}
	myfile.close();
}

// IDENTICAL TO void Monitor::printVertexSelect(Container &c,int interval){
// EXCEPT FOR int TO char* ARGUMENT SWAP
void Monitor::printVertexSelect(Container &c,char *str){
	int zero=0;
	char file[1024];
	// open log file
	sprintf(file,"%s%s.%s",OUTPUT_DATA_DIR,VERTEX_SELECTION_FILE,str);
	std::ofstream myfile (file);
	for(std::vector<Object*>::iterator i=c.o.begin();i!=c.o.end();i++){
		for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++){
			vhm_iterator t = touch_map.find(*j);
			// if vertex was touched
			if(t!=touch_map.end()){
				myfile << touch_map[*j] << endl;
			} else {
				myfile << zero << endl;
			}
		}
	}
	myfile.close();
}



double getSum(double sumSoFar, const std::pair<double,Vertex*> &p){
	return sumSoFar + p.first;
}

void Container::updateStats(double d){
	if(d){
		N++;
		md[0]=md[1];
		md[1]=md[0]*(N-1)/N+d/N;
		if(d<d_min){d_min=d;}
		if(d>d_max){d_max=d;}
	}
}

void Monitor::updateTouchMap(Vertex *v){
	// if vertex in touch_map, then increment
	if(touch_map.find(v)!=touch_map.end()){ touch_map[v]++; }
	// else init to one touch
	else { touch_map[v]=1; }
}

void Monitor::updateAvg(double e){
	// update window sum
	sum=sum-*window+e;
	// update window average
	avg_old=avg_new;
	avg_new=sum/num;
	// replace window element
	*window=e;
	// move window
	if(window!=end){window++;}
	else {window=begin;}
}

void Monitor::clearAvg(void){
	sum=avg_new=avg_old=0.0;
	window=begin;
	for(int i=0;i<num;i++){
		*(window++)=0.0;
	}
	window=begin;
}

void Monitor::initAvg(void){
	sum=avg_new=avg_old=0.0;
	num=100;
	window = new double[num];
	begin=window;
	for(int i=0;i<num;i++){
		*(window++)=0.0;
	}
	window--;
	end=window;
	window=begin;
}

void Monitor::freeAvg(void){
	window=begin;
	delete[] window;
}

void Monitor::prep(Container *c){
	loadTopN(c);
	saveOld();
	initAvg();
}

void Monitor::getBoxesPerFace(Container &c){
	bpf_min = 1E30;
	bpf_max = -1E30;
	int n=0;
	double summ=0;
    // for each object in container
	for (std::vector<Object*>::iterator i=c.o.begin();i!=c.o.end();i++) {
        // for each face in object
		for (std::vector<Face*>::iterator j=(*i)->f.begin();j!=(*i)->f.end();j++) {
			int a = (*j)->b.size();
/*			if(a==36){
				cout << endl;
				cout << endl;
				(*j)->printFace((*j)->v[0]->o->name);
				cout << endl;
				cout << endl;
				exit(0);
Face <obj>d000_FILTERED<ind>110845
[v0 97999 3943.42 5812.89 3038.56]
[v1 109604 4001.59 5668.21 3033.92]
[v2 109986 3994.99 5618.39 3034.18]
which yields face edge lengths of 156, 201, 50
			}*/
			n++;
			summ+=a;
			if(a>bpf_max){bpf_max=a;} 
			if(a<bpf_min){bpf_min=a;} 
		}
	}
	bpf_mean=summ/n;
}

void Monitor::getFacesPerBox(Space &s){
	fpb_min = 1E30;
	fpb_max = -1E30;
	int n=0;
	double summ=0;
	// for each box in space
	for (std::vector<Box*>::iterator i=s.b.begin();i!=s.b.end();i++) {
		int a = (*i)->f.size();
		if(a){
			n++;
			summ+=a;
			if(a>fpb_max){fpb_max=a;} 
			if(a<fpb_min){fpb_min=a;} 
		}
	}
	fpb_mean=summ/n;
}

void Monitor::printVerticesWithCP(void){
	cout << endl << endl;
    // for each topN
	for (tv_iterator i=topN.begin();i!=topN.end();i++) {
		(*i).second->printVertexCP();
		cout << endl;
	}
	exit(0);
}

void updateGate(Vertex *v,hashtable_v &gate,int time_out){
	// for each element in hash table
	vhm_iterator i=gate.begin();
	while(i!=gate.end()){
		vhm_iterator j = i;j++;
		// if count==1
		if ((*i).second==1){
			// remove element
			gate.erase(i);
			i=j;
		} else {
			// decrement count
			(*i).second--;
			i++;
		}
	}
	gate[v]=time_out;
}

void punishGate(Vertex *v,hashtable_v &gate,int time_out){
	gate[v]=time_out;
}

void Monitor::initRefrac(void){
	refrac_s.clear();
	refrac_l.clear();
//	next=refrac.begin();
}

bool Monitor::Refracted(Vertex *v){
	if (refrac_s.find(v)!=refrac_s.end()){return true;}
	else							{return false;}
}

void Monitor::updateRefractoryWindow(Vertex *v){
	///// create iterator to soon-to-be oldest element in set /////
/*	vs_iterator i=next;
	// this handles case where set is empty so i=next=begin==end
	if(i==refrac.end()){i=refrac.begin();}
	else {
		// all other cases
		i++;
		if(i==refrac.end()){i=refrac.begin();}
	}
*/
	///// erase oldest element /////
	if(refrac_l.size()==REFRACTORY_PERIOD){
		refrac_s.erase(refrac_l.front());
		refrac_l.pop_front();
	}
	///// store new element /////
	refrac_l.push_back(v);
	refrac_s.insert(v);
	///// update next /////
//	next=i;
}

// #####################################################
// #####################################################

