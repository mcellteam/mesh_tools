bool parse(int argc,char **argv,std::string message)
{
  bool freeze=false;
  int c;
  opterr=0;
  char *eptr=NULL;
  char filename[FILENAME_SIZE];
  char *temp;
  while((c=getopt(argc,argv,"i:o:t:f:h")) != -1)
    switch(c)
    {
      case 'i':
        // specify input directory 
        strcpy(filename,optarg);
        temp=strrchr(filename,'/');
        if(!temp) {strcat(filename,"/");}
        else if(*++temp) {strcat(filename,"/");}
        INPUT_DATA_DIR = filename;
        break;
      case 'o':
        // specify input directory 
        strcpy(filename,optarg);
        temp=strrchr(filename,'/');
        if(!temp) {strcat(filename,"/");}
        else if(*++temp) {strcat(filename,"/");}
        OUTPUT_DATA_DIR = filename;
        break;
      case 't':
        // specify target extracellular width 
        TARGET_SEPARATION = strtod(optarg,&eptr);
        break;
      case 'f':
        // specify frozen vertices file
        FROZEN_VERTICES_FILE = optarg;
        freeze = true;
        break;
      case 'h':
        cout << message << endl;
        abort ();
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

  if(optind<argc)
  {
    for (int index=optind;index<argc;index++)
    {
      cout << "\nError: extra argument found = " << argv[index] << endl;
    }
    cout << message << endl;
    exit(0);
  }
  return freeze;
}

bool checkIntSize(void){
  ///// check that assumption of 32 bit int is correct /////
  if (32==sizeof(int)*8){ return true; }
  else {
    cout << "Error. Int is not 32 bits, sizeof(int) "
          << sizeof(int) << endl;
    return false;
  }
}

// #####################################################
// #####################################################



//bool getPointEdgeDistance(double p[3],Face *f){
bool getPointEdgeDistance(double p[3],Face *f,bool inspect){
  // for each face edge
  for (int i=0;i<3;i++) {
    double Ax = f->v[i]->pN[0];
    double Ay = f->v[i]->pN[1];
    double Az = f->v[i]->pN[2];
    double Bx = f->v[(i+1)%3]->pN[0];
    double By = f->v[(i+1)%3]->pN[1];
    double Bz = f->v[(i+1)%3]->pN[2];
    double Cx = Bx-Ax;
    double Cy = By-Ay;
    double Cz = Bz-Az;
    //double uDen = (Bx-Ax)*(Bx-Ax)+(By-Ay)*(By-Ay)+(Bz-Az)*(Bz-Az);
    double uDen = Cx*Cx+Cy*Cy+Cz*Cz;

    // DEBUG
    if(inspect==true){
      cout.precision(12);
      cout << "\nContainer::getPointEdgeDistance: "
            << "face index = " << f->index
            << ", uDen " << uDen << endl;
    }
    // DEBUG

    if(uDen)
    {
      // u = AdotA-AdotB-AdotC+BdotC)/uDen
      //double u =((p[0]-Ax)*(Bx-Ax)+(p[1]-Ay)*(By-Ay)+(p[2]-Az)*(Bz-Az))/uDen;
      double u =((p[0]-Ax)*Cx+(p[1]-Ay)*Cy+(p[2]-Az)*Cz)/uDen;
      // no need to check for u ==0 and u ==1, since current 
      // vertex/face plane coincidence was checked previously.
      // Closest point on face edge line to current vertex
      // occurs on face edge between face vertices

      // DEBUG
      if(inspect==true){
        cout.precision(20);
        cout << "\nContainer::getPointEdgeDistance: "
              << "face index = " << f->index
              << ", u " << u << endl;
      }
      // DEBUG

      if (u>0 && u<1)
      {
        //double a=Ax+u*(Bx-Ax);
        //double b=Ay+u*(By-Ay);
        //double c=Az+u*(Bz-Az);

        // DEBUG
        //        if(inspect==true){ cout.precision(12);
        //          cout << "\nContainer::getPointEdgeDistance: "
        //                << "face index = " << f->index
        //                << ", must be indistinguishable: "
        //                << "[" << p[0] << "," << a << "]" << "<" << !distinguishable(p[0],a,1E-9) << "> "
        //                << "[" << p[1] << "," << b << "]" << "<" << !distinguishable(p[1],b,1E-9) << "> "
        //                << "[" << p[2] << "," << c << "]" << "<" << !distinguishable(p[2],c,1E-9) << "> "
        //                << endl;
        //        }
        // DEBUG

        double a=Ax+u*Cx;
        if(!distinguishable(p[0],a,1E-9))
        {
          double b=Ay+u*Cy;
          if(!distinguishable(p[1],b,1E-9))
          {
            double c=Az+u*Cz;
            if(!distinguishable(p[2],c,1E-9))
            {
              return true;
            }
          }
        }
      }
    }
  }
  return false;
}

// #####################################################
// #####################################################

void copyControlFile (void) {
  cout << "Copy control settings..........................";
  cout.flush();
  char line[FILENAME_SIZE],filein[FILENAME_SIZE],fileout[FILENAME_SIZE];
  // open input data file
  sprintf(filein,"controls.cc");
  std::ifstream inFile(filein);
  if (inFile.fail()) // if stream cannot be opened
  { 
    cout << "\n\ncopyControlFile: Can not open input file ["
          << filein << "]\n\n";
    exit(1); 
  }
  // open output data file
  sprintf(fileout,"%s%s",OUTPUT_DATA_DIR.c_str(),CONTROL_FILE);
  std::ofstream outFile(fileout);
  if (outFile.fail()) // if stream cannot be opened
  {
    cout << "\n\ncopyControlFile: Can not open output file ["
          << fileout << "]\n\n";
    exit(1); 
  }
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

time_t recordTime(std::ofstream& this_file, time_t before, char string[128])
{
  time_t after,diff;
  after = time (NULL);
  diff = after-before;
  this_file << string << " " << diff << " seconds\n";
  this_file.flush();
  return after;
}

// #####################################################
// #####################################################

int biggest(double *pn)
{
  const double x[3] = {fabs(pn[2]),fabs(pn[0]),fabs(pn[1])};
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
        return 2;
      } else {
        // x[0] < x[1]
        // x[2] <= x[1]
        // x[0] is smallest
        // x[1] is biggest or tied with x[2]
        return 1;
      }
    } else {
      // x[0] < x[1]
      // x[2] <= x[0]
      // x[2] is smallest or tied with x[0]
      // x[1] is biggest
      return 1;
    }
  } else {
    // x[1] <= x[0]
    if ( x[0] < x[2] ) {
      // x[1] <= x[0]
      // x[0] < x[2]
      // x1 is smallest or tied with x[0]
      // x2 is biggest
      return 2;
    } else {
      return 0;
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

void checkLineFaceIntersection(Face * const f,double lp[2][3],bool &line_flag,
                               bool &poly_flag, bool &poly_edge_flag,bool ignore_line_flag) {

  // Use line points, polygon normal and one polygon vertex
  // to compute where line intersects plane of polygon (if not parallel).
  //
  // pI = p1+u*(p2-p1) where u=dot(N,(p3-p1))/dot(N,(p2-p1))
  // pI = intersection of line and plane
  // p3 = any point on plane
  // p2,p1 = end points of line
  // N = plane normal
  // lp  = line_points

  // initialize flags
  line_flag=poly_flag=poly_edge_flag=false;
  // get face normal
  double pn[3];
  f->getNormal(pn);
  // denominator of u
  const double den = pn[0]*(lp[1][0]-lp[0][0])+pn[1]*(lp[1][1]-lp[0][1])+pn[2]*(lp[1][2]-lp[0][2]);
  // if line and polygon plane are not parallel
  if(den)
  {
    //pvc = polygon_vertex_coordinates
    Vertex * const v0=f->v[0],* const v1=f->v[1],* const v2=f->v[2]; 
    // point of intersection
    //const double u = (pn[0]*(v0->pN[0]-lp[0][0]) 
    //                + pn[1]*(v0->pN[1]-lp[0][1]) 
    //                + pn[2]*(v0->pN[2]-lp[0][2]))/den;
    const double u = (pn[0]*(v0->pN[0]-lp[0][0]) 
                      + pn[1]*(v0->pN[1]-lp[0][1]) 
                      + pn[2]*(v0->pN[2]-lp[0][2]));
    // if polygon cuts through line
    //if(u > DOUBLE_EPSILON && u < (1-DOUBLE_EPSILON)){ line_flag = true;}
    //if((u>0.0) && (u<1.0)){ line_flag = true;}
    if(den>0)
    {
      if((u>0.0) && (u<den)){ line_flag = true;}
    }
    else
    {
      if((u<0.0) && (u>den)){ line_flag = true;}
    }
    // check early exit
    if(line_flag==false && ignore_line_flag==false){return;}
    // compute polygon double-area on each of three principal planes
    //double I[3];//,d[3];
    // area = {xy,yz,zx}
    const int big = biggest(pn);
    // pI = p1+u*(p2-p1)
    // where lp[1][] is p2 and lp[0][] is p1
    const int t = (big+1)%3;
    //const int tv[2] = {big,(big+1)%3};
    //const double pI[2] = { (1-u)*lp[0][big]       + u*lp[1][big],
    //                       (1-u)*lp[0][(big+1)%3] + u*lp[1][(big+1)%3]};
    const int tv[2] = {big,t};
    //const double pI[2] = { (1-u)*lp[0][big] + u*lp[1][big],
    //                       (1-u)*lp[0][t]   + u*lp[1][t]};
    const double pI[2] = { (den-u)*lp[0][big] + u*lp[1][big],
      (den-u)*lp[0][t]   + u*lp[1][t]};
    //I[big]=pI[0];
    //I[(big+1)%3]=pI[1];
    //I[(big+2)%3] = (1-u)*lp[0][(big+2)%3] + u*lp[1][(big+2)%3];
    ////////// is point of intersection on other polygon? //////////
    //pvc = polygon_vertex_coordinates
    //bool inspect=false;
    // does point of intersection lie on polygon edge?
    //    if (getPointEdgeDistance(I,f,inspect)){poly_edge_flag = true;}
    // if point of intersection is not on polygon edge

    //		if (poly_edge_flag==false) {
    // compute three determinants
    //const double a = v0->pN[tv[0]]-pI[0];
    //const double b = v0->pN[tv[1]]-pI[1];
    //const double c = v1->pN[tv[0]]-pI[0];
    //const double d = v1->pN[tv[1]]-pI[1];
    //const double e = v2->pN[tv[0]]-pI[0];
    //const double ff = v2->pN[tv[1]]-pI[1];
    const double a  = den*v0->pN[tv[0]]-pI[0];
    const double b  = den*v0->pN[tv[1]]-pI[1];
    const double c  = den*v1->pN[tv[0]]-pI[0];
    const double d  = den*v1->pN[tv[1]]-pI[1];
    const double e  = den*v2->pN[tv[0]]-pI[0];
    const double ff = den*v2->pN[tv[1]]-pI[1];
    //const double det[3] = {
    // (v0->pN[tv[0]]-pI[0])*(v1->pN[tv[1]]-pI[1])
    //-(v1->pN[tv[0]]-pI[0])*(v0->pN[tv[1]]-pI[1]),
    // (v1->pN[tv[0]]-pI[0])*(v2->pN[tv[1]]-pI[1])
    //-(v2->pN[tv[0]]-pI[0])*(v1->pN[tv[1]]-pI[1]),
    // (v2->pN[tv[0]]-pI[0])*(v0->pN[tv[1]]-pI[1])
    //-(v0->pN[tv[0]]-pI[0])*(v2->pN[tv[1]]-pI[1])};
    const double det = a*d-c*b ;
    const double val1=det*(c*ff-e*d);
    const double val2=det*(e*b-a*ff);
    //const double det[3] = { a*d-c*b , c*ff-e*d , e*b-a*ff};
    //const double val1=det[0]*det[1];
    //const double val2=det[0]*det[2];
    // proceed if determinants are DOUBLE_EPSILON away from zero
    if (val1>0){
      if (val2>0){
        // line intersects polygon plane inside polygon
        poly_flag = true;
      }
    }
    if(fabs(val1)<DOUBLE_EPSILON*den*den*den*den && fabs(val2)<DOUBLE_EPSILON*den*den*den*den)
    {
      // does point of intersection lie on polygon edge?
      poly_edge_flag = true;
      //if (getPointEdgeDistance(I,f,inspect)){poly_edge_flag = true;}
    }
    //		}
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
      checkLineFaceIntersection(of,lp,line_flag,poly_flag,poly_edge_flag,false);
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
      checkLineFaceIntersection(cf,lp,line_flag,poly_flag,poly_edge_flag,false);
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
/*
    double Vertex::getSeparationForceEnergy(double force[3],bool flag,Container *c){
// get closest point
double pC[3];
c->computePC(cl,this,pC,false);
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
//	return (SEPARATION_WEIGHT/100.0)/2.0*se*se;
return (SEPARATION_WEIGHT/200.0)*se*se;
}*/
/*
    double Vertex::getSeparationForceEnergy(double force[3],bool flag,Container *c){
    if(cl!=NULL){
// get closest point
double pC[3];
c->computePC(cl,this,pC);
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

// set target separation
double TS = 0;
// if closest face is inside vertex neighborhood
if(c->faceInNeighborhood(cl,this)==true){ TS = LOOP_TARGET_SEPARATION*SCALE; }
else { TS = TARGET_SEPARATION*SCALE; }


double se;
if(o->vertexIsNice(this)==false){
se=sd+TS;
} else{
se=sd-TS;
}
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
return (SEPARATION_WEIGHT/200.0)*se*se;
} else {
// cl==NULL, i.e. no closest point
// therefore no contribution to force
// peg energy at maximum possible value, i.e. max separation error
// maximum separation error would occur when closest point found
// at SEARCH_RADIUS distance from current vertex AND vertex is not nice
// so max(se) = SEARCH_RADIUS + TARGET_SEPARATION
double max_se = sqrt(SEARCH_RADIUS_SQ*SCALE*SCALE) + LOOP_TARGET_SEPARATION*SCALE;
return (SEPARATION_WEIGHT/200.0)*max_se*max_se;
}
}
*/
double Vertex::getSeparationForceEnergy(double force[3],bool flag,Container *c){
  if(cl!=NULL){
    // get closest point
    double pC[3];
    c->computePC(cl,this,pC);
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

    // set target separation
    double TS = 0;
    // if closest face is inside vertex neighborhood
    if(c->faceInNeighborhood(cl,this)==true){ TS = LOOP_TARGET_SEPARATION*SCALE; }
    else { TS = TARGET_SEPARATION*SCALE; }

    double se;
    if(o->vertexIsNice(this)==false){
      se=sd+TS;
    } else{
      se=sd-TS;
    }
    if(TARGET_SEPARATION==0)
    {
      cout << "Vertex::getSeparationForceEnergy: we got problems."
            << "TARGET_SEPARATION==0.\n";
      exit(0);
    }
    se = se/TARGET_SEPARATION;
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
      //double scaled_spring_force = (SEPARATION_WEIGHT/100.0)*se/sd;
      double scaled_spring_force = SEPARATION_GAIN*(SEPARATION_WEIGHT/100.0)*se/sd;
      // force cartesian component = spring_force * unit vector
      // where unit vector = cartesian separation_distance component / separation_distance
      // force cartesian component = scaled_spring_force * cartesian separation_distance component
      force[0]+=scaled_spring_force*s[0];
      force[1]+=scaled_spring_force*s[1];
      force[2]+=scaled_spring_force*s[2];
    }
    // return energy
    return (SEPARATION_WEIGHT/200.0)*se*se;
  } else {
    // cl==NULL, i.e. no closest point
    // therefore no contribution to force
    // peg energy at maximum possible value, i.e. max separation error
    // maximum separation error would occur when closest point found
    // at SEARCH_RADIUS distance from current vertex AND vertex is not nice
    // so max(se) = SEARCH_RADIUS + TARGET_SEPARATION
    double max_se = sqrt(SEARCH_RADIUS_SQ*SCALE*SCALE) + LOOP_TARGET_SEPARATION*SCALE;
    return (SEPARATION_WEIGHT/200.0)*max_se*max_se;
  }
}

double Face::getAspectRatio(void)
{
  /* Make triangle edge vectors */
  double va[3]={v[1]->pN[0]-v[0]->pN[0],v[1]->pN[1]-v[0]->pN[1],v[1]->pN[2]-v[0]->pN[2]};
  double vb[3]={v[2]->pN[0]-v[1]->pN[0],v[2]->pN[1]-v[1]->pN[1],v[2]->pN[2]-v[1]->pN[2]};
  double vc[3]={v[0]->pN[0]-v[2]->pN[0],v[0]->pN[1]-v[2]->pN[1],v[0]->pN[2]-v[2]->pN[2]};
  double vbase[3]={0,0,0};
  double vopp[3]={0,0,0};

  /* Find length of longest edge */
  double lmax=-1e30;
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


double Edge::getStretchForceEnergy(Vertex *v,double force[3],bool flag,Container *c){
  if (0) {cout << c->num_files << endl;}
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
  // get pointer to adjacent vertex
  Vertex *av;
  if(v1==v){av=v2;}
  else {av=v1;}
  // compute separation vector
  double s[3];
  for(int i=0;i<3;i++){ s[i]=av->pN[i]-v->pN[i]; }
  // compute separation distance //////////
  double sd = sqrt(dot(s,s));
  // determine which reference length to use
  double ref1 = (f1->e[0]->l+f1->e[1]->l+f1->e[2]->l)/3.0;
  double ref2 = (f2->e[0]->l+f2->e[1]->l+f2->e[2]->l)/3.0;
  // compute separation error (signed value)
  double se1 = (sd-ref1)/ref1;
  double se2 = (sd-ref2)/ref2;
  double force_magn = ASPECT_SCALE*STRETCH_GAIN*(se1+se2)/2.0;
  if(flag){
    // force contribution
    // spring_force = spring constant * stretch
    // scaled_spring force = spring_force/edge length
    double scaled_spring_force = (EDGE_STRETCH_WEIGHT/100.0)*force_magn/sd;
    // force cartesian component = spring_force * unit vector
    // where unit vector = (adjacent vertex position - target vertex position) / edge length
    // force cartesian component = scaled_spring_force * cartesian component difference
    force[0]+=scaled_spring_force*s[0];
    force[1]+=scaled_spring_force*s[1];
    force[2]+=scaled_spring_force*s[2];
  }
  // energy contribution
  //return (EDGE_STRETCH_WEIGHT/100.0)/STRETCH_EXPONENT*force_magn*se;
  return 0.0;
}

double Vertex::getEdgeStretchForceEnergy(double force[3],bool flag,Container *c){
  double en_ergy=0;
  std::vector<Edge*> e;
  e.reserve(VECTOR_RESERVE);
  getAdjacentEdges(e);
  // for each adjacent edge of current vertex
  for (std::vector<Edge*>::iterator j=e.begin();j!=e.end();j++) {
    // compute force and energy of edge
    en_ergy+=(*j)->getStretchForceEnergy(this,force,flag,c);
  }
  return en_ergy;
}

double Vertex::getEdgeAngleFaceIntersectionForceEnergy(double force[3],bool flag){
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  double en_ergy=0;
  std::vector<Edge*> ae,oe;
  ae.reserve(VECTOR_RESERVE);
  oe.reserve(VECTOR_RESERVE);

  // for each adjacent face of current vertex
  for(std::vector<Face*>::iterator i=f.begin();i!=f.end();i++){
    // for each face edge
    for(int j=0;j<3;j++){
      (*i)->e[j]->getVertices(v1,v2,o1,o2);
      // if either edge vertex is current vertex
      if(v1==this || v2==this){
        // add edge to adjacent edge vector, ae
        ae.push_back((*i)->e[j]);
      } else {
        // add edge to non-adjacent edge vector, oe
        oe.push_back((*i)->e[j]);
      }
    }
    if(flag){
      // if adjacent face has intersection force
      if((*i)->faceInTable_iv()){
        // add intersection_force
        (*i)->getForceFromTable(force);
      }
    }
  }
  // keep unique edges
  sort(ae.begin(),ae.end());
  std::vector<Edge*>::iterator new_end = unique(ae.begin(),ae.end());
  ae.assign(ae.begin(),new_end);
  sort(oe.begin(),oe.end());
  new_end = unique(oe.begin(),oe.end());
  oe.assign(oe.begin(),new_end);
  // debug
  bool activate=false;
  //if (o->name=="a066_clipped" && index==222){activate=true;}
  // debug
  // for each non adjacent edge
  /*for (std::vector<Edge*>::iterator i=oe.begin();i!=oe.end();i++) {
    (*i)->getVertices(v1,v2,o1,o2);
  // force and energy contribution
  if (o1==this){en_ergy+=(*i)->getForceEnergy(0,force,flag,activate);}
  else		 {en_ergy+=(*i)->getForceEnergy(1,force,flag,activate);}
  }
  */
  // for each adjacent edge
  /**/	for (std::vector<Edge*>::iterator i=ae.begin();i!=ae.end();i++) {
    // force and energy contribution
    en_ergy+=(*i)->getReactionForceEnergy(force,flag,activate);
  }
  /**/
  return en_ergy;
}

void Vertex::getForceEnergy(double force[3],Container *c) {
  ////////// get force and energy contributions from membrane separation //////////
  // true -> compute force
  // DEBUG
  if (false) {cout << c->o.size() << endl;}
  bool print = false;
  double fmag = 0.0;
  double f_old[3];
  f_old[0] = force[0];
  f_old[1] = force[1];
  f_old[2] = force[2];
  /*	
        if( 
  // (o->name=="d000_clipped" && ( index==158)) ||
  (o->name=="a066_clipped" && ( index==222)) 
  )
  {
  print=true;
  cout << "Vertex::getForceEnergy: "
  << o->name << "->" << index << endl;
  }
  */      
  if(print==true){
    cout << "Vertex::getForceEnergy: force ["
          << force[0] << " "
          << force[1] << " "
          << force[2] << "\n";
  }
  // DEBUG

  getSeparationForceEnergy(force,true,c);
  // DEBUG
  if(print==true){
    double diff = sqrt(dot(force,force))-fmag;
    fmag = fmag+diff;
    f_old[0] = force[0];
    f_old[1] = force[1];
    f_old[2] = force[2];
    cout << "Vertex::getForceEnergy: sep force       ["
          << force[0] << " "
          << force[1] << " "
          << force[2] << "]"
          << ", diff = " << diff << endl;
  }

  // DEBUG
  ////////// get force and energy contributions from adjacent edge stretching //////////
  getEdgeStretchForceEnergy(force,true,c);
  // DEBUG
  if(print==true){
    double diff = sqrt(dot(force,force))-fmag;
    fmag = fmag+diff;
    cout << "Vertex::getForceEnergy: stretch force   ["
          << force[0]-f_old[0] << " "
          << force[1]-f_old[1] << " "
          << force[2]-f_old[2] << "]"
          << ", diff = " << diff << endl;
    f_old[0] = force[0];
    f_old[1] = force[1];
    f_old[2] = force[2];
  }
  // DEBUG
  ////////// get force contribution from edge angles and face intersections //////////
  getEdgeAngleFaceIntersectionForceEnergy(force,true);
  // DEBUG
  if(print==true){
    double diff = sqrt(dot(force,force))-fmag;
    fmag = fmag+diff;
    cout << "Vertex::getForceEnergy: angle force     ["
          << force[0]-f_old[0] << " "
          << force[1]-f_old[1] << " "
          << force[2]-f_old[2] << "]"
          << ", diff = " << diff << endl;
    cout << "Vertex::getForceEnergy: resultant force ["
          << force[0] << " "
          << force[1] << " "
          << force[2] << "]\n\n";
  }
  // DEBUG
}

// #####################################################
// #####################################################
void Vertex::computeNewCoords(Container *c,double pH[3],double gain) {
  // get force
  double force[3]={0,0,0};
  getForceEnergy(force,c);
  /*	// DEBUG
        if(match(2174,"d000_FILTERED_SMOOTH_SMOOTH")==true){
        cout << "Vertex::computeNewCoords: "
        << "Vertex <obj>" << o->name << "<ind>" << index << ": Force ["
        << force[0] << ", "
        << force[1] << ", "
        << force[2] << "]";
  // get closest point
  double pC[3];
  c->computePC(cl,this,pC);
  cout << ", closest point ["
  << pC[0] << ", "
  << pC[1] << ", "
  << pC[2] << "]\n";

  }
  // DEBUG*/
  // compute new vertex coordinates
  for (int k=0;k<3;k++) { pH[k] = gain*force[k]+pN[k]; }
  // update container quantities
  c->force+=sqrt(force[0]*force[0]+force[1]*force[1]+force[2]*force[2]);
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
    return 0;
  } else {
    // normal_angle = acos(normal_angle_cosine);
    // use the edge itself as a reference vector
    double refvec[3] = {v2->pN[0]-o2->pN[0],v2->pN[1]-o2->pN[1],v2->pN[2]-o2->pN[2]};
    // dot product of refvec and n1
    double d = dot(refvec,n1);
    if(!d){ return PI;}
    else  { return PI+d/fabs(d)*acos(normal_angle_cosine); }
    // I EXPECT 0 <= angle < 2*PI
  }
}

int Edge::computeFlip(void) {
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
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

/*double Edge::getForceEnergy(int i,double force[3],bool flag) {
// get polygon outward normals
double n[3];
// if i!=0, vertex requesting force is o2
if (i){f2->getNormal(n);}
// else i==0, vertex requesting force is o1
else  {f1->getNormal(n);}
// compute normal length
double L = sqrt(dot(n,n));
// compute angle between normals
double angle=getAngle();
// determine force direction
// when angle is less than PI, then angle_error is negative
// so force will be opposite direction of face normal
// when angle is greater than PI, then angle_error is positive
// so force will be in same direction as face normal
double angle_error = angle-PI;
int sign = 1;
if(angle_error<0){sign = -1;}
if(flag){
// force 
//		double force_magn = (ANGLE_STRETCH_WEIGHT/100.0)*angle_error/L;
double force_magn = (ANGLE_STRETCH_WEIGHT/100.0)*angle_error*angle_error/L;
force[0]+=sign*force_magn*n[0];
force[1]+=sign*force_magn*n[1];
force[2]+=sign*force_magn*n[2];
}
// energy
//	return (ANGLE_STRETCH_WEIGHT/100.0)/2.0*angle_error*angle_error;
return fabs((ANGLE_STRETCH_WEIGHT/100.0)/3.0*angle_error*angle_error*angle_error);
}*/

double getCurvatureLength(Vertex *o,Vertex *v1,Vertex *v2,double edge_length)
{
  double o_v1[3] = {
    o->pN[0]-v1->pN[0],
    o->pN[1]-v1->pN[1],
    o->pN[2]-v1->pN[2]
  };
  double v2_v1[3] = {
    v2->pN[0]-v1->pN[0],
    v2->pN[1]-v1->pN[1],
    v2->pN[2]-v1->pN[2]
  };
  double unum = dot(o_v1,v2_v1);
  double I[3] = {
    v1->pN[0]+unum*v2_v1[0]/edge_length/edge_length,
    v1->pN[1]+unum*v2_v1[1]/edge_length/edge_length,
    v1->pN[2]+unum*v2_v1[2]/edge_length/edge_length
  };
  double L = 
        (o->pN[0]-I[0])*(o->pN[0]-I[0])+
        (o->pN[1]-I[1])*(o->pN[1]-I[1])+
        (o->pN[2]-I[2])*(o->pN[2]-I[2]);
  return sqrt(L);
}

double computeCurvature(double L1,double L2,double angle,bool activate)
{
  double A[3] = {L1,0};
  double B[3] = {0,0};
  double C[3] = {L2*cos(angle),L2*sin(angle)};
  double a[3] = {A[0]-C[0],C[0]-B[0]+C[0]-A[0],B[0]-A[0]};
  double b[3] = {A[1]-C[1],C[1]-B[1]+C[1]-A[1],B[1]-A[1]};
  if(activate==true)
    //  if(true)
  {
    int N = 100;
    double max = -1;
    double maxt = -1;
    for(int i=0;i<N+1;i++)
    {
      double t = static_cast<double>(i)/static_cast<double>(N);
      double n1 = (3*a[0]*t*t+2*a[1]+a[2]);
      double n2 = (3*b[0]*t*t+2*b[1]+b[2]);
      double num = fabs(n1*(6*b[0]*t+2*b[1])-n2*(6*a[0]*t+2*a[1]));
      double d1 = n1*n1+n2*n2;
      double den = sqrt(d1*d1*d1);
      double doo = num/den;
      if(doo>max){max=doo;maxt=t;}
    }
    cout << "max curvature = " << max
          << ", t = " << maxt
          << ", L1 = " << L1
          << ", L2 = " << L2
          << endl;
  }

  double t = 1.0;
  //double t = L1/(L1+L2);
  double n1 = (3*a[0]*t*t+2*a[1]+a[2]);
  double n2 = (3*b[0]*t*t+2*b[1]+b[2]);
  double num = fabs(n1*(6*b[0]*t+2*b[1])-n2*(6*a[0]*t+2*a[1]));
  double d1 = n1*n1+n2*n2;
  double den = sqrt(d1*d1*d1);
  return num/den;
}

/*
    double Edge::getForceEnergy(int i,double force[3],bool flag,bool activate)
    {
// get polygon outward normals
double n[3];
// if i!=0, vertex requesting force is o2
if (i){f2->getNormal(n);}
// else i==0, vertex requesting force is o1
else  {f1->getNormal(n);}
// compute normal length
double L = sqrt(dot(n,n));
// compute angle between normals
double angle=getAngle();
// determine force direction
// when angle is less than PI, then angle_error is negative
// so force will be opposite direction of face normal
// when angle is greater than PI, then angle_error is positive
// so force will be in same direction as face normal
double angle_error = angle-PI;
int sign = 1;
if(angle_error<0){sign = -1;}
//// NEW FOR CURVATURE COMPUTATION
Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
getVertices(v1,v2,o1,o2);
// edge length
double edge_length = sqrt(getSqLength());
// get perpendicular lengths L1 and L2 from o1 and o2 vertices to edge
double L1 = getCurvatureLength(o1,v1,v2,edge_length);
double L2 = getCurvatureLength(o2,v1,v2,edge_length);
// compute curvature
double curvature = computeCurvature(L1,L2,angle,activate);
// debug
//cout << "angle_error = " << angle_error << endl;
//cout << "curvature = " << curvature << endl;
// debug
if(flag){
// force 
//  double force_magn = (ANGLE_STRETCH_WEIGHT/100.0)*angle_error*angle_error/L;
double force_magn;
if (i){ force_magn = CURVATURE_GAIN*(ANGLE_STRETCH_WEIGHT/100.0)*curvature*curvature*curvature/L2;}
else  { force_magn = CURVATURE_GAIN*(ANGLE_STRETCH_WEIGHT/100.0)*curvature*curvature*curvature/L1;}
double force_norm = force_magn/L;
//    if (activate==true)
//    {
//      cout << "curvature = " << curvature
//            << ", force_magn = " << force_magn
//            << ", angle = " << angle
//            << endl;
//      printEdge(f1->v[0]->o->name);
//    }

force[0]+=sign*force_norm*n[0];
force[1]+=sign*force_norm*n[1];
force[2]+=sign*force_norm*n[2];
}
// energy
//return fabs((ANGLE_STRETCH_WEIGHT/100.0)/3.0*angle_error*angle_error*angle_error);
//return fabs((ANGLE_STRETCH_WEIGHT/100.0)/2.0*curvature*curvature);
return fabs((ANGLE_STRETCH_WEIGHT/100.0)/4.0*curvature*curvature*curvature*curvature);
}*/


double Edge::getForceEnergy(int i,double force[3],bool flag,bool activate)
{
  // get polygon outward normals
  double n[3];
  // if i!=0, vertex requesting force is o2
  if (i){f2->getNormal(n);}
  // else i==0, vertex requesting force is o1
  else  {f1->getNormal(n);}
  // compute normal length
  double L = sqrt(dot(n,n));
  // compute angle between normals
  double angle=getAngle();
  // determine force direction
  // when angle is less than PI, then angle_error is negative
  // so force will be opposite direction of face normal
  // when angle is greater than PI, then angle_error is positive
  // so force will be in same direction as face normal
  double angle_error = angle-PI;
  int sign = 1;
  if(angle_error<0){sign = -1;}
  //// NEW FOR CURVATURE COMPUTATION
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
  // edge length
  double edge_length = sqrt(getSqLength());
  // get perpendicular lengths L1 and L2 from o1 and o2 vertices to edge
  double L1 = getCurvatureLength(o1,v1,v2,edge_length);
  double L2 = getCurvatureLength(o2,v1,v2,edge_length);
  // compute curvature
  //  double curvature = computeCurvature(L1,L2,angle,activate);
  // debug
  //cout << "angle_error = " << angle_error << endl;
  //cout << "curvature = " << curvature << endl;
  // debug
  double ae = fabs(angle_error/PI);
  //double force_metric = angle_error*angle_error*angle_error/(L1/2.0+L2/2.0);
  double force_metric = 1;
  for (int j=0;j<ANGLE_EXPONENT;j++) { force_metric *=ae; }
  force_metric = force_metric/(L1/2.0+L2/2.0);
  if(flag){
    // force 
    //  double force_magn = (ANGLE_STRETCH_WEIGHT/100.0)*angle_error*angle_error/L;
    double force_magn;
    if (i){ force_magn = ANGLE_GAIN*(ANGLE_STRETCH_WEIGHT/100.0)*force_metric/L2;}
    else  { force_magn = ANGLE_GAIN*(ANGLE_STRETCH_WEIGHT/100.0)*force_metric/L1;}
    double force_norm = force_magn/L;
    if (activate==true)
    {
      //      cout << "curvature = " << curvature
      //            << ", force_magn = " << force_magn
      //            << ", angle = " << angle
      //            << endl;
      //      printEdge(f1->v[0]->o->name);
    }

    force[0]+=sign*force_norm*n[0];
    force[1]+=sign*force_norm*n[1];
    force[2]+=sign*force_norm*n[2];
  }
  // energy
  //return fabs((ANGLE_STRETCH_WEIGHT/100.0)/3.0*angle_error*angle_error*angle_error);
  //return fabs((ANGLE_STRETCH_WEIGHT/100.0)/2.0*curvature*curvature);
  return fabs((ANGLE_STRETCH_WEIGHT/100.0)/ANGLE_EXPONENT*force_metric*ae);
}


double Edge::getReactionForceEnergy(double force[3],bool flag,bool activate) {
  // get polygon outward normals
  double n1[3],n2[3];
  f1->getNormal(n1);
  f2->getNormal(n2);
  // compute normal length
  double L1 = sqrt(dot(n1,n1));
  double L2 = sqrt(dot(n2,n2));
  // compute cosine of angle between normals
  double angle=getAngle();
  // determine force direction
  // when angle is less than PI, then angle_error is negative
  // so force will be opposite direction of face normal
  // when angle is greater than PI, then angle_error is positive
  // so force will be in same direction as face normal
  double angle_error = angle-PI;
  //// NEW FOR CURVATURE COMPUTATION
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
  // edge length
  //double edge_length = sqrt(getSqLength());
  // get perpendicular lengths L1 and L2 from o1 and o2 vertices to edge
  //double LL1 = getCurvatureLength(o1,v1,v2,edge_length);
  //double LL2 = getCurvatureLength(o2,v1,v2,edge_length);
  // compute curvature
  //  double curvature = computeCurvature(L1,L2,angle,activate);
  // debug
  //cout << "angle_error = " << angle_error << endl;
  //cout << "curvature = " << curvature << endl;
  // debug
  double ae = angle_error/PI;
  //double force_metric = angle_error*angle_error*angle_error/(L1/2.0+L2/2.0);
  double force_metric = 1;
  for (int j=0;j<ANGLE_EXPONENT;j++) { force_metric *=ae; }
  //force_metric = force_metric/(LL1/2.0+LL2/2.0);
  if(flag){
    // force 
    double force_magn1 = ANGLE_GAIN*(ANGLE_STRETCH_WEIGHT/100.0)*force_metric/L1;
    double force_magn2 = ANGLE_GAIN*(ANGLE_STRETCH_WEIGHT/100.0)*force_metric/L2;
    double force1[3]={force_magn1*n1[0],force_magn1*n1[1],force_magn1*n1[2]};
    double force2[3]={force_magn2*n2[0],force_magn2*n2[1],force_magn2*n2[2]};
    double reaction[3]={-(force1[0]+force2[0]),-(force1[1]+force2[1]),-(force1[2]+force2[2])};
    force[0]+=reaction[0]/2.0;
    force[1]+=reaction[1]/2.0;
    force[2]+=reaction[2]/2.0;
    if(activate==true)
    {
      cout << "Edge::getReactionForceEnergy: "
            << "angle_error (degrees) = " << angle_error*180/PI
            << ", ae = " << ae << endl
            << "Edge::getReactionForceEnergy: "
            << "force = ["
            << force[0] << " "
            << force[1] << " "
            << force[2] << "]";
      printEdge(f1->v[0]->o->name);
      cout << endl;
    }
  }
  // energy
  return (ANGLE_STRETCH_WEIGHT/100.0)/2.0*angle_error*angle_error;
}
/*
    double Edge::getReactionForceEnergy(double force[3],bool flag) {
// get polygon outward normals
double n1[3],n2[3];
f1->getNormal(n1);
f2->getNormal(n2);
// compute normal length
double L1 = sqrt(dot(n1,n1));
double L2 = sqrt(dot(n2,n2));
// compute cosine of angle between normals
double angle=getAngle();
// determine force direction
// when angle is less than PI, then angle_error is negative
// so force will be opposite direction of face normal
// when angle is greater than PI, then angle_error is positive
// so force will be in same direction as face normal
double angle_error = angle-PI;
if(flag){
// force 
double force_magn1 = (ANGLE_STRETCH_WEIGHT/100.0)*angle_error/L1;
double force_magn2 = (ANGLE_STRETCH_WEIGHT/100.0)*angle_error/L2;
double force1[3]={force_magn1*n1[0],force_magn1*n1[1],force_magn1*n1[2]};
double force2[3]={force_magn2*n2[0],force_magn2*n2[1],force_magn2*n2[2]};
double reaction[3]={-(force1[0]+force2[0]),-(force1[1]+force2[1]),-(force1[2]+force2[2])};
force[0]+=reaction[0]/2.0;
force[1]+=reaction[1]/2.0;
force[2]+=reaction[2]/2.0;
}
// energy
return (ANGLE_STRETCH_WEIGHT/100.0)/2.0*angle_error*angle_error;
}
*/
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

Edge* Object::findEdge(Vertex* va,Vertex* vb,hashtable_t &hm,int num_digits){
  Edge *ee=NULL;
  std::string s = keyPair(va->index,vb->index,num_digits);
  // if element exists given key, then get Edge pointer
  if (hm.count(s)>0){ ee=hm[s]; }
  return ee;
}

void Edge::getVertices(Vertex *&v1,Vertex *&v2,Vertex *&o1,Vertex *&o2){
  // find pair of vertices va and vb in common between f1 and f2
  Vertex *va=NULL,*vb=NULL;
  //if(f1->v[0]==f2->v[0] || f1->v[0]==f2->v[1] || f1->v[0]==f2->v[2]){va=f1->v[0];}
  //if(f1->v[1]==f2->v[0] || f1->v[1]==f2->v[1] || f1->v[1]==f2->v[2]){
  //  if(va==NULL){va=f1->v[1];}
  //  else {vb=f1->v[1];}
  //} 
  //if(f1->v[2]==f2->v[0] || f1->v[2]==f2->v[1] || f1->v[2]==f2->v[2]){
  //  vb=f1->v[2];
  //}

  if(f1->v[0]==f2->v[0] || f1->v[0]==f2->v[1] || f1->v[0]==f2->v[2])
  {
    va=f1->v[0];
    if(f1->v[1]==f2->v[0] || f1->v[1]==f2->v[1] || f1->v[1]==f2->v[2])
    {
      vb=f1->v[1];
    } 
    else
    {
      vb=f1->v[2];
    }
  }
  else
  {
    va=f1->v[1];
    vb=f1->v[2];
  }

  if(va==NULL || vb==NULL){
    cout << "\n\nEdge::getVertices: "
          << "common face vertices were not identified.\n\n";
    if(f1!=NULL){f1->printFace(f1->v[0]->o->name);cout << endl;}
    else{cout << "f1 is NULL.\n";}
    if(f2!=NULL){f2->printFace(f2->v[0]->o->name);cout << endl;}
    else{cout << "f2 is NULL.\n";}
    if(va!=NULL){cout << "va index = " << va->index << endl;}
    else{cout << "va is NULL.\n";}
    if(vb!=NULL){cout << "vb index = " << vb->index << endl;}
    else{cout << "vb is NULL.\n";}
    exit(0);
  }
  // identify v1 and v2 using f1
  if( (f1->v[0]==va && f1->v[1]==vb) || (f1->v[0]==vb && f1->v[1]==va)){
    v1=f1->v[0];
    v2=f1->v[1];
    o1=f1->v[2];
  } else if( (f1->v[1]==va && f1->v[2]==vb) || (f1->v[1]==vb && f1->v[2]==va)){
    v1=f1->v[1];
    v2=f1->v[2];
    o1=f1->v[0];
  } else if( (f1->v[2]==va && f1->v[0]==vb) || (f1->v[2]==vb && f1->v[0]==va)){
    v1=f1->v[2];
    v2=f1->v[0];
    o1=f1->v[1];
  } else {
    cout << "\n\nEdge::getVertices: "
          << "v1 and v2 were not successfully found.\n\n";
    exit(0);
  }
  // identify o2
  if( (f2->v[0]==va && f2->v[1]==vb) || (f2->v[0]==vb && f2->v[1]==va)){
    o2=f2->v[2];
  } else if( (f2->v[1]==va && f2->v[2]==vb) || (f2->v[1]==vb && f2->v[2]==va)){
    o2=f2->v[0];
  } else if( (f2->v[2]==va && f2->v[0]==vb) || (f2->v[2]==vb && f2->v[0]==va)){
    o2=f2->v[1];
  } else {
    cout << "\n\nEdge::getVertices: "
          << "o2 was not successfully found.\n\n";
    exit(0);
  }
}

/*
    void Edge::getVertices(Vertex *&v1,Vertex *&v2,Vertex *&o1,Vertex *&o2){
// find pair of vertices v1 and v2 in common between f1 and f2
if(f1->v[0]==f2->v[2]){
if(f1->v[1]==f2->v[1]){
v1=f1->v[0];
v2=f1->v[1];
o1=f1->v[2];
o2=f2->v[0];
} else {
v1=f1->v[2];
v2=f1->v[0];
o1=f1->v[1];
o2=f2->v[1];
}
} else if(f1->v[1]==f2->v[2]){
if(f1->v[2]==f2->v[1]){
v1=f1->v[1];
v2=f1->v[2];
o1=f1->v[0];
o2=f2->v[0];
} else {
v1=f1->v[0];
v2=f1->v[1];
o1=f1->v[2];
o2=f2->v[1];
}
} else if(f1->v[2]==f2->v[2]){
if(f1->v[0]==f2->v[1]){
v1=f1->v[2];
v2=f1->v[0];
o1=f1->v[1];
o2=f2->v[0];
} else {
v1=f1->v[1];
v2=f1->v[2];
o1=f1->v[0];
o2=f2->v[1];
}
} else if(f1->v[0]==f2->v[1]){
v1=f1->v[0];
v2=f1->v[1];
o1=f1->v[2];
o2=f2->v[2];
} else if(f1->v[1]==f2->v[1]){
v1=f1->v[1];
v2=f1->v[2];
o1=f1->v[0];
o2=f2->v[2];
} else {
v1=f1->v[2];
v2=f1->v[0];
o1=f1->v[1];
o2=f2->v[2];
}

}*/



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

void Object::createEdges(void) {
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
/*
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
*/
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
  std::vector<Edge*> e;
  getAdjacentEdges(e);
  // for each adjacent edge
  for (std::vector<Edge*>::iterator i=e.begin();i!=e.end();i++) {
    Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
    (*i)->getVertices(v1,v2,o1,o2);
    // find vertex different from self and add different vertex to vector
    if (v1!=this){a.push_back(v1);}
    else if (v2!=this) {a.push_back(v2);}
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

void printNeighborhood(std::string str,int i,std::vector<Face*> &nf,int iter){
  // for each face in neighborhood
  cout << "\n\nprintNeighborhood: iter " << iter
        << ", current vertex " << str << "->" << i << endl;
  for(std::vector<Face*>::iterator j=nf.begin();j!=nf.end();j++){
    (*j)->printFace((*j)->v[0]->o->name);
    cout << endl;
  }
}

void Vertex:: getAdjacentFaces(hashset_f &fset){
  std::vector<Edge*> e;
  getAdjacentEdges(e);
  // for each adjacent edge
  for (std::vector<Edge*>::iterator i=e.begin();i!=e.end();i++) {
    // add edge faces to set
    fset.insert((*i)->f1);
    fset.insert((*i)->f2);
  }
}

bool Object::processEdge(Edge *ee,hashtable_v_double &hood,std::vector<Edge*> &bucket,Vertex *vp){
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  ee->getVertices(v1,v2,o1,o2);
  bool empty = true;
  bool f1 =false,f2=false;
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
      hood[vp2]=sqrt( (vp2->pN[0]-vp->pN[0])*(vp2->pN[0]-vp->pN[0])+
                      (vp2->pN[1]-vp->pN[1])*(vp2->pN[1]-vp->pN[1])+
                      (vp2->pN[2]-vp->pN[2])*(vp2->pN[2]-vp->pN[2]));
    }
    if(f2){
      // if v2 found and v1 not found in hashtable
      // then add v1 to hashtable
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
    if (neighborhood[vv]<sqrt(NEIGHBORHOOD_RADIUS_SQ*SCALE*SCALE)){return false;}
    else {return true;}
  }
  return true;
}

bool Object::thawedAndAble(hashtable_v_double &hood,v_set &disabled){
  // for each vertex in hood
  for(vdhm_iterator i=hood.begin();i!=hood.end();i++){
    // if vertex is thawed and not disabled
    if(!ifFrozen(hood,(*i).first) && disabled.find((*i).first)==disabled.end()){
      // then affirm that at least one vertex is thawed and able
      return true;
    }
  }
  return false;
}
/*
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
}

///// all subsequent rounds /////

// while there are thawed vertices in neighbor list, hood,
// that are also not disabled
while(thawedAndAble(hood,disabled)){
// init collection of new faces
std::vector<Face*> new_faces_too;
collectFaces(hood,disabled,new_faces_too);

// for each face in new collection
for(std::vector<Face*>::iterator k=new_faces_too.begin();k!=new_faces_too.end();k++){
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
c.push_back(*k);
}
}
// copy local neighborhood to Object class neighborhood faces vector
(*i)->nf.assign(c.begin(),c.end());
// sort vectors
sort((*i)->nf.begin(),(*i)->nf.end());

}
}*/

bool edgeIntersectsHood(Vertex *v1,Vertex *v2,Vertex *v){
  double Ax = v1->pN[0];
  double Ay = v1->pN[1];
  double Az = v1->pN[2];
  double Bx = v2->pN[0];
  double By = v2->pN[1];
  double Bz = v2->pN[2];
  double Px = v->pN[0];
  double Py = v->pN[1];
  double Pz = v->pN[2];
  double uDen = (Bx-Ax)*(Bx-Ax)+(By-Ay)*(By-Ay)+(Bz-Az)*(Bz-Az);
  if(uDen!=0) {
    double u =((Px-Ax)*(Bx-Ax)+(Py-Ay)*(By-Ay)+(Pz-Az)*(Bz-Az))/uDen;
    if (u>0 && u<1) {
      double Ix=Ax+u*(Bx-Ax);
      double Iy=Ay+u*(By-Ay);
      double Iz=Az+u*(Bz-Az);
      // compute square of perpendicular distance from current vertex, v, to edge
      double d2 = (Px-Ix)*(Px-Ix)+(Py-Iy)*(Py-Iy)+(Pz-Iz)*(Pz-Iz);
      // if perpendicular distance is less than neighborhood radius
      if(d2<(NEIGHBORHOOD_RADIUS_SQ*SCALE*SCALE)){
        // assuming this correctly implies that the edge
        // intersects the sphere of neighborhood
        return true;
      } else {
        return false;
      }
    } else {
      // perpendicular intersection occurs outside edge vertices
      return false;
    }
  } else {
    cout << "Error. Edge vertices are indistinguishable.\n";exit(0);
  }
}

/*
    bool findNewFaces(f_set &new_faces,e_set &boundary,Vertex *v){
    new_faces.clear();
    double NR2 = NEIGHBORHOOD_RADIUS_SQ;
    bool found = false;
// for each edge in boundary
for(es_iterator i=boundary.begin();i!=boundary.end();i++){
Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
(*i)->getVertices(v1,v2,o1,o2);
// if either v1 or v2 inside NR
if( ( ( (v->pN[0]-v1->pN[0])*(v->pN[0]-v1->pN[0])+
(v->pN[1]-v1->pN[1])*(v->pN[1]-v1->pN[1])+
(v->pN[2]-v1->pN[2])*(v->pN[2]-v1->pN[2])
) < NR2) ||
( ( (v->pN[0]-v2->pN[0])*(v->pN[0]-v2->pN[0])+
(v->pN[1]-v2->pN[1])*(v->pN[1]-v2->pN[1])+
(v->pN[2]-v2->pN[2])*(v->pN[2]-v2->pN[2])
) < NR2)
){
// then add f1 and f2 to new_faces
new_faces.insert((*i)->f1);
new_faces.insert((*i)->f2);
found = true;
}
}

// if no faces were added
if(found==false){
// for each edge in boundary
for(es_iterator i=boundary.begin();i!=boundary.end();i++){
Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
(*i)->getVertices(v1,v2,o1,o2);
// if edge intersects sphere of neighborhood
if(edgeIntersectsHood(v1,v2,v)){
// then add f1 and f2 to new_faces
new_faces.insert((*i)->f1);
new_faces.insert((*i)->f2);
found = true;
}
}
}
return found;
}*/

bool findNewFaces(f_set &new_faces,e_set &boundary,Vertex *v){
  new_faces.clear();
  //	double NR2 = NEIGHBORHOOD_RADIUS_SQ;
  bool found = false;
  // for each edge in boundary
  for(es_iterator i=boundary.begin();i!=boundary.end();i++){
    Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
    (*i)->getVertices(v1,v2,o1,o2);
    double a = v->pN[0]-v1->pN[0];
    double b = v->pN[1]-v1->pN[1];
    double c = v->pN[2]-v1->pN[2];
    double d = v->pN[0]-v2->pN[0];
    double e = v->pN[1]-v2->pN[1];
    double f = v->pN[2]-v2->pN[2];
    // if either v1 or v2 inside NR
    if( ((a*a+b*b+c*c)<(NEIGHBORHOOD_RADIUS_SQ*SCALE*SCALE)) ||
        ((d*d+e*e+f*f)<(NEIGHBORHOOD_RADIUS_SQ*SCALE*SCALE))){
      // then add f1 and f2 to new_faces
      new_faces.insert((*i)->f1);
      new_faces.insert((*i)->f2);
      found = true;
    }
  }

  // if no faces were added
  if(found==false){
    // for each edge in boundary
    for(es_iterator i=boundary.begin();i!=boundary.end();i++){
      Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
      (*i)->getVertices(v1,v2,o1,o2);
      // if edge intersects sphere of neighborhood
      if(edgeIntersectsHood(v1,v2,v)){
        // then add f1 and f2 to new_faces
        new_faces.insert((*i)->f1);
        new_faces.insert((*i)->f2);
        found = true;
      }
    }
  }
  return found;
}
/*
    void updateBoundary(f_set &new_faces,e_set &boundary,f_set &hood){
// for each new_face
for(fs_iterator i=new_faces.begin();i!=new_faces.end();i++){
// for each face edge
for(int j=0;j<3;j++){
// if f1 not in hood or f2 not in hood
if( hood.find((*i)->e[j]->f1)==hood.end() ||
hood.find((*i)->e[j]->f2)==hood.end()){
// add edge to boundary
boundary.insert((*i)->e[j]);
}
}
}

// for each edge in boundary
for(es_iterator i=boundary.begin();i!=boundary.end();i++){
// if f1 and f2 in hood
if( (hood.find((*i)->f1)!=hood.end()) &&
(hood.find((*i)->f2)!=hood.end()) ){
// remove edge from boundary
boundary.erase(i);
}
}
}*/
/* OCT 2007
   void updateBoundary(f_set &new_faces,e_set &boundary,f_set &hood){
// for each edge in boundary
for(es_iterator i=boundary.begin();i!=boundary.end();i++){
// if f1 and f2 in hood
if( hood.count((*i)->f1) && hood.count((*i)->f2) ){
// remove edge from boundary
boundary.erase(i);
}
}
std::vector<Edge*> e;
// for each new_face
for(fs_iterator i=new_faces.begin();i!=new_faces.end();i++){
// grab each face edge
e.push_back((*i)->e[0]);
e.push_back((*i)->e[1]);
e.push_back((*i)->e[2]);
}
// keep unique edges
sort(e.begin(),e.end());
std::vector<Edge*>::iterator new_end = unique(e.begin(),e.end());
e.assign(e.begin(),new_end);
// for each face edge
for(std::vector<Edge*>::iterator i=e.begin();i!=e.end();i++){
// if f1 not in hood or f2 not in hood
if( !hood.count((*i)->f1) || !hood.count((*i)->f2)){
// add edge to boundary
boundary.insert(*i);
}
}
}*/

void updateBoundary(f_set &new_faces,e_set &boundary,f_set &hood){
  // for each new_face
  for(fs_iterator i=new_faces.begin();i!=new_faces.end();i++){
    // if face not in hood
    if(hood.count(*i)==0){
      // for each edge of face
      for(int j=0;j<3;j++){
        Edge *ee=(*i)->e[j];
        // look for edge in boundary
        es_iterator esi = find(boundary.begin(),boundary.end(),ee);
        // if edge found in boundary, then remove edge from boundary
        if(esi!=boundary.end()){ boundary.erase(esi); }
        // else add edge to boundary
        else { boundary.insert(ee);}
      }
    }
  }
}

/*
    void Object::newFindNeighborhoods(void){
    f_set hood,new_faces;
    e_set boundary;
// for each vertex in object
for (std::vector<Vertex*>::iterator i=v.begin();i!=v.end();i++) {
hood.clear();
new_faces.clear();
boundary.clear();
// init new faces to current vertex adjacent faces
new_faces.insert((*i)->f.begin(),(*i)->f.end());
// add new faces to hood
hood.insert(new_faces.begin(),new_faces.end());
// update hood edge boundary
updateBoundary(new_faces,boundary,hood);
// iterate until all edge vertices are outside neighborhood radius
bool found = true;
while(found){
// collect faces to add to hood
found = findNewFaces(new_faces,boundary,*i);
// add new faces to hood
hood.insert(new_faces.begin(),new_faces.end());
// update hood edge boudary
updateBoundary(new_faces,boundary,hood);
}
// copy local neighborhood to Object class neighborhood faces vector
(*i)->nf.assign(hood.begin(),hood.end());
// sort vectors
sort((*i)->nf.begin(),(*i)->nf.end());
}
}*/

void Object::buildNeighborhood(Vertex *vv){
  f_set hood,new_faces;
  e_set boundary;
  // init new faces to current vertex adjacent faces
  new_faces.insert(vv->f.begin(),vv->f.end());
  // add new faces to hood
  hood.insert(new_faces.begin(),new_faces.end());
  // update hood edge boundary
  updateBoundary(new_faces,boundary,hood);
  // iterate until all edge vertices are outside neighborhood radius
  bool found = true;
  while(found){
    // collect faces to add to hood
    found = findNewFaces(new_faces,boundary,vv);
    // add new faces to hood
    hood.insert(new_faces.begin(),new_faces.end());
    // update hood edge boudary
    updateBoundary(new_faces,boundary,hood);
  }
  // copy local neighborhood to Object class neighborhood faces vector
  vv->nf.assign(hood.begin(),hood.end());
  // sort vectors
  sort(vv->nf.begin(),vv->nf.end());
}

void Object::newFindNeighborhoods(void){
  // for each vertex in object
  for (std::vector<Vertex*>::iterator i=v.begin();i!=v.end();i++) {
    buildNeighborhood(*i);
  }
}


/*
    void Object::newFindNeighborhoods(void){
// for each vertex in object
for (std::vector<Vertex*>::iterator i=v.begin();i!=v.end();i++) {
//		buildNeighborhood(*i);
f_set hood,new_faces;
e_set boundary;
// init new faces to current vertex adjacent faces
new_faces.insert((*i)->f.begin(),(*i)->f.end());
// add new faces to hood
hood.insert(new_faces.begin(),new_faces.end());
// update hood edge boundary
updateBoundary(new_faces,boundary,hood);
// iterate until all edge vertices are outside neighborhood radius
bool found = true;
while(found){
// collect faces to add to hood
found = findNewFaces(new_faces,boundary,(*i));
// add new faces to hood
hood.insert(new_faces.begin(),new_faces.end());
// update hood edge boudary
updateBoundary(new_faces,boundary,hood);
}
// copy local neighborhood to Object class neighborhood faces vector
(*i)->nf.assign(hood.begin(),hood.end());
// sort vectors
sort((*i)->nf.begin(),(*i)->nf.end());
}
}*/

void Object::findNeighborhoods(void){
  vector<Vertex*> newverts,adj,bank;
  std::vector<Vertex*>::iterator i,k,m;
  Bit b;
  hashset_f fset;
  ////////// compute neighborhood radius in edge lengths //////////
  int radius = (int) ceil(sqrt(NEIGHBORHOOD_RADIUS_SQ*SCALE*SCALE)/getMeanEdgeLength());
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
/*
    void Object::findVertexAdjacencies(void){
// for each face, add face* to each face vertex
for (std::vector<Face*>::iterator j=f.begin();j!=f.end();j++){
((*j)->v[0])->f.push_back(*j);
((*j)->v[1])->f.push_back(*j);
((*j)->v[2])->f.push_back(*j);
}
// find neighborhoods
newFindNeighborhoods();
}*/

void Object::findVertexAdjacencies(void){
  //	int index=1;
  //	double max=f.size();
  //	double goal = 0.2;
  //	printf("0%%..");
  // for each face in object
  for (std::vector<Face*>::iterator j=f.begin();j!=f.end();j++){
    // for each vertex of face
    for(int i=0;i<3;i++){
      Vertex *vv=(*j)->v[i];
      // add face* to vertex
      vv->f.push_back(*j);
      // build neighborhood
      if(vv->nf.empty()==true){buildNeighborhood(vv);}
    }
    // track progress
    //		cout << "\nindex " << index << endl;
    //		double progress = static_cast<double>(index++)/max;
    //		if(progress>goal){
    //			printf("%d%%..",static_cast<int>(goal*100));
    //			goal+=0.2;
    //		}
    //		if(index>800){break;}
  }
  //	printf("100%%..");
  // find neighborhoods
  //	newFindNeighborhoods();
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
  char file[FILENAME_SIZE];
  // open file
  sprintf(file,"%s%s",OUTPUT_DATA_DIR.c_str(),OBJECT_LIST_FILE);
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
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
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
  char file[FILENAME_SIZE];
  sprintf(file,"%s%s",OUTPUT_DATA_DIR.c_str(),CONT_LOG_FILE);
  Cfile.open(file);
  Cfile << "\nLEGEND_______________________\n"
        << "'Iteration'\n"
        << "0th group is the original data.\n"
        << "For 1st group vertices have moved once or not at all.\n"
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
        << "Mean displacement of all vertices in all objects chosen to move during the group.\n"
        << "--\n"
        << "'Min. Edge Angle (rad.)'\n"
        << "The minimum edge angle of all edges in all objects.\n"
        << "Minimum edge angle is the angle between the two faces that share the edge.\n"
        << "___________________________________________\n\n\n\n";

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
  Cfile.width(10);
  Cfile << left << "Max Gain";
  Cfile.width(16);
  Cfile << left << "Energy";
  Cfile.width(11);
  Cfile << left << "Disp_Mean";
  Cfile.width(15);
  Cfile << right << "Disp_Min";
  Cfile.width(1);
  Cfile << ":";
  Cfile.width(15);
  Cfile << left << "Disp_Max";
  Cfile.width(19);
  Cfile << left << "MinEdgeAngle(deg)";
  Cfile.width(13);
  Cfile << left << "Duration(sec)" << endl;
  Cfile.flush();
}

void Container::sepFileInit(const int group) {
  char file[FILENAME_SIZE];
  sprintf(file,"%s%s.%d",OUTPUT_DATA_DIR.c_str(),SEP_LOG_FILE,group);
  Sfile.open(file);
  Sfile.width(22);
  Sfile << left << "#separation_distance";
  Sfile.width(22);
  Sfile << left << "object_name";
  Sfile.width(22);
  Sfile << left << "vertex_index";
  Sfile.width(22);
  Sfile << left << "virtual_displacement";
  Sfile << endl;
  Sfile.flush();
}

void Container::sepOrigLog(Monitor &stats) {
  char file[FILENAME_SIZE];
  sprintf(file,"%s%s.orig",OUTPUT_DATA_DIR.c_str(),SEP_LOG_FILE);
  std::ofstream Tfile;
  Tfile.open(file);
  Tfile.width(22);
  Tfile << left << "#separation_distance";
  Tfile.width(22);
  Tfile << left << "object_name";
  Tfile.width(22);
  Tfile << left << "vertex_index";
  Tfile.width(22);
  Tfile << left << "virtual_displacement";
  Tfile << endl;
  Tfile.flush();

  Tfile.precision(10);
  // for each object in container
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    // for each vertex in object
    for (std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
      Tfile.width(22);
      if((*j)->cl!=NULL){
        Tfile << left << sqrt((*j)->getSqSepDist(this));
      } else {
        Tfile << left << -1.0;
      }
      Tfile.width(22);
      Tfile << left << (*j)->o->name;
      Tfile.width(22);
      Tfile << left << (*j)->index;
      Tfile.width(22);
      // if vertex* found in old
      if(stats.old.find(*j)!=stats.old.end()){
        Tfile << left << sqrt(stats.old[*j]);
      } else {
        Tfile << left << -1.0;
      }
      Tfile << endl;
      Tfile.flush();
    }
  }
  Tfile.close();
}

void Container::sepSnapshotLog(Monitor &stats,const int group) {
  char file[FILENAME_SIZE];
  sprintf(file,"%s%s.snapshot.%d",OUTPUT_DATA_DIR.c_str(),SEP_LOG_FILE,group);
  std::ofstream Tfile;
  Tfile.open(file);
  Tfile.width(22);
  Tfile << left << "#separation_distance";
  Tfile.width(22);
  Tfile << left << "object_name";
  Tfile.width(22);
  Tfile << left << "vertex_index";
  Tfile.width(22);
  Tfile << left << "virtual_displacement";
  Tfile << endl;
  Tfile.flush();

  Tfile.precision(10);
  // for each object in container
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    // for each vertex in object
    for (std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
      Tfile.width(22);
      if((*j)->cl!=NULL){
        Tfile << left << sqrt((*j)->getSqSepDist(this));
      } else {
        Tfile << left << -1.0;
      }
      Tfile.width(22);
      Tfile << left << (*j)->o->name;
      Tfile.width(22);
      Tfile << left << (*j)->index;
      Tfile.width(22);
      // if vertex* found in old
      if(stats.old.find(*j)!=stats.old.end()){
        Tfile << left << sqrt(stats.old[*j]);
      } else {
        Tfile << left << -1.0;
      }
      Tfile << endl;
      Tfile.flush();
    }
  }
  Tfile.close();
}

// #####################################################
// #####################################################

void Container::fileOutit(void) {
  Cfile.close();
}

void Container::sepFileOutit(void) {
  Sfile.close();
}

// #####################################################
// #####################################################

void Container::reportVertexNiceness(int index,std::string str,Space &s) {
  // determine current niceness
  std::vector<Object*> cb;
  // for each object* in container
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    // if object names match
    if(str==(*i)->name){
      // for each vertex in object
      for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++){
        if((*j)->index==index){
          // reevaluate vertex niceness with ray trace
          collectCrossed(s,*j,cb);
          break;
        }
      }
    }
  }
  if(cb.empty()==false){
    cout << "\nContainer::reportVertexNiceness: "
          << "vertex "
          << str << "->" << index
          << ": is not nice in ray trace.\n";
  } else {
    cout << "\nContainer::reportVertexNiceness: "
          << "vertex "
          << str << "->" << index
          << ": is nice in ray trace.\n";
  }

  // report stored niceness
  // for each object* in container
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    // if object names match
    if(str==(*i)->name){
      // for each hashtable element (Vertex*->int)
      for (vhm_iterator j=(*i)->nice.begin();j!=(*i)->nice.end();j++) {
        Vertex *vv=(*j).first;
        // if vertex indices match
        if(vv->index==index){
          cout << "\nContainer::reportVertexNiceness: "
                << "vertex "
                << vv->o->name << "->" << vv->index
                << ": is not nice in hash_map(val=" << (*j).second << ")\n";
          return;
        }
      }
      cout << "\nContainer::reportVertexNiceness: "
            << "vertex "
            << str << "->" << index
            << " is nice.\n";
      return;
    }
  }
}

void Container::countNonnice(int &total,int &self) {
  total=self=0;
  // for each object* in container
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    // for each hashtable element (Vertex*->int)
    for (vhm_iterator j=(*i)->nice.begin();j!=(*i)->nice.end();j++) {
      // if vertex is self_nonnice
      if((*j).second==2){self++;}
      total++;
    }
  }
}

void Container::countIntersections(int &total,int &self) {
  total=self=0;
  // for each object* in container
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    // for each hashtable element (Face*->vector<Face*>*)
    for (htff_iterator j=(*i)->intf.begin();j!=(*i)->intf.end();j++) {
      // grab Face *key
      Face *ff=(*j).first;
      // for each Face* in vector
      std::vector<Face*> *fv=(*j).second;
      for(std::vector<Face*>::iterator k=fv->begin();k!=fv->end();k++){
        // if faces of same object
        if((*k)->v[0]->o==ff->v[0]->o){self++;}
        total++;
      }
    }
  }
}

void Container::gatherDiagnostics(void) {
  cout.precision(10);
  int v_total,v_self;
  countNonnice(v_total,v_self);
  int f_total,f_self;
  countIntersections(f_total,f_self);
  cout << "\n\nContainer::updateFile: "
        << "[mycount,container]--"
        << "total_nonnice[" << nonnice << "," << v_total << "]"
        << ", self_nonnice[" << s_nonnice << "," << v_self << "]\n"
        << "[mycount,container]--"
        << "total_intface[" << ti/2 << "," << f_total << "]"
        << ", self_intface[" << si/2 << "," << f_self << "]\n";
  cout.flush();
}

//void Container::updateFile(const int group,bool flag,double tim,double max_gain) {
void Container::updateFile(const int group,bool flag,double tim) {
  if(group==0){
    cout << "Update log files...............................";
  } else {
    cout << "Iteration " << group << ": ";
    cout << "Update log files..................";
  }
  cout.flush();
  Cfile.precision(10);
  Cfile.width(4);
  Cfile << left << group;
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
  Cfile.width(10);
  Cfile << left << "NA";
  Cfile.width(16);
  Cfile << left << energy;
  if (!flag){
    Cfile.width(11);
    Cfile << left << "NA";
    Cfile.width(15);
    Cfile << right << "NA";
    Cfile.width(1);
    Cfile << ":";
    Cfile.width(15);
    Cfile << left << "NA";
  } else {
    Cfile.width(11);
    Cfile << left << md[1];
    Cfile.width(15);
    char name[1024];
    sprintf(name,"%.6g",d_min);
    Cfile << right << name;
    Cfile.width(1);
    Cfile << ":";
    Cfile.width(15);
    sprintf(name,"%.6g",d_max);
    Cfile << left << name;
  }
  Cfile.precision(10);
  Cfile.width(19);
  Cfile << left << min_edge_angle*180.0/PI;
  Cfile.width(13);
  Cfile.precision(4);
  Cfile << left << tim << endl;
  Cfile.flush();
  cout << "complete.\n";
  cout.flush();
}

void Container::updateSepFile(double sep,std::string name,int index,double vd) {
  Sfile.precision(10);
  Sfile.width(22);
  Sfile << left << sep;
  Sfile.width(22);
  Sfile << left << name;
  Sfile.width(22);
  Sfile << left << index;
  Sfile.width(22);
  Sfile << left << vd;
  Sfile << endl;
  Sfile.flush();
}

// #####################################################
// #####################################################

void Container::createEdges(void) {
  if(!WRITE_VERBOSE_INIT){
    cout << "Creating edges.................................";
    cout.flush();
  }
  int iii=1;
  double goal = 0.2;
  printf("0%%..");
  fflush(stdout);
  // for each object, create edges
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    if(WRITE_VERBOSE_INIT){
      cout << "Creating edges for " << (*i)->name << "...";
      cout.flush();
    }
    //		(*i)->createEdges();
    // determine number of digits in largest vertex index
    int num_digits = (*i)->setNumDigits();
    // create map for finding edges
    hashtable_t hm;
    // for each face
    for (std::vector<Face*>::iterator j=(*i)->f.begin();j!=(*i)->f.end();j++) {
      (*i)->checkEdge(*j,(*j)->v[0],(*j)->v[1],hm,num_digits);
      (*i)->checkEdge(*j,(*j)->v[1],(*j)->v[2],hm,num_digits);
      (*i)->checkEdge(*j,(*j)->v[2],(*j)->v[0],hm,num_digits);
      // track progress
      double progress = static_cast<double>(iii++)/face_count;
      if(progress>goal){
        printf("%d%%..",static_cast<int>(goal*100));
        fflush(stdout);
        goal+=0.2;
      }
    }
    if(WRITE_VERBOSE_INIT){
      cout << "complete.\n";
      cout.flush();
    }
  }
  printf("100%%..");
  fflush(stdout);
  if(!WRITE_VERBOSE_INIT){
    cout << "complete.\n";
    cout.flush();
  }
}

void Container::findVertexAdjacencies(void) {
  if(!WRITE_VERBOSE_INIT){
    cout << "Finding vertex adjacencies.....................";
    cout.flush();
  }
  std::vector<Object*>::iterator i;
  // if 
  if (WRITE_NEIGHBORHOOD_TO_FILE){
    // open dat file
    char log_file[FILENAME_SIZE];
    sprintf(log_file,"%s%s",OUTPUT_DATA_DIR.c_str(),NEIGHBORHOOD_FILE);
    std::ofstream this_file (log_file);
    this_file.close();
  }
  // for each object, find vertex adjacencies
  if(WRITE_VERBOSE_INIT){ cout << endl;} 
  for (i=o.begin();i!=o.end();i++) {
    if(WRITE_VERBOSE_INIT){
      cout << "Finding adjacencies for " << (*i)->name << "...";
      cout.flush();
    }
    (*i)->findVertexAdjacencies();
    if(WRITE_VERBOSE_INIT){
      cout << "complete.\n";
      cout.flush();
    }
  }
  if(!WRITE_VERBOSE_INIT){
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
    // for each edge of face
    for(int j=0;j<3;j++){
      Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
      (*i)->e[j]->getVertices(v1,v2,o1,o2);
      // if edge does not refer to this vertex, then add edge to vector
      if(v1!=this && v2!=this){ee.push_back((*i)->e[j]);}
    }
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
  if(WRITE_VERBOSE_INIT){ cout << endl;cout.flush();}
  pdir = opendir(INPUT_DATA_DIR.c_str());
  if (!pdir) {printf("Error. Could not open %s.\n",INPUT_DATA_DIR.c_str());exit(1);}
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
      if(WRITE_VERBOSE_INIT){
        cout << "file found: " << str << "\n"; cout.flush();
      }
    }
  }
  closedir(pdir);
  if(WRITE_VERBOSE_INIT){ cout << endl;cout.flush();}

}
// #####################################################
// #####################################################

void Container::scanFiles(void) {
  std::string str;
  std::string::size_type pos1;
  Object *obj;
  char file[FILENAME_SIZE];
  //	cout << "Input Data Directory..........................."
  cout << "\nInput Data Directory\n"
        << INPUT_DATA_DIR << "\n\n"
        << "Reading Directory..............................";
  cout.flush();
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
    sprintf(file,"%s%s",INPUT_DATA_DIR.c_str(),files[count].c_str());
    if(WRITE_VERBOSE_INIT){
      cout << "loading file " << files[count] << "...";
    }
    cout.flush();
    scanFile(obj,file);
    if(WRITE_VERBOSE_INIT){
      cout << "complete.\n";
      cout.flush();
    }
  }
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
      if (PRINT_FLAG) { cout.precision(15);cout<<filename<<": Vertex "<<v->index 
        <<" "<<v->pN[0]<<" "<<v->pN[1]<<" "<<v->pN[2]<<"\n";
      }
    }
    // if first character is F for Face, add new linked list class instance
    else if (strchr("F",*str)!=NULL){
      f=new Face(str,vp);
      polygon_num++;
      obj->f.push_back(f);
      if (PRINT_FLAG) { cout.precision(15);cout<<filename<<": Face "<<f->index 
        <<" "<<f->v[0]->index<<" "<<f->v[1]->index<<" "<<f->v[2]->index<<"\n";
      }
    }
  }
  fclose(F);
}


// #####################################################
// #####################################################

Object* Container::getObjectPointer(char val[FILENAME_SIZE]){
  std::string str=val;
  // for each object* in container
  for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++){
    if(str==(*i)->name){return *i;}
  }
  return NULL;
}

void Container::loadFrozenMap(mmap_oi &frozen_map,const char *filename){
  // open file
  FILE *F;
  F = fopen(filename,"r");
  if (!F) { printf("Couldn't open input file %s\n",filename);return;}

  // for every line in file
  char val[FILENAME_SIZE];
  char *eptr;
  int i;
  Object *o1=NULL,*o2=NULL;
  char line[2048],*str;
  for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) {

    // get object name
    while (strchr(" \t,",*str)!=NULL) { str++; }
    i=0;
    while (strchr(" \t",*str)==NULL)
    {
      val[i++] = *str++;
    }
    val[i]=0;
    // val contains object name

    // get Object*
    if(o2==NULL){o1 = getObjectPointer(val);}
    else {
      // if new object name same as old object name
      if(strcmp(o2->name.c_str(),val)==false){o1=o2;}
      else {o1 = getObjectPointer(val);}
    }
    o2=o1;
    // Error
    if(o1==NULL){
      cout << "\n\nContainer::loadFrozenMap: Error.\n"
            << "No matching Object* found in container for frozen vertex:\n"
            << val << endl
            << "line from " << filename << ":\n"
            << str << endl;
      exit(0);
    }

    // get vertex index
    while (strchr(" \t,",*str)!=NULL) { str++; }
    i=0;
    while (strchr("0123456789",*str)!=NULL)
    {
      val[i++] = *str++;
    }
    val[i]=0;
    int vi = (int)strtod(val,&eptr);
    if (val==eptr)
    {
      vi=0;
      printf("Error in reading vertex index\n");
      return;
    }

    // create map entry
    frozen_map.insert(std::make_pair(o1,vi));
    //		cout << "Frozen vertex object = "
    //		<< oo->name << "->" << vi << endl;

  }

  // close file
  fclose(F);
}
/*
    void Container::readFrozen (char *filename) {
    cout << "Read frozen vertices...........................";
    cout.flush();
// load frozen map: object*->Vertex_index (int)
mmap_oi frozen_map;
oi_iterator front;
loadFrozenMap(frozen_map,filename);

// sort object* vector in container
sort(o.begin(),o.end());

// for each object* in container
std::vector<Object*>::iterator i=o.begin();
while(i!=o.end()){
// if map has elements
if(frozen_map.empty()==false){
// if object* == to front map element
front=frozen_map.begin();
if(*i==(*front).first){
// for each entry in object vertex* vector
for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++){
// if vertex*->index matches vertex index in front map element
if((*j)->index==(*front).second){
// then add vertex* to frozen vector
frozen.push_back(*j);
// erase front map element
frozen_map.erase(frozen_map.begin());
break;
}
}
}
// else increment iterator to next object* in container	
else {i++;}
}
// else all map elements have been processed	
else {break;}
}

// sort frozen vector
sort(frozen.begin(),frozen.end());

cout << "complete.\n";
cout.flush();
}
*/
void Container::readFrozen (const char *filename) {
  cout << "Read frozen vertices...........................";
  cout.flush();
  // load frozen map: object*->Vertex_index (int)
  mmap_oi frozen_map;
  oi_iterator front;
  loadFrozenMap(frozen_map,filename);

  // for each element in multimap
  for(oi_iterator i=frozen_map.begin();i!=frozen_map.end();i++){
    Object *oo=(*i).first;
    int t = (*i).second;
    int a = t-1;
    int b = oo->v[a]->index;
    do{
      if		(t==b){frozen.push_back(oo->v[a]);}
      else if (t<b){a++;}
      else if (t>b){a--;}
    } while(t!=b);
  }
  sort(frozen.begin(),frozen.end());
  cout << "complete.\n";
  cout.flush();
}

// #####################################################
// #####################################################

void Container::buildMeshAfter(const int group) {
  /*	// DEBUG
        cout << "\nContainer::buildMeshAfter: begin: group="
        << group << endl;
        cout.flush();
        */	// DEBUG
  // for each object
  for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    //          if((*i)->name!="a066_clipped"){continue;}
    /*		// DEBUG
                cout << "\nContainer::buildMeshAfter: object="
                << (*i)->name << ", group=" << group << endl;
                cout.flush();
                */		// DEBUG
    // create output filename
    char file[FILENAME_SIZE];
    if (APPEND_ITERATION_NUMBER_TO_MESH) {
      sprintf(file,"%s%s%s_%i.mesh",OUTPUT_DATA_DIR.c_str(),(*i)->name.c_str(),OUTPUT_SUFFIX,group);
    } else {
      sprintf(file,"%s%s%s.mesh",OUTPUT_DATA_DIR.c_str(),(*i)->name.c_str(),OUTPUT_SUFFIX);
    }
    // open output file
    std::ofstream newfile;
    newfile.open (file,std::ios::trunc);
    if(newfile.is_open()){
      newfile.precision(12);
      // for each vertex in object
      for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
        // print index and final coordinates
        newfile << "Vertex "
              << (*j)->index << " "
              << (*j)->pN[0] << " "
              << (*j)->pN[1] << " "
              << (*j)->pN[2] << "\n";
      }
      // for each face in object
      for(std::vector<Face*>::iterator k=(*i)->f.begin();k!=(*i)->f.end();k++) {
        newfile << "Face "
              << (*k)->index << " "
              << (*k)->v[0]->index << " "
              << (*k)->v[1]->index << " "
              << (*k)->v[2]->index << "\n";
      }
      newfile.close();
    } else {
      cout << "\nContainer::buildMeshAfter: Error:"
            << "Failed to open output file="
            << file << endl;
      cout.flush();
    }
  }
  /*	// DEBUG
        cout << "\nContainer::buildMeshAfter: after: group="
        << group << endl;
        cout.flush();
        */	// DEBUG
}

void Container::writeDistancesNOCP(const int group) {
  char file[FILENAME_SIZE];
  // create output filename
  if (APPEND_ITERATION_NUMBER_TO_DISTANCES) {
    sprintf(file,"%sclosest_point_distances_NOCP_%i.dat",OUTPUT_DATA_DIR.c_str(),group);
  } else {
    sprintf(file,"%sclosest_point_distances_NOCP.dat",OUTPUT_DATA_DIR.c_str());
  }
  // open output file
  std::ofstream newfile;
  newfile.open (file,std::ios::trunc);
  if(newfile.is_open()){
    newfile.precision(9);
    // for each object
    for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
      // for each vertex in object
      for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
        // if vertex has a closest face
        if((*j)->cl==NULL){
          // compute separation vector
          // print separation distance
          // dreamm custom points
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

void Container::writeDistances50(const int group) {
  char file[FILENAME_SIZE];
  // create output filename
  if (APPEND_ITERATION_NUMBER_TO_DISTANCES) {
    sprintf(file,"%sclosest_point_distances_%i.dat",OUTPUT_DATA_DIR.c_str(),group);
  } else {
    sprintf(file,"%sclosest_point_distances.dat",OUTPUT_DATA_DIR.c_str());
  }
  // open output file
  std::ofstream newfile;
  newfile.open (file,std::ios::trunc);
  if(newfile.is_open()){
    newfile.precision(4);
    // for each object
    for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
      // for each vertex in object
      for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
        // if vertex has a closest face
        if((*j)->cl!=NULL){
          // get closest point
          double pC[3];
          computePC((*j)->cl,*j,pC);
          // compute separation vector
          double s[3];
          for(int k=0;k<3;k++){ s[k]=pC[k]-(*j)->pN[k]; }
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

void Container::writeDistances(const int group) {
  char file[FILENAME_SIZE];
  // create output filename
  if (APPEND_ITERATION_NUMBER_TO_DISTANCES) {
    sprintf(file,"%sclosest_point_distances_%i.dat",OUTPUT_DATA_DIR.c_str(),group);
  } else {
    sprintf(file,"%sclosest_point_distances.dat",OUTPUT_DATA_DIR.c_str());
  }
  // open output file
  std::ofstream newfile;
  newfile.open (file,std::ios::trunc);
  if(newfile.is_open()){
    newfile.precision(4);
    // for each object
    for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
      // for each vertex in object
      for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
        // if vertex has a closest face
        if((*j)->cl!=NULL && (binary_search(frozen.begin(),frozen.end(),*j)==false)){
          // get closest point
          double pC[3];
          computePC((*j)->cl,*j,pC);
          // compute separation vector
          double s[3];
          for(int k=0;k<3;k++){ s[k]=pC[k]-(*j)->pN[k]; }
          // print separation distance
          if((*i)->vertexIsNice(*j)){ newfile << sqrt(dot(s,s)) << endl;}
          else {newfile << -sqrt(dot(s,s)) << endl;}
        }
      }
    }
    newfile.close();
  }
}

void Container::writeDistancesOneSidedHausdorff_noself(void) {
  cout << "Compute 1-sided Hausdorff distances (no self)..";
  cout.flush();
  char file[FILENAME_SIZE];
  // create output filename
  sprintf(file,"%sclosest_point_distances_one_sided_hausdorff_noself.dat",OUTPUT_DATA_DIR.c_str());
  // open output file
  std::ofstream newfile;
  newfile.open (file,std::ios::trunc);
  if(newfile.is_open()){
    newfile.precision(4);
    int total_vertices=0,vertices_used=0;
    std::vector<double> storage;
    // for each object
    for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
      // create map: string->double
      map_sd hausdorff;
      sd_iterator w;
      hausdorff.clear();
      // for each vertex in object
      for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
        total_vertices++;
        // if vertex has a closest face
        if((*j)->cl!=NULL && (binary_search(frozen.begin(),frozen.end(),*j)==false)){
          // get closest point
          double pC[3];
          computePC((*j)->cl,*j,pC);
          // compute separation vector
          double s[3];
          for(int k=0;k<3;k++){ s[k]=pC[k]-(*j)->pN[k]; }
          // compute separation distance
          double dd;
          if((*i)->vertexIsNice(*j)){ dd=sqrt(dot(s,s));}
          else {dd=-sqrt(dot(s,s));}
          // if object different from self
          if((*i)->name!=(*j)->cl->v[0]->o->name){
            vertices_used++;
            // if object found in map
            w=hausdorff.find((*j)->cl->v[0]->o->name);
            if(w!=hausdorff.end()){
              // if new separation distance is larger than stored value
              if(dd>(*w).second){hausdorff[(*j)->cl->v[0]->o->name]=dd;}
            } else {
              // add distance to map
              hausdorff[(*j)->cl->v[0]->o->name]=dd;
            }
          }
        }
      }
      // for one sided hausdorff, print each element of map
      for(w=hausdorff.begin();w!=hausdorff.end();w++){
        storage.push_back((*w).second);
      }
    }
    newfile << "# One-sided Hausdorff distances - no self\n"
          << "# total_vertices " << total_vertices
          << ", vertices_used " << vertices_used << endl;
    for(std::vector<double>::iterator i=storage.begin();i!=storage.end();i++){
      newfile << *i << endl;
    }
    newfile.close();
  } else {
    cout << "Error. Could not open file: " << file << endl;
  }
  cout << "complete.\n";
  cout.flush();
}

void Container::writeDistancesTwoSidedHausdorff_noself(void) {
  cout << "Compute 2-sided Hausdorff distances (no self)..";
  cout.flush();
  char file[FILENAME_SIZE];
  // create output filename
  sprintf(file,"%sclosest_point_distances_two_sided_hausdorff_noself.dat",OUTPUT_DATA_DIR.c_str());
  // open output file
  std::ofstream newfile;
  newfile.open (file,std::ios::trunc);
  if(newfile.is_open()){
    newfile.precision(4);
    // create map: string->double
    map_sd hausdorff;
    sd_iterator w,x;
    std::string ss,t;
    int total_vertices=0,vertices_used=0;
    std::vector<double> storage;
    // for each object
    for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
      // for each vertex in object
      for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
        total_vertices++;
        // if vertex has a closest face
        if((*j)->cl!=NULL && (binary_search(frozen.begin(),frozen.end(),*j)==false)){
          // get closest point
          double pC[3];
          computePC((*j)->cl,*j,pC);
          // compute separation vector
          double s[3];
          for(int k=0;k<3;k++){ s[k]=pC[k]-(*j)->pN[k]; }
          // compute separation distance
          double dd;
          if((*i)->vertexIsNice(*j)){ dd=sqrt(dot(s,s));}
          else {dd=-sqrt(dot(s,s));}
          // if object not equal to self
          if((*j)->cl->v[0]->o->name!=(*i)->name){
            vertices_used++;
            // if object found in map
            ss=(*j)->cl->v[0]->o->name+"_"+(*i)->name;
            w=hausdorff.find(ss);
            t=(*i)->name+"_"+(*j)->cl->v[0]->o->name;
            x=hausdorff.find(t);
            if(w!=hausdorff.end()){
              // if new separation distance is larger than stored value
              if(dd>(*w).second){hausdorff[ss]=dd;}
            } else if (x!=hausdorff.end()){
              // if new separation distance is larger than stored value
              if(dd>(*x).second){hausdorff[t]=dd;}
            } else {
              // add distance to map
              hausdorff[ss]=dd;
            }
          }
        }
      }
      // for two sided hausdorff, print each element of map
      for(w=hausdorff.begin();w!=hausdorff.end();w++){
        storage.push_back((*w).second);
      }
    }
    newfile << "# Two-sided Hausdorff distances - no self\n"
          << "# total_vertices " << total_vertices
          << ", vertices_used " << vertices_used << endl;
    for(std::vector<double>::iterator i=storage.begin();i!=storage.end();i++){
      newfile << *i << endl;
    }
    newfile.close();
  } else {
    cout << "Error. Could not open file: " << file << endl;
  }
  cout << "complete.\n";
  cout.flush();
}

void Container::writeDistancesOneSidedHausdorff(void) {
  cout << "Compute one-sided Hausdorff distances..........";
  cout.flush();
  char file[FILENAME_SIZE];
  // create output filename
  sprintf(file,"%sclosest_point_distances_one_sided_hausdorff.dat",OUTPUT_DATA_DIR.c_str());
  // open output file
  std::ofstream newfile;
  newfile.open (file,std::ios::trunc);
  if(newfile.is_open()){
    newfile.precision(4);
    int total_vertices=0,vertices_used=0;
    std::vector<double> storage;
    // for each object
    for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
      // create map: string->double
      map_sd hausdorff;
      sd_iterator w;
      hausdorff.clear();
      // for each vertex in object
      for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
        total_vertices++;
        // if vertex has a closest face
        if((*j)->cl!=NULL && (binary_search(frozen.begin(),frozen.end(),*j)==false)){
          vertices_used++;
          // get closest point
          double pC[3];
          computePC((*j)->cl,*j,pC);
          // compute separation vector
          double s[3];
          for(int k=0;k<3;k++){ s[k]=pC[k]-(*j)->pN[k]; }
          // compute separation distance
          double dd;
          if((*i)->vertexIsNice(*j)){ dd=sqrt(dot(s,s));}
          else {dd=-sqrt(dot(s,s));}
          // if object found in map
          w=hausdorff.find((*j)->cl->v[0]->o->name);
          if(w!=hausdorff.end()){
            // if new separation distance is larger than stored value
            if(dd>(*w).second){hausdorff[(*j)->cl->v[0]->o->name]=dd;}
          } else {
            // add distance to map
            hausdorff[(*j)->cl->v[0]->o->name]=dd;
          }
        }
      }
      // for one sided hausdorff, print each element of map
      for(w=hausdorff.begin();w!=hausdorff.end();w++){
        storage.push_back((*w).second);
      }
    }
    newfile << "# One-sided Hausdorff distances\n"
          << "# total_vertices " << total_vertices
          << ", vertices_used " << vertices_used << endl;
    for(std::vector<double>::iterator i=storage.begin();i!=storage.end();i++){
      newfile << *i << endl;
    }
    newfile.close();
  } else {
    cout << "Error. Could not open file: " << file << endl;
  }
  cout << "complete.\n";
  cout.flush();
}

void Container::writeDistancesTwoSidedHausdorff(void) {
  cout << "Compute two-sided Hausdorff distances..........";
  cout.flush();
  char file[FILENAME_SIZE];
  // create output filename
  sprintf(file,"%sclosest_point_distances_two_sided_hausdorff.dat",OUTPUT_DATA_DIR.c_str());
  // open output file
  std::ofstream newfile;
  newfile.open (file,std::ios::trunc);
  if(newfile.is_open()){
    newfile.precision(4);
    // create map: string->double
    map_sd hausdorff;
    sd_iterator w,x;
    std::string ss,t;
    int total_vertices=0,vertices_used=0;
    std::vector<double> storage;
    // for each object
    for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
      // for each vertex in object
      for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
        total_vertices++;
        // if vertex has a closest face
        if((*j)->cl!=NULL && (binary_search(frozen.begin(),frozen.end(),*j)==false)){
          vertices_used++;
          // get closest point
          double pC[3];
          computePC((*j)->cl,*j,pC);
          // compute separation vector
          double s[3];
          for(int k=0;k<3;k++){ s[k]=pC[k]-(*j)->pN[k]; }
          // compute separation distance
          double dd;
          if((*i)->vertexIsNice(*j)){ dd=sqrt(dot(s,s));}
          else {dd=-sqrt(dot(s,s));}
          // if object found in map
          ss=(*j)->cl->v[0]->o->name+"_"+(*i)->name;
          w=hausdorff.find(ss);
          t=(*i)->name+"_"+(*j)->cl->v[0]->o->name;
          x=hausdorff.find(t);
          if(w!=hausdorff.end()){
            // if new separation distance is larger than stored value
            if(dd>(*w).second){hausdorff[ss]=dd;}
          } else if (x!=hausdorff.end()){
            // if new separation distance is larger than stored value
            if(dd>(*x).second){hausdorff[t]=dd;}
          } else {
            // add distance to map
            hausdorff[ss]=dd;
          }
        }
      }
      // for two sided hausdorff, print each element of map
      for(w=hausdorff.begin();w!=hausdorff.end();w++){
        storage.push_back((*w).second);
      }
    }
    newfile << "# Two-sided Hausdorff distances\n"
          << "# total_vertices " << total_vertices
          << ", vertices_used " << vertices_used << endl;
    for(std::vector<double>::iterator i=storage.begin();i!=storage.end();i++){
      newfile << *i << endl;
    }
    newfile.close();
  } else {
    cout << "Error. Could not open file: " << file << endl;
  }
  cout << "complete.\n";
  cout.flush();
}

// #####################################################
// #####################################################
void Space::initBoxes(void) {
  // subdivide space
  num_space[0] = (int) ceil( (world[1]-world[0])/(SPACE_LENGTH*SCALE) );
  num_space[1] = (int) ceil( (world[3]-world[2])/(SPACE_LENGTH*SCALE) );
  num_space[2] = (int) ceil( (world[5]-world[4])/(SPACE_LENGTH*SCALE) );
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

  if (PRINT_FLAG) {
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
  double br[6];
  // identify face bounding box limits
  xv.clear();
  yv.clear();
  zv.clear();
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
  char file[FILENAME_SIZE];
  sprintf(file,"%s%s",OUTPUT_DATA_DIR.c_str(),SPACE_LOG_FILE);
  Sfile.open(file);
  Sfile.width(15);
  Sfile << left << "#boxes(x,y,z)";
  Sfile.width(15);
  Sfile << left << "total #boxes";	
  Sfile.width(30);
  Sfile << left << "world bounds [xmin xmax ymin ymax zmin zmax]" << endl;
  char str[15];
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

double Vertex::getSqSepDist(Container *c){
  // get closest point
  double pC[3];
  c->computePC(cl,this,pC);
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
  lp[1][0] = lp[0][0]+n[0]/L*RAY_EPSILON*SCALE;
  lp[0][1] = cy;
  lp[1][1] = lp[0][1]+n[1]/L*RAY_EPSILON*SCALE;
  lp[0][2] = cz;
  lp[1][2] = lp[0][2]+n[2]/L*RAY_EPSILON*SCALE;
}

int Container::findExtraPoint(Space &s,Vertex *v,double p[3],int index){
  double lp[2][3];
  //	std::vector<int> face_flags;
  std::vector<Object*> tmp;
  std::vector<Face*>::iterator jj;
  std::vector<Face*> crossed_faces,unique_faces,edge_faces;
  std::pair<std::vector<Face*>::iterator,std::vector<Face*>::iterator> pp;
  std::pair<std::vector<Object*>::iterator,std::vector<Object*>::iterator> ppp;
  getExtraRay(v,lp,index); // returns ray
  // look for intersected faces along ray, i.e. between face centroid and RAY ORIGIN
  collectNiceFaces(s,lp,unique_faces);
  //	findIntersectedFaces(lp,unique_faces,crossed_faces,face_flags);
  //	cout << "\nfindExtraPoint=" << unique_faces.size();
  findIntersectedFaces(lp,unique_faces,crossed_faces,edge_faces);
  sort(crossed_faces.begin(),crossed_faces.end());
  // look for adjacent face in crossed_faces
  pp=equal_range(crossed_faces.begin(),crossed_faces.end(),v->f[index]);
  // if found, then remove
  if(pp.first!=pp.second){crossed_faces.erase(pp.first);}
  // DEBUG
  if( (edge_faces.size()%2)==1 ){
    cout << "\n\nError: how can # edge_faces be odd?\n"
          << "Container::findExtraPoint: "
          << "current vertex:\n";
    v->printVertex(v->o->name);
    cout << endl;
    cout << "Container::findExtraPoint: edge_faces:\n";
    for(std::vector<Face*>::iterator i=edge_faces.begin();i!=edge_faces.end();i++){
      (*i)->printFace((*i)->v[0]->o->name);
      cout << endl;
    }
    cout << endl;
    cout << "Container::findExtraPoint: crossed_faces:\n";
    for(std::vector<Face*>::iterator i=crossed_faces.begin();i!=crossed_faces.end();i++){
      (*i)->printFace((*i)->v[0]->o->name);
      cout << endl;
    }
    cout << endl;
    exit(0);
  }
  // DEBUG
  findOddMeshes(crossed_faces,edge_faces,tmp);
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

bool Container::findCrossed1(Space &s,Vertex *v,double lp[2][3],std::vector<Object*> &c){
  // find and return crossed objects between pN and extracellular point
  //	std::vector<int> face_flags;
  std::vector<Face*>::iterator i;
  std::vector<Face*> crossed_faces,unique_faces,edge_faces;
  std::pair<std::vector<Face*>::iterator,std::vector<Face*>::iterator> pp;	
  collectNiceFaces(s,lp,unique_faces);
  //	cout << ", findCrossed1=" << unique_faces.size();
  findIntersectedFaces(lp,unique_faces,crossed_faces,edge_faces);
  /*if(v->match(10,"a002_FILTERED_SMOOTH_SMOOTH_SMOOTH")){
    if(crossed_faces.empty()==false){
    for(std::vector<Face*>::iterator k=crossed_faces.begin();k!=crossed_faces.end();k++){
    cout << "\nContainer::findCrossed1: crossed_face\n";
    (*k)->printFace((*k)->v[0]->o->name);
    cout<< endl;
    }
    } else {cout << "\nContainer::findCrossed1: crossed_faces is empty.\n";}
    if(edge_faces.empty()==false){
    for(std::vector<Face*>::iterator k=edge_faces.begin();k!=edge_faces.end();k++){
    cout << "\nContainer::findCrossed1: edge_face\n";
    (*k)->printFace((*k)->v[0]->o->name);
    cout<< endl;
    }
    } else {cout << "\nContainer::findCrossed1: edge_faces is empty.\n";}
    }
    */
  // remove current vertex adjacent faces from crossed_faces
  if(crossed_faces.empty()==false){
    sort(crossed_faces.begin(),crossed_faces.end());
    // for each adjacent face
    for(i=v->f.begin();i!=v->f.end();i++){
      pp=equal_range(crossed_faces.begin(),crossed_faces.end(),*i);
      // if adjacent face is found in crossed_faces, then remove from crossed_faces
      if(pp.first!=pp.second){crossed_faces.erase(pp.first);}
    }
  }
  /*	// DEBUG
        if( (edge_faces.size()%2)==1 ){
        cout << "\n\nError(BEFORE): how can # edge_faces be odd?\n"
        << "Container::findCrossed1: "
        << "current vertex:\n";
        v->printVertex(v->o->name);
        cout << endl;
        cout << "Container::findCrossed1: edge_faces:\n";
        for(std::vector<Face*>::iterator g=edge_faces.begin();g!=edge_faces.end();g++){
        (*g)->printFace((*g)->v[0]->o->name);
        cout << endl;
        }
        cout << endl;
        cout << "Container::findCrossed1: crossed_faces:\n";
        for(std::vector<Face*>::iterator g=crossed_faces.begin();g!=crossed_faces.end();g++){
        (*g)->printFace((*g)->v[0]->o->name);
        cout << endl;
        }
        cout << endl;
        exit(0);
        }
  // DEBUG
  */
  // remove current vertex adjacent faces from edge_faces
  if(edge_faces.empty()==false){
    sort(edge_faces.begin(),edge_faces.end());
    // for each adjacent face
    for(i=v->f.begin();i!=v->f.end();i++){
      pp=equal_range(edge_faces.begin(),edge_faces.end(),*i);
      // if adjacent face is found in edge_faces, then remove from edge_faces
      if(pp.first!=pp.second){edge_faces.erase(pp.first);}
    }
  }

  if(edge_faces.empty()==false){return false;}

  // DEBUG
  /*	if( (edge_faces.size()%2)==1 ){
        cout << "\n\nError (AFTER): how can # edge_faces be odd?\n"
        << "Container::findCrossed1: "
        << "current vertex:\n";
        v->printVertex(v->o->name);
        cout << endl;
        cout << "Container::findCrossed1: edge_faces:\n";
        for(std::vector<Face*>::iterator g=edge_faces.begin();g!=edge_faces.end();g++){
        (*g)->printFace((*g)->v[0]->o->name);
        cout << endl;
        }
        cout << endl;
        cout << "Container::findCrossed1: crossed_faces:\n";
        for(std::vector<Face*>::iterator g=crossed_faces.begin();g!=crossed_faces.end();g++){
        (*g)->printFace((*g)->v[0]->o->name);
        cout << endl;
        }
        cout << endl;
        exit(0);
        }
        */
  // DEBUG
  findOddMeshes(crossed_faces,edge_faces,c);
  return true;
}
/*
    void Container::findCrossed2(Space &s,double lp[2][3],std::vector<Object*> &c){
// find and return crossed objects between pN and extracellular point
//	std::vector<int> face_flags;
std::vector<Face*> crossed_faces,unique_faces,edge_faces;
collectNiceFaces(s,lp,unique_faces);
if(
!distinguishable(7304.70669169,lp[0][0]) &&
!distinguishable(2322.15650988,lp[0][1]) &&
!distinguishable(5593.50981227,lp[0][2])
){
for(std::vector<Face*>::iterator i=unique_faces.begin();i!=unique_faces.end();i++){
cout << "Container::findCrossed2: unique face\n";
(*i)->printFaceCP();
}
}

findIntersectedFaces(lp,unique_faces,crossed_faces,edge_faces);
if(
!distinguishable(7304.70669169,lp[0][0]) &&
!distinguishable(2322.15650988,lp[0][1]) &&
!distinguishable(5593.50981227,lp[0][2])
){
if(crossed_faces.empty()==false){
for(std::vector<Face*>::iterator i=crossed_faces.begin();i!=crossed_faces.end();i++){
cout << "Container::findCrossed2: crossed face\n";
(*i)->printFace((*i)->v[0]->o->name);
cout << endl;
}
} else {cout << "Container::findCrossed2: crossed_faces is empty.\n";}
if(edge_faces.empty()==false){
for(std::vector<Face*>::iterator i=edge_faces.begin();i!=edge_faces.end();i++){
cout << "Container::findCrossed2: edge face\n";
(*i)->printFace((*i)->v[0]->o->name);
cout << endl;
}
} else {cout << "Container::findCrossed2: edge_faces is empty.\n";}
}

// DEBUG
if( (edge_faces.size()%2)==1 ){	
cout.precision(12);
cout << "\n\nError: how can # edge_faces be odd?\n"
<< "Container::findCrossed2: "
<< "lp ["
<< lp[0][0] << " "
<< lp[0][1] << " "
<< lp[0][2] << "]["
<< lp[1][0] << " "
<< lp[1][1] << " "
<< lp[1][2] << "]\n";
cout << "Container::findCrossed2: edge_faces:\n";
for(std::vector<Face*>::iterator i=edge_faces.begin();i!=edge_faces.end();i++){
(*i)->printFace((*i)->v[0]->o->name);
cout << endl;
}
cout << endl;
cout << "Container::findCrossed2: crossed_faces:\n";
for(std::vector<Face*>::iterator i=crossed_faces.begin();i!=crossed_faces.end();i++){
(*i)->printFace((*i)->v[0]->o->name);
cout << endl;
}
cout << endl;
exit(0);
}
// DEBUG
findOddMeshes(crossed_faces,edge_faces,c);
}*/

bool Container::findCrossed2(Space &s,double lp[2][3],std::vector<Object*> &c){
  // find and return crossed objects between pN and extracellular point
  //	std::vector<int> face_flags;
  std::vector<Face*> crossed_faces,unique_faces,edge_faces;
  collectNiceFaces(s,lp,unique_faces);
  /*	if(
        !distinguishable(7304.70669169,lp[0][0]) &&
        !distinguishable(2322.15650988,lp[0][1]) &&
        !distinguishable(5593.50981227,lp[0][2])
        ){
        for(std::vector<Face*>::iterator i=unique_faces.begin();i!=unique_faces.end();i++){
        cout << "Container::findCrossed2: unique face\n";
        (*i)->printFaceCP();
        }
        }
        */

  //	cout << ", findCrossed2=" << unique_faces.size() << endl;
  findIntersectedFaces(lp,unique_faces,crossed_faces,edge_faces);
  if(edge_faces.empty()==false){return false;}
  /*if(
    !distinguishable(7304.70669169,lp[0][0]) &&
    !distinguishable(2322.15650988,lp[0][1]) &&
    !distinguishable(5593.50981227,lp[0][2])
    ){
    if(crossed_faces.empty()==false){
    for(std::vector<Face*>::iterator i=crossed_faces.begin();i!=crossed_faces.end();i++){
    cout << "Container::findCrossed2: crossed face\n";
    (*i)->printFace((*i)->v[0]->o->name);
    cout << endl;
    }
    } else {cout << "Container::findCrossed2: crossed_faces is empty.\n";}
    if(edge_faces.empty()==false){
    for(std::vector<Face*>::iterator i=edge_faces.begin();i!=edge_faces.end();i++){
    cout << "Container::findCrossed2: edge face\n";
    (*i)->printFace((*i)->v[0]->o->name);
    cout << endl;
    }
    } else {cout << "Container::findCrossed2: edge_faces is empty.\n";}
    }
    */
  // DEBUG
  if( (edge_faces.size()%2)==1 ){	
    cout.precision(12);
    cout << "\n\nError: how can # edge_faces be odd?\n"
          << "Container::findCrossed2: "
          << "lp ["
          << lp[0][0] << " "
          << lp[0][1] << " "
          << lp[0][2] << "]["
          << lp[1][0] << " "
          << lp[1][1] << " "
          << lp[1][2] << "]\n";
    cout << "Container::findCrossed2: edge_faces:\n";
    for(std::vector<Face*>::iterator i=edge_faces.begin();i!=edge_faces.end();i++){
      (*i)->printFace((*i)->v[0]->o->name);
      cout << endl;
    }
    cout << endl;
    cout << "Container::findCrossed2: crossed_faces:\n";
    for(std::vector<Face*>::iterator i=crossed_faces.begin();i!=crossed_faces.end();i++){
      (*i)->printFace((*i)->v[0]->o->name);
      cout << endl;
    }
    cout << endl;
    exit(0);
  }
  // DEBUG
  findOddMeshes(crossed_faces,edge_faces,c);
  return true;
}

void Container::collectCrossed(Space &s,Vertex *v,std::vector<Object*> &cb){
  std::vector<Object*> ca;
  std::pair<std::vector<Object*>::iterator,std::vector<Object*>::iterator> pp;
  double lp[2][3],p[3];
  // find point, p, outside of current object
  unsigned int index = 0;
  bool good = false;
  while(good==false){
    if (index>(v->f.size()-1)){
      cout << "\n\nContainer::collectCrossed: "
            << "Intersected edge along pathc from pN to ECP on every adjacent face.\n"
            << "Warning! Setting undetermined vertex to 'nice'.\n";
      v->printVertex(v->o->name);
      cout << endl;
      cb.clear();
      return;
    }
    while (!findExtraPoint(s,v,p,index)){
      index++;
      if (index>(v->f.size()-1)){
        cout << "\n\nContainer::collectCrossed: "
              << "Failed to find valid extracellular point on any adjacent face.\n"
              << "Warning! Setting undetermined vertex to 'nice'.\n";
        v->printVertex(v->o->name);
        cout << endl;
        cb.clear();
        return;
      }
    }
    // grab intersected objects between pN and RAY ORIGIN and return as ca
    lp[0][0]=v->pN[0];
    lp[0][1]=v->pN[1];
    lp[0][2]=v->pN[2];
    lp[1][0]=p[0];
    lp[1][1]=p[1];
    lp[1][2]=p[2];
    good=findCrossed1(s,v,lp,ca);
    if(good==false){index++;}
  }
  // grab intersected objects along RAY as cb
  lp[0][0]=p[0];
  lp[0][1]=p[1];
  lp[0][2]=p[2];
  lp[1][0]=p[0];
  lp[1][1]=p[1];
  lp[1][2]=p[2];
  /*	if(v->match(10,"a002_FILTERED_SMOOTH_SMOOTH_SMOOTH")){	
        cout.precision(12);
        cout << "\n\nContainer::collectCrossed: "
        << "pN=["
        << v->pN[0] << " "
        << v->pN[1] << " "
        << v->pN[2] << "]\n";
        }
        */
  // get all rays sorted from shortest to longest
  double rays[6][3];
  findClosestAxis(s,v,lp,rays); // returns ray
  // initialize axis selection index
  // applies to list of distances to world boundary sorted shortest to longest
  // 0 == shortest direction, i.e. zeroth element
  // 1 == second shortest direction, i.e. first element
  // ...
  // 5 == longest direction, i.e. fifth element
  int select=0;
  bool keep = false;
  // while an axis needs to be found and select is less than max
  while(keep==false && select<6){
    lp[0][0]=p[0];
    lp[0][1]=p[1];
    lp[0][2]=p[2];
    lp[1][0]=rays[select][0];
    lp[1][1]=rays[select][1];
    lp[1][2]=rays[select][2];
    /*		if(v->match(10,"a002_FILTERED_SMOOTH_SMOOTH_SMOOTH")){	
                cout.precision(12);
                cout << "Container::collectCrossed: "
                << "ray["
                << lp[0][0] << " "
                << lp[0][1] << " "
                << lp[0][2] << "]["
                << lp[1][0] << " "
                << lp[1][1] << " "
                << lp[1][2] << "]\n";
                }
                */
    keep=findCrossed2(s,lp,cb);
    if(keep==false){select++;}
  }
  if(select==6){
    cout << "\n\nContainer::collectCrossed: ERROR: All six principal coordinates"
          << " resulted in face edge collisions.\n";
    v->printVertex(v->o->name);
    cout << endl << endl;
    exit(0);
  }
  /*	if(v->match(10,"a002_FILTERED_SMOOTH_SMOOTH_SMOOTH")){
        if(ca.empty()==true){
        cout << "\nContainer::collectCrossed: "
        << "NO odd-intersected objects between pN and ray origin.\n";
        } else {
        cout << "\nContainer::collectCrossed: "
        << ca.size() << " intersected objects between pN and ray origin.\n";
        for(std::vector<Object*>::iterator i=ca.begin();i!=ca.end();i++){
        cout << " object " << (*i)->name << " was intersected between pN and ray origin.\n";
        }
        }
        if(cb.empty()==true){
        cout << "\nContainer::collectCrossed: "
        << "NO odd-intersected objects along ray.\n";
        } else {
        cout << "\nContainer::collectCrossed: "
        << cb.size() << " intersected objects along ray.\n";
        for(std::vector<Object*>::iterator i=cb.begin();i!=cb.end();i++){
        cout << " object " << (*i)->name << " was intersected along ray.\n";
        }
        }
        exit(0);
        }
        */
  // remove all meshes in ca from cb
  // since the object is not odd relative to pN
  sort(cb.begin(),cb.end());
  for (std::vector<Object*>::iterator ii=ca.begin();ii!=ca.end();ii++){
    pp=equal_range(cb.begin(),cb.end(),*ii);
    // if object in ca is found in cb, then remove from cb
    if(pp.first!=pp.second){cb.erase(pp.first);}
    // else add it to cb
    else {
      cb.push_back(*ii);
    }
  }
}
/*
    bool Container::updateNiceness(Vertex *v,std::vector<Object*> &cb){
    std::pair<std::vector<Object*>::iterator,std::vector<Object*>::iterator> pp;
    int old_nice = v->getVertexNiceness();
// if vertex niceness changes then set flag=true
bool flag = false;
// if cb is not empty, then vertex is not nice
if (cb.empty()==false) {
v->setVertexNiceness(1,this);
// if vertex  nice
if (!old_nice){
flag = true;
pp=equal_range(cb.begin(),cb.end(),v->o);
// if vertex is inside self object
if(pp.first!=pp.second){
v->setVertexNiceness(2,this);
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
v->setVertexNiceness(0,this);
}
return flag;
}*/

bool Container::updateNiceness(Vertex *v,std::vector<Object*> &cb){
  std::pair<std::vector<Object*>::iterator,std::vector<Object*>::iterator> pp;
  int old_nice = v->getVertexNiceness();
  int new_nice = 0;
  // if cb is not empty
  if (cb.empty()==false) {
    // vertex is not nice
    new_nice++;
    // if vertex is inside self object
    pp=equal_range(cb.begin(),cb.end(),v->o);
    if(pp.first!=pp.second){new_nice++;}
  }
  v->setVertexNiceness(new_nice,this);
  // if vertex niceness changes then return true, else return false
  return old_nice!=new_nice;
}

bool Container::checkNiceness(Space &s,Vertex *v) {
  std::vector<Object*> cb;
  // collect objects inside which vertex lies
  collectCrossed(s,v,cb);
  /*	// DEBUG
        if(v->match(10,"a002_FILTERED_SMOOTH_SMOOTH_SMOOTH")==true){
        if(cb.empty()==false){
        cout << "\n\nContainer::checkNiceness: " << "cb.size()=" << cb.size() << endl;
        for(std::vector<Object*>::iterator i=cb.begin();i!=cb.end();i++){
        cout << "Container::checkNiceness: collected object = " << (*i)->name << endl;
        }
        }
        }
  // DEBUG
  */
  // update niceness of vertex based on cb
  return updateNiceness(v,cb);
}

void Container::findNice(Space &s) {
  cout << "Find nice vertices.............................";
  cout.flush();
  int iii=1;
  double goal = 0.2;
  printf("0%%..");
  fflush(stdout);
  std::vector<Object*>::iterator i;
  std::vector<Vertex*>::iterator j;
  // for each object in container
  for (i=o.begin();i!=o.end();i++) {
    // for each vertex in object
    for (j=(*i)->v.begin();j!=(*i)->v.end();j++) {
      checkNiceness(s,*j);
      // track progress
      double progress = static_cast<double>(iii++)/vertex_count;
      if(progress>goal){
        printf("%d%%..",static_cast<int>(goal*100));
        fflush(stdout);
        goal+=0.2;
      }
    }
  }
  printf("100%%..");
  fflush(stdout);
  cout << "complete.\n";
  cout.flush();
}

// #####################################################
// #####################################################
/*
    void Container::findClosestAxis(Space &s,Vertex *v,double lp[2][3]) {
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
*/

void Container::findClosestAxis(Space &s,Vertex *v,double lp[2][3],double rays[6][3]) {
  // map: (distance to world bounds)double->(direction code)int
  map_di dist;
  dist.insert(std::make_pair(fabs(v->pN[0]-s.world[0]),0));
  dist.insert(std::make_pair(fabs(v->pN[0]-s.world[1]),1));
  dist.insert(std::make_pair(fabs(v->pN[1]-s.world[2]),2));
  dist.insert(std::make_pair(fabs(v->pN[1]-s.world[3]),3));
  dist.insert(std::make_pair(fabs(v->pN[2]-s.world[4]),4));
  dist.insert(std::make_pair(fabs(v->pN[2]-s.world[5]),5));
  int j=0;
  // for each map element
  for(di_iterator i=dist.begin();i!=dist.end();i++){
    int index=(*i).second;
    if		(index==0){
      rays[j][0] = lp[0][0]-2*(s.world[1]-s.world[0]);	// end x
      rays[j][1] = lp[0][1];							// end y
      rays[j][2] = lp[0][2];							// end z
    } else if (index==1){
      rays[j][0] = lp[0][0]+2*(s.world[1]-s.world[0]);	// end x
      rays[j][1] = lp[0][1];							// end y
      rays[j][2] = lp[0][2];							// end z
    } else if (index==2){
      rays[j][0] = lp[0][0];							// end x
      rays[j][1] = lp[0][1]-2*(s.world[3]-s.world[2]);	// end y
      rays[j][2] = lp[0][2];							// end z
    } else if (index==3){
      rays[j][0] = lp[0][0];							// end x
      rays[j][1] = lp[0][1]+2*(s.world[3]-s.world[2]);	// end y
      rays[j][2] = lp[0][2];							// end z
    } else if (index==4){
      rays[j][0] = lp[0][0];							// end x
      rays[j][1] = lp[0][1];							// end y
      rays[j][2] = lp[0][2]-2*(s.world[5]-s.world[4]);	// end z
    } else if (index==5){
      rays[j][0] = lp[0][0];							// end x
      rays[j][1] = lp[0][1];							// end y
      rays[j][2] = lp[0][2]+2*(s.world[5]-s.world[4]);	// end z
    }
    j++;
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
  if(br[1]<br[0]){int t=br[0];br[0]=br[1];br[1]=t;}
  if(br[3]<br[2]){int t=br[2];br[2]=br[3];br[3]=t;}
  if(br[5]<br[4]){int t=br[4];br[4]=br[5];br[5]=t;}
}

void Container::collectNiceFaces(Space &s,double lp[2][3],std::vector<Face*> &uf) {

  // compute box index range that contains ray
  int br[6];
  getBoxIndexRange(s,lp,br);
  ////////// collect boxes to check //////////
  std::vector<Box*> b;
  s.getBoxesFor3DIndices(br,b);
  ////////// gather faces in boxes //////////
  // for each box
  for (std::vector<Box*>::iterator i=b.begin();i!=b.end();i++) {
    // for each face in box
    for (std::vector<Face*>::iterator j=(*i)->f.begin();j!=(*i)->f.end();j++) {
      uf.push_back(*j);
    }
  }
  // keep unique faces
  sort(uf.begin(),uf.end());
  std::vector<Face*>::iterator new_end = unique(uf.begin(),uf.end());
  uf.assign(uf.begin(),new_end);
}

// #####################################################
// #####################################################

void Container::findIntersectedFaces(double lp[2][3],std::vector<Face*> &uf,
                                     std::vector<Face*> &cf,std::vector<Face*> &ef) {
  // uf = unique_faces
  // cf = crossed_faces
  // ff = face_flags
  std::vector<Face*>::iterator j;
  bool line_flag, poly_flag, poly_edge_flag;
  // for each unique polygon
  for (j=uf.begin();j!=uf.end();j++) {
    checkLineFaceIntersection(*j,lp,line_flag,poly_flag,poly_edge_flag,false);
    // does point intersect polygon
    if (line_flag==true) {
      if (poly_edge_flag==true) {
        // add polygon_index to edge array
        ef.push_back(*j);
      } else if (poly_flag==true) {
        // add polygon_index to crossed array
        cf.push_back(*j);
      }
    }
  }
}

// #####################################################
// #####################################################


void Container::findOddObjects(std::vector<Object*> &ol,std::vector<Object*> &odd){
  // sort object list
  sort(ol.begin(),ol.end());

  odd.clear();
  int count=ol.size();
  int k=0;
  // while more Object*s to process
  while (k<count) {
    // if not on last Object*
    if (k+1!=count) {
      // skip identical pairs of object indices
      if (ol[k]==ol[k+1]) {k++;k++;}
      // odd object
      else {
        odd.push_back( ol[k]);
        k++;
      }
    } else { // add remaining object to odd object list
      odd.push_back( ol[k]);
      k++;
    }
  }
}

void Container::findEvenObjects(std::vector<Object*> &ol,std::vector<Object*> &even){
  // sort object list
  sort(ol.begin(),ol.end());

  even.clear();
  int count=ol.size();
  int k=0;
  // while more Object*s to process
  while (k<count) {
    // if not on last Object*
    if (k+1!=count) {
      // skip identical pairs of object indices
      if (ol[k]==ol[k+1]) {
        even.push_back( ol[k]);
        k++;k++;
      }
      // odd object
      else {
        even.push_back( ol[k]);
        k++;
        cout << "\n\nContainer::findEvenObjects: "
              << "Weird. Was expecting only even numbers of edge_faces.\n";
      }
    } else { // add remaining object to even object list
      even.push_back( ol[k]);
      cout << "\n\nContainer::findEvenObjects: "
            << "Weird. Was expecting only even numbers of edge_faces.\n";
      k++;
    }
  }
}

void Container::findOddMeshes(std::vector<Face*> &cf,std::vector<Face*> &ef,
                              std::vector<Object*> &tmp) {
  // find mesh objects crossed an odd number of times by ray

  std::vector<Object*> ol; // object list

  // for each edge crossed face
  if(ef.empty()==false){
    std::vector<Object*> ee;
    for (std::vector<Face*>::iterator i=ef.begin();i!=ef.end();i++) {
      ee.push_back((*i)->v[0]->o);
    }
    findEvenObjects(ee,ol);
  }

  // for each crossed face
  for (std::vector<Face*>::iterator i=cf.begin();i!=cf.end();i++) {
    ol.push_back((*i)->v[0]->o);
  }
  findOddObjects(ol,tmp);

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
  //	double cn[3],on[3];
  //	cf->getNormal(cn);
  //	of->getNormal(on);
  // get face vertex coordinates
  double *opvc[3],*cpvc[3];
  cf->getVertexCoordinates(cpvc);
  of->getVertexCoordinates(opvc);
  // are faces parallel?
  //bool parallel_flag = facesParallel(cf,of);
  // are faces colinear?
  // i.e. is there a vertex on other face 
  // that lies in the plane of current face?
  //bool colinear_flag = facesColinear(cf,of);
  // are faces coplanar?
  //bool coplanar_flag=false;
  //if (facesParallel(cf,of) && facesColinear(cf,of)) {coplanar_flag=true;}
  // get number of unique vertices between current and other face
  int single_shared_vert[2]={-1,-1};
  int num_unique = numUniqueVertices(cf,of,single_shared_vert);
  //bool share_edge_flag=false;
  //bool identical_flag=false;
  //bool share_vert_flag=false;
  //if		(num_unique == 8) {share_vert_flag = true;}
  //else if (num_unique == 7) {share_edge_flag = true;}
  // else if (num_unique == 6) {identical_flag = true;}
  ////////// begin decision tree //////////
  //if (coplanar_flag) {
  if (facesParallel(cf,of) && facesColinear(cf,of)) {
    if (num_unique == 6) {
      // polygons are identical
      return true;
    } else {
      //do polygon edges intersect?	
      // if yes, intersect
      // if no, do not intersect
      if (checkEdgeEdgeIntersection(cf,of,num_unique==7)) {return true;}
      else {return false;}
    }
  } else {
    if (num_unique==8) {
      int m=0,n=1,p=0,q=1;
      // single vertex shared
      if (single_shared_vert[0]==0){m=1;n=2;}
      else if (single_shared_vert[0]==1){n=2;}
      if (single_shared_vert[1]==0){p=1;q=2;}
      else if (single_shared_vert[1]==1){q=2;}
      double lp[2][3] = {{cpvc[m][0],cpvc[m][1],cpvc[m][2]},
        {cpvc[n][0],cpvc[n][1],cpvc[n][2]}};
      bool line_flag=false, poly_flag=false, poly_edge_flag;
      checkLineFaceIntersection(of,lp,line_flag,poly_flag,poly_edge_flag,false);
      // do faces intersect?
      if (line_flag && (poly_flag || poly_edge_flag)) {return(1);}
      lp[0][0] = opvc[p][0];
      lp[0][1] = opvc[p][1];
      lp[0][2] = opvc[p][2];
      lp[1][0] = opvc[q][0];
      lp[1][1] = opvc[q][1];
      lp[1][2] = opvc[q][2];
      checkLineFaceIntersection(cf,lp,line_flag,poly_flag,poly_edge_flag,false);
      // do faces intersect?
      if (line_flag && (poly_flag || poly_edge_flag)) {return(1);}
      else {return false;}
    }else if (!(num_unique==7)) {
      // do faces intersect?
      if (checkFaceEdgeIntersection(cf,of)) {return true;}
      else {return false;}
    } else { return false;}
  }
}

// #####################################################
// #####################################################
void Space::clearBoxes(void){
  std::vector<Box*>::iterator i;
  // for each box in space, clear vector of face*
  for (i=b.begin();i!=b.end();i++) {(*i)->f.clear();}
}

void Space::recordFace(std::vector<Box*> &ptr,Face* f) {
  std::vector<Box*>::iterator i;
  // for each box, add face
  for (i=ptr.begin();i!=ptr.end();i++) {(*i)->f.push_back(f);}
}

void Container::assignFacesToBoxes(Space &s) {
  cout << "Assign faces to boxes..........................";
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
      /*			if((*j)->match(41958,"d003_FILTERED_SMOOTH_SMOOTH_SMOOTH")==true){
                                cout << endl << endl;
                                for(std::vector<Box*>::iterator k=bp.begin();k!=bp.end();k++){
                                (*k)->printBox(s.world[0],s.world[2],s.world[4]);
                                }
                                cout << endl << endl;
                                }
                                */
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

void vdErrorToRGB(double rgb[3],double vd)
{
  if(vd > MAX_VD)
  {
    rgb[0] = MAX_COLOR[0];
    rgb[1] = MAX_COLOR[1];
    rgb[2] = MAX_COLOR[2];
  }
  else if(vd < MIN_VD)
  {
    rgb[0] = MIN_COLOR[0];
    rgb[1] = MIN_COLOR[1];
    rgb[2] = MIN_COLOR[2];
  }
  else if (vd > MID_VD)
  {
    double u = (vd-MID_VD)/(MAX_VD-MID_VD);
    for(int i=0;i<3;i++)
    {
      rgb[i] = MID_COLOR[i]+u*(MAX_COLOR[i]-MID_COLOR[i]);
    }
  }
  else if (vd < MID_VD)
  {
    double u = (vd-MIN_VD)/(MID_VD-MIN_VD);
    for(int i=0;i<3;i++)
    {
      rgb[i] = MIN_COLOR[i]+u*(MID_COLOR[i]-MIN_COLOR[i]);
    }
  }
  else
  {
    rgb[0] = MID_COLOR[0];
    rgb[1] = MID_COLOR[1];
    rgb[2] = MID_COLOR[2];
  }
}

void ecwErrorToRGB(double rgb[3],double ecw_error)
{
  if(ecw_error > MAX_ERROR)
  {
    rgb[0] = MAX_COLOR[0];
    rgb[1] = MAX_COLOR[1];
    rgb[2] = MAX_COLOR[2];
  }
  else if(ecw_error < MIN_ERROR)
  {
    rgb[0] = MIN_COLOR[0];
    rgb[1] = MIN_COLOR[1];
    rgb[2] = MIN_COLOR[2];
  }
  else if (ecw_error > MID_ERROR)
  {
    double u = (ecw_error-MID_ERROR)/(MAX_ERROR-MID_ERROR);
    for(int i=0;i<3;i++)
    {
      rgb[i] = MID_COLOR[i]+u*(MAX_COLOR[i]-MID_COLOR[i]);
    }
  }
  else if (ecw_error < MID_ERROR)
  {
    double u = (ecw_error-MIN_ERROR)/(MID_ERROR-MIN_ERROR);
    for(int i=0;i<3;i++)
    {
      rgb[i] = MIN_COLOR[i]+u*(MID_COLOR[i]-MIN_COLOR[i]);
    }
  }
  else
  {
    rgb[0] = MID_COLOR[0];
    rgb[1] = MID_COLOR[1];
    rgb[2] = MID_COLOR[2];
  }
}

void Container::printVD(Monitor& stats){
  cerr << "Print virtual displacements....................";
  cerr.flush();
  double min = 1e30;
  double max = -1e30;
  std::ofstream V,F,C,D;
  std::vector<std::string> group;
  int index=1;
  double goal = 0.2;
  fprintf(stderr,"0%%..");
  fflush(stderr);
  // open DX_file
  std::string str = OUTPUT_DATA_DIR + DX_FILE;
  D.open(str.c_str());
  if(D.is_open()==false){
    fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
    exit(0);
  }
  // for each object in container
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    // open vertex_file
    str = OUTPUT_DATA_DIR + (*i)->name + VERTS_SUFFIX;
    V.open(str.c_str());
    if(V.is_open()==false){
      fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
      exit(0);
    }
    // update DX file
    D << "object \"mesh_verts_" << (*i)->name << "\" class array type float "
      << "rank 1 shape 3 items " << (*i)->v.size() << " ascii data file "
      << (*i)->name << VERTS_SUFFIX << endl
      << "  attribute \"dep\" string \"positions\"" << endl << endl;
    // open color_file
    str = OUTPUT_DATA_DIR + (*i)->name + COLORS_SUFFIX;
    C.open(str.c_str());
    if(C.is_open()==false){
      fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
      exit(0);
    }
    // update DX file
    D << "object \"mesh_colors_" << (*i)->name << "\" class array type float "
      << "rank 1 shape 3 items " << (*i)->v.size() << " ascii data file "
      << (*i)->name << COLORS_SUFFIX << endl
      << "  attribute \"dep\" string \"positions\"" << endl << endl;
    // for each vertex in object
    int vertex_num=1;
    for (std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++)
    {
      if((*j)->index!=vertex_num)
      {
        cout << "\n\nError: vertex index " << (*j)->index << " does not increase"
              << " incrementally. Must be equal to " << vertex_num << " for dx compatibility.\n\n";
        exit(0);
      }
      else
      {
        vertex_num++;
      }
      // if vertex found in topN
      td_iterator tt = stats.old.find(*j);
      if(tt!=stats.old.end())
      {
        double vd = sqrt((*tt).second);
        // write vertex location to vertex_file as "x y z\n"
        V << (*j)->pN[0] << " " << (*j)->pN[1] << " " << (*j)->pN[2] << endl;
        // map ecw error to color
        if(vd > max) {max=vd;}
        if(vd < min) {min=vd;}
        double rgb[3] = {0,0,0};
        vdErrorToRGB(rgb,vd);
        // write vertex color to color_file as "r g b\n" in range 0 to 1.0
        //C << (*j)->pN[0] << " " << (*j)->pN[1] << " " << (*j)->pN[2] << endl;
        C << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl;
      }
      else
      {
        // write vertex location to vertex_file as "x y z\n"
        V << (*j)->pN[0] << " " << (*j)->pN[1] << " " << (*j)->pN[2] << endl;
        // write vertex color to color_file as "r g b\n" in range 0 to 1.0
        C << NULL_COLOR[0] << " " << NULL_COLOR[1] << " " << NULL_COLOR[2] << endl;
      }
      // track progress
      double progress = static_cast<double>(index++)/vertex_count;
      if(progress>goal){
        fprintf(stderr,"%d%%..",static_cast<int>(goal*100));
        fflush(stderr);
        goal+=0.2;
      }
    }
    C.close();
    V.close();
    // open connections_file
    str = OUTPUT_DATA_DIR + (*i)->name + CONNECTIONS_SUFFIX;
    F.open(str.c_str());
    if(F.is_open()==false){
      fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
      exit(0);
    }
    // update DX file
    D << "object \"mesh_triangles_" << (*i)->name << "\" class array type int "
      << "rank 1 shape 3 items " << (*i)->f.size() << " ascii data file "
      << (*i)->name << CONNECTIONS_SUFFIX << endl
      << "  attribute \"ref\" string \"positions\"" << endl
      << "  attribute \"element type\" string \"triangles\"" << endl
      << endl << endl
      << "object \"my_mesh_" << (*i)->name << "\" field" << endl
      << "    component \"positions\" value \"mesh_verts_" << (*i)->name << "\"" << endl
      << "    component \"connections\" value \"mesh_triangles_" << (*i)->name << "\"" << endl
      << "    component \"colors\" value \"mesh_colors_" << (*i)->name << "\"" << endl << endl;
    std::string gs = "        member \"" + (*i)->name + "\" value \"my_mesh_" + (*i)->name + "\"";
    group.push_back(gs);
    // for each face in object
    for (std::vector<Face*>::iterator j=(*i)->f.begin();j!=(*i)->f.end();j++)
    {
      // write vertex indices to connections_files as "index1 index2 index3\n"
      F << (*j)->v[0]->index-1 << " " << (*j)->v[1]->index-1 << " " << (*j)->v[2]->index-1 << endl;
    }
    F.close();
  }
  // complete DX file
  D << "object \"meshes\" group" << endl;
  for(std::vector<std::string>::iterator i=group.begin();i!=group.end();i++)
  {
    D << *i << endl;
  }
  D.close();
  fprintf(stderr,"100%%..");
  fflush(stderr);
  cerr << "complete.\n";
  cerr.flush();
  cout << "\n\nmax virtual displacement = " << max << endl;
  cout << "min virtual displacement = " << min << endl << endl;
}


void Container::getSeparationDistancesAndECW(Space &s,Monitor& stats){
  cerr << "Get separation distances.......................";
  cerr.flush();
  double min = 1e30;
  double max = -1e30;
  std::ofstream V,F,C,D;
  std::vector<std::string> group;
  int index=1;
  double goal = 0.2;
  fprintf(stderr,"0%%..");
  fflush(stderr);
  // open DX_file
  std::string str = OUTPUT_DATA_DIR + DX_FILE;
  D.open(str.c_str());
  if(D.is_open()==false){
    fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
    exit(0);
  }
  // for each object in container
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    // open vertex_file
    str = OUTPUT_DATA_DIR + (*i)->name + VERTS_SUFFIX;
    V.open(str.c_str());
    if(V.is_open()==false){
      fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
      exit(0);
    }
    // update DX file
    D << "object \"mesh_verts_" << (*i)->name << "\" class array type float "
      << "rank 1 shape 3 items " << (*i)->v.size() << " ascii data file "
      << (*i)->name << VERTS_SUFFIX << endl
      << "  attribute \"dep\" string \"positions\"" << endl << endl;
    // open color_file
    str = OUTPUT_DATA_DIR + (*i)->name + COLORS_SUFFIX;
    C.open(str.c_str());
    if(C.is_open()==false){
      fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
      exit(0);
    }
    // update DX file
    D << "object \"mesh_colors_" << (*i)->name << "\" class array type float "
      << "rank 1 shape 3 items " << (*i)->v.size() << " ascii data file "
      << (*i)->name << COLORS_SUFFIX << endl
      << "  attribute \"dep\" string \"positions\"" << endl << endl;
    // for each vertex in object
    int vertex_num=1;
    for (std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++)
    {
      if((*j)->index!=vertex_num)
      {
        cout << "\n\nError: vertex index " << (*j)->index << " does not increase"
              << " incrementally. Must be equal to " << vertex_num << " for dx compatibility.\n\n";
        exit(0);
      }
      else
      {
        vertex_num++;
      }
      ////////// find closest point to current vertex //////////
      // false => just add value to table, do not touch existing elements
      //			there should be none anyway
      findClosest(s,*j,stats,false);
      // if vertex has a closest face and not frozen
      if((*j)->cl!=NULL && (binary_search(frozen.begin(),frozen.end(),*j)==false))
      {
        // write vertex location to vertex_file as "x y z\n"
        V << (*j)->pN[0] << " " << (*j)->pN[1] << " " << (*j)->pN[2] << endl;
        // get closest point
        double pC[3];
        //computePC((*j)->cl,*j,pC);
        computePC((*j)->cl,*j,pC);
        // compute separation vector
        double ss[3];
        for(int k=0;k<3;k++){ ss[k]=pC[k]-(*j)->pN[k]; }
        // map ecw error to color
        double ecw_error = sqrt(dot(ss,ss)) - TARGET_SEPARATION;
        if(ecw_error > max) {max=ecw_error;}
        if(ecw_error < min) {min=ecw_error;}
        double rgb[3] = {0,0,0};
        ecwErrorToRGB(rgb,ecw_error);
        // write vertex color to color_file as "r g b\n" in range 0 to 1.0
        //C << (*j)->pN[0] << " " << (*j)->pN[1] << " " << (*j)->pN[2] << endl;
        C << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl;
      }
      else
      {
        // write vertex location to vertex_file as "x y z\n"
        V << (*j)->pN[0] << " " << (*j)->pN[1] << " " << (*j)->pN[2] << endl;
        // write vertex color to color_file as "r g b\n" in range 0 to 1.0
        C << NULL_COLOR[0] << " " << NULL_COLOR[1] << " " << NULL_COLOR[2] << endl;
      }
      // track progress
      double progress = static_cast<double>(index++)/vertex_count;
      if(progress>goal){
        fprintf(stderr,"%d%%..",static_cast<int>(goal*100));
        fflush(stderr);
        goal+=0.2;
      }
    }
    C.close();
    V.close();
    // open connections_file
    str = OUTPUT_DATA_DIR + (*i)->name + CONNECTIONS_SUFFIX;
    F.open(str.c_str());
    if(F.is_open()==false){
      fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
      exit(0);
    }
    // update DX file
    D << "object \"mesh_triangles_" << (*i)->name << "\" class array type int "
      << "rank 1 shape 3 items " << (*i)->f.size() << " ascii data file "
      << (*i)->name << CONNECTIONS_SUFFIX << endl
      << "  attribute \"ref\" string \"positions\"" << endl
      << "  attribute \"element type\" string \"triangles\"" << endl
      << endl << endl
      << "object \"my_mesh_" << (*i)->name << "\" field" << endl
      << "    component \"positions\" value \"mesh_verts_" << (*i)->name << "\"" << endl
      << "    component \"connections\" value \"mesh_triangles_" << (*i)->name << "\"" << endl
      << "    component \"colors\" value \"mesh_colors_" << (*i)->name << "\"" << endl << endl;
    std::string gs = "        member \"" + (*i)->name + "\" value \"my_mesh_" + (*i)->name + "\"";
    group.push_back(gs);
    // for each face in object
    for (std::vector<Face*>::iterator j=(*i)->f.begin();j!=(*i)->f.end();j++)
    {
      // write vertex indices to connections_files as "index1 index2 index3\n"
      F << (*j)->v[0]->index-1 << " " << (*j)->v[1]->index-1 << " " << (*j)->v[2]->index-1 << endl;
    }
    F.close();
  }
  // complete DX file
  D << "object \"meshes\" group" << endl;
  for(std::vector<std::string>::iterator i=group.begin();i!=group.end();i++)
  {
    D << *i << endl;
  }
  D.close();
  fprintf(stderr,"100%%..");
  fflush(stderr);
  cerr << "complete.\n";
  cerr.flush();
  cout << "\n\nmax ecw error = " << max << endl;
  cout << "min ecw error = " << min << endl << endl;
}

void Container::getSeparationDistances(Space &s,Monitor& stats){
  cout << "Get separation distances.......................";
  cout.flush();
  int index=1;
  double goal = 0.2;
  printf("0%%..");
  fflush(stdout);
  // for each object in container
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    // for each vertex in object
    for (std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++) {
      ////////// find closest point to current vertex //////////
      // false => just add value to table, do not touch existing elements
      //			there should be none anyway
      findClosest(s,*j,stats,false);
      // track progress
      double progress = static_cast<double>(index++)/vertex_count;
      if(progress>goal){
        printf("%d%%..",static_cast<int>(goal*100));
        fflush(stdout);
        goal+=0.2;
      }
    }
  }
  printf("100%%..");
  fflush(stdout);
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
/*
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
}*/

void Container::getBoxes(std::vector<Box*> &bp,Vertex *v,Space &s){
  int cbi[3],br[6];
  // compute box index that contains current_vertex
  cbi[0] = s.location2Index(v->pN[0],"x");
  cbi[1] = s.location2Index(v->pN[1],"y");
  cbi[2] = s.location2Index(v->pN[2],"z");
  // compute bounds of box that contains current_vertex
  double xmin = s.world[0]+cbi[0]*SPACE_LENGTH*SCALE;
  double xmax = s.world[0]+(cbi[0]+1)*SPACE_LENGTH*SCALE;
  double ymin = s.world[2]+cbi[1]*SPACE_LENGTH*SCALE;
  double ymax = s.world[2]+(cbi[1]+1)*SPACE_LENGTH*SCALE;
  double zmin = s.world[4]+cbi[2]*SPACE_LENGTH*SCALE;
  double zmax = s.world[4]+(cbi[2]+1)*SPACE_LENGTH*SCALE;
  // compute minimum box_range that covers search_radius
  double sr = sqrt(SEARCH_RADIUS_SQ*SCALE*SCALE);
  br[0]=cbi[0]-static_cast<int>(ceil((sr-(v->pN[0]-xmin))/(SPACE_LENGTH*SCALE)));
  br[1]=cbi[0]+static_cast<int>(ceil((sr-(xmax-v->pN[0]))/(SPACE_LENGTH*SCALE)));
  br[2]=cbi[1]-static_cast<int>(ceil((sr-(v->pN[1]-ymin))/(SPACE_LENGTH*SCALE)));
  br[3]=cbi[1]+static_cast<int>(ceil((sr-(ymax-v->pN[1]))/(SPACE_LENGTH*SCALE)));
  br[4]=cbi[2]-static_cast<int>(ceil((sr-(v->pN[2]-zmin))/(SPACE_LENGTH*SCALE)));
  br[5]=cbi[2]+static_cast<int>(ceil((sr-(zmax-v->pN[2]))/(SPACE_LENGTH*SCALE)));
  // handle case where vertex lies on subspace boundary
  if (!(br[0]*SPACE_LENGTH*SCALE-v->pN[0])){br[0]--;}
  if (!(br[2]*SPACE_LENGTH*SCALE-v->pN[1])){br[2]--;}
  if (!(br[4]*SPACE_LENGTH*SCALE-v->pN[2])){br[4]--;}
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
/*
    bool Container::faceInNeighborhood(Face *f,Vertex *v){
// if face is in different object than vertex, then return false
//	if(f->v[0]->o!=v->o){return false;} // CACHEGRIND
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
}*/

bool Container::faceInNeighborhood(Face *f,Vertex *v){
  // if face is in different object than vertex, then return false
  //	if(f->v[0]->o!=v->o){return false;} // CACHEGRIND
  // else if face and vertex are in same object
  if(binary_search(v->nf.begin(),v->nf.end(),f)){
    // if face is in current vertex neighborhood
    return true;
  } else {
    // face is NOT in current vertex neighborhood
    // but we must check adjacent faces
    // since total neighborhood = adjacent faces + non-adjacent faces
    if(binary_search(v->f.begin(),v->f.end(),f)){
      // if face is in adjacent face vector
      return true;
    } else {
      // face is NOT in adjacent face vector
      return false;
    }
  }
}

bool Container::faceIsAdjacent(Face *f,Vertex *v){
  if(binary_search(v->f.begin(),v->f.end(),f)){
    // if face is in adjacent face vector
    return true;
  } else {
    // face is NOT in adjacent face vector
    return false;
  }
}
/*
//void Container::getCandidateFaces(vector<Box*> &bp,Vertex *v,hashset_f &cf){
void Container::getCandidateFaces(vector<Box*> &bp,Vertex *v,std::vector<Face*> &cf){
// for each box in search
for (std::vector<Box*>::iterator j=bp.begin();j!=bp.end();j++) {
// for each face in box
for (std::vector<Face*>::iterator k=(*j)->f.begin();k!=(*j)->f.end();k++) {
// if face is not in current vertex neighborhood, face is candidate
if (!faceInNeighborhood(*k,v)){
//				cf.insert(*k);
cf.push_back(*k);
}
}
}
}*/

//void Container::getCandidateFaces(vector<Box*> &bp,Vertex *v,hashset_f &cf){
void Container::getCandidateFaces(std::vector<Box*> &bp,Vertex *v,std::vector<Face*> &cf){
  // for each box in search
  for (std::vector<Box*>::iterator j=bp.begin();j!=bp.end();j++) {
    // for each face in box
    for (std::vector<Face*>::iterator k=(*j)->f.begin();k!=(*j)->f.end();k++) {
      // if face is not adjacent to current vertex, face is candidate
      if (!faceIsAdjacent(*k,v)){
        //				cf.insert(*k);
        cf.push_back(*k);
      }
    }
  }
}

bool Monitor::findTopN(Vertex *vv,tv_iterator &tt,int &rank){
  // topN is of type table_v
  //multimap<double,Vertex*> table_v;
  bool found=false;
  rank=0;
  // for each pair in topN
  for(tv_iterator i=topN.begin();i!=topN.end();i++){
    rank++;
    // if vertex* matches target
    if((*i).second==vv){
      tt=i;
      found=true;
      break;
    }
  }
  return found;
}

void Monitor::validateVertex(int index,Container &c){
  ///// check for discrepencies in separation errors /////
  // for each pair in topN
  for(tv_iterator i=topN.begin();i!=topN.end();i++){
    // if vertex index matches target
    if((*i).second->index==index){
      // if the se stored in topN does not match the computed se
      // then topN is stale, i.e. not up-to-date
      if( distinguishable((*i).first,fabs(c.computeSeparationError((*i).second,(*i).second->getSqSepDist(&c))))) {
        cout << "\n\nMonitor::validateMultimap: "
              << "Error! topN-stored se "
              << (*i).first
              << " does not match current value "
              << fabs(c.computeSeparationError((*i).second,(*i).second->getSqSepDist(&c)))
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
  bool flag = false;
  for(tv_iterator i=topN.begin();i!=topN.end();i++){
    if((*i).second->cl==NULL){
      cout << "\n\nMonitor::validateMultimap: "
            << "Error! topN contains vertex* ("
            << (*i).second->index << ") with no closest point.\n";
      flag=true;
    }
  }
  if(flag==true){exit(0);}
}

void Monitor::validateOld(void){
  tv_iterator j;
  // for each pair in old
  for(td_iterator i=old.begin();i!=old.end();i++){
    //		if((*i).second==0){cout << "Monitor::validateOld: Error. (*i).second==0.\n";exit(0);}
    if(!entryInTopN((*i).first,(*i).second,j)){
      cout << "\n\nMonitor::validateOld: "
            << "Error! old contains element (vertex->"
            << (*i).first->index << ", vd=" << (*i).second
            << ") not found in TopN.\n";
      exit(0);
    }
  }
}

bool Monitor::entryInTopN(Vertex *v,double vd_old,tv_iterator &j){
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
    cout.precision(9);
    cout << "\n\nError. Matching multimap element not found.\n";
    v->printVertex(v->o->name);
    cout << "\n\nold vd " << vd_old << endl << endl;
    // for each entry in range
    for(tv_iterator i=p.first;i!=p.second;i++){
      cout << "topN contains ["
            << (*i).first << " , "
            << (*i).second->index << "]\n";
    }
    bool booyah=true;
    for(tv_iterator i=topN.begin();i!=topN.end();i++){
      if((*i).second==v){
        cout << "current vertex found in topN: "
              << " vertex " << (*i).second->index
              << " " << (*i).first << endl << endl;
        booyah = false;
      }
    }
    if(booyah==true){
      cout << "Problem vertex was NOT found anywhere in topN.\n";
    }
    exit(0);
  }
}

//void Monitor::validateTopN(char* str){
void Monitor::validateTopN(Container *c){
  for(tv_iterator i=topN.begin();i!=topN.end();i++){
    if((*i).second->cl==NULL){
      cout << "ERROR. vertex: " << (*i).second->o->name << "->" << (*i).second->index << " in topN has no closest point.\n";
      exit(0);
    }
    // compute new vertex coords
    double ppH[3]; // new holding position coordinates (x,y,z)
    // load energy map
    (*i).second->computeNewCoords(c,ppH,gain);
    // compute virtual displacement
    double vd = (ppH[0]-(*i).second->pN[0])*(ppH[0]-(*i).second->pN[0])+
          (ppH[1]-(*i).second->pN[1])*(ppH[1]-(*i).second->pN[1])+
          (ppH[2]-(*i).second->pN[2])*(ppH[2]-(*i).second->pN[2]);
    // if vd differs from topN value
    if(distinguishable(vd,(*i).first)){
      cout << "ERROR. computed vd (" << vd << ") of vertex (" << (*i).second->o->name << "->" << (*i).second->index
            << ") differs from stored value (" << (*i).first << ") in topN.\n";
      exit(0);
    }
  }
}

void Monitor::cleanTopN(Vertex *v,double vd_old){
  tv_iterator t;
  // if find table entry with old se
  if (entryInTopN(v,vd_old,t)){
    // remove entry
    topN.erase(t);
  }
}

void Monitor::updateTopN(Vertex *v,double vd_old,double vd_new,bool flag){
  if (flag){
    tv_iterator t;
    // if find table entry with old se
    if (entryInTopN(v,vd_old,t)){
      // remove entry
      topN.erase(t);
    } else {
    }
  }
  // add new se entry to table
  topN.insert(std::make_pair(vd_new,v));
}

double Container::computeSeparationError(Vertex* v,double old_se){
  // set target separation
  double TS = 0;
  // if closest face is inside vertex neighborhood
  if(faceInNeighborhood(v->cl,v)==true){ TS = LOOP_TARGET_SEPARATION*SCALE; }
  else { TS = TARGET_SEPARATION*SCALE; }
  // compute separation error (signed value)
  if(!v->o->vertexIsNice(v)){return sqrt(old_se)+TS;}
  else{return sqrt(old_se)-TS;}
}

void Space::getSR(double sr[6],int indices[6]){
  // [xmin xmax ymin ymax zmin zmax]
  sr[0] = indices[0]*SPACE_LENGTH*SCALE+world[0];
  sr[1] = (indices[1]+1)*SPACE_LENGTH*SCALE+world[0];
  sr[2] = indices[2]*SPACE_LENGTH*SCALE+world[2];
  sr[3] = (indices[3]+1)*SPACE_LENGTH*SCALE+world[2];
  sr[4] = indices[4]*SPACE_LENGTH*SCALE+world[4];
  sr[5] = (indices[5]+1)*SPACE_LENGTH*SCALE+world[4];
}

void Container::getPD2(double pd2[6],double sr[6],Vertex *v){
  // [xmin xmax ymin ymax zmin zmax]
  pd2[0] = (sr[0]-v->pN[0])*(sr[0]-v->pN[0]);
  pd2[1] = (sr[1]-v->pN[0])*(sr[1]-v->pN[0]);
  pd2[2] = (sr[2]-v->pN[1])*(sr[2]-v->pN[1]);
  pd2[3] = (sr[3]-v->pN[1])*(sr[3]-v->pN[1]);
  pd2[4] = (sr[4]-v->pN[2])*(sr[4]-v->pN[2]);
  pd2[5] = (sr[5]-v->pN[2])*(sr[5]-v->pN[2]);
}

bool Space::expandSR(bool gate,double cd,double pd[6],int inc[6],int indices[6]){
  bool reanalyze = false;
  if(gate==true){ // if cp was found
    // expand sr in each direction if needed and allowed
    if(cd>pd[0]&&inc[0]<2){
      // expand sr in -x direction
      indices[0] = screenIndex(indices[0]-1,"x");
      inc[0]++;
      reanalyze = true;
    }
    if(cd>pd[1]&&inc[1]<2){
      // expand sr in +x direction
      indices[1] = screenIndex(indices[1]+1,"x");
      inc[1]++;
      reanalyze = true;
    }
    if(cd>pd[2]&&inc[2]<2){
      // expand sr in -y direction
      indices[2] = screenIndex(indices[2]-1,"y");
      inc[2]++;
      reanalyze = true;
    }
    if(cd>pd[3]&&inc[3]<2){
      // expand sr in +y direction
      indices[3] = screenIndex(indices[3]+1,"y");
      inc[3]++;
      reanalyze = true;
    }
    if(cd>pd[4]&&inc[4]<2){
      // expand sr in -z direction
      indices[4] = screenIndex(indices[4]-1,"z");
      inc[4]++;
      reanalyze = true;
    }
    if(cd>pd[5]&&inc[5]<2){
      // expand sr in +z direction
      indices[5] = screenIndex(indices[5]+1,"z");
      inc[5]++;
      reanalyze = true;
    }
  } else { // else if no cp was found
    // if expanding allowed in any direction then expand sr
    if(inc[0]<2){
      // expand sr in -x direction
      indices[0] = screenIndex(indices[0]-1,"x");
      inc[0]++;
      reanalyze = true;
    }
    if(inc[1]<2){
      // expand sr in +x direction
      indices[1] = screenIndex(indices[1]+1,"x");
      inc[1]++;
      reanalyze = true;
    }
    if(inc[2]<2){
      // expand sr in -y direction
      indices[2] = screenIndex(indices[2]-1,"y");
      inc[2]++;
      reanalyze = true;
    }
    if(inc[3]<2){
      // expand sr in +y direction
      indices[3] = screenIndex(indices[3]+1,"y");
      inc[3]++;
      reanalyze = true;
    }
    if(inc[4]<2){
      // expand sr in -z direction
      indices[4] = screenIndex(indices[4]-1,"z");
      inc[4]++;
      reanalyze = true;
    }
    if(inc[5]<2){
      // expand sr in +z direction
      indices[5] = screenIndex(indices[5]+1,"z");
      inc[5]++;
      reanalyze = true;
    }
  }
  return reanalyze;
}

bool Container::findClosest(Space &s,Vertex *v,Monitor& stats,bool flag) {
  if(binary_search(frozen.begin(),frozen.end(),v)==false){
    // squared distance between closest point and v
    double cd=0.0;
    // old and new Face* vectors
    std::vector<Face*> of,nf;
    of.reserve(VECTOR_RESERVE);
    nf.reserve(VECTOR_RESERVE);
    // old and new Box* vectors
    std::vector<Box*> ob,nb;
    ob.reserve(VECTOR_RESERVE);
    nb.reserve(VECTOR_RESERVE);
    // count box index changes
    int inc[6] = {0,0,0,0,0,0};
    //  sr limits [xmin xmax ymin ymax zmin zmax]
    double sr[6];
    //  6 cartesian principal distances squared (pd2) from v to sr limits
    // [xmin xmax ymin ymax zmin zmax]
    double pd2[6];
    // let search region (sr) be box containing current vertex, v
    // get box indices for sr : [xmin xmax ymin ymax zmin zmax]
    int indices[6] = { s.location2Index(v->pN[0],"x"),s.location2Index(v->pN[0],"x"),
      s.location2Index(v->pN[1],"y"),s.location2Index(v->pN[1],"y"),
      s.location2Index(v->pN[2],"z"),s.location2Index(v->pN[2],"z")};
    // get Box* for sr : Space::getBoxesFor3DLocations(v->pN)
    s.getBoxesFor3DIndices(indices,nb);
    // collect Face* in sr
    getCandidateFaces(nb,v,nf);
    // sort and keep unique Face*
    sort(nf.begin(),nf.end());
    std::vector<Face*>::iterator z = unique(nf.begin(),nf.end());
    nf.assign(nf.begin(),z);
    // compute sr limits
    s.getSR(sr,indices);
    // compute 6 cartesian principal distances (pd) from v to sr limits
    getPD2(pd2,sr,v);
    // search for closest point (cp) and closest distance (cd) among faces
    bool gate = false;
    if(searchCP(v,nf,cd)==true){gate=true;}
    // if sr is expanded, then reanalyze
    while(s.expandSR(gate,cd,pd2,inc,indices)==true){
      // copy Box* vector and get Box* for new sr
      ob = nb;
      nb.clear();
      s.getBoxesFor3DIndices(indices,nb);
      // copy Face* vector and collect Face* in new sr
      of = nf;
      nf.clear();
      getCandidateFaces(nb,v,nf);
      // sort and keep unique Face*
      sort(nf.begin(),nf.end());
      z = unique(nf.begin(),nf.end());
      nf.assign(nf.begin(),z);
      // collect Face* in sr that were NOT used before
      // for each face in nf vector
      std::vector<Face*>::iterator i=nf.begin();
      while(i!=nf.end()){
        // if face found in of vector then erase face from nf
        if(binary_search(of.begin(),of.end(),*i)){ nf.erase(i);}
        else { i++; }
      }
      // compute sr limits
      s.getSR(sr,indices);
      // compute 6 cartesian principal distances (pd) from v to sr limits
      getPD2(pd2,sr,v);
      // search for closest point (cp) and closest distance (cd) among NEW faces
      if(searchCP(v,nf,cd)==true){gate=true;}
    }	

    // if closest point was found
    if(gate){
      if(flag){
        // update sets with vertex squared virtual displacement
        stats.updateSets(v,getVertexSqD(v,stats.gain),flag);
      } 
    } else { // vertex has no closest point
      // reset pointer to closest face
      v->cl=NULL;
      if(flag){
        // if vertex* found in old
        if(stats.old.find(v)!=stats.old.end()){
          // remove vertex* from topN
          stats.cleanTopN(v,stats.old[v]);
          // remove vertex* from old
          stats.old.erase(v);
        }
      }
    }
    return gate;
  } else {
    // vertex is frozen
    // reset pointer to closest face
    v->cl=NULL;
    if(flag){
      // if vertex* found in old
      if(stats.old.find(v)!=stats.old.end()){
        // remove vertex* from topN
        stats.cleanTopN(v,stats.old[v]);
        // remove vertex* from old
        stats.old.erase(v);
      }
    }
    return false;
  }
}
/*
    bool Container::searchCP(Vertex *v,std::vector<Face*> &nf,double &squareD){
    bool gate = false;
// if candidate faces were found
if (!nf.empty()){
// get vertex normal
double n[3];
v->getNormal(n);
// define multimap
mmap_d_f df;
// collect vertices, edges, faces
std::vector<Vertex*> verts;
std::vector<Edge*> edges;
std::vector<Face*> faces;
// for each candidate face
for (std::vector<Face*>::iterator j=nf.begin();j!=nf.end();j++) {
// compute face bounding box
double c[3] = {((*j)->bb[0]+(*j)->bb[1])/2.0,
((*j)->bb[2]+(*j)->bb[3])/2.0,
((*j)->bb[4]+(*j)->bb[5])/2.0};
// compute radius^2 (sr) of smallest box-enclosing sphere
double sr = (c[0]-(*j)->bb[0])*(c[0]-(*j)->bb[0])+
(c[1]-(*j)->bb[2])*(c[1]-(*j)->bb[2])+
(c[2]-(*j)->bb[4])*(c[2]-(*j)->bb[4]);
// compute squared distance(sd) from box center to current vertex
double sd = (c[0]-v->pN[0])*(c[0]-v->pN[0])+
(c[1]-v->pN[1])*(c[1]-v->pN[1])+
(c[2]-v->pN[2])*(c[2]-v->pN[2]);
if(sd<sr){
// compute closest point on face
// if the closest point to current vertex was found on this face
verts.push_back((*j)->v[0]);
verts.push_back((*j)->v[1]);
verts.push_back((*j)->v[2]);
edges.push_back((*j)->e[0]);
edges.push_back((*j)->e[1]);
edges.push_back((*j)->e[2]);
faces.push_back(*j);
} else {
// add face to multimap: sd(double)->Face*
df.insert(std::make_pair(sd,*j));
}

}

// check faces with sd<sr
if(verts.empty()==false || edges.empty()==false || faces.empty()==false){
if(computeClosestColl(verts,edges,faces,v,squareD,n)){ gate=true;}
}

// if no closest point was found
if(gate==false){
// since multimap is sorted by sd, smallest to largest
// while no cp is found and more elements in map to check
df_iterator i=df.begin();
while(gate==false && i!=df.end()){
// check next face in list for cp
// if the closest point to current vertex was found on this face
if(computeClosest((*i).second,v,squareD,n)){ gate=true;}
// increment iterator
i++;
}
}
}
return gate;
}*/

bool Container::searchCP(Vertex *v,std::vector<Face*> &nf,double &squareD){
  bool gate = false;
  // if candidate faces were found
  if (!nf.empty()){
    // get vertex normal
    double n[3];
    v->getNormal(n);
    // define multimap
    mmap_d_f df;
    // collect vertices, edges, faces
    std::vector<Vertex*> verts;
    std::vector<Edge*> edges;
    std::vector<Face*> faces;
    // for each candidate face
    for (std::vector<Face*>::iterator j=nf.begin();j!=nf.end();j++) {
      // compute face bounding box
      double c[3] = {((*j)->bb[0]+(*j)->bb[1])/2.0,
        ((*j)->bb[2]+(*j)->bb[3])/2.0,
        ((*j)->bb[4]+(*j)->bb[5])/2.0};
      // compute radius^2 (sr) of smallest box-enclosing sphere
      double sr = (c[0]-(*j)->bb[0])*(c[0]-(*j)->bb[0])+
            (c[1]-(*j)->bb[2])*(c[1]-(*j)->bb[2])+
            (c[2]-(*j)->bb[4])*(c[2]-(*j)->bb[4]);
      // compute squared distance(sd) from box center to current vertex
      double sd = (c[0]-v->pN[0])*(c[0]-v->pN[0])+
            (c[1]-v->pN[1])*(c[1]-v->pN[1])+
            (c[2]-v->pN[2])*(c[2]-v->pN[2]);
      if(sd<sr){
        // compute closest point on face
        // if the closest point to current vertex was found on this face
        verts.push_back((*j)->v[0]);
        verts.push_back((*j)->v[1]);
        verts.push_back((*j)->v[2]);
        edges.push_back((*j)->e[0]);
        edges.push_back((*j)->e[1]);
        edges.push_back((*j)->e[2]);
        faces.push_back(*j);
      } else {
        // add face to multimap: sd(double)->Face*
        df.insert(std::make_pair(sd,*j));
      }

    }
    // check faces with sd<sr
    if(verts.empty()==false || edges.empty()==false || faces.empty()==false){
      if(computeClosestColl(verts,edges,faces,v,squareD,n)){ gate=true;}
    }
    // if no closest point was found
    if(gate==false){
      // since multimap is sorted by sd, smallest to largest
      // while no cp is found and more elements in map to check
      df_iterator i=df.begin();
      while(gate==false && i!=df.end()){
        // check next face in list for cp
        // if the closest point to current vertex was found on this face
        if(computeClosest((*i).second,v,squareD,n)){ gate=true;}
        // increment iterator
        i++;
      }
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
  checkLineFaceIntersection(f,lp,line_flag,poly_flag,poly_edge_flag,true);
  // if intersection point is on face,then save point
  if (poly_flag) {p.add(lp[1][0],lp[1][1],lp[1][2]);}
  //	return (poly_flag || poly_edge_flag);
  return (poly_flag);
}

//void Container::getEdgeIntersection(Vertex *v,double *P[3],Point &p){
//  // first pair of face vertices, P[0],P[1]
//  double e[3]={P[1][0]-P[0][0],P[1][1]-P[0][1],P[1][2]-P[0][2]};
//  double vv[3]={v->pN[0]-P[0][0],v->pN[1]-P[0][1],v->pN[2]-P[0][2]};
//  double u = dot(vv,e)/dot(e,e);
//  if (u>0 && u<1) {
//    p.add(P[0][0]+u*(e[0]),P[0][1]+u*(e[1]),P[0][2]+u*(e[2]));
//  }
//  // second pair of face vertices, P[1],P[2]
//  e[0]=P[2][0]-P[1][0];
//  e[1]=P[2][1]-P[1][1];
//  e[2]=P[2][2]-P[1][2];
//  vv[0]=v->pN[0]-P[1][0];
//  vv[1]=v->pN[1]-P[1][1];
//  vv[2]=v->pN[2]-P[1][2];
//  u = dot(vv,e)/dot(e,e);
//  if (u>0 && u<1) {
//    p.add(P[1][0]+u*(e[0]),P[1][1]+u*(e[1]),P[1][2]+u*(e[2]));
//  }
//  // third pair of face vertices, P[2],P[0]
//  e[0]=P[0][0]-P[2][0];
//  e[1]=P[0][1]-P[2][1];
//  e[2]=P[0][2]-P[2][2];
//  vv[0]=v->pN[0]-P[2][0];
//  vv[1]=v->pN[1]-P[2][1];
//  vv[2]=v->pN[2]-P[2][2];
//  u = dot(vv,e)/dot(e,e);
//  if (u>0 && u<1) {
//    p.add(P[2][0]+u*(e[0]),P[2][1]+u*(e[1]),P[2][2]+u*(e[2]));
//  }
//}

void Container::getEdgeIntersection(Vertex *v,double *P[3],Point &p){
  // first pair of face vertices, P[0],P[1]
  double e[3]={P[1][0]-P[0][0],P[1][1]-P[0][1],P[1][2]-P[0][2]};
  double vv[3]={v->pN[0]-P[0][0],v->pN[1]-P[0][1],v->pN[2]-P[0][2]};
  double num = dot(vv,e);
  double den = dot(e,e);
  //if (u>0 && u<1) {
  if ( ((num>0 && den>0)||(num<0 && den<0))
       && fabs(num)<fabs(den)) {
    double u = num/den;
    p.add(P[0][0]+u*(e[0]),P[0][1]+u*(e[1]),P[0][2]+u*(e[2]));
  }
  // second pair of face vertices, P[1],P[2]
  e[0]=P[2][0]-P[1][0];
  e[1]=P[2][1]-P[1][1];
  e[2]=P[2][2]-P[1][2];
  vv[0]=v->pN[0]-P[1][0];
  vv[1]=v->pN[1]-P[1][1];
  vv[2]=v->pN[2]-P[1][2];
  num = dot(vv,e);
  den = dot(e,e);
  //if (u>0 && u<1) {
  if ( ((num>0 && den>0)||(num<0 && den<0))
       && fabs(num)<fabs(den)) {
    double u = num/den;
    p.add(P[1][0]+u*(e[0]),P[1][1]+u*(e[1]),P[1][2]+u*(e[2]));
  }
  // third pair of face vertices, P[2],P[0]
  e[0]=P[0][0]-P[2][0];
  e[1]=P[0][1]-P[2][1];
  e[2]=P[0][2]-P[2][2];
  vv[0]=v->pN[0]-P[2][0];
  vv[1]=v->pN[1]-P[2][1];
  vv[2]=v->pN[2]-P[2][2];
  num = dot(vv,e);
  den = dot(e,e);
  //if (u>0 && u<1) {
  if ( ((num>0 && den>0)||(num<0 && den<0))
       && fabs(num)<fabs(den)) {
    double u = num/den;
    p.add(P[2][0]+u*(e[0]),P[2][1]+u*(e[1]),P[2][2]+u*(e[2]));
  }
}

/*
    void Container::computePC(Face *f,Vertex *v,double pC[3]) {
// initialize point class instance with current vertex
Point p(v->pN[0],v->pN[1],v->pN[2]);
// get face vertex coordinates
double *P[3] = {&(f->v[0])->pN[0],&(f->v[1])->pN[0],&(f->v[2])->pN[0]};
// compute vector connecting arbitrary face vertex and current vertex
double diff[3] = {P[0][0]-v->pN[0],P[0][1]-v->pN[1],P[0][2]-v->pN[2]};
// get face normal
double fn[3];
f->getNormal(fn);
// compute indicators
double num=dot(fn,diff);
//compute point on face plane that is closest to current vertex
// i.e. intersection of plane with face normal through current vertex
// add to p if intersection is on face
// return true if point lies inside or on an edge of face, false otherwise
bool inside = false;
inside = getPlaneIntersection(f,v,fn,num,dot(fn,fn),p);
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
pC[0] = p.a;
pC[1] = p.b;
pC[2] = p.c;

}*/

Face* Container::getNonAdjacentFaceOfVertex(Vertex *cv,Vertex *tv){
  // for each adjacent face of tv
  for(std::vector<Face*>::iterator i=tv->f.begin();i!=tv->f.end();i++){
    // if face is not adjacent to vertex cv
    if(binary_search(cv->f.begin(),cv->f.end(),*i)==false){
      return *i;
    }
  }
  cout << "Container::getNonAdjacentFaceOfVertex: Error:"
        << "no acceptable face was found!\n";
  exit(0);
}

Face* Container::getNonAdjacentFaceOfEdge(Vertex *cv,Edge *ee){
  // if first edge face is not adjacent to vertex cv
  if(binary_search(cv->f.begin(),cv->f.end(),ee->f1)==false){
    return ee->f1;
  }
  // if second edge face is not adjacent to vertex cv
  if(binary_search(cv->f.begin(),cv->f.end(),ee->f2)==false){
    return ee->f2;
  }
  cout << "Container::getNonAdjacentFaceOfEdge: Error:"
        << "no acceptable face was found!\n";
  exit(0);
}

bool Container::computeClosestColl(std::vector<Vertex*> &verts,std::vector<Edge*> &edges,
                                   std::vector<Face*> &faces,Vertex *v,double &squareD,double vn[3]) {
  bool signal = false;
  // initialize Face* map and triplet map
  map_df fmap;
  map_dt tmap;
  Triplet *t;
  ////////// handle vertices //////////
  // sort and keep unique vertices
  sort(verts.begin(),verts.end());
  std::vector<Vertex*>::iterator q = unique(verts.begin(),verts.end());
  verts.assign(verts.begin(),q);
  // for each vertex
  for(std::vector<Vertex*>::iterator i=verts.begin();i!=verts.end();i++){
    // squared distance between v and *i
    double sqd = ((*i)->pN[0]-v->pN[0])*((*i)->pN[0]-v->pN[0])+
          ((*i)->pN[1]-v->pN[1])*((*i)->pN[1]-v->pN[1])+
          ((*i)->pN[2]-v->pN[2])*((*i)->pN[2]-v->pN[2]);
    // create triplet
    t = new Triplet((*i)->pN[0],(*i)->pN[1],(*i)->pN[2]);
    // add squared distance to fmap
    fmap[sqd]=getNonAdjacentFaceOfVertex(v,*i);
    // add vertex to tmap
    if(!tmap.count(sqd)){tmap[sqd]=t;}
  }
  double e[3],vv[3],u;
  ////////// handle edges //////////
  // sort and keep unique edges
  sort(edges.begin(),edges.end());
  std::vector<Edge*>::iterator qq = unique(edges.begin(),edges.end());
  edges.assign(edges.begin(),qq);
  // for each edge
  for(std::vector<Edge*>::iterator i=edges.begin();i!=edges.end();i++){
    // get edge vertices
    Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
    (*i)->getVertices(v1,v2,o1,o2);
    // compute edge intersection point
    e[0]=v2->pN[0]-v1->pN[0];
    e[1]=v2->pN[1]-v1->pN[1];
    e[2]=v2->pN[2]-v1->pN[2];
    vv[0]=v->pN[0]-v1->pN[0];
    vv[1]=v->pN[1]-v1->pN[1];
    vv[2]=v->pN[2]-v1->pN[2];
    u = dot(vv,e)/dot(e,e);
    if (u>0 && u<1) {
      double ix = v1->pN[0]+u*e[0];
      double iy = v1->pN[1]+u*e[1];
      double iz = v1->pN[2]+u*e[2];
      // squared distance between v and i
      double sqd = (ix-v->pN[0])*(ix-v->pN[0])+
            (iy-v->pN[1])*(iy-v->pN[1])+
            (iz-v->pN[2])*(iz-v->pN[2]);
      // create triplet
      t = new Triplet(ix,iy,iz);
      // add squared distance to fmap
      fmap[sqd]=getNonAdjacentFaceOfEdge(v,*i);
      // add vertex to tmap
      if(!tmap.count(sqd)){tmap[sqd]=t;}
    }
  }

  ////////// handle faces //////////
  // sort and keep unique faces
  sort(faces.begin(),faces.end());
  std::vector<Face*>::iterator qqq = unique(faces.begin(),faces.end());
  faces.assign(faces.begin(),qqq);
  // initialize point class instance with current vertex
  Point p(v->pN[0],v->pN[1],v->pN[2]);
  // store face vertex coordinates
  double *P[3];
  // for each face
  for(std::vector<Face*>::iterator i=faces.begin();i!=faces.end();i++){
    // get face normal
    double fn[3];
    (*i)->getNormal(fn);
    // get face vertex coordinates
    P[0] = &((*i)->v[0])->pN[0];
    P[1] = &((*i)->v[1])->pN[0];
    P[2] = &((*i)->v[2])->pN[0];
    // compute vector connecting arbitrary face vertex and current vertex
    double diff[3];
    diff[0] = P[0][0]-v->pN[0];
    diff[1] = P[0][1]-v->pN[1];
    diff[2] = P[0][2]-v->pN[2];
    // compute indicators
    double num=dot(fn,diff);
    // if current vertex does not lie on face plane
    bool inside = false;
    p.clear();
    if (fabs(num)>DOUBLE_EPSILON) {
      //compute point on face plane that is closest to current vertex
      // i.e. intersection of plane with face normal through current vertex
      // add to p if intersection is on face
      // return true if point lies inside or on an edge of face, false otherwise
      inside = getPlaneIntersection(*i,v,fn,num,dot(fn,fn),p);
    }
    // if a point has been found, i.e. inside==true
    if(inside==true){
      // create triplet
      t = new Triplet(p.a,p.b,p.c);
      // add squared distance to fmap
      fmap[p.L]=*i;
      // add point to tmap
      if(!tmap.count(p.L)){tmap[p.L]=t;}
    } else {
    }
  }

  ////////// find closest //////////
  double sep_vec[3];
  // invert vertex normal if vertex is not nice
  double vn_copy[3];
  if (v->o->vertexIsNice(v)){vn_copy[0]=vn[0];vn_copy[1]=vn[1];vn_copy[2]=vn[2]; } 
  else 					  {vn_copy[0]=-vn[0];vn_copy[1]=-vn[1];vn_copy[2]=-vn[2]; } 
  // MAKE SURE MAP IS SORTED SMALLEST TO LARGEST
  // for each triplet in map
  mdt_iterator i = tmap.begin();
  mdf_iterator j = fmap.begin();
  mdf_iterator k = j;
  while(i!=tmap.end()){
    double a = (*i).second->x;
    double b = (*i).second->y;
    double c = (*i).second->z;
    // compute separation vector
    sep_vec[0] = a-v->pN[0];
    sep_vec[1] = b-v->pN[1];
    sep_vec[2] = c-v->pN[2];
    // compute cosine of angle between outward normal and separation vector
    // which is equal to dot product of vectors divided by vector magnitudes
    double sign = dot(sep_vec,vn_copy);
    //double cos_angle = sign/sqrt(dot(vn_copy,vn_copy)*dot(sep_vec,sep_vec));
    double cos_angle2 = sign*sign/dot(vn_copy,vn_copy)*dot(sep_vec,sep_vec);
    // is closest point located within angle window as defined in controls.cc?
    // compute square of separation distance
    // closest point must be within a specified angle of vertex normal
    // to avoid grabbing points on same object
    // alternatively, the point is allowed to be within neighborhood_radius of vertex
    // since that implies the point is not near vertex on surface
    //if ( (cos_angle > CLOSEST_POINT_COSINE)
    if ( ((sign>0.0) && (cos_angle2 > CLOSEST_POINT_COSINE*CLOSEST_POINT_COSINE))
         || (v->o!=(*j).second->v[0]->o) 
         || (*i).first<(1.0*SCALE)){
      // if square of separation distance is closer than square of SEARCH_RADIUS
      if ((*i).first<(SEARCH_RADIUS_SQ*SCALE*SCALE)) {
        // if square of separation distance is less than
        // previously saved square of separation distance
        if ((*i).first<squareD||!squareD) {
          // save 
          v->cl=(*j).second;
          squareD = (*i).first;
          signal=true;
          break;
        }
      }
    } else {
    }
    i++;
    j++;
  }
  // cleanup
  // for each tmap
  for(i=tmap.begin();i!=tmap.end();i++){
    delete (*i).second;
  }
  return signal;
}

double l2(double t[3]){
  // return length of t squared
  return t[0]*t[0]+t[1]*t[1]+t[2]*t[2];
}

void cross(double t1[3],double t2[3],double cr[3]){
  // normals cross product = t1 x t2
  cr[0] = t1[1]*t2[2]-t1[2]*t2[1];
  cr[1] = t1[2]*t2[0]-t1[0]*t2[2];
  cr[2] = t1[0]*t2[1]-t1[1]*t2[0];
}

void Container::computePC(Face *f,Vertex *v,double pC[3])
{
  // initialize point class instance with current vertex
  Point p(v->pN[0],v->pN[1],v->pN[2]);
  // get face vertex coordinates
  double *P[3] = {&(f->v[0])->pN[0],&(f->v[1])->pN[0],&(f->v[2])->pN[0]};
  // compute vector connecting arbitrary face vertex and current vertex
  double diff[3] = {P[0][0]-v->pN[0],P[0][1]-v->pN[1],P[0][2]-v->pN[2]};
  // face normal
  double fn[3];
  f->getNormal(fn);
  // compute indicators
  double num=dot(fn,diff);
  // if current vertex does not lie on face plane
  bool inside = false;
  if (fabs(num)<DOUBLE_EPSILON) {
    // vertex does lie on plane
    // compute vectors from face vertices to current vertex
    // del[0][] = {f->v[0] to v}
    // del[1][] = {f->v[1] to v}
    // del[2][] = {f->v[2] to v}
    bool onVertex=false;
    double del[3][3];
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        del[i][j]=v->pN[j]-f->v[i]->pN[j];
      }
      // if vector length is small
      if(l2(del[i])<DOUBLE_EPSILON*DOUBLE_EPSILON){
        // current vertex is on face vertex
        onVertex=true;
        p.add(f->v[i]->pN[0],f->v[i]->pN[1],f->v[i]->pN[2]);
      }
    }
    if(onVertex==false){
      // current vertex does NOT lie on face vertex
      // compute cross products
      // sqL[0] = {|del[0][] x del[1][]|^2}
      // sqL[1] = {|del[1][] x del[2][]|^2}
      // sqL[2] = {|del[2][] x del[0][]|^2}
      double sqL[3];
      bool onEdge=false;
      double c[3];
      for(int i=0;i<3;i++){
        cross(del[i],del[(i+1)%3],c);
        sqL[i]=l2(c);
        // if any squared cross product vector length is small
        if(sqL[i]<DOUBLE_EPSILON*DOUBLE_EPSILON){
          // vertex is on face edge
          // current vertex is the closest point on face to vertex 
          p.add(v->pN[0],v->pN[1],v->pN[2]);
          onEdge=true;
        }
      }
      if(onEdge==false){
        if( (sqL[0]>0 && sqL[1]>0 && sqL[2]>0 ) || 
            (sqL[0]<0 && sqL[1]<0 && sqL[2]<0 )){
          // vertex is inside face
          // current vertex is the closest point on face to vertex 
          p.add(v->pN[0],v->pN[1],v->pN[2]);
        } else {
          // vertex is outside face
          // compute closest point on face edges and vertices
          // add each face vertex
          for (int i=0;i<3;i++) {
            p.add(P[i][0],P[i][1],P[i][2]);
          }
          // add points of minimum distance between current vertex and each face edge
          getEdgeIntersection(v,P,p);
        }
      }
    }
  } else {
    // vertex does NOT lie on plane
    //compute point on face plane that is closest to current vertex
    // i.e. intersection of plane with face normal through current vertex
    // add to p if intersection is on face
    // return true if point lies inside or on an edge of face, false otherwise
    inside = getPlaneIntersection(f,v,fn,num,dot(fn,fn),p);
    // if a point has been found, i.e. inside==true
    // then it is necessarily the closet, we are done
    // but if a point has not been found, then keep looking
    // closest point must lie on an edge or vertex
    // gather and compare
    if(!inside){
      // add each face vertex
      for (int i=0;i<3;i++) {
        p.add(P[i][0],P[i][1],P[i][2]);
      }
      // add points of minimum distance between current vertex and each face edge
      getEdgeIntersection(v,P,p);
    }
  }
  // save closest point
  pC[0] = p.a;
  pC[1] = p.b;
  pC[2] = p.c;
}

bool Container::computeClosest(Face *f,Vertex *v,double &squareD,double vn[3]) {
  bool signal = false;
  // get closest point to vertex v on face f
  double pC[3];
  computePC(f,v,pC);
  // invert vertex normal if vertex is not nice
  double vn_copy[3];
  if (v->o->vertexIsNice(v)){vn_copy[0]=vn[0];vn_copy[1]=vn[1];vn_copy[2]=vn[2]; } 
  else 					  {vn_copy[0]=-vn[0];vn_copy[1]=-vn[1];vn_copy[2]=-vn[2]; } 
  // compute separation vector
  double sep_vec[3] = {pC[0]-v->pN[0],pC[1]-v->pN[1],pC[2]-v->pN[2]};
  double sign = dot(sep_vec,vn_copy);
  if (sign>0.0){
    // compute cosine of angle between outward normal and separation vector
    // which is equal to dot product of vectors divided by vector magnitudes
    //double cos_angle = sign/sqrt(dot(vn_copy,vn_copy)*dot(sep_vec,sep_vec));
    double cos_angle2 = (sign*sign)/dot(vn_copy,vn_copy)*dot(sep_vec,sep_vec);
    // closest point must be within a specified angle of vertex normal
    // to avoid grabbing points on same object
    //
    // Note all values involved are strictly positive so this optimization is valid
    if (cos_angle2>CLOSEST_POINT_COSINE*CLOSEST_POINT_COSINE){
      // compute square of separation distance
      double aa = pC[0]-v->pN[0];
      double bb = pC[1]-v->pN[1];
      double cc = pC[2]-v->pN[2];
      double temp=aa*aa+bb*bb+cc*cc;

      /*			// DEBUG
                                if(v->match(1766,"a183")==true){
                                cout << "\n\nContainer::computeClosest: current vertex:\n";
                                v->printVertex(v->o->name);
                                cout << endl;
                                cout << "Container::computeClosest: candidate closest face:\n";
                                f->printFace(f->v[0]->o->name);
                                cout << endl;
                                cout << "Container::computeClosest: "
                                << " p ["
                                << p.a << " "
                                << p.b << " "
                                << p.c << "]\n";
                                cout << "Container::computeClosest: "
                                << " temp=" << temp
                                << ", SEARCH_RADIUS_SQ=" << SEARCH_RADIUS_SQ
                                << ", SCALE=" << SCALE
                                << ", SEARCH_RADIUS_SQ*SCALE*SCALE=" << SEARCH_RADIUS_SQ*SCALE*SCALE
                                << endl;
                                }
                                */			// DEBUG

      // if square of separation distance is closer than square of SEARCH_RADIUS
      if (temp<(SEARCH_RADIUS_SQ*SCALE*SCALE)) {
        // if square of separation distance is less than
        // previously saved square of separation distance
        if (temp<squareD||!squareD) {
          // save
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
  cout << "Check edge angles..............................";
  cout.flush();
  // for each object* in container
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    // for each Edge* in object
    for (std::vector<Edge*>::iterator j=(*i)->e.begin();j!=(*i)->e.end();j++) {
      // record if angle is smallest so far
      checkAngle((*j)->getAngle()); 
    }
  }
  cout << "complete.\n";
  cout.flush();
}

void Vertex::assignHolding(double pH[3]){
  Vertex *vp;
  double d=0;
  ///// if no adjacent vertex has same position as pH /////
  std::vector<Edge*> e;
  getAdjacentEdges(e);
  // for each adjacent edge
  for(std::vector<Edge*>::iterator i=e.begin();i!=e.end();i++){
    Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
    (*i)->getVertices(v1,v2,o1,o2);
    // vp=edge vertex different from this vertex
    if (v1==this){vp=v2;}
    else {vp=v1;}
    // if non-self edge vertex is indistinguishable from pH, then displace pH
    if ( !distinguishable(vp->pN[0],pH[0])&&
         !distinguishable(vp->pN[1],pH[1])&&
         !distinguishable(vp->pN[2],pH[2])){ d=VERTEX_EPSILON*SCALE;}
  }
  pN[0]=pH[0]+d;
  pN[1]=pH[1]+d;
  pN[2]=pH[2]+d;
}

//bool Container::checkForIntersections(Vertex *v,Space &s,bool sw,hashtable_f &nb){
bool Container::checkForIntersections(Vertex *v,Space &s,hashtable_f &nb){
  // NOTE THIS FUNCTION IS NEVER CALLED WITH sw===false
  // INSTEAD getFaceIntersection IS CALLED FROM computeIntersectionForce
  //	std::vector<Face*> dummy;
  // for each adjacent face of current vertex
  for(std::vector<Face*>::iterator i=v->f.begin();i!=v->f.end();i++){
    //  if face intersects any other face, then return true
    if((*i)->getFaceIntersectionCheck(this,s,nb)==true){return true;}
    //		if(sw)	{if((*i)->getFaceIntersectionCheck(this,s,nb)==true){return true;}}
    //		else	{if((*i)->getFaceIntersection(this,false,dummy,s)==true){return true;}}
  }
  // no adjacent faces of current vertex intersect any other face
  return false;
}

void Container::collectEdgeAngles(Vertex *v,Monitor& stats){
  stats.e_angle.clear();
  ///// collect set of all unique edge angles of adjacent faces /////
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

void Face::updateBoxes(hashtable_f &ob,hashtable_f &nb){
  // this = adjacent face of current vertex
  // ob = multimap of current vertex adjacent faces and the boxes in which they perviously lay
  // nb = multimap of current vertex adjacent faces and the boxes in which they now lie

  ///// copy box*s in which this face previously lay to vector /////
  std::vector<Box*> ovec;
  ovec.reserve(VECTOR_RESERVE);
  std::pair<tf_iterator,tf_iterator> pp=ob.equal_range(this);
  // for each box* in old list (ob)
  for(tf_iterator i=pp.first;i!=pp.second;i++){
    ovec.push_back((*i).second);
  }

  ///// remove box* in common between old and new lists /////
  pp=nb.equal_range(this);
  // for each box* in new list (nb)
  tf_iterator i=pp.first;
  while(i!=pp.second){
    // if new box* found in old list (ovec)
    std::vector<Box*>::iterator qq=find(ovec.begin(),ovec.end(),(*i).second);
    // if found, then remove box* from old and new list
    if(qq!=ovec.end()){
      ovec.erase(qq);
      tf_iterator k = i;k++;
      nb.erase(i);
      i=k;
    }
    else {i++;}
  }

  ///// remove face* from remaining box* in old list /////
  if (!ovec.empty()){
    // for each remaining Box* in ovec (the old list)
    for (std::vector<Box*>::iterator j=ovec.begin();j!=ovec.end();j++){
      // sort Face* vector
      sort((*j)->f.begin(),(*j)->f.end());
      // look for this face in Box* face list
      std::pair<std::vector<Face*>::iterator,std::vector<Face*>::iterator> q;
      q=equal_range((*j)->f.begin(),(*j)->f.end(),this);
      // if found, then remove
      if(q.first!=q.second){
        (*j)->f.erase(q.first);
      }
    }
  }

  ///// add face* to remaining box* in new list /////
  pp=nb.equal_range(this);
  // for each box* in new list (nb)
  for(i=pp.first;i!=pp.second;i++){
    // add this face to Box* face vector
    (*i).second->f.push_back(this);
  }

}
/*
    void Face::getBoundingBox(double bb[6]){
// assume bb=bounding box = [xmn xmax ymin ymax zmin zmax]
// for each vertex of face
bb[0]=v[0]->pN[0];
bb[1]=v[0]->pN[0];
bb[2]=v[0]->pN[1];
bb[3]=v[0]->pN[1];
bb[4]=v[0]->pN[2];
bb[5]=v[0]->pN[2];
for(int j=1;j<3;j++){
if      (v[j]->pN[0]<bb[0]){bb[0]=v[j]->pN[0];}
else if (v[j]->pN[0]>bb[1]){bb[1]=v[j]->pN[0];}
if      (v[j]->pN[1]<bb[2]){bb[2]=v[j]->pN[1];}
else if (v[j]->pN[1]>bb[3]){bb[3]=v[j]->pN[1];}
if      (v[j]->pN[2]<bb[4]){bb[4]=v[j]->pN[2];}
else if (v[j]->pN[2]>bb[5]){bb[5]=v[j]->pN[2];}
}
}*/
/*
    void Face::getBoundingBox(double bb[6]){
// assume bb=bounding box = [xmn xmax ymin ymax zmin zmax]
// for each vertex of face
for(int j=0;j<3;j++){
if ( v[0]->pN[j] < v[1]->pN[j] ) {
if ( v[0]->pN[j] < v[2]->pN[j] ) {
if (v[1]->pN[j] < v[2]->pN[j]) {
bb[2*j]=v[0]->pN[j];
bb[2*j+1]=v[2]->pN[j];
} else {
bb[2*j]=v[0]->pN[j];
bb[2*j+1]=v[1]->pN[j];
}
} else {
bb[2*j]=v[2]->pN[j];
bb[2*j+1]=v[1]->pN[j];
}
} else {
if ( v[0]->pN[j] < v[2]->pN[j] ) {
bb[2*j]=v[1]->pN[j];
bb[2*j+1]=v[2]->pN[j];
} else {
if (v[1]->pN[j] < v[2]->pN[j]) {
bb[2*j]=v[1]->pN[j];
bb[2*j+1]=v[0]->pN[j];
} else {
bb[2*j]=v[2]->pN[j];
bb[2*j+1]=v[0]->pN[j];
}
}
}
}
}
*/
void Face::getBoundingBox(void){
  // assume bb=bounding box = [xmn xmax ymin ymax zmin zmax]
  // for each vertex of face
  for(int j=0;j<3;j++){
    if ( v[0]->pN[j] < v[1]->pN[j] ) {
      if ( v[0]->pN[j] < v[2]->pN[j] ) {
        if (v[1]->pN[j] < v[2]->pN[j]) {
          bb[2*j]=v[0]->pN[j];
          bb[2*j+1]=v[2]->pN[j];
        } else {
          bb[2*j]=v[0]->pN[j];
          bb[2*j+1]=v[1]->pN[j];
        }
      } else {
        bb[2*j]=v[2]->pN[j];
        bb[2*j+1]=v[1]->pN[j];
      }
    } else {
      if ( v[0]->pN[j] < v[2]->pN[j] ) {
        bb[2*j]=v[1]->pN[j];
        bb[2*j+1]=v[2]->pN[j];
      } else {
        if (v[1]->pN[j] < v[2]->pN[j]) {
          bb[2*j]=v[1]->pN[j];
          bb[2*j+1]=v[0]->pN[j];
        } else {
          bb[2*j]=v[2]->pN[j];
          bb[2*j+1]=v[0]->pN[j];
        }
      }
    }
  }
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
  bb[0] -= SPACE_LENGTH*SCALE;
  bb[1] += SPACE_LENGTH*SCALE;
  bb[2] -= SPACE_LENGTH*SCALE;
  bb[3] += SPACE_LENGTH*SCALE;
  bb[4] -= SPACE_LENGTH*SCALE;
  bb[5] += SPACE_LENGTH*SCALE;

  std::vector<Box*> bp;
  s.getBoxesFor3DLocations(bb,bp);
  stats.fs.clear();
  std::vector<Vertex*> t;
  t.reserve(VECTOR_RESERVE);
  // for each box
  for(std::vector<Box*>::iterator i=bp.begin();i!=bp.end();i++){
    // for each face in box
    for(std::vector<Face*>::iterator j=(*i)->f.begin();j!=(*i)->f.end();j++){
      // add face vertices to vector
      t.push_back((*j)->v[0]);
      t.push_back((*j)->v[1]);
      t.push_back((*j)->v[2]);
    }
  }
  // sort and keep unique vertices
  sort(t.begin(),t.end());
  std::vector<Vertex*>::iterator j = unique(t.begin(),t.end());
  t.assign(t.begin(),j);
  // for each vertex in vector
  for(j=t.begin();j!=t.end();j++){
    // if face vertex is not frozen, has a closest point and
    // closest point lies on an adjacent face of active vertex, v
    // then add vertex to affected vertex set
    if(	(binary_search(frozen.begin(),frozen.end(),*j)==false)
        && (*j)->cl!=NULL && binary_search(v->f.begin(),v->f.end(),(*j)->cl) ){
      stats.fs.insert(*j);
    }
  }
}


void Container::getAffectedVerticesAndEdgesAfter(Space &s,Vertex *v,double temp[3],Monitor &stats){
  double bb[6] = {temp[0],temp[0],temp[1],temp[1],temp[2],temp[2]};
  /*	METHOD 1 
        v->getBoundingBox(bb);
        bb[0] -= BIG_BOX*SCALE;
        bb[1] += BIG_BOX*SCALE;
        bb[2] -= BIG_BOX*SCALE;
        bb[3] += BIG_BOX*SCALE;
        bb[4] -= BIG_BOX*SCALE;
        bb[5] += BIG_BOX*SCALE;
        */
  /*	METHOD 2
        v->getBoundingBox(bb);
        bb[0] -= SPACE_LENGTH*SCALE;
        bb[1] += SPACE_LENGTH*SCALE;
        bb[2] -= SPACE_LENGTH*SCALE;
        bb[3] += SPACE_LENGTH*SCALE;
        bb[4] -= SPACE_LENGTH*SCALE;
        bb[5] += SPACE_LENGTH*SCALE;
        */
  /*	METHOD 3
        bb[0] = v->pN[0];
        bb[1] = v->pN[0];
        bb[2] = v->pN[1];
        bb[3] = v->pN[1];
        bb[4] = v->pN[2];
        bb[5] = v->pN[2];
        bb[0] -= BIG_BOX*SCALE;
        bb[1] += BIG_BOX*SCALE;
        bb[2] -= BIG_BOX*SCALE;
        bb[3] += BIG_BOX*SCALE;
        bb[4] -= BIG_BOX*SCALE;
        bb[5] += BIG_BOX*SCALE;
        */
  /*	METHOD 4 */
  bb[0] = v->pN[0];
  bb[1] = v->pN[0];
  bb[2] = v->pN[1];
  bb[3] = v->pN[1];
  bb[4] = v->pN[2];
  bb[5] = v->pN[2];
  bb[0] -= SPACE_LENGTH*SCALE;
  bb[1] += SPACE_LENGTH*SCALE;
  bb[2] -= SPACE_LENGTH*SCALE;
  bb[3] += SPACE_LENGTH*SCALE;
  bb[4] -= SPACE_LENGTH*SCALE;
  bb[5] += SPACE_LENGTH*SCALE;
  /**/
  std::vector<Box*> bp;
  stats.ps.clear();
  s.getBoxesFor3DLocations(bb,bp);
  //	cout << bp.size() << endl;
  // for each box
  for(std::vector<Box*>::iterator i=bp.begin();i!=bp.end();i++){
    // for each face in box
    for(std::vector<Face*>::iterator j=(*i)->f.begin();j!=(*i)->f.end();j++){
      // for each vertex of face
      for(int k=0;k<3;k++){
        // if face vertex is not frozen
        //        if(	binary_search(frozen.begin(),frozen.end(),(*j)->v[k])==false ){
        // store face vertices
        stats.ps.push_back((*j)->v[k]);
        //        }
      }
    }
  }
  // sort and keep unique 
  sort(stats.ps.begin(),stats.ps.end());
  std::vector<Vertex*>::iterator ii;
  ii = unique(stats.ps.begin(),stats.ps.end());
  stats.ps.assign(stats.ps.begin(),ii);

  // collect edges from current vertex adjacent faces
  stats.ae.clear();
  for(std::vector<Face*>::iterator i=v->f.begin();i!=v->f.end();i++){
    // store face edges
    stats.ae.push_back((*i)->e[0]);
    stats.ae.push_back((*i)->e[1]);
    stats.ae.push_back((*i)->e[2]);
  }
  // sort and keep unique 
  sort(stats.ae.begin(),stats.ae.end());
  std::vector<Edge*>::iterator j;
  j = unique(stats.ae.begin(),stats.ae.end());
  stats.ae.assign(stats.ae.begin(),j);
}

void Vertex::getAdjacentEdges(std::vector<Edge*> &vec){
  vec.clear();
  // for each adjacent face
  for(std::vector<Face*>::iterator i=f.begin();i!=f.end();i++){
    // for each face edge
    for(int j=0;j<3;j++){
      Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
      (*i)->e[j]->getVertices(v1,v2,o1,o2);
      // if either edge vertex is current vertex
      if(v1==this || v2==this){
        // add edge to vector
        vec.push_back((*i)->e[j]);
      }
    }
  }
}

void Monitor::initTable(Vertex *v,Space &s){
  std::vector<Box*> vbp;
  // get vertex adjacent faces
  //	hashset_f af;
  //	v->getAdjacentFaces(af);
  ob.clear();
  // for each adjacent face
  for(std::vector<Face*>::iterator i=v->f.begin();i!=v->f.end();i++){
    // get boxes that overlap this face's bounding box
    vbp.clear();
    s.computeBoxesToCheck(*i,vbp);
    // add box*s to multimap ob
    for(std::vector<Box*>::iterator j=vbp.begin();j!=vbp.end();j++){
      ob.insert(std::make_pair(*i,*j));
    }
  }
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

  // for each collected vertex requiring full closest point search
  for(hv_iterator kk=stats.fs.begin();kk!=stats.fs.end();kk++){
    // if vertex is not equal to current vertex and not frozen
    if(*kk!=v && (binary_search(frozen.begin(),frozen.end(),*kk)==false)){
      // if vertex currently has no closest point
      if((*kk)->cl==NULL){
        // then update neighborhood
        (*kk)->o->buildNeighborhood(*kk);
      }
      // screen vertices same as current vertex
      if(*kk==v){
        cout << "Container::updateClosest: current vertex would be processed twice!!\n";
        exit(0);
      }
      // then update collected vertex's closest point
      // true => remove existing matching element, if any
      findClosest(s,*kk,stats,true);
      // update global energy
      // false -> do not compute force, hence dummy
      double e=(*kk)->getSeparationForceEnergy(dummy,false,this);
      energy=energy-(*kk)->energy+e;
      (*kk)->energy=e;
    }
  }
  // for each collected vertex not requiring full search
  //for(hv_iterator kk=stats.ps.begin();kk!=stats.ps.end();kk++){
  for(std::vector<Vertex*>::iterator kk=stats.ps.begin();kk!=stats.ps.end();kk++){
    // if vertex not already processed, i.e. found in fs hashset,
    // not equal to current vertex and not frozen
    if(stats.fs.count(*kk)==0 && *kk!=v
       && 
       (binary_search(frozen.begin(),frozen.end(),*kk)==false)
      ){
      ///// check all adjacent faces to see /////
      ///// if affected vertex's closest point has changed /////
      double n[3];
      (*kk)->getNormal(n);
      bool gate=false;
      double squareD = 0.0;
      // compute current square of separation distance for collected vertex
      if((*kk)->cl!=NULL){ squareD = (*kk)->getSqSepDist(this); }
      // for each adjacent face of vertex *v
      for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
        // if face vertices are not in current vertex neighborhood
        if (!faceIsAdjacent(*k,*kk)){
          // if the closest point to current vertex was found on this face
          if(computeClosest(*k,*kk,squareD,n)){ gate=true; }
        }
      }
      // if new closest point was found
      if(gate){
        // update global energy
        // false -> do not compute force, hence dummy
        double e=(*kk)->getSeparationForceEnergy(dummy,false,this);
        energy=energy-(*kk)->energy+e;
        (*kk)->energy=e;
        // update sets with vertex squared virtual displacement
        stats.updateSets(*kk,getVertexSqD(*kk,stats.gain),true);
      }
    }
  }
}


void Container::removeOldIntersections(Vertex *v,hashset_v &vset2){
  // collect the vertices of all intersecting faces
  // of all adjacent faces to current vertex
  // into vset, i.e. hash_set of Vertex*s.
  //
  // for every adjacent face
  for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
    /*		// DEBUG
                cout << "\nContainer::removeOldIntersections:BEFORE "
                << "ti=" << ti 
                << ", si=" << si 
                << endl;
                */		// DEBUG
    // if adjacent face has intersecting faces
    if((*k)->faceInTable_intf()){
      // for each intersecting face of adjacent face
      std::vector<Face*> *fv=(*k)->getIntersectingFaces();
      for(std::vector<Face*>::iterator p=(*fv).begin();p!=(*fv).end();p++){
        // add intersecting face vertices to vset2
        vset2.insert((*p)->v[0]);
        vset2.insert((*p)->v[1]);
        vset2.insert((*p)->v[2]);
        // DEBUG
        /*				cout << "Container::removeOldIntersections: "
                                        << "adjface="
                                        << (*k)->v[0]->o->name << "->" << (*k)->index
                                        << ", intface="
                                        << (*p)->v[0]->o->name << "->" << (*p)->index
                                        << endl;
                                        */				// DEBUG
        // remove adjacent face from intersecting face's vector
        (*p)->removeFaceFromVector(*k,this);
        // remove intersecting face from adjacent face's vector
        //				(*k)->removeFaceFromVector(*p,this);
        /*				// DEBUG
                                        cout << "\nContainer::removeOldIntersections:MIDDLE "
                                        << "ti=" << ti 
                                        << ", si=" << si 
                                        << endl;
                                        */				// DEBUG
      }
    }
    //		// clear adjacent face's intersecting face vector
    (*k)->clearFaceFromTable_intf(this);
    // update adjacent face intersect force, i.e. set force to zero
    (*k)->clearFaceFromTable_iv();
    /*		// DEBUG
                cout << "\nContainer::removeOldIntersections:AFTER "
                << "ti=" << ti 
                << ", si=" << si 
                << endl;
                */		// DEBUG
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

void Container::updateNewIntersections(Vertex *v,hashset_v &vset2,Space &s){
  // for every adjacent face
  for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
    // if adjacent face is currently intersected
    if ((*k)->computeIntersectionForce(this,s)){
      // for each intersecting face of adjacent face
      std::vector<Face*> *fv=(*k)->getIntersectingFaces();
      for(std::vector<Face*>::iterator p=(*fv).begin();p!=(*fv).end();p++){
        // add intersecting face vertices to vset2
        vset2.insert((*p)->v[0]);
        vset2.insert((*p)->v[1]);
        vset2.insert((*p)->v[2]);
        // if intersected face was previously not intersected
        // then create new face vector in table for intersected face
        if(!(*p)->faceInTable_intf()){(*p)->addFaceToTable_intf();}
        // add adjacent face* to intersected face
        (*p)->addFaceToVector(*k,this);
        // update intersected face intersection force
        //				(*p)->calculateIntersectionForce(this);
      }
    }
  }
}

void Container::updateAdjacentFaceBoxes(Vertex *v,Monitor &stats){
  // for each adjacent face, update Box*s
  for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
    (*k)->updateBoxes(stats.ob,stats.nb);
  }
}

void Container::updateAdjacentFaceBoundingBoxes(Vertex *v){
  // for each adjacent face, update Box*s
  for(std::vector<Face*>::iterator k=v->f.begin();k!=v->f.end();k++){
    (*k)->getBoundingBox();
  }
}

void Container::getNiceCheckSet(Vertex *v,Monitor &stats,Space &s){
  // store vertices for which niceness may have changed
  stats.set_nice_check.clear();
  // add current vertex
  stats.set_nice_check.insert(v);
  //	stats.set_nice_check.insert(v);
  /*	// DEBUG
        cout << "\nContainer::getNiceCheckSet:0000 "
        << "ti=" << ti 
        << ", si=" << si 
        << endl;
        */	// DEBUG
  removeOldIntersections(v,stats.set_nice_check);
  /*	// DEBUG
        cout << "\nContainer::getNiceCheckSet:1111 "
        << "ti=" << ti 
        << ", si=" << si 
        << endl;
        */	// DEBUG
  updateNewIntersections(v,stats.set_nice_check,s);
  /*	// DEBUG
        cout << "\nContainer::getNiceCheckSet:2222 "
        << "ti=" << ti 
        << ", si=" << si 
        << endl;
        */	// DEBUG
  /*	// DEBUG
        if(v->match(20606,"healed_cut_er")==true){
        int TI = ti;
        int SI = si;
        ti=si=0;
        computeFaceIntersectionForce(s);
        exit(0);
        }
        */	// DEBUG
}

void Container::getNiceSet(Space &s,Monitor &stats){
  // for each vertex whose niceness may have changed
  for(hv_iterator i=stats.set_nice_check.begin();i!=stats.set_nice_check.end();i++){
    // update niceness
    if(checkNiceness(s,*i)){
      // if niceness changed
      // i.e. was nice and now not,
      // or was not nice and now is
      // then add vertex to set 
      // requiring full search for closest point
      stats.fs.insert(*i);
    }
  }
}

double Container::getVertexSqD(Vertex *v,double gain){
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
  cout << "Compute global energy..........................";
  cout.flush();
  energy=0;
  double dummy[3] = {0.0,0.0,0.0};
  // for each object in container
  for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++){
    // for each vertex in object
    for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++){
      // false -> do not compute force, hence dummy
      if((*j)->cl!=NULL){
        (*j)->energy=(*j)->getSeparationForceEnergy(dummy,false,this);
        energy+=(*j)->energy;
      }
    }
    // for each edge in object
    for(std::vector<Edge*>::iterator j=(*i)->e.begin();j!=(*i)->e.end();j++){
      // choice of v1 is arbitrary, v2 could have been used
      // false -> do not compute force, hence dummy

      Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
      (*j)->getVertices(v1,v2,o1,o2);
      (*j)->energy=(*j)->getStretchForceEnergy(v1,dummy,false,this)
            + (*j)->getForceEnergy(0,dummy,false,false);
      energy+=(*j)->energy;
    }
  }
  cout << "complete.\n";
  cout.flush();
}

void Container::updateMovedVertexEnergy(Vertex *v,Monitor &stats){
  double dummy[3] = {0.0,0.0,0.0};
  ///// vertex contribution /////
  // false -> do not compute force, hence dummy
  double e=(v)->getSeparationForceEnergy(dummy,false,this);
  energy=energy-(v)->energy+e;
  (v)->energy=e;
  ///// edge contribution /////
  // for each affected edge
  for(std::vector<Edge*>::iterator i=stats.ae.begin();i!=stats.ae.end();i++){
    Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
    (*i)->getVertices(v1,v2,o1,o2);
    // choice of v1 is arbitrary, v2 could have been used
    // false -> do not compute force, hence dummy
    e=(*i)->getStretchForceEnergy(v1,dummy,false,this)
          +(*i)->getForceEnergy(0,dummy,false,false);
    energy=energy-(*i)->energy+e;
    (*i)->energy=e;
  }
}

void Container::updateVertexVD(Vertex *v, Monitor &stats){
  // collect adjecnt edges
  std::vector<Edge*> e;
  v->getAdjacentEdges(e);
  // set for storing nearby vertices to update
  v_set nearby;
  // for each collected edge from current vertex adjacent faces
  for(std::vector<Edge*>::iterator i=stats.ae.begin();i!=stats.ae.end();i++){
    // if edge is not adjacent to current vertex
    if(
       find(e.begin(),e.end(),*i)==e.end()
      ){
      Vertex *v1,*v2,*o1,*o2;
      (*i)->getVertices(v1,v2,o1,o2);
      // insert all four edge vertices (v1,v2,o1,o2) into set
      nearby.insert(v1);
      nearby.insert(v2);
      nearby.insert(o1);
      nearby.insert(o2);
    }
  }
  // update collected vertices including current vertex
  // for each collected vertex
  for(vs_iterator i=nearby.begin();i!=nearby.end();i++){
    // if vertex has a closest point, then update sets
    if((*i)->cl!=NULL){
      stats.updateSets(*i,getVertexSqD(*i,stats.gain),true);
    }
  }
}

//bool Container::assignNewVertexCoords(Space &s,Vertex *v,double pH[3],Monitor &stats) {
bool Container::assignNewVertexCoords(Space &s,Vertex *v,double pH[3],Monitor &stats,bool &int_flag,bool &angle_flag) {
  /*	// DEBUG
        cout << "\nContainer::assignNewVertexCoords:0000 "
        << "ti=" << ti 
        << ", si=" << si 
        << endl;
        */	// DEBUG
  //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
  // build old face* hashtable
  stats.initTable(v,s);
  //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
  // clear new face* hashtable
  stats.nb.clear();
  // store current vertex position
  double pO[3]={v->pN[0],v->pN[1],v->pN[2]};
  // collect vertices in a local region
  // whose closest point is on adjacent face of current vertex
  getAffectedVerticesAndEdgesBefore(s,v,pO,stats);
  //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
  // collect edge angles
  collectEdgeAngles(v,stats);
  // set current position to holding position
  v->assignHolding(pH);
  // add to previous collection of vertices the collection 
  // of all vertices in the same local region as before
  // also collect edges of adjacent faces of current vertex
  getAffectedVerticesAndEdgesAfter(s,v,pO,stats);
  //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
  // store new vertex position
  //	double pR[3]={v->pN[0],v->pN[1],v->pN[2]};
  // move vertex back to original position
  //	v->pN[0]=pO[0];v->pN[1]=pO[1];v->pN[2]=pO[2];
  // move vertex back to new position
  //	v->pN[0]=pR[0];v->pN[1]=pR[1];v->pN[2]=pR[2];
  // if no faces intersect and no edges have small angles
  /*	// DEBUG
        cout << "\nContainer::assignNewVertexCoords:1111 "
        << "ti=" << ti 
        << ", si=" << si 
        << endl;
        */	// DEBUG
  //	bool int_flag = (checkForIntersections(v,s,stats.nb)==false || INTERSECTION_WEIGHT==100.0);
  int_flag = (checkForIntersections(v,s,stats.nb)==false || INTERSECTION_WEIGHT==100.0);
  //	bool int_flag = (checkForIntersections(v,s,true,stats.nb)==false || INTERSECTION_WEIGHT==100.0);
  //	bool angle_flag = !checkForSmallAngles(stats);
  angle_flag = !checkForSmallAngles(stats);
  if (int_flag && angle_flag){
    /*		// DEBUG
                cout << "\nContainer::assignNewVertexCoords:2222 "
                << "ti=" << ti 
                << ", si=" << si << endl;
                cout << "\nContainer::assignNewVertexCoords:2222 "
                << "pO ["
                << pO[0] << " "
                << pO[1] << " "
                << pO[2] << "]"
                << ", pH ["
                << pH[0] << " "
                << pH[1] << " "
                << pH[2] << "]"
                << endl;
                */		// DEBUG
    //	cout << "\nAAAA\n";
    //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
    // for each adjacent face, update Box*s
    updateAdjacentFaceBoxes(v,stats);
    //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
    // for each adjacent face, update bounding box
    updateAdjacentFaceBoundingBoxes(v);
    //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
    // scan each adjacent face and update global minimum edge angle
    updateEdgeAngles(v);
    //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
    // store vertices for which niceness may have changed
    // i.e. vertices of faces that were intersected but now aren't, and vice versa
    /*		// DEBUG
                cout << "\nContainer::assignNewVertexCoords:3333 "
                << "ti=" << ti 
                << ", si=" << si 
                << endl;
                */		// DEBUG
    getNiceCheckSet(v,stats,s);
    //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
    // detect niceness changes and collect changed vertices
    /*		// DEBUG
                cout << "\nContainer::assignNewVertexCoords:4444 "
                << "ti=" << ti 
                << ", si=" << si 
                << endl;
                */		// DEBUG
    getNiceSet(s,stats);
    //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
    // update closest point and global energy 
    // for affected vertices collected before
    updateClosest(s,v,stats);
    //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
    // update sets with squared virtual displacement of nearby vertices
    updateVertexVD(v,stats);
    //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
    // update global energy due to this vertex
    updateMovedVertexEnergy(v,stats);
    //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
    // clear face* hashtable
    /*		// DEBUG
                cout << "\nContainer::assignNewVertexCoords:5555 "
                << "ti=" << ti 
                << ", si=" << si 
                << endl;
                */		// DEBUG
    return true;
  } else { // move vertex back
    //	cout << "\nBBBB\n";
    //	stats.vertexSearchInTopN(this,103544,"d000_FILTERED_SMOOTH_SMOOTH");
    //DEBUG
    //		if(int_flag==false){ cout << "moving vertex " << v->index << " caused intersecting faces.\n"; }
    //		if(angle_flag==false){ cout << "moving vertex " << v->index << " caused small angles.\n"; }
    // DEBUG
    v->pN[0]=pO[0];
    v->pN[1]=pO[1];
    v->pN[2]=pO[2];
    /*		// DEBUG
                cout << "\nContainer::assignNewVertexCoords:6666 "
                << "ti=" << ti 
                << ", si=" << si 
                << endl;
                */		// DEBUG
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
      /*			fvec[0]=P[0]/L*INTERSECTION_WEIGHT/(SEPARATION_WEIGHT+EDGE_STRETCH_WEIGHT
                                +ANGLE_STRETCH_WEIGHT+INTERSECTION_WEIGHT) ;
                                fvec[1]=P[1]/L*INTERSECTION_WEIGHT/(SEPARATION_WEIGHT+EDGE_STRETCH_WEIGHT
                                +ANGLE_STRETCH_WEIGHT+INTERSECTION_WEIGHT) ;
                                fvec[2]=P[2]/L*INTERSECTION_WEIGHT/(SEPARATION_WEIGHT+EDGE_STRETCH_WEIGHT
                                +ANGLE_STRETCH_WEIGHT+INTERSECTION_WEIGHT) ;
                                */

      fvec[0]=P[0]/L*INTERSECTION_WEIGHT/100.0;
      fvec[1]=P[1]/L*INTERSECTION_WEIGHT/100.0;
      fvec[2]=P[2]/L*INTERSECTION_WEIGHT/100.0;
      // set intersection force
      this->addForceToFace(fvec);
    }
  }
}


bool Face::computeIntersectionForce(Container *c,Space &s) {
  //	std::vector<Face*> dummy;
  //	if(getFaceIntersection(c,false,dummy,s)){
  if(getFaceIntersection(c,s)){
    calculateIntersectionForce(c);
    return true;
  }
  return false;
}

void Container::computeFaceIntersectionForce(Space &s) {
  cout << "Compute face intersection force................";
  cout.flush();
  int ii=1;
  double max=face_count;
  double goal = 0.2;
  cout << "0%..";
  cout.flush();
  // for each object
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    // for each face in object
    for (std::vector<Face*>::iterator j=(*i)->f.begin();j!=(*i)->f.end();j++) {
      (*j)->computeIntersectionForce(this,s);
      // track progress
      double progress = static_cast<double>(ii++)/max;
      if(progress>goal){
        cout << static_cast<int>(goal*100) << "%..";
        cout.flush();
        goal+=0.2;
      }
    }
  }
  cout << "100%..complete.\n";
  cout.flush();
}

// #####################################################
// #####################################################

//bool Face::getFaceIntersectionCheck(Container *c,Space &s,hashtable_f &nb) {
//  // get boxes that overlap with this face's bounding box
//  std::vector<Box*> bp;
//  s.computeBoxesToCheck(this,bp);
//  // add box*s to multimap nb
//  for(std::vector<Box*>::iterator i=bp.begin();i!=bp.end();i++){
//    nb.insert(std::make_pair(this,*i));
//  }
//  // collect all unique faces in chosen boxes
//  f_set of;
//  for(std::vector<Box*>::iterator i=bp.begin();i!=bp.end();i++){
//    of.insert((*i)->f.begin(),(*i)->f.end());
//  }
//  // for each unique face
//  for(fs_iterator i=of.begin();i!=of.end();i++){
//    // if unique face is not same as current face
//    if(*i!=this){
//      // if faces intersect
//      if (c->checkFaceFaceIntersections(this,*i)) { 
//        if(STRICT_FACE_INTERSECTION_PREVENTION==true){
//          // if this face previously had no intersections, then reject move
//          if(!this->faceInTable_intf()){return true;}
//          // else face was previously intersected
//          else {
//            // grab previous intersecting faces of this face
//            std::vector<Face*> *fv=this->getIntersectingFaces();
//            bool same=false;
//            // for each previous intersecting face
//            for(std::vector<Face*>::iterator p=(*fv).begin();p!=(*fv).end();p++){
//              // if current intersecting face is of same object
//              //  as previous intersecting face
//              if((*p)->v[0]->o==(*i)->v[0]->o){same=true;}
//            }
//            // if intersecting objects have changed from before to now
//            // then reject move
//            if(!same){return true;}
//          }
//        } else {
//          // if current face if of same object as intersecting face, then reject move
//          if(this->v[0]->o==(*i)->v[0]->o){return true;}
//        }
//      }
//    }
//  }
//  return false;
//}

bool Face::getFaceIntersectionCheck(Container *c,Space &s,hashtable_f &nb) {
  // get boxes that overlap with this face's bounding box
  std::vector<Box*> bp;
  s.computeBoxesToCheck(this,bp);
  // add box*s to multimap nb
  for(std::vector<Box*>::iterator i=bp.begin();i!=bp.end();i++){
    nb.insert(std::make_pair(this,*i));
  }
  // collect all unique faces in chosen boxes
  std::vector<Face*> of;
  of.reserve(VECTOR_RESERVE);
  for(std::vector<Box*>::iterator i=bp.begin();i!=bp.end();i++){
    of.insert(of.end(),(*i)->f.begin(),(*i)->f.end());
  }
  // sort and keep unique faces
  sort(of.begin(),of.end());
  std::vector<Face*>::iterator j;
  j = unique(of.begin(),of.end());
  of.assign(of.begin(),j);

  // for each unique face
  for(std::vector<Face*>::iterator i=of.begin();i!=of.end();i++){
    // if unique face is not same as current face
    if(*i!=this){
      // if faces intersect
      if (c->checkFaceFaceIntersections(this,*i)) { 
        if(STRICT_FACE_INTERSECTION_PREVENTION==true){
          // if this face previously had no intersections, then reject move
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
            // if intersecting objects have changed from before to now
            // then reject move
            if(!same){return true;}
          }
        } else {
          // if current face if of same object as intersecting face, then reject move
          if(this->v[0]->o==(*i)->v[0]->o){return true;}
        }
      }
    }
  }
  return false;
}

//bool Face::getFaceIntersection(Container *c,bool just_check,std::vector<Face*> &int_f,Space &s) {
bool Face::getFaceIntersection(Container *c,Space &s) {
  std::vector<Face*>::iterator j;
  bool flag = false;
  // reset face element in table
  //	clearFaceFromTable_intf();
  // get boxes that overlap this face's bounding box
  std::vector<Box*> bp;
  bp.clear();
  s.computeBoxesToCheck(this,bp);
  // for each box in which face lies, add faces in box to vector
  std::vector<Face*> of;
  of.reserve(VECTOR_RESERVE);
  for(std::vector<Box*>::iterator i=bp.begin();i!=bp.end();i++){
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
        // save intersecting face* to this face's intersecting face vector
        addFaceToVector(*j,c);
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
        // compute new vertex coords
        double ppH[3]; // new holding position coordinates (x,y,z)
        // load energy map
        (*i)->computeNewCoords(c,ppH,gain);
        // compute virtual displacement
        double vd = (ppH[0]-(*i)->pN[0])*(ppH[0]-(*i)->pN[0])+
              (ppH[1]-(*i)->pN[1])*(ppH[1]-(*i)->pN[1])+
              (ppH[2]-(*i)->pN[2])*(ppH[2]-(*i)->pN[2]);
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

void Monitor::printVertexSelect(Container &c,const int group){
  int zero=0;
  char file[FILENAME_SIZE];
  // DEBUG
  int first=-100000,second=-100000;
  Vertex *v1=NULL,*v2=NULL;
  // DEBUG
  // open log file
  sprintf(file,"%s%s.%d",OUTPUT_DATA_DIR.c_str(),VERTEX_SELECTION_FILE,group);
  std::ofstream this_file (file);
  for(std::vector<Object*>::iterator i=c.o.begin();i!=c.o.end();i++){
    for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++){
      vhm_iterator t = touch_map.find(*j);
      // if vertex was touched
      if(t!=touch_map.end()){
        this_file << touch_map[*j] << endl;
        // DEBUG
        if (touch_map[*j]>first){
          second=first;
          v2=v1;
          first=touch_map[*j];
          v1=*j;
        }
        // DEBUG
      } else {
        this_file << zero << endl;
      }
    }
  }
  this_file.close();
  // DEBUG
  /*	cout << "\nFIRST MOST MOVED: " << first << " times\n";
        v1->printVertex(v1->o->name);
        cout << endl;
        if(v2!=NULL){
        cout << "\nSECOND MOST MOVED: " << second << " times\n";
        v2->printVertex(v2->o->name);
        cout << endl;
        }
        */	// DEBUG
}

// IDENTICAL TO void Monitor::printVertexSelect(Container &c,int interval){
// EXCEPT FOR int TO char* ARGUMENT SWAP
void Monitor::printVertexSelect(Container &c,char *str){
  int zero=0;
  char file[FILENAME_SIZE];
  // open log file
  sprintf(file,"%s%s.%s",OUTPUT_DATA_DIR.c_str(),VERTEX_SELECTION_FILE,str);
  std::ofstream this_file (file);
  for(std::vector<Object*>::iterator i=c.o.begin();i!=c.o.end();i++){
    for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++){
      vhm_iterator t = touch_map.find(*j);
      // if vertex was touched
      if(t!=touch_map.end()){
        this_file << touch_map[*j] << endl;
      } else {
        this_file << zero << endl;
      }
    }
  }
  this_file.close();
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

void Monitor::punishIfOverTouched(Vertex *v){
  // if vertex touched too much
  if(touch_map[v]>MAX_TOUCHES){punishN(v,static_cast<int>(GROUP_SIZE));}
}

bool Monitor::gainReachedSteadyState(void){
  return distinguishable(avg_old,avg_new,ENERGY_EPSILON)==false;
}

void Monitor::updateAvg(double e){
  // if energy window is NOT full
  if(num<ENERGY_WINDOW){
    // sum = old sum + new value
    sum=sum+e;
    num++;
  } else {
    // energy window is full
    // sum = old sum - oldest value in window + new value
    sum=sum-*window+e;
  }
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
  num=0;
  window=begin;
  for(int i=0;i<ENERGY_WINDOW;i++){
    *(window++)=0.0;
  }
  window=begin;
}

void Monitor::initAvg(void){
  sum=avg_new=avg_old=0.0;
  num=0;
  window = new double[ENERGY_WINDOW];
  begin=window;
  for(int i=0;i<ENERGY_WINDOW;i++){
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
/*
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
*/
//void Monitor::punishGate(Vertex *v,int time_out){
//	gate[v]=time_out;
//}

void Monitor::punishN(Vertex *v,int time_out){
  pun_n[v]=time_out;
}

void Monitor::punishInt(Vertex *v,int time_out){
  pun_int[v]=time_out;
}

void Monitor::punishCom(Vertex *v,int time_out){
  pun_com[v]=time_out;
}

void Monitor::punishAng(Vertex *v,int time_out){
  pun_ang[v]=time_out;
}

void Monitor::initRefrac(void){
  refrac_s.clear();
  refrac_l.clear();
}

bool Monitor::Refracted(Vertex *v){
  if (refrac_s.find(v)!=refrac_s.end()){return true;}
  else							{return false;}
}

void Monitor::updateRefractoryWindow(Vertex *v){
  ///// erase oldest element /////
  if(refrac_l.size()==REFRACTORY_PERIOD){
    refrac_s.erase(refrac_l.front());
    refrac_l.pop_front();
  }
  ///// store new element /////
  refrac_l.push_back(v);
  refrac_s.insert(v);
}

int Monitor::getMaxVertexMoves(void){
  std::vector<Vertex*> v;
  v.assign(refrac_l.begin(),refrac_l.end());
  sort(v.begin(),v.end());
  int c=1,max=1;
  // for each pair (i,j) of sequential vertices in vector
  std::vector<Vertex*>::iterator i=v.begin(),j;
  j=i;j++;
  while(j!=v.end()){
    // if pair are same then increment current counter, c
    if(*i==*j){ c++; }
    else { // else pair are different
      // if current counter, c, bigger than max counter, max
      // then replace max with c
      if(c>max){ max=c; }
      // reset current counter
      c=1;
    }
    i++;j++;
  }
  return max;
}

// #####################################################
// #####################################################

void Container::checkFrozenAndNoCP(void){
  // for each object
  for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++){
    // for each vertex in object
    for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++){
      // if vertex is frozen and cl
      if( (binary_search(frozen.begin(),frozen.end(),*j)==true)
          && (*j)->cl!=NULL){
        cout << "Error: vertex is frozen with cl!=NULL.\n";
        (*j)->printVertex((*j)->o->name);
        cout << endl;
        exit(0);
      }
    }
  } 
}

void Monitor::vertexSearchInTopN(Container *c,int index,std::string str){
  // get Object*
  Object *oo=NULL;
  for(std::vector<Object*>::iterator i=c->o.begin();i!=c->o.end();i++){
    // if object name matches str
    if((*i)->name==str){
      // save Object*
      oo=*i;
      break;
    }
  }
  // if matching Object* not found
  if(oo==NULL){cout << "Matching Object* not found.\n";exit(0);}

  // for each vertex in object
  Vertex *vv=NULL;
  for(std::vector<Vertex*>::iterator i=oo->v.begin();i!=oo->v.end();i++){
    // if indices match
    if((*i)->index==index){
      // save Vertex*
      vv=*i;
      break;
    }
  }
  // if matching Vertex* not found
  if(vv==NULL){cout << "Matching Vertex* not found.\n";exit(0);}

  // look for Vertex* in topN
  // for each pair in topN
  bool found=false;
  for(tv_iterator i=topN.begin();i!=topN.end();i++){
    // if Vertex*s match 
    if((*i).second==vv){
      found=true;
      cout << "\n\nMonitor::vertexSearchInTopN: "
            << "Error! The following vertex is found in topN map (se->Vertex*).\n";
      vv->printVertex(vv->o->name);cout << endl;
      cout << endl << endl;
      exit(0);
    }
  }
  if(found==false){
    cout << "Monitor::vertexSearchInTopN: "
          << str << "->" << index << " NOT found in topN.\n";
  }
}

void Container::writeObjectData(void){
  cout << "Writing object data to log files...............";
  cout.flush();
  writeObjectList();
  statusFileInit();
  cout << "complete.\n";
  cout.flush();
}

void Container::writeSeparationDistances(void){
  writeDistances(0); // ok
  // UNDO ME
  writeDistancesOneSidedHausdorff(); // ok
  writeDistancesTwoSidedHausdorff(); // ok
  writeDistancesOneSidedHausdorff_noself(); //ok
  writeDistancesTwoSidedHausdorff_noself(); // ok
  writeDistancesNOCP(0);
  // UNDO ME
}

void Container::computeMeanEdgeLengths(void){
  cout << "Compute mean edge length.......................";
  cout.flush();
  // for each object
  for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++){
    double sum=0;
    // for each edge in object
    for(std::vector<Edge*>::iterator j=(*i)->e.begin();j!=(*i)->e.end();j++){
      sum+=(*j)->l;
    }
    (*i)->mean_edge_length=sum/(*i)->e.size();
  } 
  cout << "complete.\n";
  cout.flush();
}

void Container::sortAdjacentFaces(void){
  // for each object
  for(std::vector<Object*>::iterator i=o.begin();i!=o.end();i++){
    // for each vertex in object
    for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++){
      // sort adjacent faces
      sort((*j)->f.begin(),(*j)->f.end());
    }
  } 
}

void Container::printFrozenCP(void)
{
  cout << "x_coordinate y_coordinate z_coordinate state_value x_normal y_normal z_normal" << endl;
  for(std::vector<Vertex*>::iterator i=frozen.begin();i!=frozen.end();i++)
  {
    (*i)->printVertexCP();
    cout << endl;
  }
}
