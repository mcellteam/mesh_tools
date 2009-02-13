double parse(int argc,char **argv,std::string message,std::string &filename)
{
  // if no arguments passed
  if(argc==1){
    cout << message << endl;
    exit(0);
  }

  int c;
  opterr=0;
  char *eptr=NULL;
  double threshold = -1.0;
  while((c=getopt(argc,argv,"ht:")) != -1)
    switch(c)
    {
      case 't':
        // specify max separation threshold
        threshold = strtod(optarg,&eptr);
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
    filename = argv[optind];
  }
  else
  {
    fprintf (stderr,"No input file argument found on command line.\n");
    exit(0);
  }
  return threshold;
}

void Object::scanFile(const char *filename)
{
  char line[2048],*str;
  FILE *F;
  Vertex *vv;
  Face *ff;
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

Edge* findEdge(Vertex* va,Vertex* vb,hashtable_t &hm,int num_digits){
  Edge *ee=NULL;
  std::string s = keyPair(va->index,vb->index,num_digits);
  // if element exists given key, then get Edge pointer
  if (hm.count(s)>0){ ee=hm[s]; }
  return ee;
}

void Edge::update(Face *f)
{
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

void Object::createEdge(Face *ff,Vertex* va,Vertex* vb)
{
  // new edge
  Edge *en = new Edge(ff,va,vb);
  // store edge pointer in hash table
  hm[keyPair(va->index,vb->index,num_digits)]=en;
  // add edge pointer to face
  ff->addEdge(en);
  // add edge pointer to object
  e.push_back(en);
}

void Object::buildEdge(Face *ff,Vertex *va,Vertex *vb)
{
  Edge *ee=NULL;
  ee=findEdge(va,vb,hm,num_digits);
  if(ee!=NULL){ee->update(ff);}
  else {createEdge(ff,va,vb);}
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
  num_digits = setNumDigits();
  // for each face
  for (f_iterator i=f.begin();i!=f.end();i++) {
    buildEdge(*i,(*i)->v[0],(*i)->v[1]);
    buildEdge(*i,(*i)->v[1],(*i)->v[2]);
    buildEdge(*i,(*i)->v[2],(*i)->v[0]);
  }
}

double Object::findLongestEdge(void)
{
  double max=-1.0;
  // for each edge
  for(e_iterator i=e.begin();i!=e.end();i++)
  {
    // compute squared edge length
    double d=(*i)->getLengthSq();
    if (d>max){max=d;}
  }
  return max; // return 
}


void Object::gatherFreeVertices(void)
{
  free_vertices.clear();
  // for each edge
  for(e_iterator i=e.begin();i!=e.end();i++)
  {
    // if edge is open
    if ((*i)->f2==NULL && (*i)->fvec.empty()==true)
    {
      free_vertices.push_back((*i)->v1);
      free_vertices.push_back((*i)->v2);
    }
  }
  // keep unique vertices 
  sort(free_vertices.begin(),free_vertices.end());
  v_iterator new_end = unique(free_vertices.begin(),free_vertices.end());
  free_vertices.assign(free_vertices.begin(),new_end);

  fprintf(stderr,"Border vertices found...........%d\n",free_vertices.size());
  fflush(stderr);

  // print free vertices
  if(PRINT==true)
  {
    for (v_iterator i=free_vertices.begin();i!=free_vertices.end();i++)
    {
      (*i)->printVertex();
    }
  }
}

double getSq3DDistance(double a[3],double b[3])
{
  double x = a[0]-b[0];
  double y = a[1]-b[1];
  double z = a[2]-b[2];
  return x*x+y*y+z*z;
}

bool dSortComp(const Distance lhs,const Distance rhs)
{
  return (lhs.d > rhs.d);
}

void Object::computeDistances(const double user_threshold,double max_edge_length_sq)
{
  fprintf(stderr,"Computing vertex distances......");fflush(stderr);
  // for each free vertex
  for (v_iterator i=free_vertices.begin();i!=free_vertices.end();i++)
  {
    // for each free vertex
    for (v_iterator j=free_vertices.begin();j!=free_vertices.end();j++)
    {
      // if combination is unique
      if ((*i)->index > (*j)->index)
      {
        // compute squared 3D distance between vertices
        double d = getSq3DDistance((*i)->p,(*j)->p);
        // if distance less than longest edge and user threshold, then save
        if( d < max_edge_length_sq && d < user_threshold*user_threshold )
        {
          // if vertices are not members of the same face
          // i.e. no edge defined by these two vertices
          if(findEdge(*i,*j,hm,num_digits)==NULL)
          {
            distances.push_back(Distance(d,*i,*j));
          }
        }
      }
    }
  }

  // sort distances by separation distance
  sort(distances.begin(),distances.end(),dSortComp);

  if(PRINT==true){
    for (d_iterator i=distances.begin();i!=distances.end();i++)
    {
      fprintf(stderr,"Free vertices %i and %i, sq dist = %.15g\n",
              (*i).vA->index,
              (*i).vB->index,
              (*i).d);
    }
  }
  fprintf(stderr,"complete.\n");
  fflush(stderr);
}

void Object::fixVertices(Vertex *vbad)
{
  if(PRINT==true) fprintf(stderr,"Fixing vertices.................");
  if(PRINT==true) fflush(stderr);
  std::pair<std::vector<Vertex*>::iterator,std::vector<Vertex*>::iterator> i;
  // look for bad vertex in vertex vector
  i = equal_range(v.begin(),v.end(),vbad);
  // if found, then remove
  if(i.first != i.second){v.erase(i.first);}
  if(PRINT==true) fprintf(stderr,"complete.\n");
  if(PRINT==true) fflush(stderr);
}

void Object::fixFaces(Vertex *vgood,Vertex *vbad)
{
  if(PRINT==true) fprintf(stderr,"Fixing faces....................");
  if(PRINT==true) fflush(stderr);
  // for each face
  for (f_iterator i=f.begin();i!=f.end();i++)
  {
    // replace all instances of second vertex in face list with first vertex
    for(int j=0;j<3;j++)
    {
      if(vbad==(*i)->v[j]){(*i)->v[j]=vgood;}
    }
  }
  if(PRINT==true) fprintf(stderr,"complete.\n");
  if(PRINT==true) fflush(stderr);
}

void Object::purgeEdgeFromMap(Edge *ee)
{
  // for each map element
  for(ht_iterator i=hm.begin();i!=hm.end();i++)
  {
    // if element references edge
    if((*i).second==ee)
    {
      // then erase element
      hm.erase(i);
    }
  }
}

void Object::purgeMap(void)
{
  fprintf(stderr,"Purging map elements............");
  fflush(stderr);
  // for each map element
  for(ht_iterator i=hm.begin();i!=hm.end();i++)
  {
    Edge *ee=(*i).second;
    // if neither edge vertex is free
    if(binary_search(free_vertices.begin(),free_vertices.end(),ee->v1)==false &&
       binary_search(free_vertices.begin(),free_vertices.end(),ee->v2)==false)
    {
      // then erase element
      hm.erase(i);
    }
  }
  fprintf(stderr,"complete.\n");
  fflush(stderr);
}

void Object::purgeEdges(void)
{
  fprintf(stderr,"Purging edges...................");
  fflush(stderr);
  // for each edge
  Edge *ee=NULL;
  e_iterator i=e.begin();
  while(i!=e.end())
  {
    // if neither edge vertex is free
    if(binary_search(free_vertices.begin(),free_vertices.end(),(*i)->v1)==false &&
       binary_search(free_vertices.begin(),free_vertices.end(),(*i)->v2)==false)
    {
      // then erase element
      ee=*i;
      i=e.erase(i);
      delete ee;
    }
    else{i++;}
  }
  fprintf(stderr,"complete.\n");
  fflush(stderr);
}

void Object::purgeEdgeFromFaces(Edge *ee)
{
  if(ee->f1->e[0]==ee){ee->f1->e[0]=NULL;}
  if(ee->f1->e[1]==ee){ee->f1->e[1]=NULL;}
  if(ee->f1->e[2]==ee){ee->f1->e[2]=NULL;}
  if(ee->f2!=NULL)
  {
    if(ee->f2->e[0]==ee){ee->f2->e[0]=NULL;}
    if(ee->f2->e[1]==ee){ee->f2->e[1]=NULL;}
    if(ee->f2->e[2]==ee){ee->f2->e[2]=NULL;}
  }
  for(f_iterator i=ee->fvec.begin();i!=ee->fvec.end();i++)
  {
    if((*i)->e[0]==ee){(*i)->e[0]=NULL;}
    if((*i)->e[1]==ee){(*i)->e[1]=NULL;}
    if((*i)->e[2]==ee){(*i)->e[2]=NULL;}
  }
}

void Object::updateEdges(Vertex *vgood,Vertex *vbad)
{
  if(PRINT==true) fprintf(stderr,"Update edges....................");
  if(PRINT==true) fflush(stderr);
  ///// remove edges containing good and bad vertex /////
  // for each edge
  Edge *ee=NULL;
  e_iterator i=e.begin();
  while(i!=e.end())
  {
    // if edge contains either vertex
    if((*i)->v1==vgood || (*i)->v1==vbad ||
       (*i)->v2==vgood || (*i)->v2==vbad)
    {
      // remove map elements containing this edge
      purgeEdgeFromMap(*i);
      // remove edge from adjacent faces
      purgeEdgeFromFaces(*i);
      // erase edge
      ee = *i;
      i=e.erase(i);
      delete ee;
    }
    else{i++;}
  }

  ///// add edges containing good vertex /////
  // for each face
  for (f_iterator j=f.begin();j!=f.end();j++)
  {
    // if any face vertex is good
    // then add as edge
    if ((*j)->v[0]==vgood)
    {
      buildEdge(*j,(*j)->v[0],(*j)->v[1]);
      buildEdge(*j,(*j)->v[2],(*j)->v[0]);
    }
    else if ((*j)->v[1]==vgood)
    {
      buildEdge(*j,(*j)->v[0],(*j)->v[1]);
      buildEdge(*j,(*j)->v[1],(*j)->v[2]);
    }
    else if ((*j)->v[2]==vgood)
    {
      buildEdge(*j,(*j)->v[1],(*j)->v[2]);
      buildEdge(*j,(*j)->v[2],(*j)->v[0]);
    }
  }
  if(PRINT==true) fprintf(stderr,"complete.\n");
  if(PRINT==true) fflush(stderr);
}

void Object::updateDistances(void)
{
  // for each Distance instance
  d_iterator i=distances.begin();
  while(i!=distances.end())
  {
    if(binary_search(free_vertices.begin(),free_vertices.end(),(*i).vA)==false ||
       binary_search(free_vertices.begin(),free_vertices.end(),(*i).vB)==false)
    {
      distances.erase(i);
    }
    else{i++;}
  }
}


void Object::printVerticesFaces(void)
{
  // write out vertices
  for(v_iterator i=v.begin();i!=v.end();i++)
  {
    Vertex *p=*i;
    printf("Vertex %i  %.15g %.15g %.15g\n",p->index,p->p[0],p->p[1],p->p[2]);
  }
  // write out faces
  for(f_iterator i=f.begin();i!=f.end();i++)
  {
    Face *p=*i;
    printf("Face %i  %i %i %i\n",p->index,p->v[0]->index,p->v[1]->index,p->v[2]->index);
  }
}
