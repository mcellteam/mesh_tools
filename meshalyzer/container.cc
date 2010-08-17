#include "container.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "box.h"
#include "controls.h"
#include "meshalyzer.h"
#include "object.h"
#include "space.h"

using std::cerr;
using std::cout;
using std::endl;

Container::Container (void)
:eg(NULL),o(),files(),num_files(0),area(),aspect_ratio(),edge_length(),
  edge_angle(),adjacent_face(),good_integrity(true),num_orph(0),num_mis(0),num_deg(0),num_bor(0),
  num_flip(0),num_nonman_e(0),num_nonman_v(0),num_obj(0),num_vert(0),
  num_face(0),num_edge(0),num_vol(0),num_sep(0),num_bou(0),num_indistin(0),
  num_dupl_v(0),num_dupl_f(0),num_clo_cc(0),num_clo_nn(0)
{
  // cumulative Stats
  setbb[0]=setbb[1]=setbb[2]=1E30;
  setbb[3]=setbb[4]=setbb[5]=-1E30;
  num_man[0]=num_man[1]=num_man[2]=0;
  num_cons[0]=num_cons[1]=num_cons[2]=0;
  num_out[0]=num_out[1]=num_out[2]=0;
  pairs[0][0] = 0;
  pairs[0][1] = 1;
  pairs[1][0] = 1;
  pairs[1][1] = 2;
  pairs[2][0] = 2;
  pairs[2][1] = 0;
}

int Container::countIntFace (void)
{
  int a=0;
  for (o_iterator i=o.begin();i!=o.end();i++)
  {
    a+=(*i)->getNumIntFace();
  }
  return a;
}

void Container::countOutward (int val[3])
{
  //			consistent	outward
  // val[0] 	true		true
  // val[1] 	true		false
  // val[2] 	false		true or false
  //	val[0]=0; val[1]=0; val[2]=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		if     ((*i)->consistent==true && (*i)->outward==true) {val[0]++;}
  //		else if ((*i)->consistent==true && (*i)->outward==false){val[1]++;}
  //		else {val[2]++;}
  //	}
  val[0]=num_out[0];
  val[1]=num_out[1];
  val[2]=num_out[2];
}

void Container::countConsistent (int val[3])
{
  //			manifold	consistent
  // val[0] 	true		true
  // val[1] 	true		false
  // val[2] 	false		true or false
  //	val[0]=0; val[1]=0; val[2]=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		if     ((*i)->manifold==true && (*i)->consistent==true) {val[0]++;}
  //		else if ((*i)->manifold==true && (*i)->consistent==false){val[1]++;}
  //		else {val[2]++;}
  //	}
  val[0]=num_cons[0];
  val[1]=num_cons[1];
  val[2]=num_cons[2];
}

void Container::countManifold (int val[3])
{
  //	int cc=0,nn=0; // cc=manifold, nn=not manifold
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		if ((*i)->manifold==true){cc++;}
  //		else {nn++;}
  //	}
  //	return std::make_pair(cc,nn);
  //			close		manifold
  // val[0] 	true		true			'manifold'
  // val[1] 	true		false			'nonmanifold'
  // val[2] 	false		true or false	'undefined'
  val[0]=num_man[0];
  val[1]=num_man[1];
  val[2]=num_man[2];
}

std::pair<int,int> Container::countClosed (void)
{
  //	int cc=0,nn=0; // cc=closed, nn=not closed
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		if ((*i)->closed==true){cc++;}
  //		else {nn++;}
  //	}
  //	return std::make_pair(cc,nn);
  return std::make_pair(num_clo_cc,num_clo_nn);
}

int Container::countDuplV (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->dupl_v_index.size();
  //	}
  //	return a;
  return num_dupl_v;
}

int Container::countDuplF (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->dupl_f_index.size();
  //	}
  //	return a;
  return num_dupl_f;
}

int Container::countIndistin (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->indistin_v.size();
  //	}
  //	return a;
  return num_indistin;
}

int Container::countBoundaries (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->num_bou;
  //	}
  //	return a;
  return num_bou;
}

int Container::countComponents (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->num_sep;
  //	}
  //	return a;
  return num_sep;
}

double Container::countArea (void)
{
  double a=0.0;
  for (o_iterator i=o.begin();i!=o.end();i++)
  {
    a+=(*i)->getAreaSum();
  }
  return a;
}

double Container::countVol (void)
{
  //	double a=0.0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->vol;
  //	}
  //	return a;
  return num_vol;
}

int Container::countEdge (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->e.size();
  //	}
  //	return a;
  return num_edge;
}

int Container::countFace (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->f.size();
  //	}
  //	return a;
  return num_face;
}

int Container::countVertex (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->v.size();
  //	}
  //	return a;
  return num_vert;
}

int Container::countObject (void)
{
  //	return o.size();
  return num_obj;
}

int Container::countNonmanV (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->nonman_v.size();
  //	}
  //	return a;
  return num_nonman_v;
}

int Container::countNonmanE (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->nonman_e.size();
  //	}
  //	return a;
  return num_nonman_e;
}

int Container::countFlipped (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->flipped.size();
  //	}
  //	return a;
  return num_flip;
}

int Container::countBorder (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->border.size();
  //	}
  //	return a;
  return num_bor;
}

int Container::countDegen (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->degen.size();
  //	}
  //	return a;
  return num_deg;
}

int Container::countMissing (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->missing_v.size();
  //	}
  //	return a;
  return num_mis;
}

int Container::countOrphan (void)
{
  //	int a=0;
  //	for (o_iterator i=o.begin();i!=o.end();i++){
  //		a+=(*i)->orphan.size();
  //	}
  //	return a;
  return num_orph;
}

Container::~Container (void)
{
  o_iterator i;
  for (i=o.begin();i!=o.end();i++){ delete *i; }
}

void Container::clear (void)
{
  o.clear();
}

int Container::checkEdgeEdgeIntersection (Face *cf,Face *of,bool share_edge_flag) 
{
  Controls & cs(Controls::instance());
  // cpvc = current_polygon_vertex_coordinates
  // opvc = other_polygon_vertex_coordinates
  // cv   = current_vertex
  // ov   = other_vertex
  Vertex *cv[2],*ov[2];
  // for each current face edge
  for (int i=0;i<3;i++) 
  {
    // for each other face edge
    for (int j=0;j<3;j++) 
    {
      cv[0] = cf->ptr_vertex(pairs[i][0]);
      cv[1] = cf->ptr_vertex(pairs[i][1]);
      ov[0] = of->ptr_vertex(pairs[j][0]);
      ov[1] = of->ptr_vertex(pairs[j][0]);
      // if the edges do not share a vertex
      if (cv[0]!=ov[0]&&cv[0]!=ov[1]&&cv[1]!=ov[0]&&cv[1]!=ov[1]) 
      {
        // and the edges are not parallel
        bool parallel_flag = false;
        double x[3],y[3];
        for (int k=0;k<3;k++) 
        {
          x[k] = cv[1]->getpN(k)-cv[0]->getpN(k);
          y[k] = ov[1]->getpN(k)-ov[0]->getpN(k);
        }
        double term1 = x[0]*y[0]+ x[1]*y[1]+ x[2]*y[2];
        if ( !distinguishable(term1*term1,
                              (x[0]*x[0]+ x[1]*x[1]+ x[2]*x[2])*
                              (y[0]*y[0]+ y[1]*y[1]+ y[2]*y[2])) )
        {
          parallel_flag = true;
        }
        if (!parallel_flag) 
        {
          // compute scalars
          double qDen = (cv[1]->getpN(0)-cv[0]->getpN(0))*(ov[1]->getpN(1)-ov[0]->getpN(1))
                -(cv[1]->getpN(1)-cv[0]->getpN(1))*(ov[1]->getpN(0)-ov[0]->getpN(0));
          double qNum = (cv[1]->getpN(0)-cv[0]->getpN(0))*(cv[0]->getpN(1)-ov[0]->getpN(1))
                -(cv[1]->getpN(1)-cv[0]->getpN(1))*(cv[0]->getpN(0)-ov[0]->getpN(0));
          double q = qNum/qDen;
          double uNum = ov[0]->getpN(0)-cv[0]->getpN(0)+q*(ov[1]->getpN(0)-ov[0]->getpN(0));
          double uDen = cv[1]->getpN(0)-cv[0]->getpN(0);
          if (fabs(qDen)>cs.get_double_epsilon() && fabs(uDen)>cs.get_double_epsilon()) 
          {
            double u = uNum/uDen;
            if ( 
                ((u > cs.get_double_epsilon() && u < (1.0-cs.get_double_epsilon())) && 
                 (q > cs.get_double_epsilon() && q < (1.0-cs.get_double_epsilon())) ) ||
                ( share_edge_flag && ((u > cs.get_double_epsilon() && u < (1.0-cs.get_double_epsilon())) || 
                                      (q > cs.get_double_epsilon() && q < (1.0-cs.get_double_epsilon()))) )
               ) 
            {
              return(1);
            }
          }
        }
      }
    }
  }
  return(0);
}

int Container::checkFaceEdgeIntersection (Face *cf,Face *of)
{
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
  for (int i=0;i<3;i++) 
  {
    if (!line_flag || !poly_flag) 
    {
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
  for (int i=0;i<3;i++) 
  {
    if (!line_flag || !poly_flag) 
    {
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

void Container::printBatch (Controls &cs)
{
  int num_intf = countIntFace();
  cout << "\n\n" << "/* ********************** BATCH RESULTS ********************** */\n";
  // intersecting faces
  if (num_intf==0)
  {
    cout << "# intersecting faces: none\n";
  }
  else 
  {
    cout << "# intersecting faces: " << num_intf << endl;
    //	if -p option, print offending
    if (cs.get_print_detailed_info()==1)
    {
      int j=1;
      // for each object
      for (o_iterator m=o.begin();m!=o.end();m++)
      {
        // for each intersected face
        for (ff_iterator i=(*m)->getFirstIntFace();i!=(*m)->getOnePastLastIntFace();i++)
        {
          cout << "# intersecting faces: intersected face " << j++ << endl;
          // print intersected face
          (*i).first->print(cout);
          // print intersecting faces
          for (f_iterator k=(*(*i).second).begin();k!=(*(*i).second).end();k++)
          {
            (*k)->print(cout);
            cout << endl;
          }
        }
      }
    }
  }
}

void Container::createEdges (void) 
{
  o_iterator i;
  // for each object, create edges
  for (i=o.begin();i!=o.end();i++) 
  {
    (*i)->createEdges();
  }
}

void Container::findVertexAdjacencies (void) 
{
  o_iterator i;
  // for each object, find vertex adjacencies
  for (i=o.begin();i!=o.end();i++) 
  {
    (*i)->findVertexAdjacencies();
  }
}

void Container::scanDir (const char *filename) 
{
  struct dirent *pent;			// pointer to dirent structure
  DIR *pdir = opendir(filename);	// pointer to a directory data structure
  if (!pdir) {printf("Error. Could not open %s.\n",filename);exit(1);}
  else { cerr << "\nFolder found " << filename << endl << endl;}
  while ((pent=readdir(pdir)))
  {
    // copy char array to string
    std::string str = pent->d_name;
    // if file of typ *.mesh
    std::string::size_type found = str.find(".mesh",0);
    // if found
    if (found != std::string::npos) 
    {
      // save filename
      files.push_back(str);
      // update index
      num_files++;
      // print file found to screen
      //			cout << "file found: " << str << "\n"; cout.flush();
    }
  }
  closedir(pdir);
  sort(files.begin(),files.end());
}

Object* Container::processFile (std::string filename) 
{
  // create new Object
  Object *obj = new Object(filename);
  // scan file
  scanFile(obj,filename);
  // check object contents
  if (obj->noVertices()==true || obj->noFaces()==true)
  {
    delete obj;
    cout << "\n Container::processFile: "
          << "no valid mesh object found in "
          << filename << ". Skipping file.\n";
    return NULL;
  }
  else 
  {
    // save Object* in container
    o.push_back(obj);
    // save first vertex* in object
    if (obj->noVertices()==false)
    {
      eg=obj->getFrontVertex();
    }
    else 
    {
      cout << "\nContainer::processFile: Error. "
            << "Object " << obj->getName()
            << " contains no vertices.\n";
      exit(1);
    }
    // return
    return obj;
  }
}

void Container::update (Object *oo)
{
  // update Container
  num_obj++;
  num_vert+=oo->getNumVertices();
  // DEBUG
  //	cout << "\nContainer::update: "
  //	<< "oo->f.size()=" << oo->f.size() << endl;
  // DEBUG
  num_face+=oo->getNumFaces();
  num_edge+=oo->getNumEdges();
  num_sep+=oo->getNumSep();
  //	if (oo->manifold==true){num_man_cc++;}
  //	else {num_man_nn++;}
  bool closed     = oo->isClosed();
  bool consistent = oo->isConsistent();
  bool outward    = oo->getOutward();
  bool manifold   = oo->isManifold();
  if (manifold==true){num_man[0]++;}
  else if (manifold==false && closed==true) {num_man[1]++;}
  else {num_man[2]++;}
  num_bou+=oo->getNumBoundaries();
  num_indistin+=oo->getNumIndistinVertices();
  if     (manifold==true && consistent==true) {num_cons[0]++;}
  else if (manifold==true && consistent==false){num_cons[1]++;}
  else {num_cons[2]++;}
  num_vol+=oo->getVolume();
  if (closed==true){num_clo_cc++;}
  else {num_clo_nn++;}
  //num_bor+=oo->border.size();
  num_bor+=oo->getNumBorderEdges();
  num_nonman_v+=oo->getNumNonmanVertices();
  num_nonman_e+=oo->getNumNonmanEdges();
  num_flip+=oo->getNumFlippedEdges();
  if     (consistent==true && outward==true) {num_out[0]++;}
  else if (consistent==true && outward==false){num_out[1]++;}
  else {num_out[2]++;}
  num_orph+=oo->getNumOrphanVertices();
  num_mis+=oo->getNumMissingVertices();
  num_deg+=oo->getNumDegenerateFaces();
  num_dupl_v+=oo->getNumDuplVertexIndices();
  num_dupl_f+=oo->getNumDuplFaceIndices();
}

void Container::scanFile (Object *obj,std::string filename) 
{
  char line[2048];
  // open file
  FILE *F = fopen(filename.c_str(),"r");
  if (!F) 
  {
    cout <<"Couldn't open input file "
          << filename << endl;
    exit(1);
  }
  else 
  {
    cerr << "\n\n" << "/* ********************** "
          << "OBJECT ********************** */\n";
    //	print object name 
    cerr << "name: " << obj->getName() << endl;
    cerr.flush();
    //		cout << "file found: " << filename << "\n"; cout.flush();
  }
  // for every line in file
  for (char *str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) 
  {
    // skip leading whitespace
    while (strchr(" \t,",*str)!=NULL) { str++;}
    // if first character is V for Vertex, add new linked list class instance
    if (strchr("V",*str)!=NULL)
    {
      Vertex *v=new Vertex(str,obj);
      // obj->v.push_back(v);
      obj->addVertex(v);
      //obj->vp.insert(std::make_pair(v->getIndex(),v));
      obj->addVertexIndexPair(v->getIndex(),v);
      //obj->found.insert(std::make_pair(v->getIndex(),false));
      obj->addIndexBoolPair(v->getIndex(),false);
    }
    // if first character is F for Face, add new linked list class instance
    else if (strchr("F",*str)!=NULL)
    {
      Face *f=new Face(str,obj);
      obj->addFace(f);
    }
  }
  fclose(F);
}

void Container::writeDistances (void) 
{
  char file[Controls::instance().get_filename_size()];
  // create output filename
  sprintf(file,"closest_point_distances.dat");
  // open output file
  std::ofstream newfile (file,std::ios::out);
  if (newfile.is_open())
  {
    newfile.precision(4);
    // for each object
    for (o_iterator i=o.begin();i!=o.end();i++) 
    {
      // for each vertex in object
      for (v_iterator j=(*i)->getFirstVertex();j!=(*i)->getOnePastLastVertex();j++) 
      {
        // if vertex has a closest face
        if ((*j)->ptr_closest_face()!=NULL)
        {
          // compute separation vector
          double s[3];
          for (int k=0;k<3;k++){ s[k]=(*j)->getpC(k)-(*j)->getpN(k); }
          // print separation distance
          if ((*i)->vertexIsNice(*j)){ newfile << sqrt(dot(s,s)) << endl;}
          else {newfile << -sqrt(dot(s,s)) << endl;}
        }
      }
    }
    newfile.close();
  }
}

void Container::boundWorld (Space &s) 
{
  o_iterator i;
  double xmin,xmax,ymin,ymax,zmin,zmax,range[6];
  //initialize mins and maxes
  //	xmin = o[0]->v[0]->getpN(0);
  //	xmax = o[0]->v[0]->getpN(0);
  //	ymin = o[0]->v[0]->getpN(1);
  //	ymax = o[0]->v[0]->getpN(1);
  //	zmin = o[0]->v[0]->getpN(2);
  //	zmax = o[0]->v[0]->getpN(2);
  xmin = eg->getpN(0);
  xmax = eg->getpN(0);
  ymin = eg->getpN(1);
  ymax = eg->getpN(1);
  zmin = eg->getpN(2);
  zmax = eg->getpN(2);
  ////////// loop through all objects //////////
  // for each object
  for (i=o.begin();i!=o.end();i++) 
  {
    // get range of object vertices
    // bounding box: [xmin,ymin,zmin][xmax,ymax,zmax]\n";
    (*i)->boundObject(range);
    if (range[0]<xmin) {xmin = range[0];}
    if (range[1]<ymin) {ymin = range[1];}
    if (range[2]<zmin) {zmin = range[2];}
    if (range[3]>xmax) {xmax = range[3];}
    if (range[4]>ymax) {ymax = range[4];}
    if (range[5]>zmax) {zmax = range[5];}
    if (xmin <setbb[0]) {setbb[0]=xmin;}
    if (ymin <setbb[1]) {setbb[1]=ymin;}
    if (zmin <setbb[2]) {setbb[2]=zmin;}
    if (xmax >setbb[3]) {setbb[3]=xmax;}
    if (ymax >setbb[4]) {setbb[4]=ymax;}
    if (zmax >setbb[5]) {setbb[5]=zmax;}
  }
  if (xmin<0) {s.setWorld(0,xmin*1.01);} else {s.setWorld(0,xmin*0.99);}
  if (xmax<0) {s.setWorld(1,xmax*0.99);} else {s.setWorld(1,xmax*1.01);}
  if (ymin<0) {s.setWorld(2,ymin*1.01);} else {s.setWorld(2,ymin*0.99);}
  if (ymax<0) {s.setWorld(3,ymax*0.99);} else {s.setWorld(3,ymax*1.01);}
  if (zmin<0) {s.setWorld(4,zmin*1.01);} else {s.setWorld(4,zmin*0.99);}
  if (zmax<0) {s.setWorld(5,zmax*0.99);} else {s.setWorld(5,zmax*1.01);}
}

void Container::getExtraRay (Vertex *v,double lp[2][3],int index) 
{
  Controls & cs(Controls::instance());
  // get normal info
  double n[3];
  v->f[index]->getNormal(n);
  // compute centroid of first adjacent face
  double cx = (v->f[index]->ptr_vertex(0)->getpN(0)+
               v->f[index]->ptr_vertex(1)->getpN(0)+
               v->f[index]->ptr_vertex(2)->getpN(0))/3.0;
  double cy = (v->f[index]->ptr_vertex(0)->getpN(1)+
               v->f[index]->ptr_vertex(1)->getpN(1)+
               v->f[index]->ptr_vertex(2)->getpN(1))/3.0;
  double cz = (v->f[index]->ptr_vertex(0)->getpN(2)+
               v->f[index]->ptr_vertex(1)->getpN(2)+
               v->f[index]->ptr_vertex(2)->getpN(2))/3.0;
  double L=sqrt( dot(n,n) );
  lp[0][0] = cx;
  lp[1][0] = lp[0][0]+n[0]/L*cs.get_ray_epsilon();
  lp[0][1] = cy;
  lp[1][1] = lp[0][1]+n[1]/L*cs.get_ray_epsilon();
  lp[0][2] = cz;
  lp[1][2] = lp[0][2]+n[2]/L*cs.get_ray_epsilon();
}

int Container::findExtraPoint (Space &s,Vertex *v,double p[3],int index)
{
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
  if (pp.first!=pp.second){crossed_faces.erase(pp.first);}
  findOddMeshes(crossed_faces,face_flags,num_odd_objects,tmp);
  p[0]=lp[1][0];
  p[1]=lp[1][1];
  p[2]=lp[1][2];
  if (!tmp.empty())
  {
    ppp=equal_range(tmp.begin(),tmp.end(),v->getObject());
    if (ppp.first!=ppp.second)
    {
      return 0;
    }
  }
  return 1;
}

void Container::findCrossed1 (Space &s,Vertex *v,double lp[2][3],vec_o &c)
{
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
  for (i=v->f.begin();i!=v->f.end();i++)
  {
    pp=equal_range(crossed_faces.begin(),crossed_faces.end(),*i);
    // if adjacent face is found in crossed_faces, then remove from crossed_faces
    if (pp.first!=pp.second){crossed_faces.erase(pp.first);}
  }
  findOddMeshes(crossed_faces,face_flags,num_odd_objects,c);		
}

void Container::findCrossed2 (Space &s,double lp[2][3],vec_o &c)
{
  // find and return crossed objects between pN and extracellular point
  int num_odd_objects;
  vec_i face_flags;
  vec_f crossed_faces,unique_faces;
  collectNiceFaces(s,lp,unique_faces);
  findIntersectedFaces(lp,unique_faces,crossed_faces,face_flags);
  findOddMeshes(crossed_faces,face_flags,num_odd_objects,c);		
}

void Container::collectCrossed (Space &s,Vertex *v,vec_o &cb)
{
  vec_o ca;
  std::pair<o_iterator,vec_o::iterator> pp;
  double lp[2][3],p[3];
  // find point, p, outside of current object
  unsigned int index = 0;
  while (!findExtraPoint(s,v,p,index))
  {
    index++;
    if (index>(v->f.size()-1))
    {
      cout << "\n\nEvery adjacent face failed!\n";
      cout << v->getObject()->getName() << "->" << v->getIndex() 
            << " current vertex [" << v->getpN(0) << " "
            << v->getpN(1) << " "
            << v->getpN(2) << "]"
            << endl;
      exit(1);
    }
  }
  // grab intersected objects between pN and RAY ORIGIN and return as ca
  lp[0][0]=v->getpN(0);
  lp[0][1]=v->getpN(1);
  lp[0][2]=v->getpN(2);
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
  if (v->getIndex()==531 && !strcmp(v->getObject()->getName().c_str(),"a001_FILTERED_SMOOTH_SMOOTH"))
  {
  }
  findCrossed2(s,lp,cb);
  // remove all meshes in ca from cb
  // since the object is not odd relative to pN
  sort(cb.begin(),cb.end());
  for (o_iterator ii=ca.begin();ii!=ca.end();ii++)
  {
    pp=equal_range(cb.begin(),cb.end(),*ii);
    // if object in ca is found in cb, then remove from cb
    if (pp.first!=pp.second)
    {
      cb.erase(pp.first);
    }
    // else add it to cb
    else
    {
      cb.push_back(*ii);
    }
  }
}

bool Container::updateNiceness (Vertex *v,vec_o &cb)
{
  std::pair<o_iterator,vec_o::iterator> pp;
  int old_nice = v->getVertexNiceness();
  // if vertex niceness changes then set flag=true
  bool flag = false;
  // if cb is not empty, then vertex is not nice
  if (!cb.empty()) 
  {
    v->setVertexNiceness(1);
    // if vertex was nice
    if (!old_nice)
    {
      flag = true;
      pp=equal_range(cb.begin(),cb.end(),v->getObject());
      // if vertex is inside self object
      if (pp.first!=pp.second)
      {
        v->setVertexNiceness(2);
        cout << endl << endl 
              << v->getObject()->getName() << "->" << v->getIndex() 
              << " vertex inside self [" << v->getpN(0) << " "
              << v->getpN(1) << " "
              << v->getpN(2)
              << "], cb.size() " << cb.size()
              << endl;
        for (f_iterator jj=v->nf.begin();jj!=v->nf.end();jj++)
        {
          cout << v->getObject()->getName() << "->" << (*jj)->getIndex() 
                << " adjacent face\n";
          (*jj)->print(cout);
          cout << endl;
        }
      }
    }
  }
  else
  {
    // else cb is empty, then vertex is nice
    // if vertex was nonnice, but not to self
    if (old_nice==1)
    {	
      flag=true;
    }
    // if vertex was at least nonnice to self
    else if (old_nice==2)
    {
      flag=true;
    }
    // update niceness
    v->setVertexNiceness(0);			
  }
  return flag;
}

bool Container::checkNiceness (Space &s,Vertex *v) 
{
  vec_o cb;
  // collect objects inside which vertex lies
  collectCrossed(s,v,cb);
  // update niceness of vertex based on cb
  return updateNiceness(v,cb);
}

void Container::findNice (Space &s) 
{
  cout << "Iteration 0: ";
  cout << "find nice vertices................";
  cout.flush();
  o_iterator i;
  v_iterator j;
  // for each object in container
  for (i=o.begin();i!=o.end();i++) 
  {
    // for each vertex in object
    for (j=(*i)->getFirstVertex();j!=(*i)->getOnePastLastVertex();j++) 
    {
      checkNiceness(s,*j);
    }
  }
  cout << "complete.\n";
  cout.flush();
}

void Container::findClosestAxis (Space &s,Vertex *v,double lp[2][3]) 
{
  // n=normal,r=ray
  // get normal info
  // identify nearest boundary
  double dis[6] = 
  { fabs(v->getpN(0)-s.getWorld(0)),
    fabs(v->getpN(0)-s.getWorld(1)),
    fabs(v->getpN(1)-s.getWorld(2)),
    fabs(v->getpN(1)-s.getWorld(3)),
    fabs(v->getpN(2)-s.getWorld(4)),
    fabs(v->getpN(2)-s.getWorld(5))
  };
  int i=0;
  double min = dis[i];
  for (int j=1;j<6;j++)
  {
    if (dis[j]<min){i=j;min=dis[j];}
  }
  // configure ray
  if      (i==0)
  {
    lp[1][0] = lp[0][0]-2*(s.getWorld(1)-s.getWorld(0));	// end x
    lp[1][1] = lp[0][1];				// end y
    lp[1][2] = lp[0][2];				// end z
  }
  else if (i==1)
  {
    lp[1][0] = lp[0][0]+2*(s.getWorld(1)-s.getWorld(0));	// end x
    lp[1][1] = lp[0][1];				// end y
    lp[1][2] = lp[0][2];				// end z
  }
  else if (i==2)
  {
    lp[1][0] = lp[0][0];				// end x
    lp[1][1] = lp[0][1]-2*(s.getWorld(3)-s.getWorld(2));	// end y
    lp[1][2] = lp[0][2];				// end z
  }
  else if (i==3)
  {
    lp[1][0] = lp[0][0];				// end x
    lp[1][1] = lp[0][1]+2*(s.getWorld(3)-s.getWorld(2));	// end y
    lp[1][2] = lp[0][2];				// end z
  }
  else if (i==4)
  {
    lp[1][0] = lp[0][0];				// end x
    lp[1][1] = lp[0][1];				// end y
    lp[1][2] = lp[0][2]-2*(s.getWorld(5)-s.getWorld(4));	// end z
  }
  else if (i==5)
  {
    lp[1][0] = lp[0][0];				// end x
    lp[1][1] = lp[0][1];				// end y
    lp[1][2] = lp[0][2]+2*(s.getWorld(5)-s.getWorld(4));	// end z
  }
}

void Container::getBoxIndexRange (Space &s,double lp[2][3],int br[6])
{
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

void Container::collectNiceFaces (Space &s,double lp[2][3],vec_f &uf) 
{
  int br[6];
  vec_b b;
  b_iterator i;
  f_iterator new_end;

  bool flag = false;
  if (
     !distinguishable(lp[0][0],2719.54642515,1E-8) &&
     !distinguishable(lp[0][1],4249.96388355,1E-8) &&
     !distinguishable(lp[0][2],6242.92495541,1E-8)
    ){flag=true;}

  // compute box index range that contains ray
  getBoxIndexRange(s,lp,br);

  if (flag)
  {
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
  for (i=b.begin();i!=b.end();i++) 
  {
    // for each face in box
    //for (j=(*i)->f.begin();j!=(*i)->f.end();j++) 
    for (c_f_iterator j=(*i)->first_face();j!=(*i)->one_past_last_face();j++) 
    {
      uf.push_back(*j);
    }
  }
  // keep unique faces
  sort(uf.begin(),uf.end());
  new_end = unique(uf.begin(),uf.end());
  uf.assign(uf.begin(),new_end);
}

void Container::findIntersectedFaces (double lp[2][3],std::vector<Face*> &uf,
                                      std::vector<Face*> &cf,std::vector<int> &ff) 
{
  // uf = unique_faces
  // cf = crossed_faces
  // ff = face_flags
  f_iterator j;
  bool line_flag, poly_flag, poly_edge_flag;
  // for each unique polygon
  for (j=uf.begin();j!=uf.end();j++) 
  {
    checkLineFaceIntersection(*j,lp,line_flag,poly_flag,poly_edge_flag);
    // does point intersect polygon
    if (line_flag && (poly_flag||poly_edge_flag)) 
    {
      // add polygon_index to crossed array
      cf.push_back(*j);
      if (Controls::instance().get_detect_polygon_edge_intersection()) 
      {
        // if intersection point falls on edge, signal with flag
        if ( poly_edge_flag )
        {
          ff.push_back(1);
        }
        else
        {
          ff.push_back(0);
        }
      }
    }
  }
}

void Container::findOddMeshes (vec_f & cf,vec_i & ff,int & num_odd_objects,vec_o & tmp) 
{
  // find mesh objects crossed an odd number of times by ray
  vec_o ol; // object index list
  f_iterator i,j;
  int sum,k=0,L,parity,count;
  // for each crossed face
  for (i=cf.begin();i!=cf.end();i++) 
  {
    if (Controls::instance().get_detect_polygon_edge_intersection()) 
    {
      // if ray intersected face edge
      if (ff[k++])
      {
        //look for another face in same object with edge crossed
        sum=0;
        L=k;
        for (j=i+1;j!=cf.end();j++) 
        {
          // if ray intersected face edge, and faces have same object index
          if (ff[L++] && ((*i)->getIndex()==(*j)->getIndex())){sum++;}
        }
        // compute parity of sum
        parity=sum%2;
        // if even, add instance to object list
        if (!parity) {ol.push_back((*i)->ptr_vertex(0)->getObject());}
      }
    }
    else {ol.push_back((*i)->ptr_vertex(0)->getObject());}
  }
  ///// sort object_index_list /////
  sort(ol.begin(),ol.end());
  ///// count odd objects /////
  num_odd_objects=0;
  count=ol.size();
  k=0;
  tmp.clear();
  while (k<count) 
  {
    if (k+1!=count) 
    {
      // skip identical pairs of object indices
      if (ol[k]==ol[k+1]) {k++;k++;}
      // odd object
      else 
      {
        tmp.push_back( ol[k]);
        num_odd_objects++;
        k++;
      }
    }
    else
    {
      // add remaining object to odd object list
      num_odd_objects++;
      tmp.push_back( ol[k]);
      k++;
    }
  }
}

bool Container::facesParallel (Face *cf,Face *of)
{
  // are current face and other face parallel
  // i.e. is angle between normals equal to zero?
  // i.e. is the square of the cosine of the angle equal to 1?
  double cn[3],on[3];
  // get face normals
  cf->getNormal(cn);
  of->getNormal(on);
  double term1 = cn[0]*on[0]+ cn[1]*on[1]+ cn[2]*on[2];
  if (!distinguishable(term1*term1,
                       (cn[0]*cn[0]+cn[1]*cn[1]+cn[2]*cn[2])*(on[0]*on[0]+on[1]*on[1]+on[2]*on[2]))) 
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool Container::facesColinear(Face *cf,Face *of)
{
  // cpvi = current_polygon_vertex_indices
  // opvi = other_polygon_vertex_indices
  double on[3],*opvc[3],*cpvc[3];
  int x=2,y=2;
  // get face vertex coordinates
  cf->getVertexCoordinates(cpvc);
  of->getVertexCoordinates(opvc);
  // get face vertex indices
  int cpvi[3] = {cf->ptr_vertex(0)->getIndex(),cf->ptr_vertex(1)->getIndex(),cf->ptr_vertex(2)->getIndex()};
  int opvi[3] = {of->ptr_vertex(0)->getIndex(),of->ptr_vertex(1)->getIndex(),of->ptr_vertex(2)->getIndex()};
  // get face normal
  of->getNormal(on);
  // dot product of other face normal and 
  // line connecting point on other face to point on current face
  // try to choose vertices not shared by both face
  if      ((opvi[0]!=cpvi[0])&&(opvi[0]!=cpvi[1])&&(opvi[0]!=cpvi[2])) {x=0;}
  else if ((opvi[1]!=cpvi[0])&&(opvi[1]!=cpvi[1])&&(opvi[2]!=cpvi[2])) {x=1;}
  if      ((cpvi[0]!=opvi[0])&&(cpvi[0]!=opvi[1])&&(cpvi[0]!=opvi[2])) {y=0;}
  else if ((cpvi[1]!=opvi[0])&&(cpvi[1]!=opvi[1])&&(cpvi[1]!=opvi[2])) {y=1;}
  // if polygons colinear, then normal and line are orthogonal
  // and dot product will be ~zero
  if (fabs(on[0]*(opvc[x][0]-cpvc[y][0])
           +on[1]*(opvc[x][1]-cpvc[y][1])
           +on[2]*(opvc[x][2]-cpvc[y][2]))<Controls::instance().get_double_epsilon()) {return true;}
  else {return false;}
}

int Container::numUniqueVertices (Face *cf,Face *of,int single_shared_vert[2])
{
  // how many vertices are shared between current and other face?
  int num_unique = 0;
  // for each current face vertex
  for (int j=0;j<3;j++) 
  {
    for (int k=0;k<3;k++) 
    {
      if (cf->ptr_vertex(j)!=of->ptr_vertex(k)) {num_unique++;} else {single_shared_vert[0]=j;single_shared_vert[1]=k;}
    }
  }
  return num_unique;
}

bool Container::checkFaceFaceIntersections (Face *cf,Face *of) 
{
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
  if      (num_unique == 8) {share_vert_flag = true;}
  else if (num_unique == 7) {share_edge_flag = true;}
  else if (num_unique == 6) {identical_flag = true;}
  ////////// begin decision tree //////////
  if (coplanar_flag) 
  {
    if (identical_flag) 
    {
      // polygons are identical
      return true;
    }
    else 
    {
      //do polygon edges intersect?	
      // if yes, intersect
      // if no, do not intersect
      if (checkEdgeEdgeIntersection(cf,of,share_edge_flag)) {return true;}
      else {return false;}
    }
  }
  else 
  {
    if (share_vert_flag) 
    {
      int m=0,n=1,p=0,q=1;
      // single vertex shared
      if      (single_shared_vert[0]==0){m=1;n=2;}
      else if (single_shared_vert[0]==1){n=2;}
      if      (single_shared_vert[1]==0){p=1;q=2;}
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
    }
    else if (!share_edge_flag) 
    {
      // do faces intersect?
      if (checkFaceEdgeIntersection(cf,of)) {return true;}
      else {return false;}
    }
    else
    {
      return false;
    }
  }
}

void Container::assignFacesToBoxes (Space &s) 
{
  o_iterator i;
  f_iterator j;
  vec_b bp;
  // clear boxes
  s.clearBoxes();
  ////////// identify in which boxes each face exists ////////
  // for each object in container
  for (i=o.begin();i!=o.end();i++) 
  {
    // for each face in object
    for (j=(*i)->getFirstFace();j!=(*i)->getOnePastLastFace();j++) 
    {
      // identify boxes  that overlap to any degree the face bounding box.
      // This conservative approach, i.e. always including the actual 
      // face-intersecting boxes plus others is meant to save time,
      // assuming the time to exclude the other boxes is greater than the
      // time to check intersection with other boxes later.
      bp.clear();
      s.computeBoxesToCheck(*j,bp);
      // check	
      if (bp.empty())
      { 
        cout << "ERROR: NO BOXES FOR\n" << "Face " << (*j)->getIndex() << " " 
              << ((*j)->ptr_vertex(0))->getIndex() << " " << ((*j)->ptr_vertex(1))->getIndex() << " "
              << ((*j)->ptr_vertex(2))->getIndex() << endl;
        exit(1);
      }
      // record face in boxes class
      s.recordFace(bp,*j);
    }
  }
}

void Container::getSeparationDistances (Space &s)
{
  cout << "Iteration 0: ";
  cout << "get separation distances..........";
  cout.flush();
  o_iterator i;
  v_iterator j;
  // for each object in container
  for (i=o.begin();i!=o.end();i++) 
  {
    // for each vertex in object
    for (j=(*i)->getFirstVertex();j!=(*i)->getOnePastLastVertex();j++) 
    {
      ////////// find closest point to current vertex //////////
      // false => just add value to table, do not touch existing elements
      //			there should be none anyway
      findClosest(s,*j);
    }
  }
  cout << "complete.\n";
  cout.flush();
}

void Container::getBoxes (std::vector<Box*> &bp,Vertex *v,int offset,Space &s)
{
  int cbi[3],br[6];
  // compute box index that contains current_vertex
  cbi[0] = s.location2Index(v->getpN(0),"x");
  cbi[1] = s.location2Index(v->getpN(1),"y");
  cbi[2] = s.location2Index(v->getpN(2),"z");
  // box_range
  br[0]=cbi[0]-offset;
  br[1]=cbi[0]+offset;
  br[2]=cbi[1]-offset;
  br[3]=cbi[1]+offset;
  br[4]=cbi[2]-offset;
  br[5]=cbi[2]+offset;
  // handle case where vertex lies on subspace boundary
  if (!(br[0]*s.getSpaceLength()-v->getpN(0))){br[0]--;}
  if (!(br[2]*s.getSpaceLength()-v->getpN(1))){br[2]--;}
  if (!(br[4]*s.getSpaceLength()-v->getpN(2))){br[4]--;}
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

bool Container::faceInNeighborhood (Face *f,Vertex *v)
{
  // if face is in different object than vertex, then return false
  if (f->ptr_vertex(0)->getObject()!=v->getObject()){return false;}
  // else if face and vertex are in same object
  std::pair<f_iterator,vec_f::iterator> p
        =equal_range(v->nf.begin(),v->nf.end(),f);
  // if face is in current vertex neighborhood
  if (p.first!=p.second){return true;}
  return false;
}

void Container::getCandidateFaces (std::vector<Box*> &bp,Vertex *v,hset_f &cf)
{
  // for each box in search
  for (b_iterator j=bp.begin();j!=bp.end();j++) 
  {
    // for each face in box
    for (c_f_iterator k=(*j)->first_face();k!=(*j)->one_past_last_face();k++) 
    {
      // if face vertices are not in current vertex neighborhood, face is candidate
      if (!faceInNeighborhood(*k,v))
      {
        cf.insert(*k);
      }
    }
  }
}

bool Container::findClosest (Space &s,Vertex *v) 
{
  bool gate = false;
  // declare pair
  std::pair<double,double> p;
  // get vertex normal
  double n[3];
  v->getNormal(n);
  // get Box pointers for Face* collection
  vec_b bp;
  getBoxes(bp,v,Controls::instance().get_num_adjacent_boxes(),s);
  // collect Face pointers
  hset_f cf;
  getCandidateFaces(bp,v,cf);
  // if candidate faces were found
  if (!cf.empty())
  {
    double squareD=0.0;
    // for each candidate face
    for (hf_iterator j=cf.begin();j!=cf.end();j++) 
    {
      // if the closest point to current vertex was found on this face
      if (computeClosest(*j,v,squareD,n)){ gate=true;}
    }
    // if no closest point was found
    if (!gate)
    {
      // reset pointer to closest face
      v->setClosestFacePtr(NULL);
      // set closest point to current vertex location
      for (int i=0;i<3;i++) {v->setpC(i,v->getpN(i));}
    }
  }
  return gate;
}

bool Container::getPlaneIntersection (Face *f,Vertex *v,double *n,double num,double den,Point &p)
{
  //compute point on face plane that is closest to current vertex
  // i.e. intersection of plane with face normal through current vertex
  bool line_flag, poly_flag=false,poly_edge_flag;
  double u=num/den;
  double lp[2][3] = {{v->getpN(0),v->getpN(1),v->getpN(2)},
    {v->getpN(0)+u*n[0],v->getpN(1)+u*n[1],v->getpN(2)+u*n[2]}};
  checkLineFaceIntersection(f,lp,line_flag,poly_flag,poly_edge_flag);
  // if intersection point is on face,then save point
  if (poly_flag) {p.add(lp[1][0],lp[1][1],lp[1][2]);}
  return (poly_flag || poly_edge_flag);
}

void Container::getEdgeIntersection (Vertex *v,double *P[3],Point &p)
{
  double a=dot(P[0],P[0]),
         b=dot(P[1],P[1]),
         c=dot(P[2],P[2]),
         d=dot(P[0],P[1]),
         e=dot(P[1],P[2]),
         f=dot(P[2],P[0]),
         g=dot(P[0],v->getpN_ptr()),
         h=dot(P[1],v->getpN_ptr()),
         i=dot(P[2],v->getpN_ptr());
  // first pair of face vertices, P[0],P[1]
  double uDen = a-2*d+b;
  if (uDen) 
  {
    double u = (a-d-g+h)/uDen;
    if (u>0 && u<1) 
    {
      p.add(P[0][0]+u*(P[1][0]-P[0][0]),
            P[0][1]+u*(P[1][1]-P[0][1]),
            P[0][2]+u*(P[1][2]-P[0][2]));
    }
  }
  // second pair of face vertices, P[1],P[2]
  uDen = b-2*e+c;
  if (uDen) 
  {
    double u = (b-e-h+i)/uDen;
    if (u>0 && u<1) 
    {
      p.add(P[1][0]+u*(P[2][0]-P[1][0]),
            P[1][1]+u*(P[2][1]-P[1][1]),
            P[1][2]+u*(P[2][2]-P[1][2]));
    }
  }
  // third pair of face vertices, P[2],P[0]
  uDen = c-2*f+a;
  if (uDen) 
  {
    double u = (c-f-i+g)/uDen;
    if (u>0 && u<1) 
    {
      p.add(P[2][0]+u*(P[0][0]-P[2][0]),
            P[2][1]+u*(P[0][1]-P[2][1]),
            P[2][2]+u*(P[0][2]-P[2][2]));
    }
  }
}

bool Container::computeClosest (Face *f,Vertex *v,double &squareD,double vn[3]) 
{
  Controls & cs(Controls::instance());
  bool signal = false;
  // initialize point class instance with current vertex
  Point p(v->getpN(0),v->getpN(1),v->getpN(2));
  // get face vertex coordinates
  //double *P[3] = {&(f->ptr_vertex(0)->getpN(0)),
  //                &(f->ptr_vertex(1)->getpN(0)),
  //                &(f->ptr_vertex(2)->getpN(0))};
  double *P[3] = {f->ptr_vertex(0)->getpN_ptr(),
                  f->ptr_vertex(1)->getpN_ptr(),
                  f->ptr_vertex(2)->getpN_ptr()};

  // compute vector connecting arbitrary face vertex and current vertex
  double diff[3] = {P[0][0]-v->getpN(0),P[0][1]-v->getpN(1),P[0][2]-v->getpN(2)};
  // compute indicators
  double num=dot(vn,diff);
  // if current vertex does not lie on face plane
  bool inside = false;
  if (!(num<cs.get_double_epsilon())) 
  {
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
  if (!inside)
  {
    // add each face vertex
    for (int i=0;i<3;i++) { p.add(P[i][0],P[i][1],P[i][2]); }
    // add points of minimum distance between current vertex and each face edge
    getEdgeIntersection(v,P,p);
  }
  // save closest point
  double c[3] = { p.a , p.b , p.c };
  // invert vertex normal if vertex is not nice
  double vn_copy[3];
  if (v->getObject()->vertexIsNice(v)){vn_copy[0]=vn[0];vn_copy[1]=vn[1];vn_copy[2]=vn[2]; } 
  else 					  {vn_copy[0]=-vn[0];vn_copy[1]=-vn[1];vn_copy[2]=-vn[2]; } 
  // compute separation vector
  double sep_vec[3] = {c[0]-v->getpN(0),c[1]-v->getpN(1),c[2]-v->getpN(2)};
  // compute cosine of angle between outward normal and separation vector
  // which is equal to dot product of vectors divided by vector magnitudes
  double cos_angle = dot(sep_vec,vn_copy)/sqrt(dot(vn_copy,vn_copy))/sqrt(dot(sep_vec,sep_vec));
  // is closest point located within angle window as defined in controls.cc?
  if ( cos_angle > 0)
  {
    // compute square of separation distance
    double temp=0;
    for (int i=0;i<3;i++) {temp+=(c[i]-v->getpN(i))*(c[i]-v->getpN(i));}
    // closest point must be within a specified angle of vertex normal
    // to avoid grabbing points on same object
    // alternatively, the point is allowed to be within neighborhood_radius of vertex
    // since that implies the point is not near vertex on surface
    if ( (cos_angle > cs.get_closest_point_cosine()) || temp<(cs.get_neighborhood_radius()*cs.get_neighborhood_radius()))
    {
      // if square of separation distance is closer than square of SEARCH_RADIUS
      if (temp<(cs.get_search_radius()*cs.get_search_radius())) 
      {
        // if square of separation distance is less than
        // previously saved square of separation distance
        if (temp<squareD||!squareD) 
        {
          // save 
          //for (int i=0;i<3;i++) {v->getpC(i)=c[i];}
          for (int i=0;i<3;i++) {v->setpC(i,c[i]);}
          //v->ptr_closest_face()=f;
          v->setClosestFacePtr(f);
          squareD = temp;
          signal=true;
        }
      }
    }
  }
  return signal;
}

void Container::printCumulative (void)
{
  // ISSUE WARNING IF ANY OBJECT FAILED INTEGRITY 
  // TEST OR LACK ATTRIBUTES
  // NOTE CUMULATIVE VOLUME ASSUMES ALL 
  // MESHES HAVE SAME ORIENTATION
  //	print Integrity
  Controls & cs(Controls::instance());
  printIntegrity();
  if (good_integrity==false)
  {
    cout << "\n\nWarning: Attributes and "
          << "characteristics were not evaluated,\n"
          << " since mesh file failed the integrity check.\n\n";
  }
  else 
  {
    //	print attributes
    printAttr();
    //	print characteristics
    if (cs.get_compute_attributes_only()==0) 
    {
      printChars();
    }
  }
  cout << "/* ********************** "
        << "END ********************** */\n\n";
}

void Container::printAttr (void)
{
  cout << "MESH SET ATTRIBUTES:\n\n";
  if (good_integrity==false)
  {
    cout << "    Warning: These attribute summaries may be inaccurate,\n"
          << "    since some mesh files failed the integrity check.\n";
  }
  // closed
  std::pair<int,int> pp=countClosed();
  if (pp.first==0)
  {
    cout << "    # closed mesh objects: none\n";
  }
  else 
  {
    cout << "    # closed mesh objects: " << pp.first << endl;
  }
  if (pp.second==0)
  {
    cout << "    # open mesh objects: none\n";
    cout << "    # border edges: none\n";
  }
  else 
  {
    cout << "    # open mesh objects: " << pp.second << endl;
    cout << "    # border edges: " << countBorder() << endl;
  }
  // manifold
  int val[3]={0,0,0};
  countManifold(val);
  if (val[0]==0)
  {
    cout << "    # manifold mesh objects: none\n";
  }
  else 
  {
    cout << "    # manifold mesh objects: " << val[0] << endl;
  }
  if (val[2]==0)
  {
    cout << "    # mesh objects with undefined manifoldness: none\n";
  }
  else 
  {
    cout << "    # mesh objects with undefined manifoldness: " << val[2] << endl;
  }
  if (val[1]==0)
  {
    cout << "    # nonmanifold mesh objects: none\n";
    cout << "    # nonmanifold vertices: none\n";
    cout << "    # nonmanifold edges: none\n";
  }
  else 
  {
    cout << "    # nonmanifold mesh objects: " << val[1] << endl;
    cout << "    # nonmanifold vertices: " << countNonmanV() << endl;
    cout << "    # nonmanifold edges: " << countNonmanE() << endl;
  }
  // consistent
  countConsistent(val);
  if (val[0]==0)
  {
    cout << "    # mesh objects with consistently "
          << "oriented face normals: none\n";
  }
  else 
  {
    cout << "    # mesh objects with consistently "
          << "oriented face normals: " << val[0] << endl;
  }
  if (val[1]==0)
  {
    cout << "    # mesh objects with inconsistently "
          << "oriented face normals: none\n";
  }
  else 
  {
    cout << "    # mesh objects with inconsistently "
          << "oriented face normals: " << val[1] << endl;
    cout << "       # flipped edges: " << countFlipped() << endl;
  }
  if (val[2]==0)
  {
    cout << "    # mesh objects whose face normal "
          << "orientation is undefined (because mesh is nonmanifold): none\n";
  }
  else 
  {
    cout << "    # mesh objects whose face normal "
          << "orientation is undefined (because mesh is nonmanifold): " << val[2] << endl;
  }
  // outward
  countOutward(val);
  if (val[0]==0)
  {
    cout << "    # mesh objects with outward "
          << "oriented face normals: none\n";
  }
  else 
  {
    cout << "    # mesh objects with outward "
          << "oriented face normals: " << val[0] << endl;
  }
  if (val[1]==0)
  {
    cout << "    # mesh objects with inward "
          << "oriented face normals: none\n";
  }
  else 
  {
    cout << "    # mesh objects with inward "
          << "oriented face normals: " << val[1] << endl;
  }
//  if (val[2]==0)
//  {
//    cout << "    # mesh objects whose face normal "
//          << "orientation is undefined: none\n";
//  }
//  else 
//  {
//    cout << "    # mesh objects whose face normal "
//          << "orientation is undefined: " << val[2] << endl;
//  }
  cout << endl;
}

void Container::printChars (void)
{
  cout << "MESH SET CHARACTERISTICS:\n\n";
  //	print characteristics
  if (good_integrity==false)
  {
    cout << "    Warning: These characteristics "
          << "summaries may be inaccurate,\n"
          << "    since some mesh files failed the integrity check.\n";
  }
  cout << "    # objects: " << countObject() << endl;
  cout << "    # vertices: " << countVertex() << endl;
  cout << "    # faces: " << countFace() << endl;
  cout << "    # edges: " << countEdge() << endl;
  cout << "    # components: " << countComponents() << endl;
  int val[3]={0,0,0};
  countManifold(val);
  if (val[1]>0)
  {
    cout << "    # boundaries: Since object is nonmanifold,\n"
          << "    # boundaries: the number of boundaries may be underestimated.\n";
  }
  int a=countBoundaries();
  if (a==0)
  {
    cout << "    # boundaries: none\n";
  }
  else 
  {
    cout << "    # boundaries: " << a << endl;
  }
  a=countIndistin();
  if (a==0)
  {
    cout << "    # indistinguishable vertices: none\n";
  }
  else 
  {
    cout << "    # indistinguishable vertices: " << a << endl;
  }
  // volume
  countConsistent(val);
  if (val[0]>0)
  {
    cout << "    object volume: [(data units)^3]" << endl;
    cout << "    object volume: " << countVol() << endl;
  }
  else 
  {
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
        << setbb[0] << ","
        << setbb[1] << ","
        << setbb[2] << "]["
        << setbb[3] << ","
        << setbb[4] << ","
        << setbb[5] << "]" << endl << endl;
  // vertex adjacent faces
  cout << "    Vertex adjacent face statistics [faces]:" << endl;
  adjacent_face.printStats();
  cout << "    Vertex adjacent face histogram [faces]:" << endl;
  adjacent_face.printAdjacentFaceHistogram();
  cout << endl;
  // face area
  cout << "    Face area statistics [(data units)^2]:" << endl;
  cout << "       total    " << area.getSum() << endl;
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
        for (std::vector<double>::iterator ww=edge_length.x.begin();ww!=edge_length.x.end();ww++)
        {
        double diff = *ww-edge_length.max;
        if (fabs(diff)<1E-6)
        {
        cout << "\nControls::printChars: max value=" << edge_length.max
        << " found in edge_length.x\n";
        ffound=true;
        }
        }
        if (ffound==false)
        {
        cout << "\nControls::printChars: max value=" << edge_length.max
        << " NOT found in edge_length.x\n";
        }
        */	// DEBUG
  cout << "    Edge length histogram [data units]:" << endl;
  edge_length.printHistogram();
  cout << endl;
  // if edge angles computed, i.e. if at least one mesh was conistent
  if (val[0]>0)
  {
    cout << "    Edge angle statistics [degrees]:" << endl;
    edge_angle.printStats();
    /*		// DEBUG
                ffound=false;
                for (std::vector<double>::iterator ww=edge_angle.x.begin();ww!=edge_angle.x.end();ww++)
                {
                double diff = *ww-edge_angle.max;
                if (fabs(diff)<1E-6)
                {
                cout << "\nControls::printChars: max value=" << edge_angle.max
                << " found in edge_angle.x\n";
                ffound=true;
                }
                }
                if (ffound==false)
                {
                cout << "\nControls::printChars: max value=" << edge_angle.max
                << " NOT found in edge_angle.x\n";
                }
                */		// DEBUG
    cout << "    Edge angle histogram [degrees]:" << endl;
    edge_angle.printHistogram();
    cout << endl;
  }
}

void Container::analyzeBatch (Space &s)
{
  // find intersecting faces
  if (Controls::instance().get_detect_interobject_intersections())
  {
    cout << "Find intersecting faces of all objects...";cout.flush();
    Object *oo=o.front();
    oo->findIntersectingFaces(this,s);
    cout << "complete.\n";cout.flush();
  }
}

void Container::printIntegrity (void)
{
  cout << "\nMESH SET INTEGRITY:\n\n";
  int a=countOrphan();
  if (a==0)
  {
    cout << "    # orphan vertices: none\n";
  }
  else 
  {
    cout << "    # orphan vertices: " << a << endl;
  }
  a=countMissing();
  if (a==0)
  {
    cout << "    # missing vertices: none\n";
  }
  else 
  {
    cout << "    # missin vertices: " << a << endl;
  }
  a=countDegen();
  if (a==0)
  {
    cout << "    # degenerate faces: none\n";
  }
  else 
  {
    cout << "    # degenerate faces: " << a << endl;
  }
  a=countDuplV();
  if (a==0)
  {
    cout << "    # duplicate vertex indices: none\n";
  }
  else 
  {
    cout << "    # duplicate vertex indices: " << a << endl;
  }
  a=countDuplF();
  if (a==0)
  {
    cout << "    # duplicate face indices: none\n";
  }
  else 
  {
    cout << "    # duplicate face indices: " << a << endl;
  }
  cout << endl;
}

void Container::vertexAdjacentFaces (void)
{
  //	mmap_iv af;
  // for each element in vector
  for (c_d_iterator i=adjacent_face.first_element();i!=adjacent_face.one_past_last_element();i++)
  {
    //		int c=(*i)->f.size();
    //		af.insert(std::make_pair(c,*i));
    //		adjacent_face.n++;
    adjacent_face.add2sum(*i);
    adjacent_face.add2sum2((*i)*(*i));
    adjacent_face.add2total(*i);
    if (*i<adjacent_face.getMin()) {adjacent_face.setMin(*i);}
    if (*i>adjacent_face.getMax()) {adjacent_face.setMax(*i);}
    // add to vector
    //		adjacent_face.x.push_back(c);
  }
  // build adjacent face histogram
  adjacent_face.createAdjacentFaceHistogram();
}

void Container::areaAspectRatio (void)
{
  assert(area.getSize()>0);
  // for each element in area vector
  for (c_d_iterator i=area.first_element();i!=area.one_past_last_element();i++)
  {
    //    cout << "1" << endl;cout.flush();
    area.add2sum(*i);
    area.add2sum2((*i)*(*i));
    area.add2total(*i);
    if (*i<area.getMin()) {area.setMin(*i);}
    if (*i>area.getMax()) {area.setMax(*i);}
  }

  // for each element in aspect ratio vector
  for (c_d_iterator i=aspect_ratio.first_element();i!=aspect_ratio.one_past_last_element();i++)
  {
    //    cout << "2" << endl;cout.flush();
    aspect_ratio.add2sum(*i);
    aspect_ratio.add2sum2((*i)*(*i));
    aspect_ratio.add2total(*i);
    if (*i<aspect_ratio.getMin()) aspect_ratio.setMin(*i);
    if (*i>aspect_ratio.getMax()) aspect_ratio.setMax(*i);
  }
  // build face area histogram
  //  cout << "3" << endl;cout.flush();
  area.createHistogram();
  //  cout << "4" << endl;cout.flush();
  // build aspect ratio histogram
  aspect_ratio.createAspectRatioHistogram();
  // cout << "5" << endl;cout.flush();
}

void Container::processEdgeLengths (void)
{
  // for each element in area vector
  for (c_d_iterator i=edge_length.first_element();i!=edge_length.one_past_last_element();i++)
  {
    //		double l = (*i)->l;
    //		edge_length.n++;
    edge_length.add2sum(*i);
    edge_length.add2sum2((*i)*(*i));
    edge_length.add2total(*i);
    if (*i<edge_length.getMin()) edge_length.setMin(*i);
    if (*i>edge_length.getMax()) edge_length.setMax(*i);
  }
  edge_length.createHistogram();
}

void Container::computeEdgeAngles (void)
{
  // for each element in area vector
  for (c_d_iterator i=edge_angle.first_element();i!=edge_angle.one_past_last_element();i++)
  {
    edge_angle.add2sum(*i);
    edge_angle.add2sum2((*i)*(*i));
    edge_angle.add2total(*i);
    if (*i<edge_angle.getMin()) edge_angle.setMin(*i);
    if (*i>edge_angle.getMax()) edge_angle.setMax(*i);
  }
  edge_angle.createHistogram();
}

void Container::analyzeCumulative (void)
{
  // vertices
  cerr << "Analyze vertex adjacent faces for set of all objects...................";
  cerr.flush();
  vertexAdjacentFaces();
  cerr << "complete.\n";cerr.flush();
  // faces
  cerr << "Compute face area and aspect ratio for set of all objects..............";
  cerr.flush();
  areaAspectRatio();
  cerr << "complete.\n";cerr.flush();
  // edges
  cerr << "Analyze edge lengths for set of all objects............................";
  cerr.flush();
  processEdgeLengths();
  cerr << "complete.\n";cerr.flush();
  //	if (manifold && consistent)
  //{
  // manifold and consistently oriented face normals
  cerr << "Analyze edge angles for set of all objects.............................";
  cerr.flush();
  computeEdgeAngles();
  //	}
}

