// Author: Justin Kinney
// Date: Sep 2008

#include "container.h"

#include <dirent.h>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

#include "log.h"
#include "edge.h"
#include "nice.h"
#include "controls.h"
#include "intersecting_faces.h"
#include "meshmorph.h"
#include "opttritri.h"
#include "Wm4Quaternion.h"
#include "vertex_schedule.h"
#include "octree_visitor_check_face.h"

using std::cout;
using std::endl;

Container * Container::only_one = NULL;

Container & Container::instance()
{
  // Not thread-safe.
  // -- lock mutex
  if (only_one == NULL)
    only_one = new Container();
  // -- unlock mutex
  return *only_one;
}

Container::Container (void)
:fo(NULL),files(),frozen(),world(6,0),cleft(),peri(),octree(),o()
{
  scanDir();
  scanFiles();
  createEdges();
  findVertAdj();
  computeVertexNormals();
  sortAdjacentFaces();
  boundWorld();
}

/** Get the total number of vertices in the input data.
 * \return The total number of vertices in the input data.
 */

int Container::getVertexCount (void) const
{
  int j = 0;
  // for each object
  for (o_cit i=o.begin();i!=o.end();++i)
  {
    j += i->v.size();
  }
  return j;
}

/** Get the total number of faces in the input data.
 * \return The total number of faces in the input data.
 */

int Container::getFaceCount (void) const
{
  int j = 0;
  // for each object
  for (o_cit i=o.begin();i!=o.end();++i)
  {
    j += i->f.size();
  }
  return j;
}

/** Get the total number of edges in the input data.
 * \return The total number of edges in the input data.
 */

int Container::getEdgeCount (void) const
{
  int j=0;
  // for each object
  for (o_cit i=o.begin();i!=o.end();++i)
  {
    j += i->e.size();
  }
  return j;
}

/** Find the closest face of each vertex in container
 * and compare with closest face stored in vertex class.
 *
 * \param[in] group Number of group of vertices being moved.
 * \param[in] suffix String to append to file name to identify
 * when in the program the file was written.
 */

//const int         TARGET_VERTEX_INDEX = 632;
//const std::string TARGET_VERTEX_NAME("d014");
//const int         TARGET_FACE_INDEX = 1202;
//const std::string TARGET_FACE_NAME("a170");

void Container::checkClosestFace (int const & group,std::string suffix)
{
  cout << "\nIteration " << group << ": Check closest face to vertices "
        << suffix << "...........";
  cout.flush();
  // open output file
 // Controls & cs(Controls::instance());
 // char filename[cs.get_max_filename_size()];
 // sprintf(filename,"%sclosest_face_discrepencies_%d_%s",
 //         cs.get_output_data_dir().c_str(),group,suffix.c_str());
 // std::ofstream F(filename);
 // if (F.is_open()==false)
 // {
 //   fprintf(stderr,"\nCouldn't open output file %s\n",filename);
 //   assert(F.is_open()==true);
 //   exit(1);
 // }
  vector3 p,q;
  double sqd;
  // for each object in container
  for (o_it i=o.begin();i!=o.end();++i)
  {
    //if ((*i).getName()!=TARGET_VERTEX_NAME) continue;
    // for each vertex in object
    for (v_it j=i->v.begin();j!=i->v.end();++j)
    {
      //if ((*j).getIndex()!=TARGET_VERTEX_INDEX) continue;
      if (binary_search(frozen.begin(),frozen.end(),&(*j))==false)
      {
        // DEBUG
//        cout << "Container::checkClosestFace: target vertex found.\n";
//        // for each object in container
//        for (o_it ii=o.begin();ii!=o.end();++ii)
//        {
//          // for each face in object
//          for (f_it jj=ii->f.begin();jj!=ii->f.end();++jj)
//          {
//            if (jj->isMatch(TARGET_FACE_INDEX,TARGET_FACE_NAME)==true) j->setFace(&(*jj));
//          }
//        }
        //if (ncl==NULL)
        //{
        //  cout << "Container::checkClosestFace: Error target face NOT found.\n";
        //  exit(1);
        //}
        //else
        //{
        //  ncl->print(cout);
        //}
        // DEBUG
        Face *ncl = NULL;
        //int neighbor_count = -666;
        // find closest point to current vertex
        findClosestPtToVertex(&(*j),p,sqd,ncl);
        // if closest faces are the same, then no problem
        if (j->getClosestFace()==ncl) continue;
        // closest face is incorrect
        cout << "\n\nContainer::checkClosestFace: Error. "
              << "Stored closest face to vertex does not match actual closest face.\n";
        j->print(std::cout);
        cout.precision(12);
        cout << endl;
        if (ncl != NULL)
        {
          cout << "\nActual face (good)\n"
                << "extracellular width " << sqrt(sqd) << endl
                << "closest point ["
                << p.p[0] << " "
                << p.p[1] << " "
                << p.p[2] << "]\n";
          ncl->print(std::cout);
          cout << endl;
        }
        else
        {
          cout << ", good NULL]\n\n";
        }
        if (j->getClosestFace() != NULL)
        {
          vec_fp fp;
          fp.push_back(const_cast<Face*>(j->getClosestFace()));
          ncl=NULL;
          sqd=1E30;
          int dummy;
          double search_radius = sqrt(Controls::instance().get_search_radius_sq());
          if (findClosestPtToVertexAmongFaces(&(*j),fp.begin(),fp.end(),search_radius,q,sqd,ncl,dummy)==true)
          {
            cout << "\nStored face (bad)\n"
                  << "extracellular width " << sqrt(sqd) << endl
                  << "closest point ["
                  << q.p[0] << " "
                  << q.p[1] << " "
                  << q.p[2] << "]\n";
          }
          else
          {
            cout << "\nStored face (bad)\n"
                  << "no closest point was found on stored face." << endl;
          }
          j->getClosestFace()->print(std::cout);
          cout << endl;
        }
        else
        {
          cout << "[bad NULL]\n";
        }
        writeMeshData(666);
        checkFaces("666");
        checkFacesInOctree();
        exit(0);
      }
    }
  }
  cout << "complete.\n";
  cout.flush();
}

/** For each object build instance of edge class
 * using face and vertex class data.
 */

void Container::createEdges (void)
{
  Controls & cs(Controls::instance());
  if (cs.get_write_verbose_init()==false)
  {
    cout << "Creating edges.................................";
    cout.flush();
  }
  int k=0;
  double goal = 0.2;
  double inc = 0.2;
  double a = 1.0/o.size();
  if (goal<a)
  {
    goal=a;
    inc = a;
  }
  printf("0%%..");
  fflush(stdout);
  // for each object, create edges
  for (o_it i=o.begin();i!=o.end();++i)
  {
    if (cs.get_write_verbose_init()==true)
    {
      cout << "Creating edges for " << i->getName() << "...";
      cout.flush();
    }
    i->e.reserve(4*i->v.size());
    // determine number of digits in largest vertex index
    int num_digits = i->setNumDigits();
    // create map for finding edges
    map_s_ep hm;
    // for each face
    for (f_it j=i->f.begin();j!=i->f.end();++j)
    {
      i->processEdge(&(*j),j->getVertex(0),j->getVertex(1),j->getVertex(2),hm,num_digits);
      i->processEdge(&(*j),j->getVertex(1),j->getVertex(2),j->getVertex(0),hm,num_digits);
      i->processEdge(&(*j),j->getVertex(2),j->getVertex(0),j->getVertex(1),hm,num_digits);
    }
    if (cs.get_write_verbose_init()==true)
    {
      cout << "complete.\n";
      cout.flush();
    }
    // track progress
    double progress = static_cast<double>(++k)*a;
    if (progress>goal || !distinguishable(progress,goal))
    {
      printf("%d%%..",static_cast<int>(goal*100));
      fflush(stdout);
      goal+=inc;
    }
    assert(i->getFE()==&i->e[0]);
  }
  if (cs.get_write_verbose_init()==false)
  {
    cout << "complete.\n";
    cout.flush();
  }
}

/** For each vertex in each object add pointers 
 * to adjacent faces of vertex to vertex class.
 */

void Container::findVertAdj (void)
{
  Controls & cs(Controls::instance());
  if (cs.get_write_verbose_init()==false)
  {
    cout << "Finding vertex adjacencies.....................";
    cout.flush();
  }
  // for each object, find vertex adjacencies
  if (cs.get_write_verbose_init()==true){ cout << endl;} 
  for (o_it i=o.begin();i!=o.end();++i)
  {
    if (cs.get_write_verbose_init()==true)
    {
      cout << "Finding adjacencies for " << i->getName() << "...";
      cout.flush();
    }
    i->findVertAdj();
    if (cs.get_write_verbose_init()==true)
    {
      cout << "complete.\n";
      cout.flush();
    }
  }
  if (cs.get_write_verbose_init()==false)
  {
    cout << "complete.\n";
    cout.flush();
  }
}

/** Find all .mesh files in input directory.
*/

void Container::scanDir (void)
{
  Controls & cs(Controls::instance());
  std::string::size_type found;
  DIR *pdir;		// pointer to a directory data structure
  struct dirent *pent;	// pointer to dirent structure
  if (cs.get_write_verbose_init()==true){ cout << endl;cout.flush();}
  pdir = opendir(cs.get_input_data_dir().c_str());
  if (pdir==NULL)
  {
    printf("Error. Could not open %s.\n",cs.get_input_data_dir().c_str());
    assert(pdir!=NULL);
    exit(1);
  }
  while ((pent=readdir(pdir)))
  {
    // copy char array to string
    std::string line = pent->d_name;
    // if file of typ *.mesh
    found = line.find(".mesh",0);
    // if found
    if (found != std::string::npos)
    {
      // save filename
      files.push_back(line);
      // update index
      // print file found to screen
      if (cs.get_write_verbose_init()==true)
      {
        cout << "file found: " << line << "\n"; cout.flush();
      }
    }
  }
  if (files.empty()==true)
  {
    cout << "\n\nContainer::scanDir: Error. "
          << "No .mesh files were found in input dir.\n"
          << "Input directory = " << cs.get_input_data_dir().c_str()
          << endl << endl;
    assert(files.empty()==false);
    exit(1);
  }
  closedir(pdir);
  if (cs.get_write_verbose_init()==true){ cout << endl;cout.flush();}
  o.reserve(files.size());
}

/** For each .mesh file found in input directory
 * build an instance of object class.
 *
 * \return Void.
 */

void Container::scanFiles (void)
{
  Controls & cs(Controls::instance());
  cout << "\nInput Data Directory\n"
        << cs.get_input_data_dir().c_str() << "\n\n"
        << "Reading Directory..............................";
  cout.flush();
  // for each input file
  for (uint count=0;count<files.size();++count)
  {
    // record object name
    std::string::size_type pos1 = files[count].find(".",0);
    if (pos1!=std::string::npos)
    {
      o.push_back(Object(files[count].substr(0,pos1)));
      if (count==0)
      {
        fo = &o[0];
      }
    }
    else
    {
      cout << "Error! Object name was not found in " << files[count] << "\n";
      assert(pos1==std::string::npos);
      exit(1);
    }
    // scan file
    std::string file = cs.get_input_data_dir().c_str() + files[count];
    if (cs.get_write_verbose_init()==true)
    {
      cout << "loading file " << files[count] << "...";
    }
    cout.flush();
    scanFile(&(o.back()),file.c_str());
    if (cs.get_write_verbose_init()==true)
    {
      cout << "complete.\n";
      cout.flush();
    }
  }
  // trim Object vector
  assert(fo==&o[0]);
  cout << "complete.\n";
}

/** Scan input data file once to count total number
 * of vertices and faces in the object.
 * \param[in] obj Pointer to parent Object.
 * \param[in] filename Input file name.
 */

void Container::assessFile (Object * const obj,const char *filename)
{
  char line[2048];
  int v_count=0,f_count=0;
  // open file
  FILE *F = fopen(filename,"r");
  if (!F) { printf("Couldn't open input file %s\n",filename);return;}
  // for every line in file
  for (char *str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F))
  {
    // skip leading whitespace
    while (strchr(" \t,",*str)!=NULL) { str++;}
    // if first character is V for Vertex
    if (strchr("V",*str)!=NULL){ v_count++; }
    // if first character is F for Face
    else if (strchr("F",*str)!=NULL){ f_count++; }
  }
  fclose(F);
  obj->v.reserve(v_count);
  obj->f.reserve(f_count);
}

/** Build Object from input file data.
 * \param[in] obj Pointer to parent Object.
 * \param[in] filename Input file name.
 */

//void Container::scanFile (Object * const obj,const char *filename)
//{
//  char line[2048],*str;
//  FILE *F;
//  int vertex_num=0,polygon_num=0;
//  // reserve space for vertices and faces
//  assessFile(obj,filename);
//  vec_vp vp;
//  // open file
//  F = fopen(filename,"r");
//  if (!F) { printf("Couldn't open input file %s\n",filename);return;}
//  // for every line in file
//  for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F))
//  {
//    // skip leading whitespace
//    while (strchr(" \t,",*str)!=NULL) { str++;}
//    // if first character is V for Vertex,
//    // add new linked list class instance
//    if (strchr("V",*str)!=NULL)
//    {
//      obj->v.push_back(Vertex(str,obj));
//      Vertex *vv = &obj->v.back();
//      // check that vertex indices increase sequentially from zero
//      // '+1' since vertex_num used as 0-based array index 
//      // but mesh format is 1-based vertex indexing
//      if (vertex_num+1 != vv->getIndex())
//      {
//        cout << "\n\nContainer::ScanFile: Error. "
//              << "Input vertex indices must increase "
//              << "sequentially from zero.\n"
//              << "Container::ScanFile: Vertex index must be "
//              << vertex_num+1 << ".\n"
//              << "Container::ScanFile: Vertex read:\n";
//        vv->print(cout);
//        cout << endl;
//        assert(vertex_num+1 == vv->getIndex());
//      }
//      vp.push_back(vv);
//      if (vertex_num==0)
//      {
//        obj->setFV(&obj->v[0]);
//      }
//      vertex_num++;
//      if (Controls::instance().get_write_verbose_init()==true)
//      {
//        cout<<filename<<": ";
//        vv->print(std::cout);
//      }
//    }
//    // if first character is F for Face,
//    // add new linked list class instance
//    else if (strchr("F",*str)!=NULL)
//    {
//      obj->f.push_back(Face(str,vp));
//      if (polygon_num==0)
//      {
//        obj->setFF(&obj->f[0]);
//      }
//      polygon_num++;
//      if (Controls::instance().get_write_verbose_init()==true)
//      {
//        Face *ff = &obj->f.back();
//        ff->print(std::cout);
//      }
//    }
//  }
//  assert(obj->getFV()==&obj->v[0]);
//  assert(obj->getFF()==&obj->f[0]);
//  fclose(F);
//}

void Container::scanFile (Object * const obj,const char *filename)
{
  int vertex_num=0,polygon_num=0;
  // reserve space for vertices and faces
  assessFile(obj,filename);
  vec_vp vp;
  // open file
  std::ifstream f(filename);
  if (f.is_open()==false) { printf("Couldn't open input file %s\n",filename);return;}
  // for every line in file
  std::string prefix, idx, x, y, z;
  f >> prefix;
  while(!f.eof())
  {
    //f >> prefix;
    if (prefix == "Vertex")
    {
      f >> idx >> x >> y >> z;
      obj->v.push_back(Vertex(atoi(idx.c_str()),
                              strtod(x.c_str(), NULL),
                              strtod(y.c_str(), NULL),
                              strtod(z.c_str(), NULL),obj));
      //obj->v.push_back(Vertex(str,obj));
      Vertex *vv = &obj->v.back();
      // check that vertex indices increase sequentially from zero
      // '+1' since vertex_num used as 0-based array index 
      // but mesh format is 1-based vertex indexing
      if (vertex_num+1 != vv->getIndex())
      {
        cout << "\n\nContainer::ScanFile: Error. "
              << "Input vertex indices must increase "
              << "sequentially from zero.\n"
              << "Container::ScanFile: Vertex index must be "
              << vertex_num+1 << ".\n"
              << "Container::ScanFile: Vertex read:\n";
        vv->print(cout);
        cout << endl;
        assert(vertex_num+1 == vv->getIndex());
      }
      vp.push_back(vv);
      if (vertex_num==0)
      {
        obj->setFV(&obj->v[0]);
      }
      vertex_num++;
      if (Controls::instance().get_write_verbose_init()==true)
      {
        cout<<filename<<": ";
        vv->print(std::cout);
      }
    }
    else if (prefix == "Face")
    {
      f >> idx >> x >> y >> z;
      obj->f.push_back(Face(atoi(idx.c_str()),
                            atoi(x.c_str()),
                            atoi(y.c_str()),
                            atoi(z.c_str()),vp));
      //obj->f.push_back(Face(str,vp));
      if (polygon_num==0)
      {
        obj->setFF(&obj->f[0]);
      }
      polygon_num++;
      if (Controls::instance().get_write_verbose_init()==true)
      {
        Face *ff = &obj->f.back();
        ff->print(std::cout);
      }
//      if (obj->getFF()!=&obj->f[0])
//      {
//        cout << endl << obj->getFF()
//              << " != " << &obj->f[0] << endl << endl;
//        cout << "polygon_num = " << polygon_num << endl;
//        cout << "object = " << obj->getName() << endl;
//        cout << "obj->f.size() = " << obj->f.size() << endl;
//        cout << prefix << " " << idx << " " << x << " " << y << " " << z << endl;
//      }
//      else
//      {
//        cout << obj->getFF()
//              << " == " << &obj->f[0]  << endl;
//        cout << "obj->f.size() = " << obj->f.size() << endl;
//        cout << "polygon_num = " << polygon_num << endl;
//        cout << prefix << " " << idx << " " << x << " " << y << " " << z << endl;
//        cout << "----------------------------\n";
//      }
//      assert(obj->getFF()==&obj->f[0]);
    }
    else
    {
      printf("\n\nContainer::scanFile: Hey, man, I don't know how to parse a '%s'\n", prefix.c_str());
      printf("Container::scanFile: From file '%s'\n", filename);
      exit(1);
    }
    f >> prefix;
//    if (f.eof()) break;

  }
  //  char line[2048];
  //  for (char * str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F))
  //  {
  //    // skip leading whitespace
  //    while (strchr(" \t,",*str)!=NULL) { str++;}
  //    // if first character is V for Vertex,
  //    // add new linked list class instance
  //    if (strchr("V",*str)!=NULL)
  //    {
  //    }
  //    // if first character is F for Face,
  //    // add new linked list class instance
  //    else if (strchr("F",*str)!=NULL)
  //    {
  //    }
  //  }
  assert(obj->getFV()==&obj->v[0]);
//  if (obj->getFF()!=&obj->f[0])
//  {
//    cout << endl << obj->getFF()
//          << " != " << &obj->f[0] << endl << endl;
//  }
  assert(obj->getFF()==&obj->f[0]);
  f.close();
}

/** Get pointer to object with matching name.
 * \param[in] name Name of object of interest.
 * \return Pointer to object with name that matches
 * the input name; otherwise NULL.
 */

Object* Container::getObjectPointer (std::string name) const
{
  // for each object* in container
  for (o_cit i=o.begin();i!=o.end();++i)
  {
    if (name==i->getName())
    {
      return const_cast<Object*>(&(*i));
    }
  }
  return NULL;
}

/** Parse input file of vertices (obejct name and vertex index) and store data.
 * \param[in] filename Input file name.
 * \param[out] not_found Set of object names
 * from file that were not found in input data.
 * \return Stored object data.
 */

//mmap_oi Container::loadMap(const char *filename,s_set & not_found)
//{
//  mmap_oi mymap;
//  // open file
//  FILE *F = fopen(filename,"r");
//  if (!F)
//  {
//    printf("Couldn't open input file %s\n",filename);
//    return mymap;
//  }
//  // for every line in file
//  // NOTE THE ASSUMED MAX character array length
//  char val[2048],*eptr=NULL,line[2048],*input=NULL;
//  Object *o1=NULL,*o2=NULL;
//  for (input=fgets(line,2048,F) ; input!=NULL ; input=fgets(line,2048,F))
//  {
//    // get object name
//    while (strchr(" \t,",*input)!=NULL) { input++; }
//    int i=0;
//    while (strchr(" \t",*input)==NULL)
//    {
//      val[i++] = *input++;
//    }
//    val[i]=0;
//    // val contains object name
//
//    // get Object*
//    if (o2==NULL){o1 = getObjectPointer(val);}
//    else
//    {
//      // if new object name same as old object name
//      if (strcmp(o2->getName().c_str(),val)==false){o1=o2;}
//      else {o1 = getObjectPointer(val);}
//    }
//    o2=o1;
//    if (o1==NULL)
//    {
//      not_found.insert(val);
//    }
//    else
//    {
//      // get vertex index
//      while (strchr(" \t,",*input)!=NULL) { input++; }
//      i=0;
//      while (strchr("0123456789",*input)!=NULL)
//      {
//        val[i++] = *input++;
//      }
//      val[i]=0;
//      int vi = (int)strtod(val,&eptr);
//      if (val==eptr)
//      {
//        vi=0;
//        printf("Error in reading vertex index\n");
//        return mymap;
//      }
//      // create map entry
//      mymap.insert(std::make_pair(o1,vi));
//    }
//  }
//  // close file
//  fclose(F);
//  return mymap;
//}

std::vector<Complex> Container::loadVector(const char *filename,s_set & not_found)
{
  std::vector<Complex> myvec;
  // open file
  std::ifstream f(filename);
  if (f.is_open()==false) { printf("Couldn't open input file %s\n",filename);return myvec;}
  std::string name, index;
  Object *o1=NULL,*o2=NULL;
  f >> name;
  while(!f.eof())
  {
    // get Object*
    if (o2==NULL){o1 = getObjectPointer(name.c_str());}
    else
    {
      // if new object name same as old object name
      if (strcmp(o2->getName().c_str(),name.c_str())==false){o1=o2;}
      else {o1 = getObjectPointer(name.c_str());}
    }
    o2=o1;
    if (o1==NULL)
    {
      not_found.insert(name);
    }
    else
    {
      f >> index;
      // create map entry
      Complex a(o1,atoi(index.c_str()));   
      myvec.push_back(a);
    }
    f >> name;
  }
  // close file
  f.close();
  return myvec;
}

mmap_oi Container::loadMap(const char *filename,s_set & not_found)
{
  mmap_oi mymap;
  // open file
  std::ifstream f(filename);
  if (f.is_open()==false) { printf("Couldn't open input file %s\n",filename);return mymap;}
  std::string name, index;
  Object *o1=NULL,*o2=NULL;
  f >> name;
  while(!f.eof())
  {
    f >> index;
    // get Object*
    if (o2==NULL){o1 = getObjectPointer(name.c_str());}
    else
    {
      // if new object name same as old object name
      if (strcmp(o2->getName().c_str(),name.c_str())==false){o1=o2;}
      else {o1 = getObjectPointer(name.c_str());}
    }
    o2=o1;
    if (o1==NULL)
    {
      not_found.insert(name);
    }
    else
    {
      // create map entry
      mymap.insert(std::make_pair(o1,atoi(index.c_str())));
    }
    f >> name;
  }
  // close file
  f.close();
  return mymap;
}

/** Process frozen vertex data.
 * \param[in] filename Frozen vertex file name.
 */

void Container::readCleft (const char *filename)
{
  cout << "Read cleft vertices...........................";
  cout.flush();
  // load cleft map: object*->Vertex_index (int)
  oi_it front;
  // set of object names from cleft file not found as Object name
  s_set not_found;
  // read cleft vertex file 
  mmap_oi cleft_map = loadMap(filename,not_found);

  // for each element in multimap
  for (oi_it i=cleft_map.begin();i!=cleft_map.end();++i)
  {
    Object *oo=(*i).first;
    int t = (*i).second;
    int a = t-1;
    if (static_cast<uint>(a) >= oo->v.size())
    {
      cout << "Vertex index (" << a
            << ") exceeds vector length ("
            << oo->v.size() << ") "
            << "for object " << oo->getName() << endl;
    }
    assert(static_cast<uint>(a)<oo->v.size());
    int b = oo->v[a].getIndex();
    do{
      if (t==b)
      {
        cleft.push_back(&oo->v[a]);
      }
      else if (t<b){a++;}
      else if (t>b){a--;}
    } while (t!=b);
  }
  sort(cleft.begin(),cleft.end());
  cout << "complete.\n";

  // print objects not found
  if (not_found.empty()==false)
  {
    cout << "\nContainer::readCleft: Warning.\n"
          << "No matching Object* found in container "
          << "for following cleft objects:\n";
    for (ss_it j = not_found.begin();j!=not_found.end();++j)
    {
      cout << *j << endl;
    }
    cout << endl;
  }

  cout.flush();
}

/** Process frozen vertex data.
 * \param[in] filename Frozen vertex file name.
 */

void Container::readFrozen (const char *filename)
{
  cout << "Read frozen vertices...........................";
  cout.flush();
  // load frozen map: object*->Vertex_index (int)
  oi_it front;
  // set of object names from frozen file not found as Object name
  s_set not_found;
  // read frozen vertex file 
  mmap_oi frozen_map = loadMap(filename,not_found);

  // for each element in multimap
  for (oi_it i=frozen_map.begin();i!=frozen_map.end();++i)
  {
    Object *oo=(*i).first;
    int t = (*i).second;
    int a = t-1;
    int b = oo->v[a].getIndex();
    do{
      if (t==b){frozen.push_back(&oo->v[a]);}
      else if (t<b){a++;}
      else if (t>b){a--;}
    } while (t!=b);
  }
  sort(frozen.begin(),frozen.end());
  cout << "complete.\n";

  // print objects not found
  if (not_found.empty()==false)
  {
    cout << "\nContainer::loadFrozenMap: Warning.\n"
          << "No matching Object* found in container "
          << "for following frozen objects:\n";
    for (ss_it j = not_found.begin();j!=not_found.end();++j)
    {
      cout << *j << endl;
    }
    cout << endl;
  }

  cout.flush();
}

/** Write to file the current position of all mesh objects.
 * \param[in] group The current group number from meshmorph.cc.
 */

void Container::writeMeshData (int const & group) const
{
  Controls & cs(Controls::instance());
  // for each object
  for (o_cit i=o.begin();i!=o.end();++i)
  {
    // if block wrap
    if (cs.get_recon_block_wrap()!="")
    {
      // if object is wrapper
      if (i->getName()!=cs.get_recon_block_wrap()) continue;
    }
    // create output filename
    std::string file;
    if (cs.get_append_group_number()==true)
    {
      //char sint[32];
      //sprintf(sint,"_%i.mesh",group);
      //std::string line(sint);
      //file = cs.get_output_data_dir() + i->getName() + cs.get_mesh_output_suffix() + line;
      file = Log::instance().format("%s%s%s_%i.mesh",
                                    cs.get_output_data_dir().c_str(),
                                    i->getName().c_str(),
                                    cs.get_mesh_output_suffix().c_str(),
                                    group);
    }
    else
    {
      //std::string line(".mesh");
      //file = cs.get_output_data_dir() + i->getName() + cs.get_mesh_output_suffix() + line;
      file = Log::instance().format("%s%s%s.mesh",
                                    cs.get_output_data_dir().c_str(),
                                    i->getName().c_str(),
                                    cs.get_mesh_output_suffix().c_str());
    }
    if (cs.get_write_mesh_now()>0)
    {
      if (Vertex_Schedule::instance().getNumMovedVertsGroup()==cs.get_write_mesh_now())
      {
        file = Log::instance().format("%s%s%s_iteration_%i.mesh",
                                      cs.get_output_data_dir().c_str(),
                                      i->getName().c_str(),
                                      cs.get_mesh_output_suffix().c_str(),
                                      Vertex_Schedule::instance().getNumMovedVertsGroup());
      }
    }
    // open output file
    std::ofstream newfile;
    newfile.open (file.c_str(),std::ios::trunc);
    if (newfile.is_open()==false)
    {
      cout << "\nContainer::writeMeshData: "
            << "Error. Unable to open file = "
            << newfile << endl;
      assert(newfile.is_open()==true);
      exit(1);
    }
    newfile.precision(12);
    // for each vertex in object
    for (v_cit j=i->v.begin();j!=i->v.end();++j)
    {
      // print index and final coordinates
      newfile << "Vertex "
            << j->getIndex() << " "
            << *j->getCoord(0) << " "
            << *j->getCoord(1) << " "
            << *j->getCoord(2) << "\n";
    }
    // for each face in object
    for (f_cit k=i->f.begin();k!=i->f.end();++k)
    {
      newfile << "Face "
            << k->getIndex() << " "
            << k->getVertex(0)->getIndex() << " "
            << k->getVertex(1)->getIndex() << " "
            << k->getVertex(2)->getIndex() << endl;
    }
    newfile.close();
  }
}

/** Find and record the closest face to each vertex in the model.
*/

void Container::findClosestFaceToEachVertex (void)
{
  cout << "Find closest face to each vertex...............";
  cout.flush();
  int index=0;
  double goal = 0.2;
  double inc = 0.2;
  double a = 1.0/o.size();
  if (goal<a)
  {
    goal=a;
    inc = a;
  }
  printf("0%%..");
  fflush(stdout);
  // for each object in container
  vector3 dummy_p;
  double dummy_sqd;
  for (o_it i=o.begin();i!=o.end();++i)
  {
    // for each vertex in object
    for (v_it j=i->v.begin();j!=i->v.end();++j)
    {
      if (vertexIsFrozen(&(*j))==false)
      {
        ////////// find closest point to current vertex //////////
        // true => record closest face in vertex class
        // dummy => return closest face
        Face *cl=NULL;
        //int neighbor_count = -666;
        if (findClosestPtToVertex(&(*j),dummy_p,dummy_sqd,cl)==true)
        {
          j->setFace(cl);
          //j->setNeighborCount(neighbor_count);
        }
        //assert(dummy_sqd<1E10);
        assert(dummy_sqd==dummy_sqd);
      }
    }
    // track progress
    double progress = static_cast<double>(index++)*a;
    if (progress>goal || !distinguishable(progress,goal))
    {
      printf("%d%%..",static_cast<int>(goal*100));
      fflush(stdout);
      goal+=inc;
    }
  }
//  printf("100%%..");
//  fflush(stdout);
  cout << "complete.\n";
  cout.flush();
}

void keepMin (vector3 & min,Wm4::Vector3<double> q)
{
  if (q.X() < min.p[0]) min.p[0]=q.X(); 
  if (q.Y() < min.p[1]) min.p[1]=q.Y(); 
  if (q.Z() < min.p[2]) min.p[2]=q.Z(); 
}

void keepMax (vector3 & max,Wm4::Vector3<double> q)
{
  if (q.X() > max.p[0]) max.p[0]=q.X(); 
  if (q.Y() > max.p[1]) max.p[1]=q.Y(); 
  if (q.Z() > max.p[2]) max.p[2]=q.Z(); 
}

/** Compute axis-aligned bounding box of spherical cone.
 * \param[in] pt Point of interest.
 * \param[in] n Normal vector at point of interest.
 * \param[in] cone_radius Spherical cone radius.
 * \return Bounding box along principal axes
 * [xmin,xmax,ymin,ymax,zmin,zmax].
 */

vec_d Container::getSphericalConeBoundingBox (vector3 const & pt,
                                              vector3 n,
                                              double const & cone_radius) const
{
  // Algorithm
  // First compute the bounding box
  // of a truncated cone (i.e. flat on end, not rounded)
  // by sweeping the cone normal vector
  // toward both ends of each cartesian axis.
  // To rotate the normal towards a principal axis,
  // compute a vector perpendicular to both the normal vector
  // and the principal axis via the cross product.
  // Use quaternions to rotate vector v ccw (right-hand rule)
  // around axis u a desired angle.
  // Negate the rotation angle to sweep towards the opposite end
  // of principal axis.
  //
  // Second update bounding box based on spherical cap on cone.
  // 
  // DEBUG
  bool print=false;
//  if (distinguishable(pt.p[0],6945.31147969)==false &&
//      distinguishable(pt.p[1],3532.15958884)==false &&
//      distinguishable(pt.p[2],4653.25547899)==false)
//    print = true;
  // DEBUG
  Controls & cs(Controls::instance());
  vector3 pos_x(1.0,0.0,0.0);
  vector3 pos_y(0.0,1.0,0.0);
  vector3 pos_z(0.0,0.0,1.0);
  // set normal length to 1.0 
  n *= 1.0/sqrt(n.dot(n));
//  if (distinguishable(sqrt(n.dot(n)),1.0)==true)
//  {
//    cout << "Container::getSphericalConeBoundingBox "
//          << "Error normal vector is not unit vector."
//          << " n ["
//          << n.p[0] << " "
//          << n.p[1] << " "
//          << n.p[2] << "]\n";
//    exit(1);
//  }
  // compute unit rotation axes
  // as vectors perpendicular to the plane that
  // contains the normal vector and x,y, or z axes
  vector3 axis_x = n.cross(pos_x);
  vector3 axis_y = n.cross(pos_y);
  vector3 axis_z = n.cross(pos_z);
  vector3 axes[3] = { axis_x/axis_x.length(), axis_y/axis_y.length(), axis_z/axis_z.length() };
  // cosine of angle between normal vector and each cartesian axis
  double x_val = n.dot(pos_x);
  double y_val = n.dot(pos_y);
  double z_val = n.dot(pos_z);
  // DEBUG
//  if (print==true)
//  {
//    cout << "\nContainer::getSphericalConeBoundingBox "
//          << "angle = " << cs.get_closest_point_angle()
//          << ", vertex point = ["
//          << pt.p[0] << " "
//          << pt.p[1] << " "
//          << pt.p[2] << "]"
//          << ", vertex normal = ["
//          << n.p[0] << " "
//          << n.p[1] << " "
//          << n.p[2] << "]\n";
//  }
  // DEBUG
  // set normal length equal to search cone radius
  n *= cone_radius;
  // DEBUG
//  if (print==true)
//  {
//    cout << "\nContainer::getSphericalConeBoundingBox "
//          << "scaled vertex normal = ["
//          << (n.p[0]+pt.p[0]) << " "
//          << (n.p[1]+pt.p[1]) << " "
//          << (n.p[2]+pt.p[2]) << "]\n";
//  }
  // DEBUG
  // add vertex location
  Wm4::Vector3<double> v( n.p[0], n.p[1], n.p[2]);
  Wm4::Vector3<double> p(pt.p[0],pt.p[1],pt.p[2]);
  vector3 min(1.0e30,1.0e30,1.0e30);
  vector3 max(0.0,0.0,0.0);
  keepMin(min,p);
  keepMax(max,p);
  double sign = 1.0;
  // rotate normal vector 
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<2;j++)
    {
//      if (distinguishable(sqrt(axes[i].dot(axes[i])),1.0)==true)
//      {
//        cout << "Container::getSphericalConeBoundingBox "
//              << "Error rotation axis is not unit vector."
//              << " axes[i] ["
//              << axes[i].p[0] << " "
//              << axes[i].p[1] << " "
//              << axes[i].p[2] << "]\n";
//        // << " cone_radius = 1.0"
//        exit(1);
//      }
      Wm4::Vector3<double>    u(axes[i].p[0],axes[i].p[1],axes[i].p[2]); 
      Wm4::Quaternion<double> q(u,sign*cs.get_closest_point_angle());
      q /= q.Length();
      Wm4::Vector3<double> rotated = q.Rotate(v)+p;
//      Wm4::Vector3<double> rotated_v = q.Rotate(v);
//      Wm4::Vector3<double> rotated = rotated_v+p;
      // DEBUG
//      if (print==true)
//      {
//        cout << "Container::getSphericalConeBoundingBox "
//              << "rotated normal = ["
//              << rotated.X() << " "
//              << rotated.Y() << " "
//              << rotated.Z() << "]\n";
//      }
//      vector3 diff(rotated.X()-pt.p[0],rotated.Y()-pt.p[1],rotated.Z()-pt.p[2]);
//      double length = sqrt(diff.dot(diff));
//      if (distinguishable(length,cone_radius)==true)
//      {
//        cout << "Container::getSphericalConeBoundingBox "
//              << "Error vector length changed during rotation."
//              << " cone_radius = " << cone_radius
//              // << " cone_radius = 1.0"
//              << ", length of rotated vector = " << length << endl;
//        exit(1);
//      }
//      double angle = acos(v.Dot(rotated_v)/(v.Length()*rotated_v.Length()));
//      if (distinguishable(angle,cs.get_closest_point_angle())==true)
//      {
//        cout << "Container::getSphericalConeBoundingBox "
//              << "Rotation angle does not match desired angle."
//              << " cs.get_closest_point_angle() = " << cs.get_closest_point_angle()
//              // << " cone_radius = 1.0"
//              << ", angle of rotation = " << angle << endl;
//        exit(1);
//      }
      // DEBUG
      keepMin(min,rotated);
      keepMax(max,rotated);
      // update sign
      sign *= -1.0;
    }
  }
  // determine if spherical cone cap is intersected by cartesian axis
  // if so update bounding box
  // recognizing that the limiting extent of the cone cap
  // in the principal direction a distance of
  //   cone origin + cone_radius if normal points towards + axis
  //   cone origin - cone_radius if normal points towards - axis
  //
  if      (x_val> cs.get_closest_point_cosine())
  {
    // DEBUG
    if (print==true)
    {
      cout << "Container::getSphericalConeBoundingBox +x cone cap intersect\n";
    }
    // DEBUG
    if (( pt.p[0]+cone_radius) > max.p[0])
      max.p[0] = pt.p[0]+cone_radius;
  }
  else if (x_val<-cs.get_closest_point_cosine())
  {
    // DEBUG
    if (print==true)
    {
      cout << "Container::getSphericalConeBoundingBox -x cone cap intersect\n";
    }
    // DEBUG
    if (( pt.p[0]-cone_radius) < min.p[0])
      min.p[0] = (pt.p[0]-cone_radius);
  }
  if      (y_val> cs.get_closest_point_cosine())
  {
    // DEBUG
    if (print==true)
    {
      cout << "Container::getSphericalConeBoundingBox +y cone cap intersect\n";
    }
    // DEBUG
    if (( pt.p[1]+cone_radius) > max.p[1])
      max.p[1] = (pt.p[1]+cone_radius);
  }
  else if (y_val<-cs.get_closest_point_cosine())
  {
    // DEBUG
    if (print==true)
    {
      cout << "Container::getSphericalConeBoundingBox -y cone cap intersect\n";
    }
    // DEBUG
    if (( pt.p[1]-cone_radius) < min.p[1])
      min.p[1] = (pt.p[1]-cone_radius);
  }
  if      (z_val> cs.get_closest_point_cosine())
  {
    // DEBUG
    if (print==true)
    {
      cout << "Container::getSphericalConeBoundingBox +z cone cap intersect\n";
    }
    // DEBUG
    if (( pt.p[2]+cone_radius) > max.p[2])
      max.p[2] = (pt.p[2]+cone_radius);
  }
  else if (z_val<-cs.get_closest_point_cosine())
  {
    // DEBUG
    if (print==true)
    {
      cout << "Container::getSphericalConeBoundingBox -z cone cap intersect\n";
    }
    // DEBUG
    if (( pt.p[2]-cone_radius) < min.p[2])
      min.p[2] = (pt.p[2]-cone_radius);
  }
  // build search region
  vec_d sr(6,0.0);
  sr[0] = min.p[0];
  sr[1] = max.p[0];
  sr[2] = min.p[1];
  sr[3] = max.p[1];
  sr[4] = min.p[2];
  sr[5] = max.p[2];
  return sr;
}

/** Determine if input bounding box is fully inside the input search region.
 * \param[in] sr Search region [xmin,xmax,ymin,ymax,zmin,zmax].
 * \param[in] bb Bounding box [xmin,xmax,ymin,ymax,zmin,zmax].
 * \return True if bounding box is fully inside search region; false otherwise.
 */

bool Container::boundingBoxFullyInSearchRegion(vec_d const & sr,vec_d const & bb) const
{
  Controls & cs(Controls::instance());
  // assume both sr and bb of form [xmin,xmax,ymin,ymax,zmin,zmax]
  // check x axis
  if ((bb[1]-cs.get_my_double_epsilon())<sr[1] && (bb[0]+cs.get_my_double_epsilon())>sr[0])
  {
    // check y axis
    if ((bb[3]-cs.get_my_double_epsilon())<sr[3] && (bb[2]+cs.get_my_double_epsilon())>sr[2])
    {
      // check z axis
      if ((bb[5]-cs.get_my_double_epsilon())<sr[5] && (bb[4]+cs.get_my_double_epsilon())>sr[4])
      {
        return true;
      }
    }
  }
  return false;
}

/** Find the closest point to a tile's barycenter.
 * \param[in] pt Barycenter of face tile for which closest point is searched.
 * \param[in] f Face on which tile is located.
 * \param[out] p Closest point position, if found.
 * \param[out] sqd Squared distance between barycenter
 * and closest point, if found.
 * \param[out] ncl Parent face of closest point, if found.
 * \return True if closest point found; otherwise false.
 */

bool Container::findClosestPtToBarycenter (vector3 const & pt,Face * const f,
                                           vector3 & p,
                                           double & sqd,
                                           Face *&ncl) const
{
  Controls & cs(Controls::instance());
  int num_faces=0,num_leaves = 0,total_faces_checked=0;
  bool found = false;
  double radius_step = 0.0;
  if (cs.get_number_radius_steps()>0)
  {
    radius_step = (sqrt(cs.get_search_radius_sq())-cs.get_min_search_cone_radius())
        /static_cast<double>(cs.get_number_radius_steps());
  }
  // DEBUG
  //bool flag = false;
  //if (distinguishable(pt.p[0],4741.524,1E-5)==false &&
  //    distinguishable(pt.p[1],-2525.161,1E-5)==false &&
  //    distinguishable(pt.p[2],-1302.223,1E-5)==false)
  //{
  //  flag = true;
  //}
  //if (flag==true)
  //{
  //  cout << "\n\nBEGIN Container::findClosestPtToBarycenter" << endl;
  //  cout << "Container::findClosestPtToBarycenter: point = ";
  //  pt.print(cout);
  //  cout << endl;
  //}
  //else
  //{
  //  return false;
  //}
  // DEBUG
  // for search cone radius small to max
  sqd=1E30;
  for (int i=0;i<(cs.get_number_radius_steps()+1);++i)
  {
    double search_cone_radius = cs.get_min_search_cone_radius() 
                       + radius_step*static_cast<double>(i);
    // DEBUG
    //if (flag==true)
    //{
    //  cout << "Container::findClosestPtToBarycenter: "
    //        << "i = " << i
    //        << ", search_cone_radius = " << search_cone_radius << endl;
    //  vector3 nn = *f->getNormal();
    //  cout << "Container::findClosestPtToBarycenter: vertex normal = ["
    //        << nn.p[0] << " "
    //        << nn.p[1] << " "
    //        << nn.p[2] << "]\n";
    //}
    // DEBUG
    // calculate search region (i.e. bounding box of search cone)
    //vec_d sr = getTruncatedConeBoundingBox(pt,*f->getNormal(),search_cone_radius);
    vec_d sr = getSphericalConeBoundingBox(pt,*f->getNormal(),search_cone_radius);
    // DEBUG
    //if (flag==true)
    //{
    //  cout << "Container::findClosestPtToBarycenter: "
    //        << " sr ["
    //        << sr[0] << " "
    //        << sr[1] << " "
    //        << sr[2] << " "
    //        << sr[3] << " "
    //        << sr[4] << " "
    //        << sr[5] << "]\n";
    //}
    // DEBUG
    // collect faces from search region 
    // (possible optimization is not check same face twice)
    Octree_Visitor_Face visitor(Vector3r(sr[0],sr[2],sr[4]),
                                Vector3r(sr[1],sr[3],sr[5]));
    octree->visit( visitor );
    num_faces += visitor.num_faces();
    num_leaves += visitor.num_leaves();
    // DEBUG
    //if (flag==true)
    //{
    //  cout << "Container::findClosestPtToBarycenter: "
    //        << "# faces to check = " << visitor.num_faces()
    //        << ", # leaves checked = " << visitor.num_leaves() << endl;
    //}
    // DEBUG
    // find closest point to vertex among faces
    sqd=1E30;
    int faces_checked=0;
    found = findClosestPtToBarycenterAmongFaces(pt,f,visitor.mybegin(),
                                                visitor.myend(),search_cone_radius,p,sqd,
                                                ncl,faces_checked);
    total_faces_checked += faces_checked;
    // if found closest among collected faces
    if (found==true)
    {
      // DEBUG
      //if (flag==true)
      //{
      //  cout << "Container::findClosestPtToBarycenter: "
      //        << "sep_dis = " << sqrt(sqd) << endl;
      //}
      // DEBUG
      // already exhaustively searched all faces
      // in the search region defined by sd<search_cone_radius
      // So if the separating distance found is less than
      // search_cone_radius then it must be the closest
      // because we checked ALL faces a distance search_cone_radius
      // away from pt and closer. In other words, if there
      // was a closer one then we would have already found it.
      if (sqd < (search_cone_radius*search_cone_radius)) break;
    }
    // DEBUG
    //else
    //{
    //  if (flag==true)
    //  {
    //    cout << "Container::findClosestPtToBarycenter: "
    //          << "closest point not found.\n";
    //  }
    //}
    // DEBUG
  }
  // DEBUG
  //if (flag==true)
  //{
  //  if (ncl!=NULL)
  //  {
  //    cout << "\n closest_point ["
  //          << p.p[0] << " "
  //          << p.p[1] << " "
  //          << p.p[2] << "], d = " << sqrt(sqd) << endl;
  //    ncl->print(cout);
  //    cout << endl;
  //    //exit(1);
  //  }
  //  else
  //  {
  //    cout << "\n closest_point not found, d = ???\n";
  //    cout << " found = " << found << "\n";
  //  }
  //  //exit(1);
  //}
  // DEBUG
  // update statistics of count of face and boxes used in search
  Log::instance().updateClosestPtStats(num_faces,num_leaves,total_faces_checked);
  return found;
}

///** Ray trace from vertex to discover the number of neighboring objects.
// * \param[in] v Vertex for which nieghboring objects are counted.
// * \return Number of objects around vertex.
// */
//
//int Container::calculateVertexNeighborCount (Vertex const * const v) const
//{
//  Controls & cs(Controls::instance());
//  o_set neighbors;
//  double search_cone_radius = sqrt(cs.get_search_radius_sq());
//  // collect rays defined by search cone
//  std::vector<vector3> rays;
//  getSphericalConeRays(*v->getPos(),v->getNormal(),search_cone_radius,rays);
//  // for each ray
//  for (std::vector<vector3>::iterator j=rays.begin();j!=rays.end();j++)
//  {
//    // find intersected faces along ray
//    face_grp fg = Nice::instance().findIntFacesAlongRay(v,*v->getPos(),*j,false);
//    Object const * closest_object = NULL;
//    double min_sqd = 1E30;
//    // for each intersected face
//    for (fp_it k=fg.crossed_faces.begin();k!=fg.crossed_faces.end();++k)
//    {
//      // compute squared distance to face
//      vector3 closest_point;
//      double sqd;
//      findClosestPtInFaceToLocation(*v->getPos(),*k,closest_point,sqd);
//      // if distance is least
//      if (sqd<min_sqd)
//      {
//        // record distance and object
//        min_sqd = sqd;
//        closest_object = (*k)->getObject();
//      }
//    }
//    for (fp_it k=fg.edge_faces.begin();k!=fg.edge_faces.end();++k)
//    {
//      // compute squared distance to face
//      vector3 closest_point;
//      double sqd;
//      findClosestPtInFaceToLocation(*v->getPos(),*k,closest_point,sqd);
//      // if distance is least
//      if (sqd<min_sqd)
//      {
//        // record distance and object
//        min_sqd = sqd;
//        closest_object = (*k)->getObject();
//      }
//    }
//    if (closest_object!=NULL) neighbors.insert(closest_object);
//  }
//  // return number of unique objects
//  return neighbors.size();
//}

/** Ray trace from vertex to discover the number of neighboring objects.
 * \param[in] v Vertex for which nieghboring objects are counted.
 * \return Number of objects around vertex.
 */

int Container::calculateVertexNeighborCount (Vertex const * const v,vector3 const * pos,vector3 const & normal) const
{
  Controls & cs(Controls::instance());
  o_set neighbors;
  double search_cone_radius = sqrt(cs.get_search_radius_sq());
  // collect rays defined by search cone
  std::vector<vector3> rays;
  //getSphericalConeRays(*v->getPos(),v->getNormal(),search_cone_radius,rays);
  getSphericalConeRays(*pos,normal,search_cone_radius,rays);
  // for each ray
  for (std::vector<vector3>::iterator j=rays.begin();j!=rays.end();j++)
  {
    // find intersected faces along ray
    //face_grp fg = Nice::instance().findIntFacesAlongRay(v,*v->getPos(),*j,false);
    face_grp fg = Nice::instance().findIntFacesAlongRay(v,*pos,*j,false);
    Object const * closest_object = NULL;
    double min_sqd = 1E30;
    // for each intersected face
    for (fp_it k=fg.crossed_faces.begin();k!=fg.crossed_faces.end();++k)
    {
      // compute squared distance to face
      vector3 closest_point;
      double sqd;
      //findClosestPtInFaceToLocation(*v->getPos(),*k,closest_point,sqd);
      findClosestPtInFaceToLocation(*pos,*k,closest_point,sqd);
      // if distance is least
      if (sqd<min_sqd)
      {
        // record distance and object
        min_sqd = sqd;
        closest_object = (*k)->getObject();
      }
    }
    for (fp_it k=fg.edge_faces.begin();k!=fg.edge_faces.end();++k)
    {
      // compute squared distance to face
      vector3 closest_point;
      double sqd;
      //findClosestPtInFaceToLocation(*v->getPos(),*k,closest_point,sqd);
      findClosestPtInFaceToLocation(*pos,*k,closest_point,sqd);
      // if distance is least
      if (sqd<min_sqd)
      {
        // record distance and object
        min_sqd = sqd;
        closest_object = (*k)->getObject();
      }
    }
    if (closest_object!=NULL) neighbors.insert(closest_object);
  }
  // return number of unique objects
  return neighbors.size();
}

/** Compute rays along perimeter of cone originating at vertex.
 * \param[in] pt Point of interest.
 * \param[in] n Normal vector at point of interest.
 * \param[in] cone_radius Spherical cone radius.
 * \param[out] rays Endpoints of rays.
 */

void Container::getSphericalConeRays (vector3 const & pt,
                                      vector3 n,
                                      double const & cone_radius,
                                      std::vector<vector3> & rays) const
{
  // Algorithm
  // Sweeping the cone normal vector
  // toward both ends of each cartesian axis.
  // To rotate the normal towards a principal axis,
  // compute a vector perpendicular to both the normal vector
  // and the principal axis via the cross product.
  // Use quaternions to rotate vector v ccw (right-hand rule)
  // around axis u a desired angle.
  // Negate the rotation angle to sweep towards the opposite end
  // of principal axis.
  Controls & cs(Controls::instance());
  vector3 pos_x(1.0,0.0,0.0);
  vector3 pos_y(0.0,1.0,0.0);
  vector3 pos_z(0.0,0.0,1.0);
  // set normal length to 1.0 
  n *= 1.0/sqrt(n.dot(n));
  // compute unit rotation axes
  // as vectors perpendicular to the plane that
  // contains the normal vector and x,y, or z axes
  vector3 axis_x = n.cross(pos_x);
  vector3 axis_y = n.cross(pos_y);
  vector3 axis_z = n.cross(pos_z);
  vector3 axes[3] = { axis_x/axis_x.length(), axis_y/axis_y.length(), axis_z/axis_z.length() };
  //// cosine of angle between normal vector and each cartesian axis
  //double x_val = n.dot(pos_x);
  //double y_val = n.dot(pos_y);
  //double z_val = n.dot(pos_z);
  // set normal length equal to search cone radius
  n *= cone_radius;
  // add vertex location
  Wm4::Vector3<double> v( n.p[0], n.p[1], n.p[2]);
  Wm4::Vector3<double> p(pt.p[0],pt.p[1],pt.p[2]);
  rays.push_back(n+pt);
  //vector3 min(1.0e30,1.0e30,1.0e30);
  //vector3 max(0.0,0.0,0.0);
  //keepMin(min,p);
  //keepMax(max,p);
  double sign = 1.0;
  // rotate normal vector 
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<2;j++)
    {
      Wm4::Vector3<double>    u(axes[i].p[0],axes[i].p[1],axes[i].p[2]); 
      Wm4::Quaternion<double> q(u,sign*cs.get_closest_point_angle());
      q /= q.Length();
      Wm4::Vector3<double> rotated = q.Rotate(v)+p;
      rays.push_back(vector3(rotated[0],rotated[1],rotated[2]));
      //keepMin(min,rotated);
      //keepMax(max,rotated);
      // update sign
      sign *= -1.0;
    }
  }
}


/** Find the closest point to a vertex.
 * \param[in] v Vertex for which closest point is searched.
 * \param[out] p Closest point position, if found.
 * \param[out] sqd Squared distance between vertex
 * and closest point, if found.
 * \param[out] ncl Parent face of closest point, if found.
 * \param[out] neighbor_count Number of different objects encountered in search for closest face.
 * \return True if closest point found; otherwise false.
 */

bool Container::findClosestPtToVertex (Vertex const * const v,
                                       vector3 & p,
                                       double & sqd,
                                       Face *&ncl) const
{
  Controls & cs(Controls::instance());
  int num_faces=0,num_leaves = 0,total_faces_checked=0;
  bool found = false;
  double radius_step = 0.0;
  if (cs.get_number_radius_steps()>0)
  {
    radius_step = (sqrt(cs.get_search_radius_sq())-cs.get_min_search_cone_radius())
        /static_cast<double>(cs.get_number_radius_steps());
  }
  // for search cone radius small to max
  //// DEBUG
  //if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
  //{
  //  cout << "\n\nBEGIN Container::findClosestPtToVertex" << endl;
  //  v->print(cout);
  //}
  //else
  //{
  //  return false;
  //}
  //// DEBUG
  // DEBUG
//  if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
//  {
//    // for each object
//    for (o_cit i=o.begin();i!=o.end();++i)
//    {
//      if ((*i).getName()!=TARGET_FACE_NAME) continue;
//      // for each face in object
//      for (f_cit j=i->f.begin();j!=i->f.end();++j)
//      {
//        if ((*j).getIndex()!=TARGET_FACE_INDEX) continue;
//          cout << "\nContainer::findClosestPtToVertex: target face\n";
//          j->print(cout);
//          // calculate bounding box of face
//          Vector3r lower,upper;
//          j->getBoundingBox(lower,upper);
//          cout << "Container::findClosestPtToVertex: target face bounding box\n"
//                << " bb ["
//                << lower.getX() << " "
//                << upper.getX() << " "
//                << lower.getY() << " "
//                << upper.getY() << " "
//                << lower.getZ() << " "
//                << upper.getZ() << "]\n";
//      }
//    }
//  }
  // DEBUG
  sqd=1E30;
  for (int i=0;i<(cs.get_number_radius_steps()+1);++i)
  {
    double search_cone_radius = cs.get_min_search_cone_radius() 
          + radius_step*static_cast<double>(i);
    //// DEBUG
    //if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
    //{
    //  cout << "Container::findClosestPtToVertex: "
    //        << "i = " << i
    //        << ", search_cone_radius = " << search_cone_radius << endl;
    //  vector3 nn = v->getNormal();
    //  cout << "Container::findClosestPtToVertex: vertex normal = ["
    //        << nn.p[0] << " "
    //        << nn.p[1] << " "
    //        << nn.p[2] << "]\n";
    //}
    //// DEBUG
    // calculate search region (i.e. bounding box of search cone)
    //vec_d sr = getTruncatedConeBoundingBox(*v->getPos(),v->getNormal(),search_cone_radius);
    vec_d sr = getSphericalConeBoundingBox(*v->getPos(),v->getNormal(),search_cone_radius);
    //// DEBUG
    //if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
    //{
    //  cout << "Container::findClosestPtToVertex: "
    //        << " sr ["
    //        << sr[0] << " "
    //        << sr[1] << " "
    //        << sr[2] << " "
    //        << sr[3] << " "
    //        << sr[4] << " "
    //        << sr[5] << "]\n";
    //}
    //// DEBUG
    // collect faces from search region 
    // (possible optimization is not check same face twice)
    int faces_checked=0;
    Octree_Visitor_Face visitor(Vector3r(sr[0],sr[2],sr[4]),
                                Vector3r(sr[1],sr[3],sr[5]));
    octree->visit( visitor );
    num_faces += visitor.num_faces();
    num_leaves += visitor.num_leaves();
      //// DEBUG
      //if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
      //{
      //  cout << "Container::findClosestPtToVertex: "
      //        << "# faces to check = " << visitor.num_faces()
      //        << ", # leaves checked = " << visitor.num_leaves() << endl;
      //}
      //// DEBUG
    // find closest point to vertex among faces
    if (findClosestPtToVertexAmongFaces(v,visitor.mybegin(),
                                        visitor.myend(),search_cone_radius,p,sqd,
                                        ncl,faces_checked)==true) found=true;
    total_faces_checked += faces_checked;
    // if found closest among collected faces
    if (found==true)
    {
      //// DEBUG
      //if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
      //{
      //  cout << "Container::findClosestPtToVertex: "
      //        << "sep_dis = " << sqrt(sqd) << endl;
      //}
      //// DEBUG
      //if (cs.get_count_neighbors()==true)
      //{
      //  neighbor_count = countObjects(visitor.mybegin(),visitor.myend());
      //}
      if (sqd < (search_cone_radius*search_cone_radius)) break;
    }
    //// DEBUG
    //else
    //{
    //  if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
    //  {
    //    cout << "Container::findClosestPtToVertex: "
    //          << "closest point not found.\n";
    //  }
    //}
    //// DEBUG
  }
  //// DEBUG
  //if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
  //{
  //  if (ncl!=NULL)
  //  {
  //    cout << "\n closest_point ["
  //          << p.p[0] << " "
  //          << p.p[1] << " "
  //          << p.p[2] << "], d = " << sqrt(sqd) << endl;
  //    ncl->print(cout);
  //    cout << endl;
  //    v->print(cout);
  //    //exit(1);
  //  }
  //  else
  //  {
  //    cout << "\n closest_point not found, d = ???\n";
  //  }
  //  exit(1);
  //}
  //// DEBUG
  // update statistics of count of face and boxes used in search
  Log::instance().updateClosestPtStats(num_faces,num_leaves,total_faces_checked);
  return found;
}

/** Determine if closest point is inside search cone. 
 * \param[in] pt Point of interest.
 * \param[in] p Closest point of interest.
 * \param[in] n Vertex normal vector.
 * \param[in] sqd_sep_dist Squared distance between point and vertex.
 * \return True if point is inside vertex search cone; false otherwise;
 */

bool Container::closestPtIsInSearchCone (vector3 const & pt,
                                         vector3 const & p,
                                         vector3 const & n,
                                         double const & sqd_sep_dist) const
{
  Controls & cs(Controls::instance());
  double a = n.dot(n);
  // compute separation vector
  vector3 sep_vec(p-pt);
  // compute cosine of angle between outward normal and separation vector
  // which is equal to dot product of vectors divided by vector magnitudes
  double sign = sep_vec.dot(n);
  if (sign>0.0)
  {
    double cos_angle2 = sign*sign/(a*sep_vec.dot(sep_vec));
    // is closest point located within search cone
    // with angle window as defined in controls.cc?
    // OR located very close to current vertex?
    if ( (cos_angle2 > cs.get_closest_point_cosine()*cs.get_closest_point_cosine())
         || (sqd_sep_dist<cs.get_small_ecw_threshold()))
    {
      // if square of extracellular width is less than max allowed
      if (sqd_sep_dist<(cs.get_search_radius_sq()))
      {
        return true;
      }
    }
  }
  return false;
}

/** Find closest point to tile's barycenter among input faces.
 * \param[in] pt Barycenter of face tile for which closest point is searched.
 * \param[in] f Face on which tile is located.
 * \param[in] begin Iterator pointing to first face in vector of faces to check.
 * \param[in] end Iterator pointing to one past the last face in vector of faces to check.
 * \param[out] p Closest point position, if found.
 * \param[out] squareD Squared distance between vertex
 * and closest point, if found.
 * \param[out] ncl Parent face of closest point, if found.
 * \param[out] faces_checked Cumulative number of faces checked.
 * \return True if closest point found; otherwise false.
 */

bool Container::findClosestPtToBarycenterAmongFaces (vector3 const & pt,
                                                     Face * const f,
                                                     fp_cit begin,fp_cit end,
                                                     double const & cone_radius,
                                                     vector3 & p,
                                                     double & squareD,
                                                     Face * & ncl,
                                                     int & faces_checked) const
{
  // DEBUG
  //bool flag = false;
  //if (distinguishable(pt.p[0],4741.524,1E-5)==false &&
  //    distinguishable(pt.p[1],-2525.161,1E-5)==false &&
  //    distinguishable(pt.p[2],-1302.223,1E-5)==false)
  //{
  //  flag = true;
  //  cout << "Container::findClosestPtToBarycenterAmongFaces: This face:\n";
  //  f->print(cout);
  //}
  //if (flag==false)
  //{
  //  return false;
  //}
  // DEBUG
  // signal if closest point was found
  bool closest_point_found = false;
  // if candidate faces were NOT found
  if (begin==end) return closest_point_found;
  // get vertex normal
  vector3 n(*f->getNormal());
  vector3 closest_point;
  double min_d2;
  // NOTE the niceness of barycenters should be calculated.
  // For now, I assume no nonnice vertices or barycenters.
  // invert vertex normal if vertex is not nice
  //  if (Nice::instance().vertexIsNice(v)==false)
  //  {
  //    for (int i=0;i<3;++i) n.p[i]=-n.p[i];
  //  }
  faces_checked=0;
  // for any number of candidate faces
  // (including the case of a single input candidate face)
  // for each candidate face
  for (fp_cit j=begin;j!=end;++j)
  {
    // DEBUG
    //if (flag==true)
    //{
    //  cout << "Container::findClosestPtToBarycenterAmongFaces: test face:\n";
    //  (*j)->print(cout);
    //}
    // DEBUG
    (*j)->clearFlag();
    // reject face of tile
    if (*j==f)
    {
      //if (flag==true)
      //{
      //  cout << "Container::findClosestPtToBarycenterAmongFaces: face rejected as same.\n";
      //}
      continue;
    }
    // reject faces on wrong side of orthogonal plane to vertex normal
    if (faceLiesOppositeToNormal(pt,*j,n)==true)
    {
      //if (flag==true)
      //{
      //  cout << "Container::findClosestPtToBarycenterAmongFaces: face rejected as wrong side.\n";
      //}
      continue;
    }
    // reject coplanar faces
    if (Intersecting_Faces::instance().facesCoplanar(f,*j)==true)
    {
      //if (flag==true)
      //{
      //  cout << "Container::findClosestPtToBarycenterAmongFaces: reject coplanar.\n";
      //}
      continue;
    }
    // if current vertex, v, lies inside sphere
    faces_checked++;
    // then process face
    findClosestPtInFaceToLocation(pt,*j,closest_point,min_d2);
    // is point closest overall
    if (distinguishable(min_d2,squareD)==false)
    {
      // keep face with lower face index
      if ((*j)->getIndex()<ncl->getIndex())
      {
        ncl=*j;
        closest_point_found=true;
        // DEBUG
        //if (flag==true)
        //{
        //  cout << "Container::findClosestPtToBarycenterAmongFaces: aaa face updated.\n"; 
        //}
        // DEBUG
      }
    }
    else if (min_d2 < squareD)
    {
      // if closest point is within search cone
      if (closestPtIsInSearchCone(pt,closest_point,n,min_d2)==true)
      {
        // save face
        ncl=*j;
        // save squared extracellular width
        squareD = min_d2;
        // save point
        p = closest_point;
        closest_point_found=true;
        // DEBUG
        //if (flag==true)
        //{
        //  cout << "Container::findClosestPtToBarycenterAmongFaces: bbb face updated.\n"; 
        //}
        // DEBUG
      }
    }
    // ADDED 10 Nov 2008
    //vector3 nn = *(*j)->getNormal();
    findPtInFaceWhereNormalInt(pt,n,cone_radius,*j,closest_point,min_d2);
    //findClosestPtInFaceToLocation(pt,*j,closest_point,min_d2);
    // is point closest overall
    if (distinguishable(min_d2,squareD)==false)
    {
      // keep face with lower face index
      if ((*j)->getIndex()<ncl->getIndex())
      {
        ncl=*j;
        closest_point_found=true;
        // DEBUG
        //if (flag==true)
        //{
        //  cout << "Container::findClosestPtToBarycenterAmongFaces: ccc face updated.\n"; 
        //}
        // DEBUG
      }
    }
    else if (min_d2 < squareD)
    {
      // if closest point is within search cone
      if (closestPtIsInSearchCone(pt,closest_point,n,min_d2)==true)
      {
        // save face
        ncl=*j;
        // save squared extracellular width
        squareD = min_d2;
        // save point
        p = closest_point;
        closest_point_found=true;
        // DEBUG
        //if (flag==true)
        //{
        //  cout << "Container::findClosestPtToBarycenterAmongFaces: ddd face updated.\n"; 
        //}
        // DEBUG
      }
    }
  }
  return closest_point_found;
}

/** Count number of different objects contained in collection of input faces.
 *
 * \param[in] begin Iterator pointing to first face in vector of faces to check.
 * \param[in] end Iterator pointing to one past the last face in vector of faces to check.
 * \return Number of different objects.
 */

int Container::countObjects (fp_cit begin,fp_cit end) const
{
  // if candidate faces were NOT found
  if (begin==end)
  {
    cout << "\n\nContainer::countObject: ERROR. Zero candidate faces but closest face found.!\n";
    exit(1);
    //return 0;
  }
  // store pointer to face parent object in set
  o_set unique_objects;
  // for each candidate face
  for (fp_cit j=begin;j!=end;++j)
  {
    unique_objects.insert((*j)->getObject());
  }
  return unique_objects.size();
}

/** Find closest point to vertex among input faces.
 *
 * \param[in] v Vertex of interest.
 * \param[in] begin Iterator pointing to first face in vector of faces to check.
 * \param[in] end Iterator pointing to one past the last face in vector of faces to check.
 * \param[out] p Closest point position, if found.
 * \param[out] squareD Squared distance between vertex
 * and closest point, if found.
 * \param[out] ncl Parent face of closest point, if found.
 * \param[out] faces_checked Cumulative number of faces checked.
 * \return True if closest point found; otherwise false.
 */

bool Container::findClosestPtToVertexAmongFaces (Vertex const * const v,
                                                 fp_cit begin,fp_cit end,
                                              double const & cone_radius,
                                                 vector3 & p,
                                                 double & squareD,
                                                 Face * & ncl,
                                                 int & faces_checked) const
{
  // signal if closest point was found
  bool closest_point_found = false;
  // if candidate faces were NOT found
  if (begin==end) return closest_point_found;
  // get vertex normal
  vector3 n(v->getNormal());
  vector3 closest_point;
  double min_d2;
  // invert vertex normal if vertex is not nice
  if (Nice::instance().vertexIsNice(v)==false)
  {
    for (int i=0;i<3;++i) n.p[i]=-n.p[i];
  }
  faces_checked=0;
  // for any number of candidate faces
  // (including the case of a single input candidate face)
  // for each candidate face
  for (fp_cit j=begin;j!=end;++j)
  {
    //// DEBUG
    //if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
    //{
    //  (*j)->print(cout);
    //}
    //// DEBUG
    (*j)->clearFlag();
    // reject adjacent faces of vertex
    if (v->faceIsAdjacent(*j)==true) continue;
    // reject faces on wrong side of orthogonal plane to vertex normal
    if (faceLiesOppositeToNormal(*v->getPos(),*j,n)==true) continue;
    faces_checked++;
    findClosestPtInFaceToLocation(*v->getPos(),*j,closest_point,min_d2);
    // if new separation distance is the same as old separation distance
    if (distinguishable(min_d2,squareD)==false)
    {
      // keep face with lower face index
      if ((*j)->getIndex()<ncl->getIndex())
      {
        ncl=*j;
        closest_point_found=true;
        //// DEBUG
        //if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
        //{
        //  cout << "Container::findClosestPtToVertexAmongFaces: aaa face updated.\n"; 
        //}
        //// DEBUG
      }
    }
    // else if new separation distance is less than old separation distance
    else if (min_d2 < squareD)
    {
      // if closest point is within search cone
      if (closestPtIsInSearchCone(*v->getPos(),closest_point,n,min_d2)==true)
      {
        // save face
        ncl=*j;
        // save squared extracellular width
        squareD = min_d2;
        // save point
        p = closest_point;
        closest_point_found=true;
        //// DEBUG
        //if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
        //{
        //  cout << "Container::findClosestPtToVertexAmongFaces: bbb face updated.\n"; 
        //}
        //// DEBUG
      }
    }
    // ADDED 10 Nov 2008
    findPtInFaceWhereNormalInt(*v->getPos(),n,cone_radius,*j,closest_point,min_d2);
    // if new separation distance is the same as old separation distance
    if (distinguishable(min_d2,squareD)==false)
    {
      // keep face with lower face index
      if ((*j)->getIndex()<ncl->getIndex())
      {
        ncl=*j;
        closest_point_found=true;
        // DEBUG
        //if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
        //{
        //  cout << "Container::findClosestPtToVertexAmongFaces: ccc face updated.\n"; 
        //}
        // DEBUG
      }
    }
    // else if new separation distance is less than old separation distance
    else if (min_d2 < squareD)
    {
      // if closest point is within search cone
      if (closestPtIsInSearchCone(*v->getPos(),closest_point,n,min_d2)==true)
      {
        // save face
        ncl=*j;
        // save squared extracellular width
        squareD = min_d2;
        // save point
        p = closest_point;
        closest_point_found=true;
        // DEBUG
        //if (v->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
        //{
        //  cout << "Container::findClosestPtToVertexAmongFaces: ddd face updated.\n"; 
        //}
        // DEBUG
      }
    }
  }
  return closest_point_found;
}

/** Find closest point in face to input location.
 *
 * \param[in] pt Location of interest.
 * \param[in] face Face on which closest point to input location will be identified.
 * \param[out] closest_point Closest point position, if found.
 * \param[out] fSqrDistance Squared distance between input location
 * and closest point, if found.
 * \return True if closest point found; otherwise false.
 */

void Container::findClosestPtInFaceToLocation (vector3 const & pt,
                                               Face const * const face,
                                               vector3 & closest_point,
                                               double & fSqrDistance) const
{
  // Algorithm from ./doc/Wm4DistVector3Triangle3.cc
  //
  // Wild Magic Source Code
  // David Eberly
  // http://www.geometrictools.com
  // Copyright (c) 1998-2008
  //
  // This library is free software; you can redistribute it and/or modify it
  // under the terms of the GNU Lesser General Public License as published by
  // the Free Software Foundation; either version 2.1 of the License, or (at
  // your option) any later version.  The license is available for reading at
  // either of the locations:
  //     http://www.gnu.org/copyleft/lgpl.html
  //     http://www.geometrictools.com/License/WildMagicLicense.pdf
  //
  // Version: 4.0.1 (2007/05/06)

  vector3 kDiff(*(face->getVertex(0)->getPos()) - pt);
  vector3 kEdge0(*(face->getVertex(1)->getPos()) - *(face->getVertex(0)->getPos()));
  vector3 kEdge1(*(face->getVertex(2)->getPos()) - *(face->getVertex(0)->getPos()));
  double fA00 = kEdge0.dot(kEdge0);
  double fA01 = kEdge0.dot(kEdge1);
  double fA11 = kEdge1.dot(kEdge1);
  double fB0 = kDiff.dot(kEdge0);
  double fB1 = kDiff.dot(kEdge1);
  double fC = kDiff.dot(kDiff);
  double fDet = fabs(fA00*fA11-fA01*fA01);
  double fS = fA01*fB1-fA11*fB0;
  double fT = fA01*fB0-fA00*fB1;

  if (fS + fT <= fDet)
  {
    if (fS < 0.0)
    {
      if (fT < 0.0)  // region 4
      {
        if (fB0 < 0.0)
        {
          fT = 0.0;
          if (-fB0 >= fA00)
          {
            fS = 1.0;
            fSqrDistance = fA00+(2.0)*fB0+fC;
          }
          else
          {
            fS = -fB0/fA00;
            fSqrDistance = fB0*fS+fC;
          }
        }
        else
        {
          fS = 0.0;
          if (fB1 >= 0.0)
          {
            fT = 0.0;
            fSqrDistance = fC;
          }
          else if (-fB1 >= fA11)
          {
            fT = 1.0;
            fSqrDistance = fA11+(2.0)*fB1+fC;
          }
          else
          {
            fT = -fB1/fA11;
            fSqrDistance = fB1*fT+fC;
          }
        }
      }
      else  // region 3
      {
        fS = 0.0;
        if (fB1 >= 0.0)
        {
          fT = 0.0;
          fSqrDistance = fC;
        }
        else if (-fB1 >= fA11)
        {
          fT = 1.0;
          fSqrDistance = fA11+(2.0)*fB1+fC;
        }
        else
        {
          fT = -fB1/fA11;
          fSqrDistance = fB1*fT+fC;
        }
      }
    }
    else if (fT < 0.0)  // region 5
    {
      fT = 0.0;
      if (fB0 >= 0.0)
      {
        fS = 0.0;
        fSqrDistance = fC;
      }
      else if (-fB0 >= fA00)
      {
        fS = 1.0;
        fSqrDistance = fA00+(2.0)*fB0+fC;
      }
      else
      {
        fS = -fB0/fA00;
        fSqrDistance = fB0*fS+fC;
      }
    }
    else  // region 0
    {
      // minimum at interior point
      double fInvDet = (1.0)/fDet;
      fS *= fInvDet;
      fT *= fInvDet;
      fSqrDistance = fS*(fA00*fS+fA01*fT+(2.0)*fB0) +
            fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
    }
  }
  else
  {
    double fTmp0, fTmp1, fNumer, fDenom;

    if (fS < 0.0)  // region 2
    {
      fTmp0 = fA01 + fB0;
      fTmp1 = fA11 + fB1;
      if (fTmp1 > fTmp0)
      {
        fNumer = fTmp1 - fTmp0;
        fDenom = fA00-2.0f*fA01+fA11;
        if (fNumer >= fDenom)
        {
          fS = 1.0;
          fT = 0.0;
          fSqrDistance = fA00+(2.0)*fB0+fC;
        }
        else
        {
          fS = fNumer/fDenom;
          fT = 1.0 - fS;
          fSqrDistance = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
        }
      }
      else
      {
        fS = 0.0;
        if (fTmp1 <= 0.0)
        {
          fT = 1.0;
          fSqrDistance = fA11+(2.0)*fB1+fC;
        }
        else if (fB1 >= 0.0)
        {
          fT = 0.0;
          fSqrDistance = fC;
        }
        else
        {
          fT = -fB1/fA11;
          fSqrDistance = fB1*fT+fC;
        }
      }
    }
    else if (fT < 0.0)  // region 6
    {
      fTmp0 = fA01 + fB1;
      fTmp1 = fA00 + fB0;
      if (fTmp1 > fTmp0)
      {
        fNumer = fTmp1 - fTmp0;
        fDenom = fA00-(2.0)*fA01+fA11;
        if (fNumer >= fDenom)
        {
          fT = 1.0;
          fS = 0.0;
          fSqrDistance = fA11+(2.0)*fB1+fC;
        }
        else
        {
          fT = fNumer/fDenom;
          fS = 1.0 - fT;
          fSqrDistance = fS*(fA00*fS+fA01*fT+(2.0)*fB0) +
                fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
        }
      }
      else
      {
        fT = 0.0;
        if (fTmp1 <= 0.0)
        {
          fS = 1.0;
          fSqrDistance = fA00+(2.0)*fB0+fC;
        }
        else if (fB0 >= 0.0)
        {
          fS = 0.0;
          fSqrDistance = fC;
        }
        else
        {
          fS = -fB0/fA00;
          fSqrDistance = fB0*fS+fC;
        }
      }
    }
    else  // region 1
    {
      fNumer = fA11 + fB1 - fA01 - fB0;
      if (fNumer <= 0.0)
      {
        fS = 0.0;
        fT = 1.0;
        fSqrDistance = fA11+(2.0)*fB1+fC;
      }
      else
      {
        fDenom = fA00-2.0f*fA01+fA11;
        if (fNumer >= fDenom)
        {
          fS = 1.0;
          fT = 0.0;
          fSqrDistance = fA00+(2.0)*fB0+fC;
        }
        else
        {
          fS = fNumer/fDenom;
          fT = 1.0 - fS;
          fSqrDistance = fS*(fA00*fS+fA01*fT+(2.0)*fB0) +
                fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
        }
      }
    }
  }

  // account for numerical round-off error
  if (fSqrDistance < 0.0)
  {
    fSqrDistance = 0.0;
  }

  kEdge0 *= fS;
  kEdge1 *= fT;
  closest_point = *(face->getVertex(0)->getPos()) + kEdge0 + kEdge1;
}

/** Find intersection point of vertex normal vector and input face.
 *
 * \param[in] pt Location of interest.
 * \param[in] face Face on which intersection point of normal vector will be identified.
 * \param[out] intersection_point Intersection point position, if found.
 * \param[out] fSqrDistance Squared distance between input location
 * and intersection point, if found.
 * \return True if intersection point found; otherwise false.
 */

void Container::findPtInFaceWhereNormalInt (vector3 const & pt,
                                            vector3 n,
                                            double const & cone_radius,
                                            Face const * const face,
                                            vector3 & intersection_point,
                                            double & fSqrDistance) const
{
  // DEBUG
  //if (face->isMatch(TARGET_FACE_INDEX,TARGET_FACE_NAME)==true)
  //{
  //  cout << "Container::findPtInFaceWhereNormalInt: "
  //       << "normal ["
  //       << n.p[0] << " "
  //       << n.p[1] << " "
  //       << n.p[2] << "]\n";
  //}
  // DEBUG
  // set normal length to cone radius 
  //double scale = 1.0/sqrt(n.dot(n));
  //n *= scale;
  // DEBUG
//  if (face->isMatch(TARGET_FACE_INDEX,TARGET_FACE_NAME)==true)
//  {
//    cout << "Container::findPtInFaceWhereNormalInt: "
//         << "normal length = " << sqrt(n.dot(n)) << endl;
//  }
  // DEBUG
  n *= cone_radius/sqrt(n.dot(n));
  // DEBUG
  //if (face->isMatch(TARGET_FACE_INDEX,TARGET_FACE_NAME)==true)
  //{
  //  cout << "Container::findPtInFaceWhereNormalInt: "
  //       << "normal length = " << sqrt(n.dot(n))
  //       << ", normal ["
  //       << n.p[0] << " "
  //       << n.p[1] << " "
  //       << n.p[2] << "]\n";
  //}
  // DEBUG
  n += pt;
  // DEBUG
//  if (face->isMatch(TARGET_FACE_INDEX,TARGET_FACE_NAME)==true)
//  {
//    cout << "Container::findPtInFaceWhereNormalInt: "
//         << "segment end ["
//         << n.p[0] << " "
//         << n.p[1] << " "
//         << n.p[2] << "]\n";
//    double test1 = sqrt(n.dot(n));
//    vector3 boo(n);
//    boo *= 1.0/test1;
//    double yah = sqrt(boo.dot(boo));
//    cout << "Container::findPtInFaceWhereNormalInt: "
//         << "cone radius = " << cone_radius
//         << ", scaled length = " << test1
//         << ", descaled length = " << yah << endl; 
//  }
  // DEBUG
  // check for intersection
  result r;
  intersect_triangle3(&pt,&n,
                      face->getNormal(),
                      face->getVertex(0)->getPos(),
                      face->getVertex(1)->getPos(),
                      face->getVertex(2)->getPos(),r);
  // does point intersect polygon
  if (r.line_flag==true  && (r.poly_edge_flag==true || r.poly_flag==true) )
  {
    // compute point of intersection
    vector3 diff(n-pt);
    double d = face->getVertex(0)->getPos()->dot(*face->getNormal());
    double t = (d-pt.dot(*face->getNormal()))/(diff.dot(*face->getNormal()));
    diff *= t;
    intersection_point = pt+diff;
    vector3 seg(intersection_point-pt);
    fSqrDistance = seg.dot(seg); 
    // DEBUG
    //if (face->isMatch(TARGET_FACE_INDEX,TARGET_FACE_NAME)==true)
    //{
    //  cout << "Container::findPtInFaceWhereNormalInt: "
    //        << "Target face was intersected."
    //        << "intersection point = ["
    //        << intersection_point.p[0] << " "
    //        << intersection_point.p[1] << " "
    //        << intersection_point.p[2] << "]\n";
    //}
    // DEBUG
  }
  //// DEBUG
  //else
  //{
  //  if (face->isMatch(TARGET_FACE_INDEX,TARGET_FACE_NAME)==true)
  //  {
  //    cout << "Container::findPtInFaceWhereNormalInt: "
  //          << "Target face was NOT intersected.\n";
  //    cout << "Container::findPtInFaceWhereNormalInt: "
  //          << "segment origin ["
  //          << pt.p[0] << " "
  //          << pt.p[1] << " "
  //          << pt.p[2] << "], "
  //          << "segment end ["
  //          << n.p[0] << " "
  //          << n.p[1] << " "
  //          << n.p[2] << "]\n";
  //    cout << "Container::findPtInFaceWhereNormalInt: "
  //          << "line_flag = " << r.line_flag
  //          << ", poly_edge_flag = " << r.poly_edge_flag
  //          << ", ploy_flag = " << r.poly_flag << endl;
  //  }
  //}
  // DEBUG
}

/** Find the edge in model with smallest angle and return angle.
 * \return Smallest edge angle in radians.
 */

double Container::getMinEdgeAngle (void) const
{
  double min = 1E30;
  // for each object* in container
  for (o_cit i=o.begin();i!=o.end();++i)
  {
    // for each Edge* in object
    for (e_cit j=i->e.begin();j!=i->e.end();++j)
    {
      // record if angle is smallest so far
      double a = j->getAngle(); 
      if (a<min) min=a;
    }
  }
  return min;
}

/** Calculate and record limits of model along each prinicpal axis.
*/

void Container::boundWorld (void)
{
  //initialize mins and maxes
  Vertex * const v = &o[0].v[0];
  double xmin = *v->getCoord(0);
  double xmax = *v->getCoord(0);
  double ymin = *v->getCoord(1);
  double ymax = *v->getCoord(1);
  double zmin = *v->getCoord(2);
  double zmax = *v->getCoord(2);
  // calculate model limits
  // for each object
  for (o_cit i=o.begin();i!=o.end();++i)
  {
    // get range of object vertices
    double range[6];
    i->boundObject(&range[0]);
    if (range[1]>xmax) {xmax = range[1];}
    if (range[0]<xmin) {xmin = range[0];}
    if (range[3]>ymax) {ymax = range[3];}
    if (range[2]<ymin) {ymin = range[2];}
    if (range[5]>zmax) {zmax = range[5];}
    if (range[4]<zmin) {zmin = range[4];}
  }
  // Extend limits a little so model
  // not tagent to bounding box
  if (xmin<0) {world[0]=xmin*1.01;} else {world[0]=xmin*0.99;}
  if (xmax<0) {world[1]=xmax*0.99;} else {world[1]=xmax*1.01;}
  if (ymin<0) {world[2]=ymin*1.01;} else {world[2]=ymin*0.99;}
  if (ymax<0) {world[3]=ymax*0.99;} else {world[3]=ymax*1.01;}
  if (zmin<0) {world[4]=zmin*1.01;} else {world[4]=zmin*0.99;}
  if (zmax<0) {world[5]=zmax*0.99;} else {world[5]=zmax*1.01;}
}

/** Write description and summary of this class to output stream.
 * \param[in] target Pre-initialized output stream.
 */

void Container::writeSummary (std::ostream & target)
{
  target.width(30);
  target << std::left << "world bounds [xmin xmax ymin ymax zmin zmax]" << endl;
  target << "[" << world[0] << " " << world[1] << " " << world[2] << " "
        << world[3] << " " << world[4] << " " << world[5] << "]\n";
}

/** Return the world limit along input axis.
 * \param[in] axis Requested direction
 * 0,1,2,3,4,5 == -x,+x,-y,+y,-z,+z.
 * \return World limit along ith axis.
 */

double Container::getWorld (int const & axis) const
{
  assert(axis > -1 && axis < 6);
  return world[axis];
}

/** Calculate and record the normal vector of each vertex in the model.
*/

void Container::computeVertexNormals (void)
{
  // for each object in container
  for (o_it i=o.begin();i!=o.end();++i)
  {
    // for each vertex in object
    for (v_it j=i->v.begin();j!=i->v.end();++j)
    {
      j->setNormal(); 
    }
  }
}

/** Sort adjacent faces of each vertex in model.
*/

void Container::sortAdjacentFaces (void)
{
  // for each object
  for (o_it i=o.begin();i!=o.end();++i)
  {
    // for each vertex in object
    for (v_it j=i->v.begin();j!=i->v.end();++j)
    {
      // sort adjacent faces
      j->sortAdjacentFaces();
    }
  }
}

/** Check that all faces have flag set to false.
 * \param[in] s Error message if face found with flag set to true.
 */

void Container::checkFaces (std::string s)
{
  // for each object
  for (o_it i=o.begin();i!=o.end();++i)
  {
    // for each face in object
    for (f_it j=i->f.begin();j!=i->f.end();++j)
    {
      if (j->getFlag()==true)
      {
        cout << "\n" << s
              << " Error. Face flag is true.\n";
        j->print(cout);
        assert(j->getFlag()==false);
        exit(1);
      }
    }
  }
}

/** Determine if entire input face lies in hemispace opposite to normal vector at input location.
 * \param[in] pt Location of interest.
 * \param[in] face Face of interest.
 * \param[in] nn Normal vector of this face. 
 * \return True if entire face lies opposite to normal vector; false otherwise.
 */

bool Container::faceLiesOppositeToNormal(vector3 const & pt,Face const * const face,vector3 const & nn) const
{
  vector3 t0(*face->getVertex(0)->getPos()-pt); 
  vector3 t1(*face->getVertex(1)->getPos()-pt); 
  vector3 t2(*face->getVertex(2)->getPos()-pt); 
  return (nn.dot(t0)<0.0 && nn.dot(t1)<0.0 && nn.dot(t2)<0.0);
}

/** Print octree overlapping region.
*/

void Container::printRegionInOctree (void)
{
  Vector3r lower(7908.28430532,3708.83766659,7477.33967687);
  Vector3r upper(7968.06669069,3768.61532283,7532.54764615);
  // check that face is found in all octree cells overlapping bounding box
  Octree_Visitor_Check_Face visitor(NULL,lower,upper);
  octree->visit( visitor );
}
/** Check that all faces are present in appropriate octree cells.
*/

void Container::checkFacesInOctree (void)
{
  // for each object
  for (o_it i=o.begin();i!=o.end();++i)
  {
    // for each face in object
    for (f_it j=i->f.begin();j!=i->f.end();++j)
    {
      // calculate bounding box of face
      Vector3r lower,upper;
      j->getBoundingBox(lower,upper);
      // check that face is found in all octree cells overlapping bounding box
      Octree_Visitor_Check_Face visitor(&(*j),lower,upper);
      octree->visit( visitor );
    }
  }
}

/** Check if moved vertex has breached octree boundary. 
 *
 * \param[in] new_pos The new position of current vertex.
 * \return True if vertex has moved to a location outside of octree bounday,
 * false otherwise.
 */

bool Container::vertexOutsideOctreeBounds (vector3 const * const new_pos)
{
  Controls const & cs(Controls::instance());
  vector3  low(cs.get_octree_min_x(),
               cs.get_octree_min_y(),
               cs.get_octree_min_z());
  vector3 high(cs.get_octree_min_x()+cs.get_octree_width(),
               cs.get_octree_min_y()+cs.get_octree_width(),
               cs.get_octree_min_z()+cs.get_octree_width());
  return new_pos->p[0]<low.p[0]  ||  
         new_pos->p[0]>high.p[0] ||  
         new_pos->p[1]<low.p[1]  ||  
         new_pos->p[1]>high.p[1] ||  
         new_pos->p[2]<low.p[2]  ||  
         new_pos->p[2]>high.p[2];  
}

void Container::writeOrthVerts2Dreamm (std::vector<const Vertex *> & orthogonal_vertices) 
{
  char mystr[256];
  sprintf(mystr,"%sorth_verts.mdl",Controls::instance().get_output_data_dir().c_str());
  //std::string orth_verts_filename = "orth_verts.mdl";
  std::ofstream orth_verts_handle (mystr);
  sprintf(mystr,"%sreleases.mdl",Controls::instance().get_output_data_dir().c_str());
  //std::string releases_filename   = "releases.mdl";
  std::ofstream releases_handle (mystr);

  orth_verts_handle << "orth_verts OBJECT\n{\n";
  std::string body = "\n  NUMBER_TO_RELEASE = 1\n  SITE_DIAMETER = 0\n  RELEASE_PATTERN = rp\n}\n\n";
  int total_count = 0;
  std::set<std::string> molecule_set;
  std::set<std::string>::iterator i;

  // instantiate molecule releases (t=0,zero diffusion,single molecule)
  // centered at each orthogonal vertex
  //
  // for each orthogonal vertex
  for (cvp_cit j=orthogonal_vertices.begin();j!=orthogonal_vertices.end();++j)
  {

    vector3 const * p = (*j)->getPos();
    double sB = p->p[0]/1000.0;
    double sC = p->p[1]/1000.0;
    double sD = p->p[2]/1000.0;
    char name[256];
    sprintf(name,"rs_%d",total_count);
    char str[1024];
    sprintf(str,"%s SPHERICAL_RELEASE_SITE\n{\n",name);
    releases_handle << str;
    sprintf(str,"  LOCATION = [%g,%g,%g]\n",sB,sC,sD);
    releases_handle <<  str;
    sprintf(str,"  LIGAND = %s",(*j)->getObjectName().c_str());
    molecule_set.insert((*j)->getObjectName());
    releases_handle <<  str;
    releases_handle <<  body;
    sprintf(str,"  %s OBJECT %s {}\n",name,name);
    orth_verts_handle <<  str;
    total_count++;
  }
  orth_verts_handle << "}";
  orth_verts_handle.close();
  releases_handle.close();

  cout << "num different molecules (i.e. objects) found  = " << molecule_set.size() << endl;
  sprintf(mystr,"%smolecules.mdl",Controls::instance().get_output_data_dir().c_str());
  //std::string molecule_filename = "molecules.mdl";
  std::ofstream molecule_handle (mystr);
  for (i=molecule_set.begin();i!=molecule_set.end();i++)
  {
    molecule_handle << "DEFINE_MOLECULE "
          << *i << " {DIFFUSION_CONSTANT_3D=0}" << endl;
  }
  molecule_handle.close();
}

void Container::writeContoursSingleSection (std::vector<const Vertex*> & single_section, double z) 
{
  std::ofstream contour_handle;
  openContourFile(contour_handle,z);

  // organize Vertex * by object name
  std::multimap<std::string,const Vertex*>  map_s_v;
  std::multimap<std::string,const Vertex*>::iterator j,k;

  // for each vertex in section
  for (cvp_cit i=single_section.begin();i!=single_section.end();i++)
  {
    std::string name = (*i)->getObjectName();
    map_s_v.insert(std::make_pair(name,*i));
  }

  k = map_s_v.begin();
  std::string last = (*k).first;

  //double r = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
  //double g = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
  //double b = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
  double r = 1.0;
  double g = 0;
  double b = 1.0;

  //initializeContour(contour_handle,last);

  // for each orthogonal vertex in section
  // presorted by object name
  for (j=map_s_v.begin();j!=map_s_v.end();j++)
  {
    if ((*j).first==last)
    {
      initializeContour(contour_handle,last,r,g,b);
      addVertexToFile(contour_handle,(*j).second);
      finalizeContour(contour_handle);
    }
    else
    {
      //finalizeContour(contour_handle);
      last = (*j).first;
      //r = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
      //g = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
      //b = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
      r = 1.0;
      g = 0;
      b = 1.0;
      initializeContour(contour_handle,last,r,g,b);
      addVertexToFile(contour_handle,(*j).second);
      finalizeContour(contour_handle);
    }
  }
  //finalizeContour(contour_handle);
  contour_handle << "</Section>";
  contour_handle.close();
}

void Container::writeOrthVerts2Recon3D (std::vector<const Vertex*> & orthogonal_vertices) 
{
  std::multimap<double,const Vertex*>  map_d_v;
  std::multimap<double,const Vertex*>::iterator i,k;

  // for each orthogonal vertex
  for (cvp_cit j=orthogonal_vertices.begin();j!=orthogonal_vertices.end();++j)
  {
    const double z = *((*j)->getCoord(2));
    map_d_v.insert(std::make_pair(z,*j));
  }

  i = map_d_v.begin();
  double last = (*i).first;
  int myvert_count = 0;
  const double spacing = 50.0;
  std::vector<const Vertex*> single_section;

  // for each orthogonal vertex
  // presorted numerically, descreasing by z value
  for (i=map_d_v.begin();i!=map_d_v.end();i++)
  {
    double a = (*i).first / spacing;
    float rem = a - floor(a);
    if (fabs(rem) > Controls::instance().get_epsilon()) continue;

    ++myvert_count;
    if (distinguishable((*i).first,last))
    {
      if (!single_section.empty())
      {
        writeContoursSingleSection(single_section,last);
        single_section.clear();
      }
      last = (*i).first;
    }
    single_section.push_back((*i).second);
  }

  cout << "num aligned orth verts (written to file) = " << myvert_count << endl;
}

void Container::finalizeContour (std::ofstream & contour_handle)
{
  contour_handle << "	  \"/>\n</Transform>\n\n";
}

void Container::initializeContour (std::ofstream & contour_handle,const std::string & name,
                                   const double & r,const double & g,const double & b)
{
  contour_handle
    << "<Transform dim=\"0\"\n"
    << " xcoef=\" 0 1 0 0 0 0\"\n"
    << " ycoef=\" 0 0 1 0 0 0\">\n"
    << "<Contour name=\"" << name
    //<< "\" hidden=\"false\" closed=\"true\" simplified=\"true\" border=\"0.501961 0 1\" fill=\"1 1 0.72549\" mode=\"9\"\n"
    << "\" hidden=\"false\" closed=\"true\" simplified=\"true\" border=\""
    << r << " "
    << g << " "
    << b << "\" fill=\""
    << r << " "
    << g << " "
    << b << "\" mode=\"9\"\n"
    << " points=\"";
}

void Container::addVertexToFile (std::ofstream & contour_handle,const Vertex * vp)
{
  const double offset = 0.002;
  double const * a = vp->getCoord(0);
  double const * b = vp->getCoord(1);
  char str[256];
  //sprintf(str,"        %g %g,\n",(*a)/1000.0,(*b)/1000.0);
  sprintf(str,"        %g %g,\n",(*a)/1000.0-offset,(*b)/1000.0-offset);
  contour_handle << str;
  sprintf(str,"        %g %g,\n",(*a)/1000.0-offset,(*b)/1000.0+offset);
  contour_handle << str;
  sprintf(str,"        %g %g,\n",(*a)/1000.0+offset,(*b)/1000.0-offset);
  contour_handle << str;
  sprintf(str,"        %g %g,\n",(*a)/1000.0+offset,(*b)/1000.0+offset);
  contour_handle << str;
}

//const char img_data[3][256] = {
//{"\"A1-B_S12.160.jpg\" />"},
//{"\"A1-B_S12.160.jpg\" />"},
//{"\"A1-B_S12.160.jpg\" />"}}

const char img_data[101][256] = {
{"<Transform dim=\"3\"\n xcoef=\" -0.567106 0.908516  0.279141 0 0 0\"\n ycoef=\"  1.87787  -0.329688 0.8712   0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.160.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.865917 0.923752  0.283536 0 0 0\"\n ycoef=\"  1.49572  -0.314659 0.897645 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.161.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.831868 0.929041  0.260561 0 0 0\"\n ycoef=\"  1.27392  -0.293663 0.895238 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.162.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.768474 0.943801  0.247365 0 0 0\"\n ycoef=\"  1.13071  -0.285576 0.905262 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.163.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.805632 0.945825  0.228786 0 0 0\"\n ycoef=\"  1.07788  -0.26981  0.908425 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.164.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.757899 0.954317  0.220822 0 0 0\"\n ycoef=\"  0.918611 -0.259237 0.930146 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.165.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.645357 0.956967  0.18769  0 0 0\"\n ycoef=\"  0.861434 -0.233107 0.92383  0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.166.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.575172 0.954218  0.160956 0 0 0\"\n ycoef=\"  0.876024 -0.216296 0.881098 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.167.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.672134 0.965437  0.167847 0 0 0\"\n ycoef=\"  0.801169 -0.206858 0.884693 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.168.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.64109  0.945975  0.167473 0 0 0\"\n ycoef=\"  0.775758 -0.195263 0.869563 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.169.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.757725 0.975328  0.176162 0 0 0\"\n ycoef=\"  0.467626 -0.19215  0.939955 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.170.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.574634 0.967270  0.15724  0 0 0\"\n ycoef=\"  0.387368 -0.176538 0.931817 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.171.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.345244 0.955863  0.120641 0 0 0\"\n ycoef=\"  0.354572 -0.156131 0.911487 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.172.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.232312 0.957871  0.099471 0 0 0\"\n ycoef=\"  0.34869  -0.129712 0.857674 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.173.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.263422 0.990501  0.083953 0 0 0\"\n ycoef=\"  0.289201 -0.123626 0.855877 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.174.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.092615 0.982441  0.041032 0 0 0\"\n ycoef=\"  0.323915 -0.089680 0.820583 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.175.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.221774 0.951882  0.085049 0 0 0\"\n ycoef=\"  0.726552 -0.131599 0.806829 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.176.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\"  0.086517 0.946791  0.057829 0 0 0\"\n ycoef=\"  0.327542 -0.122283 0.846494 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.177.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\"  0.362542 0.93642   0.035024 0 0 0\"\n ycoef=\"  0.061469 -0.101223 0.879123 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.178.jpg\" />\n"},
{"<Transform dim=\"5\"\n xcoef=\"  0.275379 0.91144   0.082182 0 0 0\"\n ycoef=\" -0.186991 -0.047717 0.876127 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.179.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\"  0.327895 0.889087  0.061304 0 0 0\"\n ycoef=\" -0.162636 -0.011406 0.857243 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.180.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\"  0.555568 0.90288   0.036076 0 0 0\"\n ycoef=\"  0.280108 -0.051341 0.857085 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.181.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\"  0.373342 0.908789  0.052171 0 0 0\"\n ycoef=\"  0.612486 -0.102106 0.864458 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.182.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\"  0.166105 0.912235  0.093130 0 0 0\"\n ycoef=\"  1.03847  -0.172595 0.877695 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.183.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.128142 0.912321  0.121956 0 0 0\"\n ycoef=\"  0.857257 -0.205699 0.917056 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.184.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.299415 0.938929  0.116941 0 0 0\"\n ycoef=\"  0.822884 -0.232976 0.941297 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.185.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.263943 0.975146  0.135819 0 0 0\"\n ycoef=\"  1.16998  -0.222324 0.92108  0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.186.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.144616 0.981653  0.119271 0 0 0\"\n ycoef=\"  1.25191  -0.207466 0.925867 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.187.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.435291 1.01591   0.10602  0 0 0\"\n ycoef=\"  0.840253 -0.149542 0.960153 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.188.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.448744 1.02806   0.085348 0 0 0\"\n ycoef=\"  0.772701 -0.118098 0.959777 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.189.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.372011 1.0242    0.064173 0 0 0\"\n ycoef=\"  0.491639 -0.087172 0.966604 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.190.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.303778 1.01996   0.053958 0 0 0\"\n ycoef=\"  0.484059 -0.077729 0.968171 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.191.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.208073 1.01628   0.034365 0 0 0\"\n ycoef=\"  0.420249 -0.054410 0.970022 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.192.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.055000 1.00581   0.018570 0 0 0\"\n ycoef=\"  0.34236  -0.052159 0.975393 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.193.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\"  0.132279 0.990768  0.002700 0 0 0\"\n ycoef=\"  0.286624 -0.035613 0.979478 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.194.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\"  0.204218 1.00125  -0.026736 0 0 0\"\n ycoef=\"  0.587619 -0.026077 0.919795 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.195.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\"  0.073291 0.992519  0.002385 0 0 0\"\n ycoef=\"  0.506158 -0.010611 0.892593 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.196.jpg\" />\n"},
{"<Transform dim=\"1\"\n xcoef=\"  0        1         0        0 0 0\"\n ycoef=\"  0         0        1        0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.197.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\"  0.27253  0.989669 -0.039327 0 0 0\"\n ycoef=\"  0.348675 -0.021699 0.968097 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.198.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\"  0.325072 0.984014 -0.039448 0 0 0\"\n ycoef=\"  0.416136 -0.023148 0.953598 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.199.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\"  0.330289 0.97511  -0.031893 0 0 0\"\n ycoef=\"  0.404986 -0.032236 0.9648   0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.200.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\"  0.136082 0.995151 -0.032178 0 0 0\"\n ycoef=\"  0.121593 -0.000853 0.978634 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.201.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.09065  1.00328  -0.036164 0 0 0\"\n ycoef=\" -0.097670  0.035650 0.952329 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.202.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.059485 1.0137   -0.038432 0 0 0\"\n ycoef=\" -0.409962  0.074925 0.975278 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.203.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.104324 1.01133  -0.007118 0 0 0\"\n ycoef=\"  0.047449  0.024072 0.945988 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.204.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.192996 1.02209  -0.007660 0 0 0\"\n ycoef=\"  0.024973  0.038426 0.909541 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.205.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\"  0.155243 0.981077 -0.047352 0 0 0\"\n ycoef=\" -0.135798  0.081202 0.897562 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.206.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.159184 1.01125   0.001690 0 0 0\"\n ycoef=\" -0.016627  0.022632 0.932037 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.207.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.27954  1.02071   0.020250 0 0 0\"\n ycoef=\" -0.050841  0.023077 0.93891  0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.208.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.311175 1.00913   0.052492 0 0 0\"\n ycoef=\" -0.157851  0.054966 0.964588 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.209.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.292632 1.00852   0.055014 0 0 0\"\n ycoef=\" -0.247812  0.097687 0.952239 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.210.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.352196 1.0113    0.051685 0 0 0\"\n ycoef=\" -0.320865  0.072630 0.972166 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.211.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.517736 1.01814   0.057519 0 0 0\"\n ycoef=\"  0.088489  0.028498 0.961353 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.212.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.618524 1.03833   0.046930 0 0 0\"\n ycoef=\" -0.095403  0.034972 0.968201 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.213.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.640915 1.0417    0.024168 0 0 0\"\n ycoef=\"  0.124188  0.033604 0.910679 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.214.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.659717 1.01871   0.022580 0 0 0\"\n ycoef=\"  0.419468  0.028311 0.85005  0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.215.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.495591 1.004     0.017757 0 0 0\"\n ycoef=\"  0.533418  0.048065 0.842181 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.216.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.620478 1.01627   0.016763 0 0 0\"\n ycoef=\"  0.122618  0.056544 0.903929 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.536\" brightness=\"-0.525\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.217.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.376811 1.00842   0.023513 0 0 0\"\n ycoef=\" -0.543506  0.044956 0.996859 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.218.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.347052 1.00961   0.026795 0 0 0\"\n ycoef=\" -0.410749  0.038119 0.981788 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.219.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.663176 1.05941   0.041255 0 0 0\"\n ycoef=\" -0.314064  0.005371 0.993747 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.220.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.763934 1.06496   0.058509 0 0 0\"\n ycoef=\" -0.412361  0.010700 0.974029 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.221.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.784827 1.05772   0.069075 0 0 0\"\n ycoef=\"  0.031003 -0.005328 0.969049 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.222.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.64403  1.05804   0.049062 0 0 0\"\n ycoef=\" -0.135916  0.023965 0.962666 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.223.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.596945 1.05616   0.046746 0 0 0\"\n ycoef=\" -0.236727  0.028652 0.96731  0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.224.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.481917 1.06206   0.034466 0 0 0\"\n ycoef=\" -0.219585  0.022361 0.965894 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.225.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.514526 1.07158   0.013272 0 0 0\"\n ycoef=\" -0.177502  0.037844 0.962736 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.226.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.755211 1.05395   0.034348 0 0 0\"\n ycoef=\" -0.063100  0.033946 0.956009 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.227.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -0.696786 1.04958   0.034329 0 0 0\"\n ycoef=\" -0.181886  0.047434 0.972457 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.228.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -0.715407 1.04153   0.060353 0 0 0\"\n ycoef=\"  0.318194 -0.005487 0.946767 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.229.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -0.710064 1.01696   0.085684 0 0 0\"\n ycoef=\"  0.643598 -0.045513 0.879763 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.230.jpg\" />\n"},
{"<Transform dim=\"5\"\n xcoef=\" -1.16379  1.08583   0.075300 0 0 0\"\n ycoef=\"  0.377913 -0.009395 0.897096 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.231.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -1.25433  1.07844   0.068434 0 0 0\"\n ycoef=\"  0.360629 -0.039616 0.991683 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.232.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -1.24871  1.06816   0.075640 0 0 0\"\n ycoef=\"  0.396387 -0.037202 0.997071 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.233.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -1.24117  1.05192   0.085490 0 0 0\"\n ycoef=\"  0.509605 -0.064906 0.961517 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.234.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -1.45834  1.0646    0.106936 0 0 0\"\n ycoef=\"  0.326807 -0.055998 0.975617 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.235.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -1.4302   1.06469   0.127169 0 0 0\"\n ycoef=\"  0.532403 -0.070365 0.96537  0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.236.jpg\" />\n"},
{"<Transform dim=\"5\"\n xcoef=\" -1.52792  1.07429   0.143841 0 0 0\"\n ycoef=\"  0.649129 -0.103984 0.971546 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.237.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -1.72439  1.07799   0.165071 0 0 0\"\n ycoef=\"  0.788722 -0.104344 0.951925 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.238.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -1.56433  1.06795   0.154614 0 0 0\"\n ycoef=\"  0.896773 -0.104341 0.952253 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.239.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -1.70677  1.08287   0.165337 0 0 0\"\n ycoef=\"  0.873286 -0.105163 0.981151 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.240.jpg\" />\n"},
{"<Transform dim=\"5\"\n xcoef=\" -1.93115  1.09273   0.166158 0 0 0\"\n ycoef=\"  1.03756  -0.109699 0.960169 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.241.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -2.19938  1.09911   0.196347 0 0 0\"\n ycoef=\"  1.32725  -0.143896 0.954941 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.242.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -2.11207  1.08446   0.191679 0 0 0\"\n ycoef=\"  1.29066  -0.123522 0.94935  0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.243.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -2.26635  1.1152    0.143889 0 0 0\"\n ycoef=\"  0.735342 -0.072433 0.966669 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.244.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -2.13363  1.11236   0.127081 0 0 0\"\n ycoef=\"  0.700579 -0.058386 0.957397 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.245.jpg\" />\n"},
{"<Transform dim=\"3\"\n xcoef=\" -2.3825   1.12435   0.134263 0 0 0\"\n ycoef=\"  0.92583  -0.065302 0.978261 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.246.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -2.54387  1.12134   0.175902 0 0 0\"\n ycoef=\"  1.28503  -0.113561 0.964423 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.247.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -1.99091  1.06984   0.192905 0 0 0\"\n ycoef=\"  1.13652  -0.131673 0.952262 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.248.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -2.45699  1.09211   0.209045 0 0 0\"\n ycoef=\"  0.844755 -0.105703 0.985923 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.249.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -2.59333  1.10207   0.207521 0 0 0\"\n ycoef=\"  1.63528  -0.180552 0.923165 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.102\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.250.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -2.61576  1.11149   0.171091 0 0 0\"\n ycoef=\"  1.09958  -0.12532  0.929613 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.251.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -2.39431  1.12068   0.091140 0 0 0\"\n ycoef=\"  0.729446 -0.029549 0.876252 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.252.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -1.54375  0.995467  0.100091 0 0 0\"\n ycoef=\"  1.17229  -0.044727 0.834971 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.253.jpg\" />\n"},
{"<Transform dim=\"5\"\n xcoef=\" -1.71362  0.996894  0.128574 0 0 0\"\n ycoef=\"  1.03975  -0.064026 0.887317 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.254.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -1.38566  0.936509  0.136393 0 0 0\"\n ycoef=\"  1.12851  -0.110456 0.927976 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.255.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -2.08588  1.01388   0.144965 0 0 0\"\n ycoef=\"  1.39503  -0.152661 0.97403  0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.256.jpg\" />\n"},
{"<Transform dim=\"5\"\n xcoef=\" -2.01715  1.02005   0.14098  0 0 0\"\n ycoef=\"  1.12674  -0.105629 0.965139 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.257.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -2.17679  1.03498   0.107633 0 0 0\"\n ycoef=\"  0.803071 -0.065977 0.965064 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.258.jpg\" />\n"},
{"<Transform dim=\"6\"\n xcoef=\" -2.22523  1.03516   0.121136 0 0 0\"\n ycoef=\"  0.974713 -0.092563 0.961757 0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.259.jpg\" />\n"},
{"<Transform dim=\"4\"\n xcoef=\" -2.24839  1.08261   0.028368 0 0 0\"\n ycoef=\"  0.922613 -0.061297 0.93662  0 0 0\">\n<Image mag=\"0.0022159\" contrast=\"1.103\" brightness=\"-0.345\" red=\"true\" green=\"true\" blue=\"true\"\n src=\"R34CA1-B_S12.260.jpg\" />\n"}};

void Container::openContourFile (std::ofstream & contour_handle,double z)
{
  int index = z/50.0;
  char str[256];
  sprintf(str,"%sorthogonal_vertices.%d",Controls::instance().get_output_data_dir().c_str(),
                                         index);
  contour_handle.open(str);
  assert (contour_handle.is_open());
  //contour_handle << "<?xml version=\"1.0\"?>\n"
  //  << "<!DOCTYPE Section SYSTEM \"section.dtd\">\n\n"
  //  << "<Section index=\"" << index << "\" thickness=\"0.05\" alignLocked=\"true\">\n"
  //  << "<Transform dim=\"3\"\n"
  //  << " xcoef=\" -0.567106 0.908516 0.279141 0 0 0\"\n"
  //  << " ycoef=\" 1.87787 -0.329688 0.8712 0 0 0\">\n"
  //  << "<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n"
  //  << " src=\"R34CA1-B_S12." << index+100 << ".jpg\" />\n"
  //  << "<Contour name=\"domain1\" hidden=\"false\" closed=\"true\" simplified=\"false\" border=\"1 0 1\" fill=\"1 0 1\" mode=\"11\"\n"
  //  << " points=\"0 0,\n"
  //  << "	4096 0,\n"
  //  << "	4096 4096,\n"
  //  << "	0 4096,\n"
  //  << "	\"/>\n"
  //  << "</Transform>\n\n";
  contour_handle << "<?xml version=\"1.0\"?>\n"
    << "<!DOCTYPE Section SYSTEM \"section.dtd\">\n\n"
    << "<Section index=\"" << index << "\" thickness=\"0.05\" alignLocked=\"true\">\n"
    //<< "<Transform dim=\"3\"\n"
    //<< " xcoef=\" -0.567106 0.908516 0.279141 0 0 0\"\n"
    //<< " ycoef=\" 1.87787 -0.329688 0.8712 0 0 0\">\n"
    //<< "<Image mag=\"0.0022159\" contrast=\"0.882\" brightness=\"-0.075\" red=\"true\" green=\"true\" blue=\"true\"\n"
    //<< " src=\"R34CA1-B_S12." << index+100 << ".jpg\" />\n"
    << img_data[index-60]
    << "<Contour name=\"domain1\" hidden=\"false\" closed=\"true\" simplified=\"false\" border=\"1 0 1\" fill=\"1 0 1\" mode=\"11\"\n"
    << " points=\"0 0,\n"
    << "	4096 0,\n"
    << "	4096 4096,\n"
    << "	0 4096,\n"
    << "	\"/>\n"
    << "</Transform>\n\n";
}

void Container::binVertexNormalAngle (double * bins,int * count,const double & cos_angle)
{
  int a = 12;
  unsigned int start=0;
  unsigned int end=a-1;
  while (end-start>1)
  {
    unsigned int mid=(start+end)/2;
    if (cos_angle <= bins[start] && cos_angle >= bins[mid]) { end=mid;}
    else { start=mid; }
  }
  if (cos_angle <= bins[start] && cos_angle >= bins[end]) { count[start]++;}
  else { count[end]++; }
}

void Container::findOrthogonalRegions (void)
{
  double pi = 3.14159;
  double threshold = cos(pi/2.0*Controls::instance().get_region_orthogonality_threshold());
  // instantiate vector of vertex pointers (or vertex vector index)
  std::vector<const Vertex*> orthogonal_vertices; 
  // initialize vertex normal angle histogram
  double bins[13];
  bins[0]  = 1.0;
  bins[1]  = 1.0;
  for (int i=2;i<11;i++) bins[i]  = cos(pi/2.0*(static_cast<double>(i-1)/10.0));
  bins[11] = 0.0;
  bins[12] = 0.0;
  int   count[12];
  for (int i=0;i<12;i++) count[i]=0;

  cout << "# ecw num_neighbors\n";
  // for each object in container
  for (o_it i=o.begin();i!=o.end();++i)
  {
    // for each vertex in object
    for (v_it j=i->v.begin();j!=i->v.end();++j)
    {
      // get vertex normal
      // (previously calculated and stored during initialization)
      vector3 n = j->getNormal(); 
      // compute cosine of angle between vertex normal and section plane
      double a = n.p[0]*n.p[0];
      double b = n.p[1]*n.p[1];
      double c = n.p[2]*n.p[2];
      double d = a+b;
      double e = d+c;
      double cos_angle = d/sqrt(d*e);
      // bin vertex normal angle
      binVertexNormalAngle (&bins[0],&count[0],cos_angle);
      // if normal satisfies orthogonality threshold then keep
      if (cos_angle>= threshold)
      {
        orthogonal_vertices.push_back(&(*j));
        // if vertex has a closest face
        if (j->getClosestFace()!=NULL)
        {
          // get closest point
          vector3 closest_point;
          double sqd;
          findClosestPtInFaceToLocation(*j->getPos(),j->getClosestFace(),closest_point,sqd);
          double dd = sqrt(sqd);
          cout << dd << " " << calculateVertexNeighborCount(&(*j),j->getPos(),j->getNormal())
                << endl;
        }
      }
      //orthogonal_vertices.push_back(&(*j));
    }
  }
  cout << "            total num vertices = " << getVertexCount() << endl;
  cout << "       num orthogonal vertices = " << orthogonal_vertices.size() << endl;

  // 1) output to contour files
  writeOrthVerts2Recon3D(orthogonal_vertices);
  // 2) output to visualize in dreamm
  writeOrthVerts2Dreamm(orthogonal_vertices);
  // write histogram to stdout
  cout << "\n\n"
        << "Histogram of vertex normal angle relative to plane of section (degrees):\n\n";
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         acos(bins[0])*180.0/pi, acos(bins[1])*180.0/pi, count[0], acos(bins[6])*180.0/pi, acos(bins[7])*180.0/pi, count[6]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         acos(bins[1])*180.0/pi, acos(bins[2])*180.0/pi, count[1], acos(bins[7])*180.0/pi, acos(bins[8])*180.0/pi, count[7]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         acos(bins[2])*180.0/pi, acos(bins[3])*180.0/pi, count[2], acos(bins[8])*180.0/pi, acos(bins[9])*180.0/pi, count[8]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         acos(bins[3])*180.0/pi, acos(bins[4])*180.0/pi, count[3], acos(bins[9])*180.0/pi, acos(bins[10])*180.0/pi, count[9]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         acos(bins[4])*180.0/pi, acos(bins[5])*180.0/pi, count[4], acos(bins[10])*180.0/pi, acos(bins[11])*180.0/pi, count[10]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         acos(bins[5])*180.0/pi, acos(bins[6])*180.0/pi, count[5], acos(bins[11])*180.0/pi, acos(bins[12])*180.0/pi, count[11]);
  cout << "\n\n";
  //printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
  //       bins[6], bins[7], count[6], bins[13], bins[15], count[14]);
  //printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
  //       bins[7], bins[8], count[7], bins[15], bins[16], count[15]);
}

void Container::findVertexNeighbors (const int & group)
{
  cout << "Iteration " << group << ": ";
  cout << "Counting vertex neighbor objects...............";
  cout.flush();
  // for each object
  for (o_it i=o.begin();i!=o.end();++i)
  {
    // for each vertex in object
    for (v_it j=i->v.begin();j!=i->v.end();++j)
    {
      j->setNeighborCount(calculateVertexNeighborCount(&(*j),j->getPos(),j->getNormal()));
    }
  }
  cout << "complete.\n";
  cout.flush();
}

double Container::calculateCleftArea (void) const
{
  double area = 0.0;
  // for each cleft vertex
  for (vp_cit i = cleft.begin();i!=cleft.end();i++)
  {
    area += (*i)->getArea();
  }
  return area;
}

void Container::findPerisynaptic (void)
{
  // since clefts do NOT cover entire synapse junction,
  // then perisynaptic regions should be twice (not equal)
  // as wide as clefts to adequately sample ECS around synapses.
  //
  // in doc/perisynaptic.m determine that a perisynaptic region twice
  // as wide as cleft will have area 8 times larger than cleft area.

  cout << "Finding perisynaptic vertices..................";
  cout.flush();
  v_set newest,keepers;
  v_set next_to_check(cleft.begin(),cleft.end());
  double cleft_area = calculateCleftArea();
  double peri_area  = 0.0;

  while (true)
  {
    for (vs_it i=next_to_check.begin();i!=next_to_check.end();i++)
    {
      vec_vp adjacent_vertices;
      (*i)->getAdjVertices(adjacent_vertices);
      for (vp_it j=adjacent_vertices.begin();j!=adjacent_vertices.end();j++)
      {
        // if vertex is cleft or keeper, then continue
        if (vertexIsCleft(*j) || (keepers.count(*j)>0)) continue;
        // vertex is perisynaptic
        newest.insert(*j);
        peri_area += (*j)->getArea();
      }
    }

    keepers.insert(newest.begin(),newest.end());
    if (peri_area >= 8.0*cleft_area) break;
    //next_to_check.clear();
    //next_to_check;.insert(newest.begin(),newest.end());
    next_to_check = newest;
    newest.clear();
  }
  peri.assign(keepers.begin(),keepers.end());
  cout << "complete.\n";
  cout.flush();
}
