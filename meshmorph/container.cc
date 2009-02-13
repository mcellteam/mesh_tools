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
:fo(NULL),files(),frozen(),world(6,0),octree(),o()
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
//  F.close();
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
  int k=1;
  double goal = 0.2;
  printf("0%%..");
  fflush(stdout);
  // for each object, create edges
  double a = 1.0/o.size();
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
    if (progress>goal)
    {
      printf("%d%%..",static_cast<int>(goal*100));
      fflush(stdout);
      goal+=0.2;
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
  int index=1;
  double goal = 0.2;
  printf("0%%..");
  fflush(stdout);
  // for each object in container
  double a = 1.0/o.size();
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
        if (findClosestPtToVertex(&(*j),dummy_p,dummy_sqd,cl)==true)
        {
          j->setFace(cl);
        }
        //assert(dummy_sqd<1E10);
        assert(dummy_sqd==dummy_sqd);
      }
    }
    // track progress
    double progress = static_cast<double>(index++)*a;
    if (progress>goal)
    {
      printf("%d%%..",static_cast<int>(goal*100));
      fflush(stdout);
      goal+=0.2;
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
  double radius_step = (sqrt(cs.get_search_radius_sq())-cs.get_min_search_cone_radius())
        /static_cast<double>(cs.get_number_radius_steps());
  // DEBUG
  //bool flag = false;
  //if (distinguishable(pt.p[0],4431.686,1E-7)==false &&
  //    distinguishable(pt.p[1],5570.356,1E-7)==false &&
  //    distinguishable(pt.p[2],6260.725,1E-7)==false)
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
  //    exit(1);
  //  }
  //  else
  //  {
  //    cout << "\n closest_point not found, d = ???\n";
  //  }
  //  exit(1);
  //}
  // DEBUG
  // update statistics of count of face and boxes used in search
  Log::instance().updateClosestPtStats(num_faces,num_leaves,total_faces_checked);
  return found;
}

/** Find the closest point to a vertex.
 * \param[in] v Vertex for which closest point is searched.
 * \param[out] p Closest point position, if found.
 * \param[out] sqd Squared distance between vertex
 * and closest point, if found.
 * \param[out] ncl Parent face of closest point, if found.
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
  double radius_step = (sqrt(cs.get_search_radius_sq())-cs.get_min_search_cone_radius())
        /static_cast<double>(cs.get_number_radius_steps());
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
  //if (distinguishable(pt.p[0],4431.686,1E-7)==false &&
  //    distinguishable(pt.p[1],5570.356,1E-7)==false &&
  //    distinguishable(pt.p[2],6260.725,1E-7)==false)
  //{
  //  flag = true;
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
    //  (*j)->print(cout);
    //}
    // DEBUG
    (*j)->clearFlag();
    // reject face of tile
    if (*j==f) continue;
    // reject faces on wrong side of orthogonal plane to vertex normal
    if (faceLiesOppositeToNormal(pt,*j,n)==true) continue;
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


