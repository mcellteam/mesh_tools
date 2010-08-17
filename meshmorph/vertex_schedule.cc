// Author: Justin Kinney
// Date: Sep 2008

#include "vertex_schedule.h"

#include <cassert>
#include <iostream>

#include "log.h"
#include "nice.h"
#include "vertex.h"
#include "controls.h"
#include "container.h"
#include "refractory.h"
#include "virtual_disp.h"
#include "gain_schedule.h"
#include "intersecting_faces.h"

using std::cout;
using std::endl;

Vertex_Schedule * Vertex_Schedule::only_one = NULL;

Vertex_Schedule & Vertex_Schedule::instance(void)
{
  // Not thread-safe.
  // -- lock mutex
  if (only_one == NULL)
    only_one = new Vertex_Schedule();
  // -- unlock mutex
  return *only_one;
}

Vertex_Schedule::Vertex_Schedule (void)
  :cv(NULL),seed(NULL),vset(),pH(0,0,0),count(0),count_ref(0)
{
}

/** Process vertex move sequence data.
 * \param[in] filename File containing a list of vertices to be moved in order.
 */

void Vertex_Schedule::readVertexSequence (const char *filename)
{
  cout << "Read sequence of vertex moves..................";
  cout.flush();
  // load frozen map: object*->Vertex_index (int)
  oi_it front;
  // set of object names from frozen file not found as Object name
  s_set not_found;
  // read frozen vertex file 
  std::vector<Complex> myvec = Container::instance().loadVector(filename,not_found);

  // for each element in multimap
  // DEBUG
//  int mofo = 0;
  // DEBUG
  for (std::vector<Complex>::iterator i=myvec.begin();i!=myvec.end();++i)
  {
    Object *oo=(*i).o;
    int t = (*i).index;
  // DEBUG
//    cout << "Vertex_Schedule::readVertexSequence: mofo = " << mofo++ << endl;
//    cout << "Vertex_Schedule::readVertexSequence: object = " << oo->getName()
//          << ", index = " << t << endl;
  // DEBUG
    int a = t-1;
    int b = oo->v[a].getIndex();
    do{
      if (t==b){vset.push_back(&oo->v[a]);}
      else if (t<b){a++;}
      else if (t>b){a--;}
    } while (t!=b);
  }
  cout << "complete.\n";

  // print objects not found
//  if (not_found.empty()==false)
//  {
//    cout << "\nVertex_Schedule::ReadVertexSequence: Warning.\n"
//          << "No matching Object* found in container "
//          << "for following objects:\n";
//    for (ss_it j = not_found.begin();j!=not_found.end();++j)
//    {
//      cout << *j << endl;
//    }
//    cout << endl;
//  }
  cout.flush();
}

/** Calculate and record move location of input vertex.
 * \param[in] v Vertex attempting a move.
 */

void Vertex_Schedule::calculateMoveLocation (Vertex * const v)
{
  cv=v;
  // compute new holding position coordinates (x,y,z)
  //cout << "222";cout.flush();
  pH = cv->getNewPos(Gain_Schedule::instance().getGain());
  // impose max displacement policy
  Refractory::instance().enforceMaxdisplacement(cv,pH);
  // displacement along x,y,z
  double a=pH.p[0]-*cv->getCoord(0);
  double b=pH.p[1]-*cv->getCoord(1);
  double d=pH.p[2]-*cv->getCoord(2);
  // compute and store squared displacement 
  State::instance().setVD2(a*a+b*b+d*d);
}

/** Identify and record a collection of vertices to move next.
 */

void Vertex_Schedule::collectVerticesToMoveNext (const int & group)
{
  Controls & cs(Controls::instance());
  if (cs.get_vertex_sequence_file().empty()==true)
  {
    seed=NULL;
    // try to pick nonnice vertex as set seed
    // for each element of nice hashmap
    // Note sequence of elements in hash_map is not necessarily repeatable.
    for (vhm_cit j=Nice::instance().beginNice();j!=Nice::instance().endNice();++j)
    {
      Vertex *vv=(*j).first;
      // if vertex is not refracted and not refracted
      // and not frozen and has a closest face
//      if ( Refractory::instance().isRefracted(vv)==false &&
//           vv->getClosestFace()!=NULL &&
//           (Container::instance().vertexIsFrozen(vv)==false))
      if (Refractory::instance().vertexIsMoveCandidate(vv)==true)
      {
        seed = vv;break;
      }
    }

    // else try to pick intersected face vertex as set seed
    if (seed==NULL)
    {
      Intersecting_Faces & i_f(Intersecting_Faces::instance());
      // for each element of intf hashmap
      // Note sequence of elements in hash_map is not necessarily repeatable.
      for (htff_cit j=i_f.begin();j!=i_f.end();++j)
      {
        Face *ff=(*j).first;
        // for each vertex of face
        for (int k=0;k<3;++k)
        {
          Vertex *vv=ff->getVertex(k);
          // if vertex is not refracted and not refracted
          // and not frozen and has a closest face
//          if (Refractory::instance().isRefracted(vv)==false &&
//              vv->getClosestFace()!=NULL && 
//              (Container::instance().vertexIsFrozen(vv)==false))
          if (Refractory::instance().vertexIsMoveCandidate(vv)==true)
          {
            seed = vv;break;
          }
        }
        if (seed!=NULL){break;}
      }
    }
//    // if block wrap
//    if (cs.get_recon_block_wrap()!="")
//    {
//      Virtual_Disp & vd(Virtual_Disp::instance());
//      seed = (*(vd.getSeed())).second;
////             (Container::instance().vertexIsFrozen(seed)==true) ||
//      vset.clear();
//      while (vd.getSeed()!=vd.rendVirtDispMap())
//      {
//        if (seed->getObjectName()==cs.get_recon_block_wrap())
//        {
//          if (Refractory::instance().vertexIsMoveCandidate(seed)==true)
//          {
//            vset.push_back(seed);
//          }
//        }
//        vd.advanceSeedToNextLargestVert();
//        seed = (*(vd.getSeed())).second;
//      }
//      assert(vset.empty()==false);
//      cout << "Vertex_Schedule:: Aquired " << vset.size()
//           << " vertices from " << cs.get_recon_block_wrap() << "object.\n";
//      cout.flush();
//      return;
//    }

    // if block wrap
    if (cs.get_recon_block_wrap()!="")
    {
      Container & c(Container::instance());
      vset.clear();
      for (o_it i=c.o.begin();i!=c.o.end();++i)
      {
        // if object is wrapper
        if (i->getName()==cs.get_recon_block_wrap())
        {
          // for each vertex in object
          for (v_it j=i->v.begin();j!=i->v.end();++j)
          {
            // if vertex is move candidate
            if (Refractory::instance().vertexIsMoveCandidate(&(*j))==true)
            {
              vset.push_back(&(*j));
            }
            //else
            //{
            //  cout << "wrapper vertex is NOT move candidate: ";
            //  cout.flush();
            //  if (j->getClosestFace()==NULL)
            //  {
            //    cout << "No closest face.\n";
            //    cout.flush();
            //  }
            //  j->print(cout);
            //}
          }
        }
      }
      assert(vset.empty()==false);
      cout << "Vertex_Schedule:: Aquired " << vset.size()
           << " vertices from " << cs.get_recon_block_wrap() << "object.\n";
      cout.flush();
      return;
    }

    // else pick vertex from sorted virtual displacement list
    Virtual_Disp & vd(Virtual_Disp::instance());
    if (seed==NULL)
    {
      seed = (*(vd.getSeed())).second;
      // starting at current location of tvi in topN
      // find vertex with largest virtual displacement that is not refracted
      // and not frozen and has a closest face
//      while ((Refractory::instance().isRefracted(seed)==true)   ||
//             (Container::instance().vertexIsFrozen(seed)==true) ||
//             (seed->getClosestFace()==NULL))
      while (Refractory::instance().vertexIsMoveCandidate(seed)==false)
      {
        vd.advanceSeedToNextLargestVert();
        seed = (*(vd.getSeed())).second;
        if (vd.getSeed()==vd.rendVirtDispMap())
        {vset.clear();return;}
      }
    }
    if (cs.get_freeze_sheets()==false)
    {
      // create vertex update region
      vector3 origin,end;
      origin.p[0] = *seed->getCoord(0)-cs.get_seed_region_size();
      origin.p[1] = *seed->getCoord(1)-cs.get_seed_region_size();
      origin.p[2] = *seed->getCoord(2)-cs.get_seed_region_size();
      end.p[0]    = *seed->getCoord(0)+cs.get_seed_region_size();
      end.p[1]    = *seed->getCoord(1)+cs.get_seed_region_size();
      end.p[2]    = *seed->getCoord(2)+cs.get_seed_region_size();
      // make a visitor
      Octree_Visitor_Face visitor(Vector3r(origin.p[0],origin.p[1],origin.p[2]),
                                  Vector3r(end.p[0],end.p[1],end.p[2]));
      // execute visitor
      Container::instance().octree->visit( visitor );
      //visitor.sort();
      vec_vp bin;
      bin.reserve(cs.get_vector_reserve());
      // for each face in box
      for (fp_cit j=visitor.mybegin();j!=visitor.myend();++j)
      {
        (*j)->clearFlag();
        // add face vertices to vector
        bin.push_back((*j)->getVertex(0));
        bin.push_back((*j)->getVertex(1));
        bin.push_back((*j)->getVertex(2));
      }
      // sort and keep unique Vertex*
      sort(bin.begin(),bin.end(),my_ltv());
      bin.assign(bin.begin(),unique(bin.begin(),bin.end()));
      // build vset
      vset.clear();	
      vset.push_back(seed);
      for (vp_it i=bin.begin();i!=bin.end();++i)
      {
        vset.push_back(*i);
      }
    }
    else
    {
      //seed->getAdjVertices(vset);
      vset.clear();	
      vset.reserve(100000);
      vset.push_back(seed);
      vd.advanceSeedToNextLargestVert();
      while (vd.getSeed()!=vd.rendVirtDispMap())
      {
        seed = (*(vd.getSeed())).second;
        if (Refractory::instance().vertexIsMoveCandidate(seed)==true)
        {
          vset.push_back(seed);
        }
        vd.advanceSeedToNextLargestVert();
      }
    }
  }
  // error checking
  if (vset.empty()==true)
  {
    cout << "\nVertex_Schedule::collectVerticesToMoveNext: "
          << "No vertices were found to move.\n";
    assert(vset.empty()==false);
  }
  // to detect if any vertices from set were moved
  count_ref = count;
  // record seed virtual displacement in for steady-state analysis
  Virtual_Disp::instance().addVdToSeedVd(count+static_cast<int>(cs.get_group_size()*(group-1.0)));
  // collect statistics on number of vertices per set
  Log::instance().updateVertexSchedulingStats(vset.size());
}
