// Author: Justin Kinney
// Date: Sep 2008

#include "virtual_disp.h"

#include <iostream>
#include <cmath>
#include <cassert>

#include "vertex.h"
#include "controls.h"
#include "container.h"
#include "gain_schedule.h"

using std::cout;
using std::endl;

Virtual_Disp * Virtual_Disp::only_one = NULL;

Virtual_Disp & Virtual_Disp::instance(void)
{
  // Not thread-safe.
  // -- lock mutex
  if (only_one == NULL)
    only_one = new Virtual_Disp();
  // -- unlock mutex
  return *only_one;
}

Virtual_Disp::Virtual_Disp (void)
  :vd2_to_v(),v_to_vd2(),seed(),seed_vd(),seed_ad()
{
}

/** Build v_to_vd2, the complement map of vd2_to_v.
 */

void Virtual_Disp::buildVirtDispMapComplement (void)
{
  v_to_vd2.clear();
  // find vertex in vd2_to_v
  for (tv_it i=vd2_to_v.begin();i!=vd2_to_v.end();++i)
  {
    assert((*i).second!=NULL);
    // add to v_to_vd2
    v_to_vd2[(*i).second]=(*i).first;
  }
}

/** Find and return element in vd2_to_v that matches input
 * vertex and squared virtual displacement.
 * \param[in] v Vertex of interest.
 * \param[in] sqd_vd Squared virtual displacement of vertex.
 * \param[out] match Iterator to matching vd2_to_v element, if found;
 * \return True if match found; false otherwise.
 */

bool Virtual_Disp::findVirtDispToVertAssoc (Vertex * const v,
                                            float const & sqd_vd,
                                            tv_it & match)
{
  // vd2_to_v has descending order
  // kby key (double, squared virtual displacement)
  //
  // check for elements in vd2_to_v
  // with matching squared virtual displacement
  std::pair<tv_it,tv_it> p = vd2_to_v.equal_range(sqd_vd);
  // if no match found then return
  if (p.first==p.second) return false;
  // for each match 
  for (tv_it i=p.first;i!=p.second;++i)
  { 
    // if the displacment match is associate with input vertex
    if (v==(*i).second)
    {
      // record iterator
      match = i;
      return true;
    }
  }
  // print error message
  // no match is associated with v
  cout.precision(9);
  cout << "\n\nVirtual_Disp::findVirtDispToVertAssoc: "
        << "Error. No element in vd2_to_v matches both input "
        << "squared virtual displacement and associated vertex.\n"
        << "input squared virtual displacement = " << sqd_vd << endl
        << "input vertex:\n";
  v->print(std::cout);
  cout << "\nVirtual_Disp::findVirtDispToVertAssoc: "
        << "print vd2_to_v elements that match input displacement:\n";
  // for each entry in range
  for (tv_it i=p.first;i!=p.second;++i)
  {
    assert((*i).second!=NULL);
    if ((*i).second==NULL)
    {
      cout << "vd2_to_v contains ["
            << (*i).first << " , NULL]\n";
    }
    else
    {
      cout << "vd2_to_v contains ["
            << (*i).first << " , "
            << (*i).second->getIndex() << "]\n";
    }
  }
  bool match_not_found=true;
  cout << "\nVirtual_Disp::findVirtDispToVertAssoc: "
        << "print vd2_to_v elements that match input vertex:\n";
  for (tv_it i=vd2_to_v.begin();i!=vd2_to_v.end();++i)
  {
    if ((*i).second==v)
    {
      cout << "current vertex found in vd2_to_v: "
            << " vertex " << (*i).second->getIndex()
            << " " << (*i).first << endl << endl;
      match_not_found = false;
    }
  }
  if (match_not_found==true)
  {
    cout << "Input vertex was NOT found anywhere in vd2_to_v.\n";
  }
  cout << endl;
  assert(true==false);
  exit(1);
}

/** Build vd2_to_v (which maps squared virtual dispacement to vertex pointer).
 */

void Virtual_Disp::buildVirtDispMap (const int & group)
{
  cout << "\n\nIteration " << group << ": ";
  cout << "Computing vertex virtual displacements.........";
  cout.flush();
  vd2_to_v.clear();
  Container & c(Container::instance());
  // DEBUG
//  double maxd = 0.0;
  // DEBUG
  // for each object in model
  for (o_it j=c.o.begin();j!=c.o.end();++j)
  {
    // for each vertex in object
    for (v_it i=j->v.begin();i!=j->v.end();++i)
    {
      // if vertex has a closest point
      if (i->getClosestFace()!=NULL)
      {
        // compute squared virtual displacement
        // DEBUG
//        if (Gain_Schedule::instance().getRefGain()<1e-5)
//        {
//          cout.precision(12);
//          cout << "Virtual_Disp::buildVirtDispMap: gain = "
//                << Gain_Schedule::instance().getRefGain() << endl;
//          exit(1);
//        }
        // DEBUG
        //double sq_vd = i->getSqVirtualDisp(Gain_Schedule::instance().getRefGain());
        double sq_vd = i->getSqVirtualDisp(Gain_Schedule::instance().getGain());
        // DEBUG
//        if (sq_vd>maxd) maxd = sq_vd;
//        cout << "max = " << sqrt(maxd) << endl;
//        if (sq_vd>(12000*12000))
//        {
//          cout << "\nVirtual_Disp::buildVirtDispMap: "
//                << "huge vd = " << sqrt(sq_vd) << endl; 
//          //exit(1);
//        }
        // DEBUG
        if (i->getClosestFace()==NULL)
        {
          cout << "\nVirtual_Disp::buildVirtDispMap: "
                << "Error. Vertex had a closest face"
                << "and now it doesn't. What changed?\n";
          assert(i->getClosestFace()!=NULL);
          exit(1);
        }
        // Record association between squared virtual displacement
        // (as key) and input vertex.
        vd2_to_v.insert(std::make_pair(sq_vd,&(*i)));

      }
    }
  }
  cout << "complete.\n";
  cout.flush();
}

/** Calculate and record the virtual displacement and input vertex association.
 * \param[in] v Vertex of interest.
 * \param[in] vertex_has_closest_pt If true then add vertex to maps;
 * \param[in] gain Proportionality constant between force and position.
 * otherwise remove associations referencing input vertex.
 */

void Virtual_Disp::updateVirtualDisp (Vertex * const v,
                                      bool const vertex_has_closest_pt,
                                      float const & gain)
{
  // if closest point was found
  if (vertex_has_closest_pt==true)
  {
    // update sets with vertex squared virtual displacement
    setVirtualDisp(v,v->getSqVirtualDisp(gain));
  }
  else
  {
    // vertex has no closest point
    removeVertFromAllMaps(v);
  }
}

/** Record association between input vertex and squared virtual displacement,
 * removing any previous associations.
 * \param[in] v Vertex of interest.
 * \param[in] sqd_vd Squared virtual displacement of vertex.
 */

void Virtual_Disp::setVirtualDisp (Vertex * const v,
                               float const & sqd_vd)
{
  // use v_to_vd2 to update vd2_to_v
  // if vertex* found in v_to_vd2 remove first
  if (v_to_vd2.find(v)!=v_to_vd2.end())
  {
    tv_it t;
    // if find table entry with old se
    if (findVirtDispToVertAssoc(v,v_to_vd2[v],t))
    {
      // Remove association between squared virtual displacement
      // (as key) and input vertex.
      vd2_to_v.erase(t);
    }
    // Record association between squared virtual displacement
    // (as key) and input vertex.
    vd2_to_v.insert(std::make_pair(sqd_vd,v));
  }
  // else just add to vd2_to_v
  else
  {
    // Record association between squared virtual displacement
    // (as key) and input vertex.
    vd2_to_v.insert(std::make_pair(sqd_vd,v));
  }
  // update v_to_vd2
  v_to_vd2[v]=sqd_vd;
}

/** Find and remove associations between input squared virtual displacement
 * and input vertex in vd2_to_v map.
 * \param[in] v Vertex of interest.
 * \param[in] sqd_vd Squared virtual displacement associated with vertex.
 */

void Virtual_Disp::removeVirtDispToVertAssoc (Vertex * const v,
                                              float const & sqd_vd)
{
  tv_it t;
  // find assoc in vd2_to_v that matches virtual disp and vertex
  if (findVirtDispToVertAssoc(v,sqd_vd,t))
  {
    // remove entry
    vd2_to_v.erase(t);
  }
}

/** Find vertex in map, record rank in ordering by squared virtual displacment,
 * and return iterator to element.
 * \param[in] v Vertex of interest.
 * \param[out] match Iterator to first element that matches input vertex.
 * \param[out] rank Numbering of match in sort by displacement (larget to smallest).
 * \return True if match found; false otherwise.
 */

bool Virtual_Disp::getVertAndRank (Vertex * const v,tv_it & match,int & rank)
{
  bool match_found=false;
  rank=0;
  // for each pair in vd2_to_v
  for (tv_it i=vd2_to_v.begin();i!=vd2_to_v.end();++i)
  {
    // if this vertex matches input vertex
    if ((*i).second==v)
    {
      match = i;
      match_found = true;
      break;
    }
    rank++;
  }
  return match_found;
}

/** Perform data integrity check on vertex to virtual displacement map.
 */

void Virtual_Disp::validateVirtDispMapComplement (void)
{
  tv_it j;
  // for each pair in v_to_vd2
  for (td_it i=v_to_vd2.begin();i!=v_to_vd2.end();++i)
  {
    assert((*i).first!=NULL);
    // if no match is found in vd2_to_v
    bool match_found = findVirtDispToVertAssoc((*i).first,(*i).second,j);
    if (match_found==false)
    {
      cout << "\n\nVirtual_Disp::validateVirtDispMapComplement: "
            << "Error! v_to_vd2 contains element (vertex->"
            << (*i).first->getIndex() << ", vd=" << (*i).second
            << ") with no match found in vd2_to_v.\n";
      assert(match_found==true);
      exit(1);
    }
  }
}

/** Perform data integrity check on virtual displacement to vertex map
 * by checking for duplicate vertices.
 */

void Virtual_Disp::validateVirtDispMap (void) const
{
  ///// check for duplicate vertex* /////
  // instantiate set of vertex*
  v_set uniq_set;
  // for each pair in vd2_to_v
  for (tv_cit i=vd2_to_v.begin();i!=vd2_to_v.end();++i)
  {
    // load v_set
    uniq_set.insert((*i).second);
  }
  // if the two sets are not the same size
  // then vd2_to_v likely contains duplicate vertex* entries
  if (uniq_set.size()!=vd2_to_v.size())
  {
    cout << "\n\nVirtual_Disp::validateVirtDispMap: "
          << "Error! vd2_to_v likely contains duplicate vertex* entries.\n";
    assert(uniq_set.size()==vd2_to_v.size());
    exit(1);
  }
}

/** Perform data integrity check on virtual displacement to vertex map
 * by checking that displacements are not NaN, that the vertex is not NULL,
 * that the displacement is less than threshold value,
 * that closest face to vertex is not NULL, and comparing
 * a calculated virtual displacement to stored displacement.
 */

void Virtual_Disp::validateVirtDispMap2 (void)
{
  // for each pair in vd2_to_v
  for (tv_it i=vd2_to_v.begin();i!=vd2_to_v.end();++i)
  {
    // check if displacement is Nan
    assert(isnan((*i).first)==false);
    // check if vertex is NULL
    assert((*i).second!=NULL);
    // check if squared virtual displacement is large
    if ((*i).first > 1E100)
    {
      cout << "\n\n Virtual_Disp::validateVirtDispMap2: Error. "
            << "virtual displacement > 1e100.\n";
      if ((*i).second!=NULL)
      {
        (*i).second->print(cout);
        cout << endl;
      }
      else
      {
        cout << "Virtual_Disp::validateVirtDispMap2: ERROR. "
              << "Vertex is NULL.\n\n";
      }
      exit(1);
    }
    // vertex closest face is NULL
    if ((*i).second->getClosestFace()==NULL)
    {
      cout << "Virtual_Disp::validateVirtDispMap2: "
            << "ERROR. vertex: " << (*i).second->getObject()->getName()
           << "->" << (*i).second->getIndex()
           << " in vd2_to_v has no closest point.\n";
      assert((*i).second->getClosestFace()!=NULL);
      exit(1);
    }
    // THE FOLLOWING CHECK DOES NOT WORK SINCE DIFFERENT GAINS
    // WERE USED FOR THE STORED VIRTUAL DISPLACEMENTS
    //// calculate the squared virtual displacment of this vertex
    //double sq_vd = (*i).second->getSqVirtualDisp(Gain_Schedule::instance().getRefGain());
    //// if calculated squared virtual displacement
    //// differs from stored vd2_to_v value
    //if (distinguishable(sq_vd,(*i).first)==true)
    //{
    //  cout << "\n\nVirtual_Disp::validateVirtDispMap2: "
    //        << "Error. computed squared vd (" << sq_vd
    //        << ") differ from stored value (" << (*i).first
    //        << ") for the following vertex in vd2_to_v.\n";
    //  (*i).second->print(cout);
    //  cout << endl;
    //  assert(distinguishable(sq_vd,(*i).first)==false);
    //  exit(1);
    //}
  }
}

/** Re-initialize class before beginning a new group.
 */

void Virtual_Disp::resetForNewGroup (const int & group)
{
  // clear and rebuild association maps
  // of squared virtual displacement and vertex pointer
  buildVirtDispMap(group);
  buildVirtDispMapComplement();
  // point iterator to pair with largest separation error
  seed = vd2_to_v.rbegin();
  // clear seed virtual displacement history
  seed_vd.clear();
  seed_ad.clear();
}

/** Remove associations of input vertex from all virtual displacement maps.
 * \param[in] v Vertex of to be removed.
 */

void Virtual_Disp::removeVertFromAllMaps (Vertex * v)
{
  // if vertex* found in v_to_vd2
  if (v_to_vd2.find(v)!=v_to_vd2.end())
  {
    if (Controls::instance().get_disable_messages()==false)
    {
      cout << "\nVirtual_Disp::removeVertFromAllMaps: "
            << "removing vertex from virtual displacement maps.\n";
      v->print(cout);
      cout << endl;
    }
    // remove vertex* from vd2_to_v
    removeVirtDispToVertAssoc(v,v_to_vd2[v]);
    // remove vertex* from v_to_vd2
    v_to_vd2.erase(v);
  }
}

