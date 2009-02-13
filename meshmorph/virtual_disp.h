// Author: Justin Kinney
// Date: Sep 2008

#ifndef VIRTUAL_DISP_H
#define VIRTUAL_DISP_H 1

#include <cmath>

#include "log.h"
#include "controls.h"
#include "meshmorph.h"

typedef std::map<Vertex*,float,ltv>                                 map_d;
typedef std::map<Vertex*,float,ltv>::iterator                       td_it;
typedef std::multimap<float,Vertex*,ltd>                            map_v;
typedef std::multimap<float,Vertex*,ltd>::iterator                  tv_it;
typedef std::multimap<float,Vertex*,ltd>::reverse_iterator          tv_rit;
typedef std::multimap<float,Vertex*,ltd>::const_reverse_iterator    tv_crit;
typedef std::multimap<float,Vertex*,ltd>::const_iterator            tv_cit;

typedef std::vector<std::string>::const_iterator                    s_cit;

class Virtual_Disp
{
private:
  static Virtual_Disp * only_one;
  map_v   vd2_to_v; // map of virtual displacement squared to vertex*
  map_d   v_to_vd2; // map of vertex* to virtual displacement squared 
  tv_rit      seed; // virt disp map seed of vertex move collection
  //vec_d    seed_vd; // store virtual displacements of seed vertices
  vec_s    seed_vd; // store virtual displacements of seed vertices
                    // along with iteration number
  vec_s    seed_ad; // store actual displacements of seed vertices
                    // along with iteration number
  Virtual_Disp                          (void);
  Virtual_Disp                          (Virtual_Disp const &);
  Virtual_Disp & operator =             (Virtual_Disp const &);
  ~Virtual_Disp                         (void);
public:
  static  Virtual_Disp & instance       (void);
  void    buildVirtDispMapComplement    (void);
  bool    findVirtDispToVertAssoc       (Vertex * const,float const &,tv_it &);
  void    buildVirtDispMap              (const int &);
  void    setVirtualDisp                (Vertex * const,float const &);
  void    updateVirtualDisp             (Vertex * const,bool const,float const &);
  void    removeVirtDispToVertAssoc     (Vertex * const,float const &);
  bool    getVertAndRank                (Vertex * const,tv_it &,int &);
  void    resetForNewGroup              (const int &);
  void    removeVertFromAllMaps         (Vertex *);
  void    validateVirtDispMapComplement (void);
  void    validateVirtDispMap2          (void);
  void    validateVirtDispMap           (void) const;

  /** Store seed actual displacement to detect steady state of model.
   */

  void addAdToSeedAd (const int & count,double ad)
  {
    seed_ad.push_back(Log::instance().format("%d %g",count,ad));
  }

  /** Store seed virtual displacement to detect steady state of model.
   */

  void addVdToSeedVd (const int & count)
  {
    //double myvd = sqrt((*seed).first)/sqrt(Controls::instance().get_search_radius_sq());
    double myvd = sqrt((*seed).first);
    seed_vd.push_back(Log::instance().format("%d %g",count,myvd));
    //seed_vd.push_back((*seed).first);
  }

  /** Get an iterator pointing to the first element in
   * seed actual displacement vector.
   * \return Iterator pointing to the first seed actual displacement element.
   */

  //d_cit beginSeedActDisp (void)
  s_cit beginSeedActDisp (void)
  {
    return seed_ad.begin();
  }

  /** Get an iterator pointing to one past the last element in
   * seed actual displacement vector.
   * \return Iterator pointing to one past the last seed actual displacement element.
   */

  //d_cit endSeedActDisp (void)
  s_cit endSeedActDisp (void)
  {
    return seed_ad.end();
  }

  /** Get an iterator pointing to the first element in
   * seed virtual displacement vector.
   * \return Iterator pointing to the first seed virtual displacement element.
   */

  //d_cit beginSeedVirtDisp (void)
  s_cit beginSeedVirtDisp (void)
  {
    return seed_vd.begin();
  }

  /** Get an iterator pointing to one past the last element in
   * seed virtual displacement vector.
   * \return Iterator pointing to one past the last seed virtual displacement element.
   */

  //d_cit endSeedVirtDisp (void)
  s_cit endSeedVirtDisp (void)
  {
    return seed_vd.end();
  }

  /** Retrieve an iterator pointing to seed of collection of vertices to move.
   * \return Iterator pointing to element of virtual displacment map.
   */

  tv_crit getSeed (void)
  {
    return seed;
  }

  /** Retrieve number of vertices in virtual displacement map.
   * \return Number of elements in virtual displacement map.
   */

  int getNumVertsInVirtDispMap (void)
  {
    return vd2_to_v.size();
  }

  /** Advance iterator to point to next element in map
   * sorted by virtual displacement.
   */

  void advanceSeedToNextLargestVert (void)
  {
    seed++;
  }

  /** Reset iterator to point to first element in map
   * sorted by virtual displacement.
   */

  void resetSeedToLargestVert (void)
  {
    seed = vd2_to_v.rbegin();
  }

  /** Get an iterator pointing to the first element in
   * virtual displacement map.
   * \return Iterator pointing to the first virtual displacement element.
   */

  tv_cit beginVirtDispMap (void)
  {
    return vd2_to_v.begin();
  }

  /** Get a reverse iterator pointing to one past the last element in
   * virtual displacement map.
   * \return Reverse iterator pointing to one past
   * the last virtual displacement element.
   */

  tv_crit rendVirtDispMap (void)
  {
    return vd2_to_v.rend();
  }

  /** Get an iterator pointing to one past the last element in
   * virtual displacement map.
   * \return Iterator pointing to one past
   * the last virtual displacement element.
   */

  tv_cit endVirtDispMap (void)
  {
    return vd2_to_v.end();
  }

};

#endif
