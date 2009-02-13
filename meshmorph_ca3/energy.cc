// Author: Justin Kinney
// Date: Sep 2008

#include "energy.h"

#include <iostream>

#include "meshmorph.h"

#include "controls.h"
#include "container.h"
#include "edge.h"

using std::cout;
using std::cerr;
using std::endl;
using std::left;

Energy * Energy::only_one = NULL;

Energy & Energy::instance(void)
{
  // Not thread-safe.
  // -- lock mutex
  if (only_one == NULL)
    only_one = new Energy();
  // -- unlock mutex
  return *only_one;
}

Energy::Energy (void)
:total_energy()
{
}

///** Calculate cumulative energy in entire model,
// * and store in class member.
// */
//
//void Energy::computeGlobalEnergy (void)
//{
//  // Note face intersection force is not included in energy calculation.
//  cout << "Compute global energy..........................";
//  cout.flush();
//  double energy=0;
//  vector3 dummy;
//  // for each object in container
//  Container & c(Container::instance());
//  for (o_it i=c.o.begin();i!=c.o.end();++i)
//  {
//    // for each vertex in object
//    for (v_it j=i->v.begin();j!=i->v.end();++j)
//    {
//      // Note that if vertex does not have a closest face,
//      // then the max ecw energy is returned.
//      // add vertex ecw energy
//      // 'false' means 'do not compute force', hence dummy
//      energy += j->getEcwForceEnergy(dummy,false);
//    }
//    // for each edge in object
//    for (e_it j=i->e.begin();j!=i->e.end();++j)
//    {
//      // choice of v1 is arbitrary, v2 could have been used
//      // 'false' means 'do not compute force', hence dummy
//      energy += (j->getStretchForceEnergy(j->getV1(),dummy,false)
//                +j->getAngleForceEnergy(0,dummy,false));
//    }
//  }
//  total_energy.push_back(energy);
//  cout << "complete.\n";
//  cout.flush();
//}

/** Print total model energy statistics.
 *
 * \param[in] target Pre-initialized output stream.
 */

void Energy::writeStats (std::ostream & target)
{
  target << "\nTotal model energy (sampled every "
        << Controls::instance().get_energy_sample_period() << " vertex moves) = "; 

  for (d_it i=total_energy.begin();i!=total_energy.end();++i)
  {
    target << *i << endl;
  }
  target << endl;
}

