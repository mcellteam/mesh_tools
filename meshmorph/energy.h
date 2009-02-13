// Author: Justin Kinney
// Date: Sep 2008

#ifndef ENERGY_H
#define ENERGY_H 1

#include "meshmorph.h"

class Energy
{
private:
  vec_d   total_energy; // samples of total energy of model 
  static Energy * only_one;
  Energy                      (void);
  Energy                      (Energy const &);
  Energy & operator =         (Energy const &);
  ~Energy                     (void);
public:
  static Energy & instance    (void);
  void   computeGlobalEnergy  (void);
  void   writeStats           (std::ostream &);
};

#endif
