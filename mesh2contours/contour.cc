// Author: Justin Kinney
// Date: Apr 2010

#include <iostream>
#include <string.h>

#include "contour.h"
#include "vertex.h"

using std::cout;
using std::endl;

Contour::Contour (void)
  :p()
{
}

void Contour::print (std::string name)
{
  //if (strcmp(cs.style,"cp")==false)
  //{
  //  // for each vertex * in contour
  //  for (v_iterator i=p.begin();i!=p.end();i++)
  //  {
  //    cout << (*i)->pN[0] << " "
  //          << (*i)->pN[1] << " "
  //          << (*i)->pN[2] << " 1 0 0 1\n";
  //  }
  //}
  //else
  //{
    cout << "<Contour name=\"" << name
          << "\" hidden=\"false\" closed=\"true\" "
          << "simplified=\"true\" border=\"0 0 1\" "
          //<< "fill=\"1 1 0.501961\" mode=\"9\"\n"
          << "fill=\"0 0 1\" mode=\"-13\"\n"
          << " points=\"";
    // for each vertex* in contour
    //for (v_iterator i=p.begin();i!=p.end();i++)
    cout.precision(20);
    for (std::vector<Vertex*>::iterator i=p.begin();i!=p.end();i++)
    {
      //cout << (*i)->pN[0] << " "
      //      << (*i)->pN[1] << ",\n";
      cout << (*i)->getpN(0) << " "
            << (*i)->getpN(1) << ",\n";
    }
    cout << "\"/>\n";
  //}
}
	
