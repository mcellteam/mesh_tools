#include "box.h"

#include <iostream>

#include "container.h"
#include "object.h"
#include "space.h"
#include "vertex.h"

using std::cout;
using std::endl;

void Box::printBox (Space *s)
{
  cout << "Box indices ["
        << x << " "
        << y << " "
        << z << "]\n"
        << "Box range ["
        << xmin(s->world[0],s->space_length) << " "
        << xmax(s->world[0],s->space_length) << " "
        << ymin(s->world[2],s->space_length) << " "
        << ymax(s->world[2],s->space_length) << " "
        << zmin(s->world[4],s->space_length) << " "
        << zmax(s->world[4],s->space_length) << "]";
}

bool Box::faceIntersectionAlreadyKnown (Face *a,Face *b)
{
  // Face a is the intersected face
  // Face b is the intersecting face
  // 
  // get iterator pointing to location of Face a in it's object list
  ff_iterator i= a->findFaceInTable_intf();
  // if Face a is in it's object list
  if (i!=a->v[0]->o->intf.end())
  {
    // if Face b is in Face a's vector of intersecting faces
    if (find((*(*i).second).begin(),(*(*i).second).end(),b)!=(*(*i).second).end())
    {
      return true;
    }
    else 
    {
      return false;
    }
  }
  else 
  {
    return false;
  }
}


void Box::getFaceIntersection (Container *c) 
{
  // for each face in box
  for (f_iterator i=f.begin();i!=f.end();i++)
  {
    // for each face in box
    for (f_iterator j=f.begin();j!=f.end();j++)
    {
      // if faces are different
      if (*i!=*j)
      {
        // if faces intersect
        if (c->checkFaceFaceIntersections(*i,*j)) 
        {
          // if (*i) face is not already in object hashtable of intersected faces
          if ((*i)->faceInTable_intf()==false)
          {
            (*i)->addFaceToTable_intf();
          }
          // add (*j) face to vector of faces intersecting (*i) face
          // NOTE THE RISK OF DUPLICATE FACES IN THESE VECTORS
          (*i)->addFaceToVector(*j);
        }
      }
    }
  }
}

