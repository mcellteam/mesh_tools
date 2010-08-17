#include "boundary.h"

#include <algorithm>
#include <iostream>

#include "edge.h"
#include "object.h"
#include "vertex.h"

using std::cout;
using std::endl;

Boundary::Boundary (void)
:e(),end(NULL),begin(NULL),open(false)
{
}

bool Boundary::closed (void)
{
  if (begin==end)
  {
    open=false;
    return true;
  }
  else
  {
    return false;
  }
}

void Boundary::print (void)
{
  cout << "\nBEGIN VERTEX:\n";
  begin->print(cout);
  cout << endl;
  cout << "END VERTEX:\n";
  end->print(cout);
  cout << endl;
  for (e_iterator i=e.begin();i!=e.end();i++)
  {
    (*i)->print(cout);
    cout << endl;
  }
}

void Boundary::init (Edge *ee)
{
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  ee->getVertices(v1,v2,o1,o2);
  end=v1;
  begin=v2;
  open=true;
  e.clear();
}

bool Boundary::add (Edge *ee)
{
  // if edge not found in vector
  if (find(e.begin(),e.end(),ee)==e.end())
  {
    e.push_back(ee);
    return true;
  }
  else
  {
    return false;
  }
}

bool Boundary::edgeExtendsBoundary (Edge *ee)
{
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  ee->getVertices(v1,v2,o1,o2);
  if (v1==begin)
  {
    if (add(ee))
    {
      begin=v2;
    }
    else
    {
      cout << "Error: tried to add edge to boundary twice!\n";
      exit(1);
    }
    return true;
  }
  else if (v2==begin)
  {
    if (add(ee))
    {
      begin=v1;
    }
    else
    {
      cout << "Error: tried to add edge to boundary twice!\n";
      exit(1);
    }
    return true;
  }
  else if (v1==end)
  {
    if (add(ee))
    {
      end=v2;
    }
    else
    {
      cout << "Error: tried to add edge to boundary twice!\n";
      exit(1);
    }
    return true;
  }
  else if (v2==end)
  {
    if (add(ee))
    {
      end=v1;
    }
    else
    {
      cout << "Error: tried to add edge to boundary twice!\n";
      exit(1);
    }
    return true;
  }
  return false;
}

