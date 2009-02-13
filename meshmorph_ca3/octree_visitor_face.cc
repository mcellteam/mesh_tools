// Author: Justin Kinney
// Date: Sep 2008

#include "octree_visitor_face.h"

#include <iostream>

#include "face.h"
#include "controls.h"
#include "Primitives.h"

using std::cout;
using std::endl;

/// standard object services ---------------------------------------------------
Octree_Visitor_Face::Octree_Visitor_Face(Vector3r low,Vector3r up)
:mylower(&low),myupper(&up),matches(),leaves_visited(0)
{
  matches.reserve(Controls::instance().get_vector_reserve());
}

//virtual void Octree_Visitor_Face::dump(std::string const &indent) = 0;

Octree_Visitor_Face::~Octree_Visitor_Face()
{
}


/// commands -------------------------------------------------------------------

void Octree_Visitor_Face::visitRoot
(
   const OctreeCell* pRootCell,
   const OctreeData& octreeData
)
{
   // continue visit traversal
   OctreeRoot::continueVisit( pRootCell, octreeData, *this );
}


void Octree_Visitor_Face::visitBranch
(
   const OctreeCell* subCells[8],
   const OctreeData& octreeData
)
{
   // subcell numbering:
   //    y z       6 7
   //    |/   2 3  4 5
   //     -x  0 1
   //
   // (in binary:)
   //    y z           110 111
   //    |/   010 011  100 101
   //     -x  000 001
   //

   const Vector3r& lower( octreeData.getBound().getLowerCorner() );
   const Vector3r& upper( octreeData.getBound().getUpperCorner() );

   // if subcell does not overlap region of interest 
   // then do not visit subcells
   if (myupper->getX()<lower.getX() ||
       myupper->getY()<lower.getY() ||
       myupper->getZ()<lower.getZ() ||
       mylower->getX()>upper.getX() ||
       mylower->getY()>upper.getY() ||
       mylower->getZ()>upper.getZ())
     return;

   // step through subcells (can be in any order...)
   for( dword i = 0;  i < 8;  ++i )
   {
      // continue visit traversal
      OctreeBranch::continueVisit( subCells, octreeData, i, *this );
   }
}


void Octree_Visitor_Face::visitLeaf
(
   const Array<const Face*>& items,
   const OctreeData&          octreeData
)
{
   // make some short aliases
   const Vector3r& lower( octreeData.getBound().getLowerCorner() );
   const Vector3r& upper( octreeData.getBound().getUpperCorner() );

   // if subcell does not overlap region of interest 
   // then do not visit subcells
   if (myupper->getX()<lower.getX() ||
       myupper->getY()<lower.getY() ||
       myupper->getZ()<lower.getZ() ||
       mylower->getX()>upper.getX() ||
       mylower->getY()>upper.getY() ||
       mylower->getZ()>upper.getZ())
     return;

   // write items
   for( dword i = 0, end = items.getLength();  i < end;  ++i )
   {
     if (items[i]->getFlag()==false)
     {
       matches.push_back(const_cast<Face*>(items[i]));
       const_cast<Face*>(items[i])->setFlag();
     }
   }
   // update stats
   leaves_visited++;
}

