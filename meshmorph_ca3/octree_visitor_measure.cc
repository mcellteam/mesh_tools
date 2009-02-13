// Author: Justin Kinney
// Date: Sep 2008

#include "octree_visitor_measure.h"

#include <iostream>

#include "face.h"
#include "controls.h"
#include "Primitives.h"

using std::cout;
using std::endl;

/// standard object services ---------------------------------------------------
Octree_Visitor_Measure::Octree_Visitor_Measure(void)
  :max_depth(0),num_leaves(0)
{
}

//virtual void Octree_Visitor_Face::dump(std::string const &indent) = 0;

Octree_Visitor_Measure::~Octree_Visitor_Measure()
{
}


/// commands -------------------------------------------------------------------

void Octree_Visitor_Measure::visitRoot
(
   const OctreeCell* pRootCell,
   const OctreeData& octreeData
)
{
   // continue visit traversal
   OctreeRoot::continueVisit( pRootCell, octreeData, *this );
}


void Octree_Visitor_Measure::visitBranch
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

   // step through subcells (can be in any order...)
   for( dword i = 0;  i < 8;  ++i )
   {
      // continue visit traversal
      OctreeBranch::continueVisit( subCells, octreeData, i, *this );

      // maybe do something here...
   }
}


void Octree_Visitor_Measure::visitLeaf
(
   const Array<const Face*>& items,
   const OctreeData&          octreeData
)
{
  // avoid warning message
  if (items.getLength()<-100) {}
  //
  num_leaves++;

  if (octreeData.getLevel()>max_depth)
    max_depth = octreeData.getLevel();
}

