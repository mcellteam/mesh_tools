// Author: Justin Kinney
// Date: Sep 2008

#include "octree_visitor_update.h"

#include <iostream>

#include "face.h"
#include "Primitives.h"

using std::cout;
using std::endl;

/// standard object services ---------------------------------------------------
Octree_Visitor_Update::Octree_Visitor_Update(Face * f,Vector3r oldlow,Vector3r oldup,
                                                       Vector3r newlow,Vector3r newup)
  :face(f),oldlower(&oldlow),oldupper(&oldup),newlower(&newlow),newupper(&newup)
{
}


Octree_Visitor_Update::~Octree_Visitor_Update()
{
}


/// commands -------------------------------------------------------------------
void Octree_Visitor_Update::visitRoot
(
   const OctreeCell* pRootCell,
   const OctreeData& octreeData
)
{
   // continue visit traversal
   OctreeRoot::continueVisit( pRootCell, octreeData, *this );
}


void Octree_Visitor_Update::visitBranch
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

   // if subcell does not overlap old region of interest 
   // and does not overlap new region of interest 
   // then do not visit subcells
   if ((oldupper->getX()<lower.getX() ||
        oldlower->getX()>upper.getX()) &&
       (newupper->getX()<lower.getX() ||
        newlower->getX()>upper.getX())) return;
   if ((oldupper->getY()<lower.getY() ||
        oldlower->getY()>upper.getY()) &&
       (newupper->getY()<lower.getY() ||
        newlower->getY()>upper.getY())) return;
   if ((oldupper->getZ()<lower.getZ() ||
        oldlower->getZ()>upper.getZ()) &&
       (newupper->getZ()<lower.getZ() ||
        newlower->getZ()>upper.getZ())) return;

   // step through subcells (can be in any order...)
   for( dword i = 0;  i < 8;  ++i )
   {
      // continue visit traversal
      OctreeBranch::continueVisit( subCells, octreeData, i, *this );

   }
}


void Octree_Visitor_Update::visitLeaf
(
   const Array<const Face*>& items,
   const OctreeData&          octreeData
)
{
   // make some short aliases
   const Vector3r& lower( octreeData.getBound().getLowerCorner() );
   const Vector3r& upper( octreeData.getBound().getUpperCorner() );

   // does subcell overlap old region of interest?
   const bool old_overlap = oldupper->getX()>=lower.getX() && 
                            oldlower->getX()<=upper.getX() &&
                            oldupper->getY()>=lower.getY() &&
                            oldlower->getY()<=upper.getY() &&
                            oldupper->getZ()>=lower.getZ() &&
                            oldlower->getZ()<=upper.getZ();

   // does subcell overlap new region of interest?
   const bool new_overlap = newupper->getX()>=lower.getX() && 
                            newlower->getX()<=upper.getX() &&
                            newupper->getY()>=lower.getY() &&
                            newlower->getY()<=upper.getY() &&
                            newupper->getZ()>=lower.getZ() &&
                            newlower->getZ()<=upper.getZ();

   if (old_overlap==true && new_overlap==false)
   {
     //   the face should be in this cell
     //   remove face from this cell
     // for each face in leaf
     for( dword i = 0, end = items.getLength();  i < end;  ++i )
     {
       // find face
       if (items[i]==face)
       {
         // remove face
         const_cast< Array<const Face*>& >(items).remove(i);
         break;
       }
     }
   }
   else if (old_overlap==false && new_overlap==true)
   {
     //   face should NOT be in this cell
     //   add face to this cell
     const_cast< Array<const Face*>& >(items).append(face);
   }
   // cell overlaps new AND old
   //   face should already be in this cell
   // cell overlaps NEITHER new NOR old
   //   face should not be in this cell

}

