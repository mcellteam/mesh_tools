// Author: Justin Kinney
// Date: Sep 2008

#include "octree_visitor_check_face.h"

#include <iostream>

#include "face.h"
#include "Primitives.h"

using std::cout;
using std::endl;

// standard object services ---------------------------------------------------
Octree_Visitor_Check_Face::Octree_Visitor_Check_Face(Face * f,Vector3r low,Vector3r up)
  :face(f),face_lower(&low),face_upper(&up)
{
}


Octree_Visitor_Check_Face::~Octree_Visitor_Check_Face()
{
}


// commands -------------------------------------------------------------------
void Octree_Visitor_Check_Face::visitRoot
(
   const OctreeCell* pRootCell,
   const OctreeData& octreeData
)
{
   // continue visit traversal
   OctreeRoot::continueVisit( pRootCell, octreeData, *this );
}

void dump_this(const Vector3r& lower,const Vector3r& upper)
{
  printf("[%.15g, %.15g, %.15g]-[%.15g, %.15g, %.15g]\n",
          lower.getX(), lower.getY(), lower.getZ(), upper.getX(), upper.getY(), upper.getZ());
}

void Octree_Visitor_Check_Face::visitBranch
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
   if (face_upper->getX()<lower.getX() ||
       face_upper->getY()<lower.getY() ||
       face_upper->getZ()<lower.getZ() ||
       face_lower->getX()>upper.getX() ||
       face_lower->getY()>upper.getY() ||
       face_lower->getZ()>upper.getZ())
     return;

   // step through subcells (can be in any order...)
   for( dword i = 0;  i < 8;  ++i )
   {
      // continue visit traversal
      OctreeBranch::continueVisit( subCells, octreeData, i, *this );
   }
}

void Octree_Visitor_Check_Face::visitLeaf
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
   if (face_upper->getX()<lower.getX() ||
       face_upper->getY()<lower.getY() ||
       face_upper->getZ()<lower.getZ() ||
       face_lower->getX()>upper.getX() ||
       face_lower->getY()>upper.getY() ||
       face_lower->getZ()>upper.getZ())
     return;

   // look for face in item list
   for( dword i = 0, end = items.getLength();  i < end;  ++i )
   {
     if (face!=NULL)
     {
       if (items[i]==face) return;
     }
     else
       return;
   }
   // face not found, signal error
   if (face!=NULL)
   {
     cout.precision(12);
     cout << "\nOctree_Visitor_Check_Face::visitLeaf: "
           << "Error face overlaps octree cell but is not contained therein.\n";
     face->print(cout);
     cout << "face bounding box ["
           << face_lower->getX() << " "
           << face_upper->getX() << " "
           << face_lower->getY() << " "
           << face_upper->getY() << " "
           << face_lower->getZ() << " "
           << face_upper->getZ() << "]\n";
     cout << "octree cell ["
           << lower.getX() << " "
           << upper.getX() << " "
           << lower.getY() << " "
           << upper.getY() << " "
           << lower.getZ() << " "
           << upper.getZ() << "]\n";
     exit(1); 
   }
}

