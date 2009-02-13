// Author: Justin Kinney
// Date: Sep 2008

#ifndef OCTREE_VISITOR_UPDATE_H
#define OCTREE_VISITOR_UPDATE_H

#include "meshmorph.h"

#include "face.h"
#include "Octree.h" 
#include <iostream>

using hxa7241_graphics::OctreeVisitor;
using hxa7241_graphics::Vector3r;
using hxa7241_graphics::OctreeCell;
using hxa7241_graphics::OctreeData;
using hxa7241_graphics::Array;
using hxa7241_graphics::OctreeRoot;
using hxa7241_graphics::dword;
using hxa7241_graphics::OctreeBranch;


class Octree_Visitor_Update
   : public OctreeVisitor<Face>
{
/// standard object services ---------------------------------------------------
public:
            Octree_Visitor_Update(Face *,Vector3r,Vector3r,Vector3r,Vector3r);

   virtual ~Octree_Visitor_Update();

private:
            Octree_Visitor_Update( const Octree_Visitor_Update& );
   Octree_Visitor_Update& operator=( const Octree_Visitor_Update& );


/// commands -------------------------------------------------------------------
/// octree visitor overrides
protected:
   virtual void  visitRoot  ( const OctreeCell* pRootCell,
                              const OctreeData& octreeData );
   virtual void  visitBranch( const OctreeCell* subCells[8],
                              const OctreeData& octreeData );
   virtual void  visitLeaf  ( const Array<const Face*>& items,
                              const OctreeData& octreeData );

   // any other commands ...


/// queries --------------------------------------------------------------------
   // any queries ...


/// fields ---------------------------------------------------------------------
private:
   Face * face;
   Vector3r * oldlower;
   Vector3r * oldupper;
   Vector3r * newlower;
   Vector3r * newupper;
};

#endif
