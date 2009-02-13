// Author: Justin Kinney
// Date: Sep 2008

#ifndef OCTREE_VISITOR_MEASURE_H
#define OCTREE_VISITOR_MEASURE_H

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


class Octree_Visitor_Measure
   : public OctreeVisitor<Face>
{
/// standard object services ---------------------------------------------------
public:
            Octree_Visitor_Measure(void);

   virtual ~Octree_Visitor_Measure();

            int get_max_depth (void)
            {
              return max_depth;
            }
            int get_num_leaves (void)
            {
              return num_leaves;
            }
private:
            Octree_Visitor_Measure( const Octree_Visitor_Measure& );
   Octree_Visitor_Measure& operator=( const Octree_Visitor_Measure& );


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
   int max_depth;
   int num_leaves; 
};

#endif
