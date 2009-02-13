// Author: Justin Kinney
// Date: Sep 2008

#ifndef OCTREE_VISITOR_FACE_H
#define OCTREE_VISITOR_FACE_H

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


class Octree_Visitor_Face
   : public OctreeVisitor<Face>
{
/// standard object services ---------------------------------------------------
public:
            Octree_Visitor_Face(Vector3r,Vector3r);

   virtual ~Octree_Visitor_Face();

            //void print (void)
            //{
            //  for (fp_it i=matches.begin();i!=matches.end();i++)
            //  {
            //    (*i)->print(std::cout);
            //  }
            //}

            fp_it mybegin (void)
            {
              return matches.begin(); 
            }

            fp_it myend (void)
            {
              return matches.end(); 
            }
//            void sort (void)
//            {
//              std::sort(matches.begin(),matches.end());
//              matches.assign(matches.begin(),unique(matches.begin(),matches.end()));
//            }
            int num_faces (void)
            {
              return matches.size();
            }
            int num_leaves (void)
            {
              return leaves_visited;
            }
//            virtual void dump_(std::string const &indent);
//            void dump_node(std::string const &indent);
//            void dump_leaf(std::string const &indent);
private:
            Octree_Visitor_Face( const Octree_Visitor_Face& );
   Octree_Visitor_Face& operator=( const Octree_Visitor_Face& );


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
   Vector3r * mylower;
   Vector3r * myupper;
//   hashset_fp matches;
   vec_fp matches;
   int leaves_visited;
};

#endif
