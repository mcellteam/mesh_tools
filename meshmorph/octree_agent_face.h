// Author: Justin Kinney
// Date: Sep 2008

#ifndef OCTREE_AGENT_FACE_H
#define OCTREE_AGENT_FACE_H

#include "meshmorph.h"

#include "face.h"
#include "Octree.h" 

using namespace hxa7241_graphics;

class Octree_Agent_Face
: public OctreeAgent<Face>
{
  /// standard object services ---------------------------------------------------
public:
  Octree_Agent_Face() {};

  virtual ~Octree_Agent_Face() {};
private:
  Octree_Agent_Face( const Octree_Agent_Face& );
  Octree_Agent_Face& operator=( const Octree_Agent_Face& );


  /// queries --------------------------------------------------------------------
  /// octree agent overrides
protected:
  virtual bool  isOverlappingCell ( const Face&    item,
                                    const Vector3r& lowerCorner,
                                    const Vector3r& upperCorner )        const;

  // could also override getSubcellOverlaps to provide more efficent
  // calculation (boundary testing can be shared).
};

#endif
