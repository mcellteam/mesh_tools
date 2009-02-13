// Author: Justin Kinney
// Date: Sep 2008

#include "octree_agent_face.h"

/// queries --------------------------------------------------------------------
bool Octree_Agent_Face::isOverlappingCell
(
 const Face&    item,
 const Vector3r& lowerCorner,
 const Vector3r& upperCorner
 ) const
{
  //const Vector3r& itemLower( item.getLower() );
  //const Vector3r& itemUpper( item.getUpper() );
  Vector3r itemLower,itemUpper;
  item.getBoundingBox(itemLower,itemUpper);

  // check the two ranges overlap in every dimension
  bool isOverlap = true;
  for( int i = 3;  i-- > 0; )
  {
    //isOverlap &= (itemLower[i] < upperCorner[i]) &
    //      (itemUpper[i] > lowerCorner[i]);
    isOverlap &= (itemLower[i] <= upperCorner[i]) &
          (itemUpper[i] >= lowerCorner[i]);
  }
   
  return isOverlap;
}

