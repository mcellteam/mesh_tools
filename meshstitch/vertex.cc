#include "vertex.h"

Vertex::Vertex(char *triplet)
:x(0),y(0),z(0),index(0)
{
  char val[80];
  char *eptr;

  // get past 'Vertex'
  while (strchr("Vertx",*triplet)!=NULL) {triplet++;}

  // grab vertex index
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  int i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  index = (int) strtod(val,&eptr);
  if (val==eptr)
  {
    index=0;
    x=y=z=0;
    printf("Error in reading vertex index\n");
    return;
  }

  // grab x coord
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  x = strtod(val,&eptr);
  if (val==eptr)
  {
    x=y=z=0;
    printf("Error in reading vertex\n");
    return;
  }

  // grab y coord
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  y = strtod(val,&eptr);
  if (val==eptr)
  {
    x=y=z=0;
    printf("Error in reading vertex\n");
    return;
  }

  // grab z coord
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  z = strtod(val,&eptr);
  if (val==eptr)
  {
    x=y=z=0;
    printf("Error in reading vertex\n");
    return;
  }
}

