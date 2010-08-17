#include "vertex.h"

#include <stdlib.h>
#include <string.h>

using std::cout;
using std::endl;

Vertex::Vertex (char* triplet)
  :index(0),p(),projection(0)
{

  char val[80];
  char *eptr;
  int i;

  char *cp=triplet;

  // get past 'Vertex'
  while (strchr("Vertx",*triplet)!=NULL) {triplet++;}

  // grab vertex index
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  index = atoi(val);
  if (val==eptr)
  {
    index=0;
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
  p.p[0]=strtod(val,&eptr);
  if (val==eptr)
  {
    printf("Error in reading vertex\n");
    printf("Error in reading vertex: string %s\n",cp);
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
  p.p[1]=strtod(val,&eptr);
  if (val==eptr)
  {
    printf("Error in reading vertex\n");
    printf("Error in reading vertex: string %s\n",cp);
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
  p.p[2]=strtod(val,&eptr);
  if (val==eptr)
  {
    printf("Error in reading vertex\n");
    printf("Error in reading vertex: string %s\n",cp);
    return;
  }
}
