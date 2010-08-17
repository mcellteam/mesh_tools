#include "face.h"

Face::Face (char const * triplet)
:index(0),v1(0),v2(0),v3(0)
{
  char val[80];
  char *eptr;
  int i;

  // get past 'Face'
  while (strchr("Face",*triplet)!=NULL) {triplet++;}

  // grab Face index
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  index = (int) strtod(val,&eptr);
  if (val==eptr)
  {
    index=0;
    v1=v2=v3=0;
    printf("Error in reading face index\n");
    return;
  }

  // grab first vertex index
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  v1 = (int) strtod(val,&eptr);
  if (val==eptr)
  {
    v1=v2=v3=0;
    printf("Error in reading vertex index\n");
    return;
  }

  // grab second vertex index
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  v2 = (int) strtod(val,&eptr);
  if (val==eptr)
  {
    v1=v2=v3=0;
    printf("Error in reading vertex index\n");
    return;
  }

  // grab third vertex index
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  v3 = (int) strtod(val,&eptr);
  if (val==eptr)
  {
    v1=v2=v3=0;
    printf("Error in reading vertex index\n");
    return;
  }
}

