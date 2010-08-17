#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "contour.h"
#include "extra_vertex.h"
#include "face.h"
#include "group.h"
#include "meshstitch.h"
#include "site.h"
#include "vertex.h"

int main(int argc,char *argv[])
{

  if (argc != 5)
  {
    printf("\nSyntax: meshstitch input_file1 input_file2 z_value orientation\n\n");
    printf("Detail: The z_value refers to the z axis coordinate of\n");
    printf("        the common contour between the two input meshes. \n");
    printf("        Orientation is '1' if the polygons to be removed\n");
    printf("        from input_file2 lie on the +z side of the z_value, and\n");
    printf("        and orientation is '-1' if the polygons to be removed\n");
    printf("        from input_file2 lie on the -z side of the z_value.\n");
    printf("        The polygons to be removed from input_file1 are assumed\n");
    printf("        to lie on the opposite of the z_value specified by orientation.\n");
    printf("Output: input_file1.clipped and input_file2.clipped.\n\n");
    return 1;
  }

  fprintf(stderr,"Warning. This algorithm may assumme that face vertices v1 v2 v3 are listed\n");
  fprintf(stderr,"in ccw order and that the vector formed by the cross product of\n");
  fprintf(stderr,"vectors v1v2 and v2v3 points out of the mesh.\n");

  double epsilon = 1e-14;
  // vertex linked lists for first and second file
  std::list<Vertex> vlh1,vlh2;

  // face linked lists for first and second file
  std::list<Face> flh1,flh2;

  // vertex with target z_value linked lists for first and second file
  std::list<Vertex> v1_match_h,v2_match_h;

  // shared vertex linked lists for first and second file
  std::list<Vertex> v1_shared_h,v2_shared_h;

  // candidate bad vertex linked lists for first and second file
  std::list<Vertex> cbv1h,cbv2h;

  // candidate bad face linked lists for first and second file
  std::list<Face> cbf1h,cbf2h;

  // bad vertex linked lists for first and second file
  std::list<Vertex> bv1h,bv2h;

  // bad face linked lists for first and second file
  std::list<Face> bf1h,bf2h;

  // extra vertex face search pool
  std::list<Face> evfsp1,evfsp2;

  // get z_value
  char *eptr;
  int z_value = (int) (strtod(argv[3],&eptr));

  // get target z_value
  int orientation = (int) strtod(argv[4],&eptr);

  ////////// get data /////////

  // open first file
  char *infile1 = argv[1];
  FILE *F = fopen(infile1,"r");
  if (!F)
  {
    printf("Couldn't open input file %s\n",infile1);
    return 1;
  }

  // for every line in first file
  char line[2048];
  for (char *str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F))
  {

    // skip leading whitespace
    while (strchr(" \t,",*str)!=NULL)
    { str++;}

    // if first character is V for Vertex, add new linked list class instance
    if (strchr("V",*str)!=NULL)
    {
      vlh1.push_back(Vertex(str));
    } 
    // if first character is F for Face, add new linked list class instance
    else if (strchr("F",*str)!=NULL)
    {
      flh1.push_back(Face(str));
    }
  }
  fclose(F);

  // opend second file
  char *infile2 = argv[2];
  F = fopen(infile2,"r");
  if (!F)
  {
    printf("Couldn't open input file %s\n",infile2);
    return 1;
  }

  // for every line in second file
  for (char *str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F))
  {

    // skip leading whitespace
    while (strchr(" \t,",*str)!=NULL) { str++;}

    // if first character is V for Vertex, add new linked list class instance
    if (strchr("V",*str)!=NULL)
    {
      vlh2.push_back(Vertex(str));
    } 
    // if first character is F for Face, add new linked list class instance
    else if (strchr("F",*str)!=NULL)
    {
      flh2.push_back(Face(str));
    }
  }
  fclose(F);

  // find maximum vertex index
  int max_vert1=maxVert(vlh1);
  int max_vert2=maxVert(vlh2);
  // find max face index
  int max_face1=maxFace(flh1);
  int max_face2=maxFace(flh2);

  // create array of pointers to vertex list
  Vertex* vert_array1[max_vert1+1];
  Vertex* vert_array2[max_vert2+1];

  // fill arrays
  for (iter_list_v p=vlh1.begin();p!=vlh1.end();p++)
  {
    vert_array1[p->index]=&(*p);
  }
  for (iter_list_v p=vlh2.begin();p!=vlh2.end();p++)
  {
    vert_array2[p->index]=&(*p);
  }

  ////////// find shared vertices /////////

  // build linked list of first and second vertex elements with z_value matching target
  matchZ(vlh1,v1_match_h,z_value);
  matchZ(vlh2,v2_match_h,z_value);

  char name[32];
  if(0){strcpy(name,"v1_match_h");printVertices(v1_match_h,name);}
  if(0){strcpy(name,"v2_match_h");printVertices(v2_match_h,name);}

  // find shared vertices
  findSharedVertices(v1_match_h,v2_match_h,v1_shared_h,v2_shared_h,epsilon);

  if(0){strcpy(name,"v1_shared_h");printVertices(v1_shared_h,name);}
  if(0){strcpy(name,"v2_shared_h");printVertices(v2_shared_h,name);}

  // get extra vertex face search pool
  getFaceSearchPool(flh1,evfsp1,z_value,vert_array1);
  getFaceSearchPool(flh2,evfsp2,z_value,vert_array2);

  if(0){strcpy(name,"evfsp1");printFaces(evfsp1,name);}
  if(0){strcpy(name,"evfsp2");printFaces(evfsp2,name);}

  // group vertices
  Group g1,g2;
  g1.getContours(evfsp1,vert_array1,max_vert1,z_value,v1_shared_h);
  g2.getContours(evfsp2,vert_array2,max_vert2,z_value,v2_shared_h);

  if(0){strcpy(name,"FILE 1");g1.printContourFacesVertices(name);}
  if(0){strcpy(name,"FILE 2");g2.printContourFacesVertices(name);}

  ///// STRATEGY /////
  // (1) collect all vertices with z value half section
  //		thickness away from target z_value in orientation direction
  // (2) collect all polygons connected to candidate bad vertices from step 1
  // (3) for each collected polygon, if at least one vertex is a shared vertex from above
  // 		then any vertex on that polygon with z value half section
  //		thickness away from target z_value in orientation direction is bad
  // (4) for each collected polygon, if at least one vertex is a shared vertex
  // 		then polygon is bad

  // gather candidate bad vertices from first and second vertex linked list
  int val = z_value+orientation*25*(-1);
  candidateBadVertices(vlh1,cbv1h,val);
  val = z_value+orientation*25;
  candidateBadVertices(vlh2,cbv2h,val);

  // contour_tiler occassionally creates an extra vertex, so need to find them
  // gather extra vertices from first and second match vertex linked lists
  g1.extraVertices(v1_shared_h);
  g2.extraVertices(v2_shared_h);

  ///// find faces that have at least one vertex from cbv lists /////
  candidateBadFaces(flh1,cbf1h,cbv1h);
  candidateBadFaces(flh2,cbf2h,cbv2h);

  ///// find bad vertices and faces from first vertex and first face linked lists
  identifyBadVerticesAndFaces(cbf1h,cbv1h,v1_shared_h,g1,bf1h,bv1h);

  // rescan candidate bad faces using known bad vertices
  // looking for bad faces that do not have face vertex that is shared
  // Apparently, all bad vertices are members of a bad faces with shared vertices,
  // so no further searching for bad vertices is performed in findBadFaces.
  findBadFaces(cbf1h,bv1h,bf1h);

  ///// find bad vertices and faces from second vertex and second face linked lists
  identifyBadVerticesAndFaces(cbf2h,cbv2h,v2_shared_h,g2,bf2h,bv2h);

  // rescan candidate bad faces using known bad vertices
  // looking for bad faces that do not have face vertex that is shared 
  // Apparently, all bad vertices are members of a bad faces with shared vertices,
  // so no further searching for bad vertices is performed in findBadFaces.
  findBadFaces(cbf2h,bv2h,bf2h);

  if(0){strcpy(name,"bv1h");printVertices(bv1h,name);}
  if(0){strcpy(name,"bv2h");printVertices(bv2h,name);}

  if(0){strcpy(name,"vlh1");printVertices(vlh1,name);}
  if(0){strcpy(name,"vlh2");printVertices(vlh2,name);}

  if(0){strcpy(name,"bf1h");printFaces(bf1h,name);}
  if(0){strcpy(name,"bf2h");printFaces(bf2h,name);}

  if(0){strcpy(name,"flh1");printFaces(flh1,name);}
  if(0){strcpy(name,"flh2");printFaces(flh2,name);}

  // need to identify the vertices on either side of the extra vertices
  g1.lateralVertices(bf1h,z_value,vert_array1);
  g2.lateralVertices(bf2h,z_value,vert_array2);

  g1.convertLateral(vlh2,epsilon);
  g2.convertLateral(vlh1,epsilon);

  val = z_value+orientation*25;
  g1.gatherThirdDeleteFace(flh2,val,vert_array2);
  val = z_value+orientation*25*(-1);
  g2.gatherThirdDeleteFace(flh1,val,vert_array1);

  if(0){strcpy(name,"FILE 1");g1.printOtherContourData(name);}
  if(0){strcpy(name,"FILE 2");g2.printOtherContourData(name);}

  if(0){printf("ADDED TO FILE 1\n");}
  g2.addFacesVertices(flh1,vlh1,max_vert1,max_face1);
  if(0){printf("ADDED TO FILE 2\n");}
  g1.addFacesVertices(flh2,vlh2,max_vert2,max_face2);

  ///// remove bad vertices from first and second vertex list /////
  removeBadVertices(vlh1,bv1h);
  removeBadVertices(vlh2,bv2h);

  ///// remove bad faces from first and second face list /////
  removeBadFaces(flh1,bf1h);
  removeBadFaces(flh2,bf2h);

  if(0){strcpy(name,"flh1");printFaces(flh1,name);}
  if(0){strcpy(name,"flh2");printFaces(flh2,name);}

  ////////// write output to files //////////

  // open first file
  char outfile1[128];
  sprintf(outfile1,"%s.clipped",infile1);
  printf("\n%s written.\n",outfile1);
  F = fopen(outfile1,"w");
  if (!F)
  {
    printf("Couldn't open output file %s\n",outfile1);
    return 1;
  }

  // for every vertex in first file
  for (iter_list_v p=vlh1.begin();p!=vlh1.end();p++) 
  {
    sprintf(line,"Vertex %i  %.15g %.15g %.15g\n",p->index,p->x,p->y,p->z);
    fputs(line,F);
  }

  // for every face in first file
  for (iter_list_f p=flh1.begin();p!=flh1.end();p++)
  {
    sprintf(line,"Face %i  %i %i %i\n",p->index,p->v1,p->v2,p->v3);
    fputs(line,F);
  }
  fclose(F);

  // open second file
  char outfile2[128];
  sprintf(outfile2,"%s.clipped",infile2);
  printf("%s written.\n\n",outfile2);
  F = fopen(outfile2,"w");
  if (!F)
  {
    printf("Couldn't open output file %s\n",outfile2);
    return 1;
  }

  // for every vertex in second file
  for (iter_list_v p=vlh2.begin();p!=vlh2.end();p++) 
  {
    sprintf(line,"Vertex %i  %.15g %.15g %.15g\n",p->index,p->x,p->y,p->z);
    fputs(line,F);
  }

  // for every face in second file
  for (iter_list_f p=flh2.begin();p!=flh2.end();p++)
  {
    sprintf(line,"Face %i  %i %i %i\n",p->index,p->v1,p->v2,p->v3);
    fputs(line,F);
  }
  fclose(F);
}

int distinguishable(double a,double b,double eps)
{
  double c;

  c=a-b;

  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;

  if (b<a) eps*=a;
  else eps*=b;
  return (c>eps);
}

int maxVert (std::list<Vertex> & v)
{
  int max=-1;
  for (iter_list_v i=v.begin();i!=v.end();i++)
  {
    if (max<i->index)
    {
      max=i->index;
    }
  }
  return max;
}

void matchZ (std::list<Vertex> & v,std::list<Vertex> & th,int z_value)
{
  // build linked list of first vertex elements with z_value matching target
  for (iter_list_v q=v.begin();q!=v.end();q++)
  {
    if (q->z == z_value)
    {
      th.push_back(*q);
    }
  }
}

void candidateBadVertices (std::list<Vertex> & v, std::list<Vertex> & th, int val)
{
  // gather candidate bad vertices from first vertex linked list
  // for each vertex
  for (iter_list_v q=v.begin();q!=v.end();q++)
  {
    // if z value matches target
    if (q->z == val)
    {
      // store vertex 
      th.push_back(*q);
    }
  }
}

void printVertices (std::list<Vertex> & v,char *str)
{
  for (iter_list_v q=v.begin();q!=v.end();q++)
  {
    printf("%s %i %.15g %.15g %.15g\n",str,q->index,q->x,q->y,q->z);
  }
}

void removeBadVertices (std::list<Vertex> & V,std::list<Vertex> & B)
{
  // for each vertex
  iter_list_v p = V.begin();
  while (p!=V.end())
  {
    // search bad vertex list
    bool flag = true;
    iter_list_v q=B.begin();
    while (q!=B.end() && flag)
    {
      // if vertex is bad
      if ( p->index == q->index )
      {
        V.erase(p);
        B.erase(q);
        flag = false;
      }
      else
      {
        q++;
      }
    }
    if (flag) { p++;}
  }
}

void findSharedVertices (std::list<Vertex> & v1_match_h,
                         std::list<Vertex> & v2_match_h,
                         std::list<Vertex> & v1_shared_h,
                         std::list<Vertex> & v2_shared_h,
                         double epsilon)
{
  // find shared vertices
  iter_list_v p = v1_match_h.begin();
  while (p!=v1_match_h.end())
  {
    bool flag = true;
    iter_list_v q = v2_match_h.begin();
    while (q!=v2_match_h.end() && flag)
    {
      if ( !distinguishable(p->x,q->x,epsilon) &&
           !distinguishable(p->y,q->y,epsilon) &&
           !distinguishable(p->z,q->z,epsilon))
      {
        // store vertex from first file
        v1_shared_h.push_back(*p);
        v1_match_h.erase(p);
        // store vertex from second file
        v2_shared_h.push_back(*q);
        v2_match_h.erase(q);
        flag = false;
      }
      else
      {
        q++;
      }
    }
    if(flag){p++;}
  }
}

void identifyBadVerticesAndFaces (std::list<Face> & cbf,
                                  std::list<Vertex> & cbv,
                                  std::list<Vertex> & v_shared,
                                  Group &g,
                                  std::list<Face> & bfh,
                                  std::list<Vertex> & bvh)
{
  ///// find bad vertices and faces from first vertex and first face linked lists
  // for each candidate bad face
  iter_list_f p=cbf.begin();
  while (p!=cbf.end())
  {
    // search shared vertex list v_shared
    bool flag = true;
    iter_list_v q=v_shared.begin();
    while(q!=v_shared.end() && flag)
    {
      // if shared vertex is at least one of three face vertices
      if (p->vertexInFace(&(*q)))
      {
        // add face to bad face list
        bfh.push_back(*p);
        flag = false;
        // store pointer
        Face *temp = &(*p);
        cbf.erase(p);
        // search candidate bad vertex list for each face vertex
        iter_list_v pp=cbv.begin();
        while(pp!=cbv.end())
        {
          // if bad vertex is one of three face vertices
          if (temp->vertexInFace(&(*pp)))
          {
            // add vertex to bad vertex list
            bvh.push_back(*pp);
            cbv.erase(pp);
          }
          else
          {
            pp++;
          }
        }
      }
      if (q==v_shared.end()) break;
      q++;
    }
    if(flag)
    {
      // search extra vertex list
      // for each group
      for(int i=0;i<g.count;i++)
      {
        // for each extra vertex in group
        iter_vec_ev qq=g.c[i].extra.begin();
        while(qq!=g.c[i].extra.end() && flag)
        {
          // if extra vertex is at least one of three face vertices
          if (p->vertexInFace(qq->v))
          {
            // add face to bad face list
            bfh.push_back(*p);
            flag = false;
            // store pointer
            Face *temp = &(*p);
            cbf.erase(p);
            // search candidate bad vertex list for each face vertex
            iter_list_v pp=cbv.begin();
            while(pp!=cbv.end())
            {
              // if bad vertex is one of three face vertices
              if (temp->vertexInFace(&(*pp)))
              {
                // add vertex to bad vertex list
                bvh.push_back(*pp);
                cbv.erase(pp);
              }
              else
              {
                pp++;
              }
            }
          }
          qq++;
        }
      }
    }
    if(flag){ p++;}
  }
}

void candidateBadFaces (std::list<Face> & F,
                        std::list<Face> & th,
                        std::list<Vertex> & V)
{
  ///// find faces that have at least one vertex from cbv lists /////
  // for each face from first face linked list
  for (iter_list_f q=F.begin();q!=F.end();q++)
  {
    // search cbv1 list 
    bool flag = true;
    iter_list_v p=V.begin();
    while(p!=V.end() && flag)
    {
      // if candidate bad vertex is at least one of three face vertices
      if (q->vertexInFace(&(*p)))
      {
        // add face to candidate bad face list
        th.push_back(*q);
        flag = false;
      }
      if (p==V.end()) break;
      p++;
    }
  }
}

void findBadFaces (std::list<Face> & CBF,
                         std::list<Vertex> & V,
                         std::list<Face> & BF)
{
  // scan candidate bad faces searching for known bad vertices
  // for each candidate bad face
  iter_list_f p=CBF.begin();
  while (p!=CBF.end())
  {
    // search bad vertex list
    bool flag = true;
    iter_list_v q=V.begin();
    while(q!=V.end() && flag)
    {
      // if bad vertex is at least one of three face vertices
      if (p->vertexInFace(&(*q)))
      {
        // add face to bad face list
        BF.push_back(*p);
        flag = false;
        CBF.erase(p);
      }
      q++;
    }
    if(flag){ p++;}
  }
}

int maxFace (std::list<Face> & L)
{
  int max=0;
  for (iter_list_f p=L.begin();p!=L.end();p++)
  {
    if (max<p->index) {max=p->index;}
  }
  return max;
}

void getFaceSearchPool (std::list<Face> & F,
                        std::list<Face> & qh,
                        int z_value,
                        Vertex** vert_array)
{
  ///// gather faces with at least two vertices from match list
  // for each face
  for (iter_list_f p=F.begin();p!=F.end();p++)
  {
    // if at least two vertices match z_value
    if((( vert_array[p->v1]->z==z_value)&&
        ( vert_array[p->v2]->z==z_value))||
        ((vert_array[p->v1]->z==z_value)&&
         (vert_array[p->v3]->z==z_value))||
        ((vert_array[p->v2]->z==z_value)&&
         (vert_array[p->v3]->z==z_value)))
    {
      // add face pointer to list
      qh.push_back(*p);
    }
  }
}

void printFaces (std::list<Face> & L,char *str)
{
  for(iter_list_f q=L.begin();q!=L.end();q++)
  {
    printf("%s %i %i %i %i\n",str,q->index,q->v1,q->v2,q->v3);
  }
}

void removeBadFaces (std::list<Face> & F,std::list<Face> & B)
{
  // for each face
  iter_list_f p=F.begin();
  while (p!=F.end())
  {
    // search bad face list
    bool flag = true;
    iter_list_f q=B.begin();
    while (q!=B.end() && flag)
    {
      // if face is bad
      if ( p->index == q->index )
      {
        F.erase(p);
        B.erase(q);
        flag = false;
      }
      else
      {
        q++;
      }
    }
    if (flag) {p++;}
  }
}

int compare (const void* a, const void* b)
{
  return ( *(int*)a - *(int*)b );
}

