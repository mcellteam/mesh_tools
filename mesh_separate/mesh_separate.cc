#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

class Vertex
{
public:
  double x,y,z;
  int index;
  Vertex(char *triplet);
};

Vertex::Vertex(char *triplet)
:x(0),y(0),z(0),index(0)
{
  char val[80];
  char *eptr;
  int i;

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

class Face
{
public:
  int index;	// Face index
  int v1,v2,v3;	// vertex indices
  Face(char *triplet);
};

Face::Face(char *triplet)
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

int compare (const void* a, const void* b )
{
  return ( *((int*)b+1) - *((int*)a+1) );
}

int main(int argc,char *argv[])
{
  if (argc != 2)
  {
    printf("\nSyntax: mesh_separate input_file\n\n");
    printf("Detail: The input mesh file is scanned for separate closed objects.\n");
    printf("Output: One file per separate object with the form input_file.#.mesh.\n\n");
    return 1;
  }

  fprintf(stderr,"\nInput mesh file is assumed to have vertex and\n");
  fprintf(stderr,"faces with sequentially increasing indeces. Run\n");
  fprintf(stderr,"the mesh input file through mesh_renumber first if\n");
  fprintf(stderr,"this criterion is not met.\n");

  std::vector<Vertex> vlh;
  std::vector<Face> flh;

  // open first file
  char *infile = argv[1];
  FILE *F = fopen(infile,"r");
  if (!F)
  {
    printf("Couldn't open input file %s\n",infile);
    return 1;
  }

  // for every line in first file
  char line[2048];
  for (char *str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F))
  {
    // skip leading whitespace
    while (strchr(" \t,",*str)!=NULL) { str++;}

    // if first character is V for Vertex, add new linked list class instance
    if (strchr("V",*str)!=NULL){
      vlh.push_back(Vertex(str));
    } 
    // if first character is F for Face, add new linked list class instance
    else if (strchr("F",*str)!=NULL){
      flh.push_back(Face(str));
    }
  }
  fclose(F);

  int num_vertices = vlh.size()+1;
  int num_faces = flh.size()+1;

  std::vector<int> faces(num_faces,0);
  std::vector<int> vertex(num_vertices,0);
  std::vector<std::vector<int> > vertices(num_vertices);

  ///// load vertices array /////
  // for each face
  for (std::vector<Face>::iterator q=flh.begin();q!=flh.end();q++)
  {

    // load face index into v1
    // load face index into v2
    // load face index into v3
    vertices[q->v1].push_back(q->index);
    vertices[q->v2].push_back(q->index);
    vertices[q->v3].push_back(q->index);
  }


  ///// coarsely separate closed surfaces /////
  int next_group = 1;
  int current_group;
  // for each vertex
  for (int i=1;i<num_vertices;i++)
  {
    current_group = 0;
    // for each associated face
    //for (q=vertices[i];q!=NULL;q=q->next) {
    for (std::vector<int>::iterator q=vertices[i].begin();q!=vertices[i].end();q++)
    {
      // look up face group
      // if not last link in list
      assert(*q<static_cast<int>(faces.size()));
      if (faces[*q])
      {
        current_group = faces[*q];
      }
    }

    // if a group was found
    if (current_group)
    {
      // set all faces to current_group
      // for each associated face
      //for (q=vertices[i];q!=NULL;q=q->next) {
      for (std::vector<int>::iterator q=vertices[i].begin();q!=vertices[i].end();q++)
      {
        //assert(*q<faces.size());
        assert(*q<static_cast<int>(faces.size()));
        faces[*q] = current_group;
      }
    } else
    {
      // no group was found
      // set all faces to next_group
      // for each associated face
      //for (q=vertices[i];q!=NULL;q=q->next) {
      for (std::vector<int>::iterator q=vertices[i].begin();q!=vertices[i].end();q++)
      {
        //assert(*q<faces.size());
        assert(*q<static_cast<int>(faces.size()));
        faces[*q] = next_group;
      }
      // increment next_group
      next_group++;
    }

    }

    ///// amalgamate groups /////
    bool amalgamate = true;
    int smallest_group;
    bool swap;
    while(amalgamate)
    {
      amalgamate = false;
      // for each vertex
      for (int i=1;i<num_vertices;i++)
      {

        // for each associated face
        swap = false;
        //q=vertices[i];
        assert(vertices[i].front()<static_cast<int>(faces.size()));
        smallest_group = faces[vertices[i].front()];
        //for (q=vertices[i];q!=NULL;q=q->next) {
        for (std::vector<int>::iterator q=vertices[i].begin();q!=vertices[i].end();q++)
        {
          // look up face group
          //assert(*q<faces.size());
          assert(*q<static_cast<int>(faces.size()));
          if (faces[*q]!=smallest_group){ swap=true;amalgamate = true;}
          if (faces[*q]<smallest_group){smallest_group=faces[*q];}
        }

        // if two or more groups were found
        if (swap)
        {
          // set all faces to smallest_group
          // for each associated face
          //for (q=vertices[i];q!=NULL;q=q->next) {
          for (std::vector<int>::iterator q=vertices[i].begin();q!=vertices[i].end();q++)
          {
            // if not smallest_group
            //assert(*q<faces.size());
            assert(*q<static_cast<int>(faces.size()));
            if (faces[*q]!=smallest_group){
              //swap all occurences of group in faces array with smallest_group
              for (int j=1;j<num_faces;j++)
              {
                //assert(*q<faces.size());
                assert(*q<static_cast<int>(faces.size()));
                //assert(j<faces.size());
                assert(j<static_cast<int>(faces.size()));
                if(faces[j]==faces[*q]) { faces[j]=smallest_group; }
              }
            }
          }
        }

        }
      }


      // for each vertex
      for (int i=1;i<num_vertices;i++)
      {
        // for first associated face
        //q=vertices[i];
        // set vertex group to face group
        assert(vertices[i].front()<static_cast<int>(faces.size()));
        vertex[i] = faces[vertices[i].front()];
      }

      ///// count number of faces in each group /////
      int max_groups = 65536;
      int face_count[max_groups][2];
      int num_groups = 0;
      //for each face
      for (int i=1;i<num_faces;i++)
      {
        int loc=-1;
        //search for group in face_count
        bool flag = false;
        for (int j=0;j<num_groups;j++)
        {
          // if face group found in face_count record location
          assert(i<static_cast<int>(faces.size()));
          if(face_count[j][0]==faces[i]){loc = j;flag = true;}
        }
        // if face group found in face_count, increment count
        if (flag) {assert(loc>=0);face_count[loc][1]++;}
        else
        {
          //add face group to face_count
          //assert(i<faces.size());
          assert(i<static_cast<int>(faces.size()));
          face_count[num_groups][0] = faces[i];
          face_count[num_groups][1] = 1;
          //check if room for more groups
          if (num_groups+1 < max_groups) {num_groups++;}
          else {printf("Increase size of max_groups (currently = %i)\n",max_groups);}
        }
      }

      ///// sort face_count /////
      qsort(face_count,num_groups,2*sizeof(int),compare);

      ///// output groups /////
      // for each group
      for (int j=0;j<num_groups;j++)
      {

        // open file
        char prefix[1024];
        int L = strlen(infile)-5;
        strncpy (prefix,infile,L);
        prefix[L]='\0';
        //sprintf(outfile,"%s.%i",infile,j);
        char outfile[128];
        sprintf(outfile,"%s_%i.mesh",prefix,j);
        printf("\n%s written (%i faces).\n\n",outfile,face_count[j][1]);
        F = fopen(outfile,"w");
        if (!F)
        {
          printf("Couldn't open output file %s\n",outfile);
          return 1;
        }

        // for each vertex
        //for (p=vend;p!=NULL;p=p->previous) {
        for (std::vector<Vertex>::iterator p=vlh.begin();p!=vlh.end();p++)
        {
          // if vertex is in current group
          if(vertex[p->index]==face_count[j][0])
          {
            sprintf(line,"Vertex %i  %.15g %.15g %.15g\n",p->index,p->x,p->y,p->z);
            fputs(line,F);
          }
        }

        // for each face
        //for (p=fend;p!=NULL;p=p->previous) {
        for (std::vector<Face>::iterator p=flh.begin();p!=flh.end();p++)
        {
          // if face is in current group
          assert(p->index<static_cast<int>(faces.size()));
          if(faces[p->index]==face_count[j][0])
          {
            sprintf(line,"Face %i  %i %i %i\n",p->index,p->v1,p->v2,p->v3);
            fputs(line,F);
          }
        }
        fclose(F);

      }
      }
