#include "space.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "box.h"
#include "controls.h"
#include "vertex.h"

using std::cout;
using std::endl;

Space::~Space (void)
{
  b_iterator i;
  for (i=b.begin();i!=b.end();i++)
  {
    delete *i;
  }
}

int Space::screenIndex (int i,const char *c)
{
  int a=0;
  if      (!strcmp(c,"x")){a=num_space[0]-1;}
  else if (!strcmp(c,"y")){a=num_space[1]-1;}
  else if (!strcmp(c,"z")){a=num_space[2]-1;}
  else {cout << "Received unexpected string.\n";exit(1);}
  //
  if (i<0){i=0;}
  if (i>a){i=a;}
  return i;
}

void Space::getBoxesFor3DLocations (double dr[6],vec_b& bp)
{
  // assume p = [xmin,xmax,ymin,ymax,zmin,zmax]
  int ir[6]; // index range
  // compute 3D index range that includes 3D location range
  ir[0] = location2Index(dr[0],"x");   // -x
  ir[1] = location2Index(dr[1],"x");	//  x
  ir[2] = location2Index(dr[2],"y");	// -y
  ir[3] = location2Index(dr[3],"y");	//  y
  ir[4] = location2Index(dr[4],"z");   // -z
  ir[5] = location2Index(dr[5],"z");	//  z
  /*	// DEBUG
        cout << "\nSpace::getBoxesFor3DLocations: "
        << "data range ["
        << dr[0] << ","
        << dr[1] << ","
        << dr[2] << ","
        << dr[3] << ","
        << dr[4] << ","
        << dr[5] << "]\n";
        cout << "\nSpace::getBoxesFor3DLocations: "
        << "num_space ["
        << num_space[0] << ","
        << num_space[1] << ","
        << num_space[2] << "]\n";
        cout << "\nSpace::getBoxesFor3DLocations: "
        << "world ["
        << world[0] << ","
        << world[1] << ","
        << world[2] << ","
        << world[3] << ","
        << world[4] << ","
        << world[5] << "]\n";
        cout << "\nSpace::getBoxesFor3DLocations: "
        << "index range ["
        << ir[0] << ","
        << ir[1] << ","
        << ir[2] << ","
        << ir[3] << ","
        << ir[4] << ","
        << ir[5] << "]\n";
        */	// DEBUG
  getBoxesFor3DIndices(ir,bp);
}

void Space::getBoxesFor3DIndices (int br[6],vec_b& bp)
{
  // assume br = [xmin,xmax,ymin,ymax,zmin,zmax]
  int x0=br[0],x1=br[1],y0=br[2],y1=br[3],z0=br[4];
  b_iterator i;
  i = b.begin()+num_space[0]*((z0-1)*num_space[1]+(y0-1));
  for (int z = br[4];z<br[5]+1;z++)
  {
    i += num_space[0]*num_space[1];
    for (int y = br[2];y<br[3]+1;y++)
    {
      i += num_space[0];
      if (x1==x0)
      {
        bp.push_back(*(i+x0));
      }
      else
      {
        bp.insert(bp.end(),i+x0,i+x1+1);
      }
    }
    i-=(y1-y0+1)*num_space[0];
  }
  if (bp.empty()==true)
  {
    cout << "br ["
          << br[0] << ","
          << br[1] << ","
          << br[2] << ","
          << br[3] << ","
          << br[4] << ","
          << br[5] << "]\n";
    exit(1);
  }
}

int Space::indices2Index (int x,int y,int z)
{
  //	return z*num_space[0]*num_space[1]+y*num_space[0]+x;
  return num_space[0]*(z*num_space[1]+y)+x;
}

void Space::index2Range (int i,double r[2],char *c)
{
  r[0] = i*space_length;
  r[1] = (i+1)*space_length;
  if (!strcmp(c,"x"))
  {
    r[0] += world[0];
    r[1] += world[0];
  }
  else if (!strcmp(c,"y"))
  {
    r[0] += world[2];
    r[1] += world[2];
  }
  else if (!strcmp(c,"z"))
  {
    r[0] += world[4];
    r[1] += world[4];
  }
}

int Space::location2Index (double ss,const char * c)
{
  // since world[] is calculated only once, face limits may be less than world[0,2,4]
  // and greater than world[1,3,5]. Consequently, box_range as computed above may
  // yield negative indices or indices beyond declared range, respectively.
  // Therefore, constrain box_range to known index range.
  // Effectively, border boxes are understood to extend to infinity
  int a=0;
  if (!strcmp(c,"x"))
  {
    a = (int) floor((ss-world[0])/space_length);
    //		cout << "\nSpace::location2Index: "
    //		<< "x index (proposed,accepted) = " << a;
    a = screenIndex(a,"x");
    //		cout << "," << a << endl;
  }
  else if (!strcmp(c,"y"))
  {
    a = (int) floor( (ss-world[2])/space_length );
    //		cout << "\nSpace::location2Index: "
    //		<< "y index (proposed,accepted) = " << a;
    a = screenIndex(a,"y");
    //		cout << "," << a << endl;
  }
  else if (!strcmp(c,"z"))
  {
    a = (int) floor( (ss-world[4])/space_length );
    //		cout << "\nSpace::location2Index: "
    //		<< "z index (proposed,accepted) = " << a;
    a = screenIndex(a,"z");
    //		cout << "," << a << endl;
  }
  else {cout << "Received unexpected string.\n"; exit(1);}
  return a;
}

void Space::initBoxes (double num_faces) 
{
  // subdivide space
  double world_volume = (world[1]-world[0])*(world[3]-world[2])*(world[5]-world[4]);
  double box_volume = Controls::instance().get_faces_per_box() * world_volume / num_faces;
  	cout << "\nSpace::initBoxes: "
  	<< "FACES_PER_BOX = " << Controls::instance().get_faces_per_box()
  	<< ", world_volume = " << world_volume
  	<< ", num_faces = " << num_faces << "\n";
  space_length = pow ( fabs ( box_volume ), 1.0 / 3.0 );
  	cout << "Space::initBoxes: "
  	<< "space_length = " << space_length << "\n";
  	cout << "Space::initBoxes: "
        << "world ["
        << world[0] << ","
        << world[1] << ","
        << world[2] << ","
        << world[3] << ","
        << world[4] << ","
        << world[5] << "]\n";
        
  num_space[0] = (int) ceil( (world[1]-world[0])/space_length );
  num_space[1] = (int) ceil( (world[3]-world[2])/space_length );
  num_space[2] = (int) ceil( (world[5]-world[4])/space_length );
  	cout << "Space::initBoxes: "
        << "num_space ["
        << num_space[0] << ","
        << num_space[1] << ","
        << num_space[2] << "]\n\n";
        
  num_boxes = num_space[0]*num_space[1]*num_space[2];
  // allocate memory for boxes
  b.reserve(num_boxes);
  // store box limits in boxes class
  // for each box
  for (int z =0;z<num_space[2];z++) 
  {
    for (int y =0;y<num_space[1];y++) 
    {
      for (int x =0;x<num_space[0];x++) 
      {
        b.push_back(new Box(x,y,z));
      }
    }
  }
}

void Space::computeBoxesToCheck (Face *f,vec_b &bp) 
{
  vec_d xv,yv,zv;
  // identify face bounding box limits
  xv.push_back(f->ptr_vertex(0)->getpN(0));
  xv.push_back(f->ptr_vertex(1)->getpN(0));
  xv.push_back(f->ptr_vertex(2)->getpN(0));
  yv.push_back(f->ptr_vertex(0)->getpN(1));
  yv.push_back(f->ptr_vertex(1)->getpN(1));
  yv.push_back(f->ptr_vertex(2)->getpN(1));
  zv.push_back(f->ptr_vertex(0)->getpN(2));
  zv.push_back(f->ptr_vertex(1)->getpN(2));
  zv.push_back(f->ptr_vertex(2)->getpN(2));
  sort(xv.begin(),xv.end());
  sort(yv.begin(),yv.end());
  sort(zv.begin(),zv.end());
  // grab face 3D location range
  double dr[6];
  dr[0] = xv[0];  // -x
  dr[1] = xv[2];	//  x
  dr[2] = yv[0];	// -y
  dr[3] = yv[2];	//  y
  dr[4] = zv[0];  // -z
  dr[5] = zv[2];	//  z
  /*	// DEBUG
        cout << "\nSpace::computeBoxesToCheck: "
        << "data range ["
        << dr[0] << ","
        << dr[1] << ","
        << dr[2] << ","
        << dr[3] << ","
        << dr[4] << ","
        << dr[5] << "]\n";
        */	// DEBUG
  // collect boxes to check
  getBoxesFor3DLocations(dr,bp);
}

void Space::deleteBoxes (void)
{
  b_iterator i;
  // for each box in space, clear vector of face*
  for (i=b.begin();i!=b.end();i++) 
  {
    delete (*i);
  }
  b.clear();
}

void Space::clearBoxes (void)
{
  b_iterator i;
  // for each box in space, clear vector of face*
  for (i=b.begin();i!=b.end();i++) {(*i)->clearFaces();}
}

void Space::recordFace (std::vector<Box*> &ptr,Face* f) 
{
  b_iterator i;
  // for each box, add face
  for (i=ptr.begin();i!=ptr.end();i++) 
  {
    (*i)->addFace(f);
  }
}

