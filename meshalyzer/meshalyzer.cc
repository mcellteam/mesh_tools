#include <cmath>
#include <iostream>

#include "container.h"
#include "controls.h"
#include "edge.h"
#include "meshalyzer.h"
#include "object.h"
#include "space.h"

using std::cout;
using std::cerr;
using std::endl;

int main(int argc,char **argv)
{
  // instantiate controls class
  Controls & cs(Controls::instance());

  // parse command line 
  cs.parse(argc,argv,cs.getUsageMessage());

  // create container, objects, vertices, faces, edges, and find adjacencies
  Container c;

  // initialize space data structure
  Space s;

  // if single input file was found
  if (cs.get_folder_passed()==0)
  {
    // save filename
    c.addFile(cs.get_inpath());
    // update index
    //c.num_files++;
    c.incrementNumFiles();
    // build data structure of mesh
    Object *obj=c.processFile(c.getFile(0));
    // if either no vertices or no faces were found
    if (obj!=NULL)
    {
      // check mesh integrity
      obj->checkIntegrity();
      // if mesh file is not fatally flawed
      // i.e. no missing or orphaned vertices
      // and must be contiguous vertex and face indexing.
      // Note score may be set to 1 for other reasons in analyze()
      if (obj->goodIntegrity()==true)
      {
        obj->createEdges();
        obj->findVertexAdjacencies();
        // partition space
        c.boundWorld(s);
        s.initBoxes(obj->getNumFaces());
        c.assignFacesToBoxes(s);
        // analyze and respond
        obj->analyze(&c,cs,s);
      }
      obj->print(cs);
    }
  }
  else
  { 
    // else scan folder
    c.scanDir(cs.get_inpath().c_str());
    // for each file in folder
    for (int count=0;count<c.getNumFiles();count++)
    {
      // build data structure of mesh
      Object *obj=c.processFile(cs.get_inpath()+c.getFile(count));
      // if neither zero vertices nor zero faces were found
      if (obj!=NULL)
      {
        obj->checkIntegrity();
        if (obj->goodIntegrity()==true)
        {
          obj->createEdges();
          obj->findVertexAdjacencies();
          // partition space
          c.boundWorld(s);
          s.deleteBoxes();
          s.initBoxes(obj->getNumFaces());
          c.assignFacesToBoxes(s);
          // analyze and respond
          obj->analyze(&c,cs,s);
          c.update(obj);
          obj->store(c);
        }
        else
        {
          c.setIntegrity(false);
        }
        // print stats
        if (cs.get_print_set_volume_only()==0) { obj->print(cs); }
        // if no request for inter-object
        // face intersection detection
        if (cs.get_detect_interobject_intersections()==0)
        {
          // then save space, delete object
          delete obj;
          c.clearObjects();
        }
      }
    }
    c.boundWorld(s);
    // analyze data as set
    cerr << "\n" << "/* ********************** "
         << "SET OF OBJECTS ********************** */\n\n";
    if (cs.get_compute_attributes_only()==0) c.analyzeCumulative();
    // print cumulative surface area, volume, 
    if (cs.get_print_set_volume_only()==0) { c.printCumulative(); }
    else               { cout << c.countVol() << endl; }
    // detect intersecting faces between objects
    if (cs.get_detect_interobject_intersections()==1)
    {
      // partition space
      s.deleteBoxes();
      s.initBoxes(c.countFace());
      c.assignFacesToBoxes(s);
      // analyze and respond
      c.analyzeBatch(s);
      c.printBatch(cs);
    }
  }
}

int distinguishable (double a,double b,double epsilon)
{
  double c;
  c=a-b;
  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;
  if (b<a) return (c>a*epsilon);
  else return (c>b*epsilon);
}

int distinguishable (double a,double b)
{
  double c;
  c=a-b;
  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;
  if (b<a) return (c>a*Controls::instance().get_double_epsilon());
  else     return (c>b*Controls::instance().get_double_epsilon());
}

double dot (double a[3],double b[3])
{
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

bool getPointEdgeDistance (double p[3],Face *f)
{
  // for each face edge
  for (int i=0;i<3;i++) 
  {
    double Ax = f->ptr_vertex(i)->getpN(0);
    double Ay = f->ptr_vertex(i)->getpN(1);
    double Az = f->ptr_vertex(i)->getpN(2);
    double Bx = f->ptr_vertex((i+1)%3)->getpN(0);
    double By = f->ptr_vertex((i+1)%3)->getpN(1);
    double Bz = f->ptr_vertex((i+1)%3)->getpN(2);
    double uDen = (Bx-Ax)*(Bx-Ax)+(By-Ay)*(By-Ay)+(Bz-Az)*(Bz-Az);
    if (uDen) 
    {
      // u = AdotA-AdotB-AdotC+BdotC)/uDen
      double u =((p[0]-Ax)*(Bx-Ax)+(p[1]-Ay)*(By-Ay)+(p[2]-Az)*(Bz-Az))/uDen;
      // no need to check for u ==0 and u ==1, since current 
      // vertex/face plane coincidence was checked previously.
      // Closest point on face edge line to current vertex
      // occurs on face edge between face vertices
      if (u>0 && u<1) 
      {
        double a=Ax+u*(Bx-Ax);
        double b=Ay+u*(By-Ay);
        double c=Az+u*(Bz-Az);
        if (!distinguishable(p[0],a) && 
            !distinguishable(p[1],b) && 
            !distinguishable(p[2],c)) { return true; }
      }
    }
  }
  return false;
}

void biggest (double *x, int& big) 
{
  if ( x[0] < x[1] ) 
  {
    // x[0] < x[1]
    if ( x[0] < x[2] ) 
    {
      // x[0] < x[1]
      // x[0] < x[2]
      if (x[1] < x[2]) 
      {
        // x[0] < x[1]
        // x[1] < x[2]
        // x[0] is smallest
        // x[2] is biggest
        big = 2;
      }
      else 
      {
        // x[0] < x[1]
        // x[2] <= x[1]
        // x[0] is smallest
        // x[1] is biggest or tied with x[2]
        big = 1;
      }
    }
    else 
    {
      // x[0] < x[1]
      // x[2] <= x[0]
      // x[2] is smallest or tied with x[0]
      // x[1] is biggest
      big = 1;
    }
  }
  else 
  {
    // x[1] <= x[0]
    if ( x[0] < x[2] ) 
    {
      // x[1] <= x[0]
      // x[0] < x[2]
      // x1 is smallest or tied with x[0]
      // x2 is biggest
      big = 2;
    }
    else 
    {
      big = 0;
    }
  }
}

void threeValueSort (double p1,double p2,double p3, double &biggest, double &smallest) 
{
  if ( p1 < p2 ) 
  {
    // x[0] < x[1]
    if ( p2 < p3 ) 
    {
      // x[0] < x[1]
      // x[1] < x[2]
      smallest = p1;
      biggest = p3;
    }
    else 
    {
      // x[0] < x[1]
      // x[2] <= x[1]
      if (p1 < p3) 
      {
        // x[0] < x[1]
        // x[0] < x[2]
        // x[2] <= x[1]
        smallest = p1;
        biggest = p2;
      }
      else 
      {
        // x[0] < x[1]
        // x[2] <= x[1]
        // x[2] <= x[0]
        smallest = p3;
        biggest = p2;
      }
    }
  }
  else 
  {
    // x[1] <= x[0]
    if ( p1 < p3 ) 
    {
      // x[1] <= x[0]
      // x[0] < x[2]
      smallest = p2;
      biggest = p3;
    }
    else 
    {
      // x[1] <= x[0]
      // x[2] <= x[0]
      if (p2 < p3) 
      {
        // x[1] <= x[0]
        // x[2] <= x[0]
        // x[1] < x[2]
        smallest = p2;
        biggest = p1;
      }
      else 
      {
        // x[1] <= x[0]
        // x[2] <= x[0]
        // x[2] <= x[1]
        smallest = p3;
        biggest = p1;
      }
    }
  }
}

void checkLineFaceIntersection (Face *f,double lp[2][3],bool &line_flag,
                                bool &poly_flag, bool &poly_edge_flag) 
{
  //lp  = line_points
  // initialize flags
  line_flag=poly_flag=poly_edge_flag=false;
  // get face normal
  double pn[3];
  f->getNormal(pn);
  // use line points, polygon normal and one polygon vertex
  // to compute where line intersects plane of polygon (if not parallel)
  // pI=p1+u*(p2-p1) where u=dot(N,(p3-p1))/dot(N,(p2-p1))
  // pI = intersection of line and plane
  // p3 = any point on plane
  // p2,p1 = end points of line
  // N = plane normal
  // denominator of u
  double den = pn[0]*(lp[1][0]-lp[0][0])+pn[1]*(lp[1][1]-lp[0][1])+pn[2]*(lp[1][2]-lp[0][2]);
  // if line and polygon plane are not parallel
  if (den) 
  {
    //pvc = polygon_vertex_coordinates
    Vertex *v0=f->ptr_vertex(0); 
    // point of intersection
    double u = (pn[0]*(v0->getpN(0)-lp[0][0]) 
                + pn[1]*(v0->getpN(1)-lp[0][1]) 
                + pn[2]*(v0->getpN(2)-lp[0][2]))/den;
    // if polygon cuts through line
    if (u > Controls::instance().get_double_epsilon() && u < (1.0-Controls::instance().get_double_epsilon())) 
    {
      line_flag = true;
    }
    // compute polygon double-area on each of three principal planes
    int tv[2],big;
    double pI[2],I[3];//,d[3];
    double area[3] = { fabs(pn[2]),// xy
      fabs(pn[0]),// yz
      fabs(pn[1])// zx
    };
    biggest(area,big);
    // pI = p1+u*(p2-p1)
    // where lp[1][] is p2 and lp[0][] is p1
    tv[0] = big;
    tv[1] = (big+1)%3;
    pI[0] = (1-u)*lp[0][big] + u*lp[1][big];
    pI[1] = (1-u)*lp[0][(big+1)%3] + u*lp[1][(big+1)%3];
    I[big]=pI[0];
    I[(big+1)%3]=pI[1];
    I[(big+2)%3] = (1-u)*lp[0][(big+2)%3] + u*lp[1][(big+2)%3];
    ////////// is point of intersection on other polygon? //////////
    //pvc = polygon_vertex_coordinates
    Vertex *v1=f->ptr_vertex(1),*v2=f->ptr_vertex(2); 
    // does point of intersection lie on polygon edge?
    if (getPointEdgeDistance(I,f)){poly_edge_flag = true;}
    // if point of intersection is not on polygon edge
    if (poly_edge_flag==false) 
    {
      // compute three determinants
      double det[3];
      det[0] = (v0->getpN(tv[0])-pI[0])*(v1->getpN(tv[1])-pI[1])
              -(v1->getpN(tv[0])-pI[0])*(v0->getpN(tv[1])-pI[1]);
      det[1] = (v1->getpN(tv[0])-pI[0])*(v2->getpN(tv[1])-pI[1])
              -(v2->getpN(tv[0])-pI[0])*(v1->getpN(tv[1])-pI[1]);
      // proceed if determinants are DOUBLE_EPSILON away from zero
      if ((det[0]*det[1])>0)
      {
        det[2]=(v2->getpN(tv[0])-pI[0])*(v0->getpN(tv[1])-pI[1])
              -(v0->getpN(tv[0])-pI[0])*(v2->getpN(tv[1])-pI[1]);
        if ((det[0]*det[2])>0)
        {
          // line intersects polygon plane inside polygon
          poly_flag = true;
        }
      }
    }
  }
}

std::string keyPair (int a,int b,int num_digits)
{
  char str[128],format[32];
  sprintf(format,"%%0%dd%%0%dd",num_digits,num_digits);
  if (a<b) { sprintf(str,format,a,b); }
  else     { sprintf(str,format,b,a); }
  return str;
}

bool edgeMatch (Edge *e,int va,int vb) 
{
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  e->getVertices(v1,v2,o1,o2);
  if ( (v1->getIndex()==va && v2->getIndex()==vb) ||
       (v1->getIndex()==vb && v2->getIndex()==va) ){return true;}
  else {return false;}
}

bool overlap (double a1,double a2,double b1,double b2)
{
  // if no overlap of f1 bb (a1,a2) and f2 bb (b1,b2)
  if (a2<b1 || b2<a1) { return false; }
  else                { return true;  }
}

bool faceBBsOverlap (double bb[6],Face *f)
{
  double bb2[6];

  double *p1=f->ptr_vertex(0)->getpN_ptr();
  double *p2=f->ptr_vertex(1)->getpN_ptr();
  double *p3=f->ptr_vertex(2)->getpN_ptr();

  threeValueSort(*p1,*p2,*p3,bb2[1],bb2[0]);
  // if overlap in x
  if (overlap(bb[0],bb[1],bb2[0],bb2[1]))
  {
    p1++;p2++;p3++;
    threeValueSort(*p1,*p2,*p3,bb2[3],bb2[2]);
    // if overlap in y
    if (overlap(bb[2],bb[3],bb2[2],bb2[3]))
    {
      p1++;p2++;p3++;
      threeValueSort(*p1,*p2,*p3,bb2[5],bb2[4]);
      // if overlap in z
      if (overlap(bb[4],bb[5],bb2[4],bb2[5]))
      {
        return true;
      } else {return false;}
    } else {return false;}
  } else {return false;}
}

