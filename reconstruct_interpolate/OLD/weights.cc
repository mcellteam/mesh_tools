#include <cmath>
#include <string>

#include "contour.h"
#include "controls.h"
#include "point.h"
#include "weights.h"

Weights::Weights (void)
:bx(),by()
{
}

bool Weights::checkDeviation (Contour const * c,
                                  int limit)
{
  Controls & cs(Controls::instance());
  ///// check deviation at each point in contour /////
  bool flag=false;
  // for each deviation
  for (int i=0;i<limit;i++)
  {
    // if deviation exceeds threshold
    if (c->getDeviation(i)>cs.getDeviationThreshold())
    {
      flag = true;
      // edit weight array
      if (!i)
      {
        bx[4*(limit-1)+3] = bx[4*(limit-1)+2];
        by[4*(limit-1)+3] = by[4*(limit-1)+2];
        bx[4*i+0] = bx[4*i+1];
        by[4*i+0] = by[4*i+1];
      }
      else
      {
        bx[4*(i-1)+3] = bx[4*(i-1)+2];
        by[4*(i-1)+3] = by[4*(i-1)+2];
        bx[4*i+0] = bx[4*i+1];
        by[4*i+0] = by[4*i+1];
      }
    }
  }
  return flag;
}

//SplinePoint* Weights::computeSplines (SplinePoint* sp,double dt,int limit)
void Weights::computeSplines (v_sp & sp,double dt,int limit)
{
  Controls & cs (Controls::instance()); 
  //cout << "Weights::computeSplines : limit = " << limit
  //      << ", getNum = " << cs.getNum()
  //      << ", begin sp.size = " << sp.size() << "\n";
  // for each point in contour
  for (int i=0;i<limit;i++)
  {
    // for each parameter increment
    for (int j=0;j<cs.getSplineSamplesPerPoint();j++)
    {
      int index = i*cs.getSplineSamplesPerPoint()+j;
      //cout << "Weights::computeSplines : index = " << index
      //      << "\n";
      double inc1 = j*dt;
      double inc2 = inc1*inc1;
      double inc3 = inc2*inc1;
      sp[index].t   = i+inc1;
      // DEBUG
      //cout << "\nWeights::computeSplines: i = " << i
      //      << ", limit = " << limit
      //      << ", 4*i+3 = " << 4*i+3
      //      << ", bx.size = " << bx.size() << endl;
      // DEBUG
      assert((4*i+3)<static_cast<int>(bx.size()));
      assert((4*i+3)<static_cast<int>(by.size()));
      sp[index].x   = 1.0/6.0*(bx[4*i+0]*(-inc3+3.0*inc2-3.0*inc1+1.0) 
                               + bx[4*i+1]*(3.0*inc3-6.0*inc2+4.0) 
                               + bx[4*i+2]*(-3.0*inc3+3.0*inc2+3.0*inc1+1.0) 
                               + bx[4*i+3]*inc3);
      sp[index].y   = 1.0/6.0*(by[4*i+0]*(-inc3+3.0*inc2-3.0*inc1+1.0) 
                               + by[4*i+1]*(3.0*inc3-6.0*inc2+4.0) 
                               + by[4*i+2]*(-3.0*inc3+3.0*inc2+3.0*inc1+1.0) 
                               + by[4*i+3]*inc3);
      double xdot    = 1.0/6.0*(bx[4*i+0]*(-3.0*inc2+6.0*inc1-3.0)
                                + bx[4*i+1]*(9.0*inc2-12.0*inc1)
                                + bx[4*i+2]*(-9.0*inc2+6.0*inc1+3.0)
                                + bx[4*i+3]*3.0*inc2);
      double ydot    = 1.0/6.0*(by[4*i+0]*(-3.0*inc2+6.0*inc1-3.0)
                                + by[4*i+1]*(9.0*inc2-12.0*inc1)
                                + by[4*i+2]*(-9.0*inc2+6.0*inc1+3.0)
                                + by[4*i+3]*3.0*inc2);
      double xdotdot = 1.0/6.0*(bx[4*i+0]*(-6.0*inc1+6.0)
                                + bx[4*i+1]*(18.0*inc1-12.0)
                                + bx[4*i+2]*(-18.0*inc1+6.0)
                                + bx[4*i+3]*6.0*inc1);
      double ydotdot = 1.0/6.0*(by[4*i+0]*(-6.0*inc1+6.0)
                                + by[4*i+1]*(18.0*inc1-12.0)
                                + by[4*i+2]*(-18.0*inc1+6.0)
                                + by[4*i+3]*6.0*inc1);
      double den     = fabs(xdot*ydotdot-ydot*xdotdot);
      double num_part = sqrt(xdot*xdot+ydot*ydot);
      sp[index].intfac = num_part;
      if (den) { sp[index].r = num_part*num_part*num_part/den; }
      else { sp[index].r = cs.getMaxRad(); }
      if (sp[index].r>cs.getMaxRad()) { sp[index].r=cs.getMaxRad(); }
    }
  }
  //cout << "Weights::computeSplines : end sp.size = " << sp.size()
  //      << "\n";
  // for each point in contour
  //return sp;
}

void Weights::loadWeightArrays (p_l_iterator raw_points_begin,
                                p_l_iterator raw_points_end,
                                int const & limit)
{
  ///// load arrays of weights /////
  // END OF RAW POINTS ARRAY WAS PASSED AS P 
  // AND LIST IS TRAVERSED IN 'PREVIOUS' DIRECTION
  // spline 0
  p_l_iterator p = raw_points_end;
  p--;
  bx.push_back(p->getX());
  by.push_back(p->getY());
  //i++;
  p = raw_points_begin;
  bx.push_back(p->getX());
  by.push_back(p->getY());
  //i++;
  p++;
  bx.push_back(p->getX());
  by.push_back(p->getY());
  //i++;
  p++;
  bx.push_back(p->getX());
  by.push_back(p->getY());
  //i++;
  // splines 1 through m-2
  p = raw_points_begin;
  p_l_iterator q = p;q++;
  p_l_iterator r = q;r++;
  p_l_iterator s = r;s++;
  int myadded = 0;
  for (int j=1;j<(limit-2);j++)
  { // note '-2' because m = contour_array[i]-1
    bx.push_back(p->getX());
    by.push_back(p->getY());
    myadded++;
    bx.push_back(q->getX());
    by.push_back(q->getY());
    myadded++;
    bx.push_back(r->getX());
    by.push_back(r->getY());
    myadded++;
    bx.push_back(s->getX());
    by.push_back(s->getY());
    myadded++;
    p++;
    q++;
    r++;
    s++;
  }
  // spline m-1
  int k=0;
  for (p=raw_points_begin;p!=raw_points_end;p++)
  {
    if (k==(limit-3)) { break; }
    k++;
  }
  bx.push_back(p->getX());
  by.push_back(p->getY());
  //i++;
  p++;
  bx.push_back(p->getX());
  by.push_back(p->getY());
  //i++;
  p++;
  bx.push_back(p->getX());
  by.push_back(p->getY());
  //i++;
  p=raw_points_begin;
  bx.push_back(p->getX());
  by.push_back(p->getY());
  //i++;
  // spline m
  k=0;
  for (p=raw_points_begin;p!=raw_points_end;p++)
  {
    if (k==(limit-2)) { break; }
    k++;
  }
  bx.push_back(p->getX());
  by.push_back(p->getY());
  //i++;
  p++;
  bx.push_back(p->getX());
  by.push_back(p->getY());
  //i++;
  p=raw_points_begin;
  bx.push_back(p->getX());
  by.push_back(p->getY());
  //i++;
  p++;
  bx.push_back(p->getX());
  by.push_back(p->getY());
  assert(bx.size()>0);
  assert(by.size()>0);
}

