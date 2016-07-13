#include <cassert>
#include <cmath>
#include <iostream>
#include <stdlib.h>

#include "contour.h"
#include "controls.h"
#include "histogram.h"
#include "point.h"
#include "weights.h"

Contour::Contour (Contour const & rhs)
:section(rhs.section),name(rhs.name),
  raw_points(rhs.raw_points),spline_samples(rhs.spline_samples),deviations(rhs.deviations)
{
}

Contour & Contour::operator = (const Contour& rhs)
{
  std::cout << "Assignment operator prohibited on instances of contour class.\n";
  std::cout << "Object " << rhs.name << std::endl;
  exit(1);
}

Contour::Contour (char *str, int sec)
:section(sec),name(),raw_points(),
  spline_samples(),deviations()
{
  char *ptr = str;
  while(ptr!=NULL){
    ptr=strpbrk(str,"#()");
    if(ptr!=NULL){
      *ptr='_';
    }
  }
  name = str;
  assert(name!="");
}

void Contour::interpPoints (p_l_iterator & lhs,p_l_iterator & rhs)
{
  Controls & cs (Controls::instance()); 
  // compute distance between points
  double distx = rhs->getX()-lhs->getX();
  double disty = rhs->getY()-lhs->getY();
  double dist = sqrt(distx*distx+disty*disty);
  int num_interpolated_points = static_cast<int>(dist/cs.getLinearThreshold());
  //if (num_interpolated_points > 1000)
  //{
  //  cout << "name = " << name
  //        << ", num_interpolated_points = " << num_interpolated_points << endl;
  //} 
  //cout << "lhs (" << lhs->getX() << " " << lhs->getY() << ")"
  //     << ", rhs (" << rhs->getX() << " " << rhs->getY() << ")"
  //     << ", num_interpolated_points = " << num_interpolated_points << endl;
  if (num_interpolated_points)
  {
    // linearly interpolate num_interpolated_points evenly spaced points
    double num_segments = num_interpolated_points+1;
    int i = num_interpolated_points;
    while (i)
    {
      // insert point
      double x = lhs->getX()+(1.0-static_cast<double>(i)/num_segments)*distx;
      double y = lhs->getY()+(1.0-static_cast<double>(i)/num_segments)*disty;
      raw_points.insert(rhs,Point(x,y,lhs->getZ()));
      i--;
    }
  }
}

void Contour::removeDuplicates (void)
{
  p_l_iterator i=raw_points.begin();
  p_l_iterator j=i;
  j++;
  while (j!=raw_points.end())
  {
    if (pointsAreIdentical(*i,*j))
    {
      raw_points.erase(i);
      j=i;
      j++;
    }
    else
    {
      i++;
      j++;
    }
  }
  // compare first and last link in list
  i = raw_points.begin();
  j = raw_points.end();
  j--;
  if (pointsAreIdentical(*i,*j))
  {
    raw_points.erase(i);
  }
}

void Contour::linearlyInterp (void)
{
  p_l_iterator i=raw_points.begin();
  p_l_iterator j=i;
  j++;
  while (j!=raw_points.end())
  {
    interpPoints(i,j);
    i++;j++;
  }
  // compare first and last link in list
  i=raw_points.begin();
  j=raw_points.end();
  j--;
  interpPoints(i,j);
}

void Contour::fitSplines (Histogram & h)
{
  Controls & cs (Controls::instance()); 
  Weights w;
  int i=0;
  int current_section = -1;
  if (getSection()!=current_section)
  {
    current_section = getSection();
  }
  if (raw_points.size()<3)
  {
    printf("contour has less than 3 points and ");
    printf("was skipped: contour %s, section %d,",getName(),getSection());
    printf(" num_points %d\n",static_cast<int>(raw_points.size()));
  }
  else
  {
    ///// load arrays of weights /////
    w.allocateX(4*raw_points.size());
    w.allocateY(4*raw_points.size());
    w.loadWeightArrays(raw_points.begin(),raw_points.end(),raw_points.size());
    // W ARRAY IS CONSTRUCTED WITH FIRST WEIGHT CORRESPONDING TO FIRST RAW POINT IN INPUT FILE.
    ///// compute splines for current contour /////
    double dt = 1.0/static_cast<double>(cs.getSplineSamplesPerPoint());
    std::vector<SplinePoint> sp(raw_points.size()*cs.getSplineSamplesPerPoint(),SplinePoint());
    w.computeSplines(sp,dt,raw_points.size());
    ///// sample splines /////
    sampleSplines(sp,dt,raw_points.size(),i);
    ///// compute deviation distance between raw points and splines /////
    computeDeviation(sp,raw_points.size(),cs.getSplineSamplesPerPoint());
    ///// check deviations /////
    int cnt = cs.getMaxDeviationAdjustments();
    if (cs.getDeviationThreshold())
    {
      while (w.checkDeviation(this,raw_points.size()) && cnt)
      {
        // recompute splines
        w.computeSplines(sp,dt,raw_points.size());
        // clear inter_points in contour
        clearSpline();
        sampleSplines(sp,dt,raw_points.size(),i);
        computeDeviation(sp,raw_points.size(),cs.getSplineSamplesPerPoint());
        cnt--;
      }
      if (!cnt)
      {
        cout << "Deviation adjustment of contour (" << name << ") ended by threshold ("
              << cs.getMaxDeviationAdjustments() << ")\n";
      }
    }
    ///// update min and max deviation distances /////
    for (uint ii=0;ii<raw_points.size();ii++)
    {
      // DEBUG
      //cout << "\n";cout.flush();
      //cout << "Contour::fitSplines: getDeviationThreshold = " << cs.getDeviationThreshold() << "\n";cout.flush();
      //cout << "Contour::fitSplines: ii = " << ii << "\n";cout.flush();
      //cout << "Contour::fitSplines: raw_points.size = " << raw_points.size() << "\n";cout.flush();
      //cout << "Contour::fitSplines: getDeviation(ii) = " << getDeviation(ii) << "\n";cout.flush();
      // DEBUG
      h.update(getDeviation(ii));
    }
    if (cs.getDiag())
    {
      printDiagnostics(sp,raw_points.size()*cs.getSplineSamplesPerPoint(),cs.getOutputDir(),i);
    }
    // clean up	
    w.freeX();
    w.freeY();
  }
  //m++;
  i++;
}

void Contour::sampleSplines (v_sp & sp,double dt,int limit,int tag)
{
  Controls & cs (Controls::instance()); 
  ///// sample splines /////
  // dmin = minimum spline sampling distance
  // dmax = maximum spline sampling distance
  // decrease tau to increase steepness
  // inflection point of sampling function = tau*rI
  // decrease rI to sample finer
  //int k=0;
  double intfac_array[3];
  double r_array[2];
  const double vmax=cs.getDmax()/cs.getT();
  double t_accum=0.0;
  std::vector<double> rvec;
  //int    kvec[limit*cs.getSplineSamplesPerPoint()/2];
  //std::vector<int > kvec;
  //double vvec[limit*cs.getSplineSamplesPerPoint()/2];
  std::vector<double> vvec;
  //double tvec[limit*cs.getSplineSamplesPerPoint()/2];
  std::vector<double> tvec;

  for (int i=0;i<limit;i++)
  {
    for (int j=0;j<cs.getSplineSamplesPerPoint();j++)
    {
      double inc=static_cast<double>(i)+static_cast<double>(j)*dt;
      int myswitch = j%2;
      if (inc && !myswitch)
      {
        intfac_array[0] = sp[(int)((inc-2.0*dt)/dt)].intfac;
        intfac_array[1] = sp[(int)((inc-dt)/dt)].intfac;
        intfac_array[2] = sp[(int)(inc/dt)].intfac;
        r_array[0]      = sp[(int)((inc-2.0*dt)/dt)].r;
        r_array[1]      = sp[(int)(inc/dt)].r;
        double xval     = sp[(int)(inc/dt)].x;
        double yval     = sp[(int)(inc/dt)].y;
        // increment along spline length
        double dl = 2*dt/3.0*(intfac_array[0]+4.0*intfac_array[1]+intfac_array[2]);
        // mean radius of curvature
        double r_mean = (r_array[0]+r_array[1])/2.0;
        // mean velocity
        double v_mean = sqrt(cs.getAmax()*r_mean);
        if (v_mean>vmax) {v_mean=vmax;}
        // time increment
        double delt = dl/v_mean;
        // accumulated time
        t_accum += delt;
        // sample data
        if (t_accum >= cs.getT())
        {
          // add interpolated point to contour
          spline_samples.push_back(Point(xval,yval,getSection()*cs.getSectionThickness()));
          // clear variables
          t_accum = 0.0;
        }
        if (cs.getDiag())
       // if (name=="Post" && section==488)
        {
          //rvec[k]=r_mean;
          rvec.push_back(r_mean);
          //vvec[k]=v_mean;
          vvec.push_back(v_mean);
          //tvec[k]=delt;
          tvec.push_back(delt);
          //kvec[k]=k;
          //kvec.push_back(k);
          //k++;
        }
      }
    }
  }

  // store number of sampled points
//  setNumInterpPoints(num_sampled);

  // diagnostics
  if (cs.getDiag())
  //if (name=="Post" && section==488)
  {
    char filename[256],line[2048];
    FILE *F;
    ///// print radius of curvature /////
    sprintf(filename,"%s%s_%i_%i.rad",cs.getOutputDir(),getName(),getSection(),tag);
    F = fopen(filename,"w");
    if (!F) { printf("Could not open output file %s\n",filename); return; }
    // for each point
    for(int i=0;i<limit*cs.getSplineSamplesPerPoint()/2-1;i++){
      //sprintf(line,"%i %.15g\n",kvec[i],rvec[i]);
      sprintf(line,"%i %.15g\n",i,rvec[i]);
      fputs(line,F);
    }
    fclose(F);
    ///// print velocity /////
    sprintf(filename,"%s%s_%i_%i.vel",cs.getOutputDir(),getName(),getSection(),tag);
    F = fopen(filename,"w");
    if (!F) { printf("Could not open output file %s\n",filename); return; }
    // for each point
    for(int i=0;i<limit*cs.getSplineSamplesPerPoint()/2-1;i++){
      //sprintf(line,"%i %.15g\n",kvec[i],vvec[i]);
      sprintf(line,"%i %.15g\n",i,vvec[i]);
      fputs(line,F);
    }
    fclose(F);
    ///// print incremental time /////
    sprintf(filename,"%s%s_%i_%i.tim",cs.getOutputDir(),getName(),getSection(),tag);
    F = fopen(filename,"w");
    if (!F) { printf("Could not open output file %s\n",filename); return; }
    // for each point
    for(int i=0;i<limit*cs.getSplineSamplesPerPoint()/2-1;i++){
      //sprintf(line,"%i %.15g\n",kvec[i],tvec[i]);
      sprintf(line,"%i %.15g\n",i,tvec[i]);
      fputs(line,F);
    }
    fclose(F);
  }
}

//void Contour::computeDeviation (SplinePoint* sp,int count,int num)
void Contour::computeDeviation (v_sp & sp,const int count,const int num)
{
  ///// compute deviation distance between raw points and splines /////
  allocateDeviations(count);
  double diffx,diffy,dist0,distneg,distneg0,distpos,distpos0;
  int j=0,i=0,m,n;
  for (p_l_iterator p = raw_points.begin();p!=raw_points.end();p++)
  {
    //cout << "\nname = " << name;cout.flush();
    //cout << ", p->getX = " << p->getX();cout.flush();
    //cout << ", p->getY = " << p->getY();cout.flush();
    //cout << ", count = " << count;cout.flush();
    //cout << ", num = " << num;cout.flush();
    //cout << ", j = " << j;cout.flush();
    //cout << ", sp.size = " << sp.size();cout.flush();
    //cout << ", sp.[j].x = " << sp[j].x;cout.flush();
    //cout << ", sp.[j].y = " << sp[j].y;cout.flush();
    assert(count>0);
    assert(j<static_cast<int>(sp.size()));
    //assert(j<count);
    //cout << endl;cout.flush();
    // compare raw point to first SplinePoint
    diffx = p->getX()-sp[j].x;
    diffy = p->getY()-sp[j].y;
    dist0 = sqrt(diffx*diffx+diffy*diffy);
    // check negative dir
    // search SplinePoints from largest to smallest index
    // stop when distance to raw point ceases to decrease
    m=j;
    distneg=dist0;
    do {
      distneg0=distneg;
      if(!m){m=num*count-1;}
      else {m--;}
      diffx = p->getX()-sp[m].x;
      diffy = p->getY()-sp[m].y;
      distneg = sqrt(diffx*diffx+diffy*diffy);
    } while (distneg<distneg0);
    // check positive dir
    // search SplinePoints from smallest to largest index
    // stop when distance to raw point ceases to decrease
    n=j;
    distpos=dist0;
    do {
      distpos0=distpos;
      if(n==num*count){n=0;}
      else {n++;}
      diffx = p->getX()-sp[n].x;
      diffy = p->getY()-sp[n].y;
      distpos = sqrt(diffx*diffx+diffy*diffy);
    } while (distpos<distpos0);

    //if (distpos0==distneg0)
    //{
    //  cout << "\nContour::computeDeviation: dist0    = " << dist0 << endl;
    //  cout << "Contour::computeDeviation: "
    //        << "distneg = " << distneg
    //        << ", distneg0 = " << distneg0 << endl;
    //  cout << "Contour::computeDeviation: "
    //        << "distpos = " << distpos
    //        << ", distpos0 = " << distpos0 << endl;
    //  cout << "Contour::computeDeviation: " 
    //        << "m = " << m << " of " << num*count-1 << endl;
    //  cout << "Contour::computeDeviation: " 
    //        << "n = " << n << " of " << "0" << endl;
    //  assert(distpos0!=distneg0);
    //}

    // if original distance was lowest
    if (dist0<=distneg0 && dist0<=distpos0)
    {
      setDeviation(i,dist0);
    }
    else if(distneg0<=dist0 && distneg0<=distpos0)
    {
      setDeviation(i,distneg0);
    }
    else if (distpos0<=dist0 && distpos0<=distneg0)
    {
      setDeviation(i,distpos0);
    }
    else
    { 
      printf("\n\nweird error!\n");
      printf("Contour %s, section %d, #raw points %d, num %d\n",getName(),getSection(),count,num);
      printf("rawx %.15g, rawy %.15g\n",p->getX(),p->getY());
      printf("j %i, dist0 %.15g, dist0x %.15g, dist0y %.15g\n",j,dist0,sp[j].x,sp[j].y);
      printf("m %i, distneg0 %.15g, distneg %.15g," ,m,distneg0,distneg);
      printf(" distnegx %.15g, distnegy %.15g\n" ,sp[m].x,sp[m].y);
      printf("n %i, distpos0 %.15g, distpos %.15g," ,n,distpos0,distpos);
      printf(" distposx %.15g, distposy %.15g\n" ,sp[n].x,sp[n].y);
      cout << "case1: dist0(" << dist0
            << ")<distneg0(" << distneg0
            << ") && dist0(" << dist0
            << ")<distpos0(" << distpos0
            << ")\n";
      cout << "case2: distneg0(" << distneg0
            << ")<dist0(" << dist0
            << ") && distneg0(" << distneg0
            << ")<distpos0(" << distpos0
            << ")\n";
      cout << "case3: distpos0(" << distpos0
            << ")<dist0(" << dist0
            << ") && distpos0(" << distpos0
            << ")<distneg0(" << distneg0
            << ")\n";
    }
    j+=num;
    if (j==count) { break; }
    i++;
  }
}

//void Contour::printDiagnostics (SplinePoint* sp,int num,char const * const outdir,int tag)
void Contour::printDiagnostics (v_sp & sp,int num,char const * const outdir,int tag)
{
  int j;
  char filename[256],line[2048];
  FILE *F;
  ///// print raw points /////
  sprintf(filename,"%s%s_%i_%i.raw",outdir,getName(),getSection(),tag);
  F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  // for each point
  for (p_l_iterator qq=raw_points.begin();qq!=raw_points.end();qq++)
  {
    sprintf(line,"%.15g %.15g\n",qq->getX(),qq->getY());
    fputs(line,F);
  }
  fclose(F);
  ///// print spline points /////
  sprintf(filename,"%s%s_%i_%i.spline",outdir,getName(),getSection(),tag);
  F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  // for each spline point
  for(j=0;j<num;j++){
    sprintf(line,"%.15g %.15g\n",sp[j].x,sp[j].y);
    fputs(line,F);
  }
  fclose(F);
  ///// print interpolated points /////
  sprintf(filename,"%s%s_%i_%i.interp",outdir,getName(),getSection(),tag);
  F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  // for each point
  for (c_p_iterator i = spline_samples.begin();i!=spline_samples.end();i++)
  {
    sprintf(line,"%.15g %.15g\n",i->getX(),i->getY());
    fputs(line,F);
  }
  fclose(F);
  ///// print grace script /////
  sprintf(filename,"%sbfile_%s_%i_%i",outdir,getName(),getSection(),tag);
  F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  sprintf(line,"#Obligatory descriptive comment\n");
  fputs(line,F);
  sprintf(line,"READ xy \"%s_%i_%i.raw\"\n",getName(),getSection(),tag);
  fputs(line,F);
  sprintf(line,"READ xy \"%s_%i_%i.spline\"\n",getName(),getSection(),tag);
  fputs(line,F);
  sprintf(line,"READ xy \"%s_%i_%i.interp\"\n",getName(),getSection(),tag);
  fputs(line,F);
  sprintf(line,"legend on\ns0 legend ");
  sprintf(line,"%s\"Raw\"\ns1 legend \"Spline\"\ns2 legend \"Interpolated\"\n",line);
  fputs(line,F);
  sprintf(line,"s0 symbol 1\ns0 symbol color 2\ns0 symbol fill color 2\n");
  sprintf(line,"%ss0 line type 0\ns0 line color 2\ns0 errorbar color 2\n",line);
  fputs(line,F);
  sprintf(line,"s1 symbol color 1\ns1 symbol fill color 1\n");
  sprintf(line,"%ss1 line color 1\ns1 errorbar color 1\n",line);
  fputs(line,F);
  sprintf(line,"s2 symbol 2\ns2 symbol color 4\ns2 symbol fill color 4\n");
  sprintf(line,"%ss2 line type 0\ns2 line color 4\ns2 errorbar color 4\n",line);
  fputs(line,F);
  fclose(F);
  ///// print shell script /////
  sprintf(filename,"%s%s_%i_%i.sh",outdir,getName(),getSection(),tag);
  F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  sprintf(line,"#!/bin/bash\n");
  fputs(line,F);
  sprintf(line,"xmgrace -noask -nosafe -batch bfile_%s_%i_%i\n",getName(),getSection(),tag);
  fputs(line,F);
  fclose(F);
}

void Contour::computeHistogram (Histogram & h)
{
  double foo=0;
  int m=0;
  h.setMean(h.getSum()/h.getNum());
  // for each contour
  if (raw_points.size()<3)
  {
    printf("contour has less than 3 points and was skipped:");
    printf(" contour %s, section %d, num_points %d\n",getName(),getSection(),static_cast<int>(raw_points.size()));
  }
  else
  {
    // for each raw point in contour
    for (uint i=0;i<raw_points.size();i++)
    {
      // compute stddev scratch work
      foo+=(getDeviation(i)-h.getMean())*(getDeviation(i)-h.getMean());
      // bin deviation
      h.incrementCount(static_cast<int>(getDeviation(i)/(h.getMax()/14)));
    }
  }
  m++;
  h.setStddev(sqrt(foo/(h.getNum()-1)));
}

