#include <cassert>
#include <iostream>
#include <stdlib.h>

#include "container.h"
#include "contour.h"
#include "controls.h"
#include "point.h"

using std::cout;
using std::endl;

Container::Container (void)
:o(),points_per_contour(),num_files_read(0)
{
  o.reserve(Controls::instance().getNumReserve());
}

int Container::getMatchingObject (char const * myname)
{
  int a = 0;
  for (std::vector<Object>::iterator i = o.begin();i!=o.end();i++)
  {
    if( !strcmp(i->getName(),myname) )
    {
      return a;
    }
    a++;
  }
  return -1;
}

int Container::newObject (char * myname,int i)
{
  assert(o.size()<o.capacity());
  o.push_back(Object(myname,i));
  assert(o.size()>0);
  assert(o.size()>0);
  return o.size()-1;
}

void Container::addContour2Object (char * myname,int section)
{
  // has object been created with same name?
  int index = getMatchingObject(myname);
  if (index < 0)
  {
    // no, create new object
    index = newObject(const_cast<char*>(myname),section);
  }
  // add contour to pointer object
  assert(index>=0);
  //assert(static_cast<uint>(index)<o.size());
  assert(index<static_cast<int>(o.size()));
  o[index].addContour(Contour(myname,section));
  // update pointer object min and max
  if (section<o[index].getMinSection()){o[index].setMinSection(section);}
  if (section>o[index].getMaxSection()){o[index].setMaxSection(section);}
}

void Container::initTransform (double *t)
{
  t[0]=0.0; t[1]=1.0; t[2]=0.0; t[3]=0.0; t[4]=0.0; t[5]=0.0;
  t[6]=0.0; t[7]=0.0; t[8]=1.0; t[9]=0.0; t[10]=0.0; t[11]=0.0;
}

void Container::setTransform (double *t,char *p)
{
  char val[80],*eptr;
  int i;
  // grab '1'
  while (strchr(" \t,",*p)!=NULL) { p++; }
  i=0;
  while (strchr("0123456789+-eE.",*p)!=NULL){val[i++]=*p++;}
  val[i]=0;
  t[0]=strtod(val,&eptr);
  if (val==eptr) {
    t[0]=0.0; t[1]=1.0; t[2]=0.0; t[3]=0.0; t[4]=0.0; t[5]=0.0;
    printf("Error in reading '1' coefficient\n"); printf("str =%s\n",p);
    return;
  }
  // grab 'x'
  while (strchr(" \t,",*p)!=NULL) { p++; }
  i=0;
  while (strchr("0123456789+-eE.",*p)!=NULL){val[i++]=*p++;}
  val[i]=0;
  t[1]=strtod(val,&eptr);
  if (val==eptr) {
    t[0]=0.0; t[1]=1.0; t[2]=0.0; t[3]=0.0; t[4]=0.0; t[5]=0.0;
    printf("Error in reading 'x' coefficient\n"); printf("str =%s\n",p);
    return;
  }
  // grab 'y'
  while (strchr(" \t,",*p)!=NULL) { p++; }
  i=0;
  while (strchr("0123456789+-eE.",*p)!=NULL){val[i++]=*p++;}
  val[i]=0;
  t[2]=strtod(val,&eptr);
  if (val==eptr) {
    t[0]=0.0; t[1]=1.0; t[2]=0.0; t[3]=0.0; t[4]=0.0; t[5]=0.0;
    printf("Error in reading 'y' coefficient\n"); printf("str =%s\n",p);
    return;
  }
  // grab 'xy'
  while (strchr(" \t,",*p)!=NULL) { p++; }
  i=0;
  while (strchr("0123456789+-eE.",*p)!=NULL){val[i++]=*p++;}
  val[i]=0;
  t[3]=strtod(val,&eptr);
  if (val==eptr) {
    t[0]=0.0; t[1]=1.0; t[2]=0.0; t[3]=0.0; t[4]=0.0; t[5]=0.0;
    printf("Error in reading 'xy' coefficient\n"); printf("str =%s\n",p);
    return;
  }
  // grab 'x*x'
  while (strchr(" \t,",*p)!=NULL) { p++; }
  i=0;
  while (strchr("0123456789+-eE.",*p)!=NULL){val[i++]=*p++;}
  val[i]=0;
  t[4]=strtod(val,&eptr);
  if (val==eptr) {
    t[0]=0.0; t[1]=1.0; t[2]=0.0; t[3]=0.0; t[4]=0.0; t[5]=0.0;
    printf("Error in reading 'x*x' coefficient\n"); printf("str =%s\n",p);
    return;
  }
  // grab 'y*y'	
  while (strchr(" \t,",*p)!=NULL) { p++; }
  i=0;
  while (strchr("0123456789+-eE.",*p)!=NULL){val[i++]=*p++;}
  val[i]=0;
  t[5]=strtod(val,&eptr);
  if (val==eptr) {
    t[0]=0.0; t[1]=1.0; t[2]=0.0; t[3]=0.0; t[4]=0.0; t[5]=0.0;
    printf("Error in reading 'y*y' coefficient\n"); printf("str =%s\n",p);
    return;
  }
}

c_l_iterator Container::getLastContourFromObject (char * name)
{
  int index = getMatchingObject(name);
  if (index < 0)
  {
    cout << "Error. Matching object not found with name = " << name << endl;
    assert(index > 0);
  }
  return o[index].getLastContour();
}

void Container::getContours (void)
{
  Controls & cs(Controls::instance());
  char name[80];
  char filename[1024];
  sprintf(filename,"%s%s",cs.getInputDir(),cs.getPrefix());
  // for each reconstruct input file
  int num_contours=0;
  for (int i=cs.getMinSection();i<cs.getMaxSection()+1;i++)
  {
    // open file
    char infile[1024];
    sprintf(infile,"%s.%d",filename,i);
    FILE *F = fopen(infile,"r");
    if (!F) { printf("Could not open input file %s\n",infile); exit(1);}
    bool contour_flag = false;
    // initialize Transform
    double transform[12];
    initTransform(transform);
    // for every line in file
    char line[2048];
    for (char *str=fgets(line,2048,F);str!=NULL;str=fgets(line,2048,F))
    {
      if (strstr(str,"Transform dim")!=NULL)
      {
        // get next line
        str=fgets(line,2048,F);
        if(str==NULL)
        {
          printf("Nothing after Transform Dim.");
          printf(" Reconstruct file may be corrupted: %s.\n",infile);
          exit(1);
        }
        // get xcoeff
        if (strstr(str,"xcoef=")==NULL)
        {
          printf("No xcoef. Reconstruct file may be corrupted.: %s.\n",infile);
          exit(1);
        }
        char *coef = strstr(str,"xcoef=");
        coef += 8; // advance pointer to start of coefficients
        // 8, because 'xcoeff="' is full string
        setTransform(transform,coef);
        // get next line
        str=fgets(line,2048,F);
        if(str==NULL)
        {
          printf("Nothing after xcoef.");
          printf(" Reconstruct file may be corrupted: %s.\n",infile);
          exit(1);
        }
        // get ycoeff
        if (strstr(str,"ycoef=")==NULL)
        {
          printf("No ycoef. Reconstruct file may be corrupted: %s.\n",infile);
          exit(1);
        }
        coef = strstr(str,"ycoef=");
        coef += 8; // advance pointer to start of coefficients
        // 8, because 'ycoeff="' is full string
        setTransform(transform+6,coef);
      }
      // if start of contour
      else if (strstr(str,"<Contour")!=NULL)
      {
        // find name
        char *ptr = strstr(str,"name=");
        ptr += 6; // advance pointer to start of name
        int j=0;
        while (strchr("\"",*ptr)==NULL){name[j++]=*ptr++;}
        if (j>79) {printf("Error. Ran off end of array.\n");exit(1);}
        name[j]=0;
        // convert bad characters to underscore
        char *myptr = &name[0];
        while(myptr!=NULL){
          myptr=strpbrk(name,"#()");
          if(myptr!=NULL){
            *myptr='_';
          }
        }
        // skip ignored_contours
        if (cs.contourIsIgnored(name)==true)
        {
          printf("contour %s ignored.\n",name);
          continue;
        }
        addContour2Object(name,i);
        num_contours++;
        // set contour flag
        contour_flag = true;
      } 
      else if (strstr(str,"/>")!=NULL)
      {
        contour_flag = false;
      }
      else if (contour_flag)
      {
        // add point to contour
        c_l_iterator c_it = getLastContourFromObject(name); 
        c_it->addRawPoint(Point(str,i,cs.getSectionThickness(),&transform[0]));
      }
    }
    fclose(F);
    num_files_read++;
  }
  printf("Read %d contours from %d input files.\n",num_contours,num_files_read);
}

void Container::purgeBadContours (void)
{
  // for each object
  for (o_iterator i=o.begin();i!=o.end();i++)
  {
    i->purgeBadContours();
  }
}

void Container::removeDuplicates (void)
{
  // for each object
  for (o_iterator i=o.begin();i!=o.end();i++)
  {
    i->removeDuplicates();
  }
}

void Container::interpolateRawPoints (void)
{
  for (o_iterator i = o.begin();i!=o.end();i++)
  {
    i->interpolateRawPoints();
  }
}

void Container::countPtsPerContour (void)
{
  points_per_contour.clear();
  // for each object
  for (o_iterator i = o.begin();i!=o.end();i++)
  {
    i->countPtsPerContour(points_per_contour);
  }
}

void Container::fitSplines (Histogram & h)
{
  for (o_iterator i = o.begin();i!=o.end();i++)
  {
    //i->fitSplines(h,points_per_contour);
    i->fitSplines(h);
  }
}

void Container::computeHistogram (Histogram & h)
{
  for (o_iterator i = o.begin();i!=o.end();i++)
  {
    i->computeHistogram(h);
  }
}

void Container::clearPtsFiles (void)
{
  // clear any existing pts files
  Controls & cs(Controls::instance());
  char filename[256],line[2048];
  FILE *F;
//  for (c_o_iterator i = firstObject();i!=onePastLastObject();i++)
//  {
//    for (c_c_l_iterator j = i->firstContour();j!=i->onePastLastContour();j++)
//    {
//      if (strcmp(i->getName(),j->getName()))
//      {
//        cout << "Object name (" << i->getName() << ") not equal to contour name (" << j->getName() << ")\n";
//      }
//      else
//      {
//        cout << "Object name (" << i->getName() << ") IS equal to contour name (" << j->getName() << ")\n";
//      }
//      sprintf(filename,"%s%s%d.pts",cs.getOutputDir(),j->getName(),j->getSection());
//      F = fopen(filename,"w");
//      if (!F) { printf("Couldn't open output file %s\n",filename); exit(0); }
//      fclose(F);
//    }
//  }
  for (c_o_iterator i = firstObject();i!=onePastLastObject();i++)
  {
    for (int j = i->getMinSection();j<=i->getMaxSection();j++)
    {
      sprintf(filename,"%s%s%d.pts",cs.getOutputDir(),i->getName(),j);
      F = fopen(filename,"w");
      if (!F) { printf("Couldn't open output file %s\n",filename); exit(0); }
      fclose(F);
    }
  }

  ////////// clear script files //////////
  sprintf(filename,"%smesh.sh",cs.getOutputDir());
  F = fopen(filename,"w");
  if (!F) { printf("Couldn't open output file %s\n",filename); exit(0); }
  sprintf(line,"#!/bin/bash\n\n");
  fputs(line,F);
  fclose(F);
  sprintf(filename,"%sconvert.sh",cs.getOutputDir());
  F = fopen(filename,"w");
  if (!F) { printf("Couldn't open output file %s\n",filename); exit(0); }
  sprintf(line,"#!/bin/bash\n\n");
  fputs(line,F);
  fclose(F);
}

int Container::getNumRawPoints (void)
{
  int sum = 0;
  for (o_iterator i = o.begin();i!=o.end();i++)
  {
    sum += i->getNumRawPoints();
  }
  return sum;
}


