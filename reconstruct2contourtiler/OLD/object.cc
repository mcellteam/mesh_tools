#include <iostream>
#include <string>
#include <stdlib.h>

#include "contour.h"
#include "controls.h"
#include "object.h"

Object::Object (char * str,int sec)
:name(str),min_section(sec),max_section(sec),contours()
{
}

Object::Object (Object const & rhs)
:name(rhs.name),min_section(rhs.min_section),
  max_section(rhs.max_section),contours(rhs.contours)
{
}

Object& Object::operator = (const Object& rhs)
{
  std::cout << "Assignment operator prohibited on instances of object class.\n";
  std::cout << "Object " << rhs.name << std::endl;
  exit(1);
}

void Object::addContour (Contour mycontour)
{
  contours.push_back(mycontour);
}

void Object::removeDuplicates (void)
{
  // for each contour 
  for (c_l_iterator i=contours.begin();i!=contours.end();i++)
  {
    i->removeDuplicates();
  }
}

void Object::purgeBadContours (void)
{
  // for each contour 
  c_l_iterator i=contours.begin();
  while (i!=contours.end())
  {
    int j = i->getNumRawPoints();
    if (j<3)
    {
      cout << "Contour " << i->getName();
      cout << " on section " << i->getSection();
      cout << " not processed. Number of contour points (" << j;
      cout << ") < 3.\n";
      i = contours.erase(i); 
    }
    else
    {
      i++;
    }
  }
}

//bool Object::degenerateObject (void) const
//{
//  // for each contour
//  for (c_c_l_iterator i=contours.begin();i!=contours.end();i++)
//  {
//    if(i->getNumSplineSamples()<3){return true;}
//    if(i->getNumSplineSamples()<3){return true;}
//  }
//  return false;
//}

bool Object::degenerateObject (void) const
{
  // for each contour
  bool flag = false;
  for (c_c_l_iterator i=contours.begin();i!=contours.end();i++)
  {
    if(i->getNumSplineSamples()<3)
    {
      cout << "Warning. Contour <" << i->getName() << "," << i->getSection() << "> is degenerate (has only "
            << i->getNumSplineSamples() << " spline samples and "
            << i->getNumRawPoints() << " raw points).\n";
      flag=true;
    }
  }
  return flag;
}

void Object::printCaps (void) const
{
  Controls & cs(Controls::instance());
  ///// print min capping pts file /////
  char filename[256],line[2048];
  FILE *F;
  sprintf(filename,"%s%s%d.pts",cs.getOutputDir(),name.c_str(),min_section-1);
  // open file
  F = fopen(filename,"w");
  if (!F) { printf("Couldn't open output file %s\n",filename);exit(1);}
  // print file contents
  // sprintf(line,"1\n0.0 0.0 %d\n",(int)((min_section-1)*cs.getSectionThickness()*cs.getScale()));
  sprintf(line,"1\n0.0 0.0 %g\n",static_cast<double>(min_section-1)*cs.getSectionThickness()*cs.getScale());
  fputs(line,F);
  // close pts file
  fclose(F);
  ///// print max capping pts file /////
  sprintf(filename,"%s%s%d.pts",cs.getOutputDir(),name.c_str(),max_section+1);
  // open file
  F = fopen(filename,"w");
  if (!F) { printf("Couldn't open output file %s\n",filename);exit(1);}
  // print file contents 
  sprintf(line,"1\n0.0 0.0 %d\n",(int)((max_section+1)*cs.getSectionThickness()*cs.getScale()));
  fputs(line,F);
  // close pts file
  fclose(F);
}

void Object::interpolateRawPoints (void)
{
  for (c_l_iterator i = contours.begin();i!=contours.end();i++)
  {
    i->linearlyInterp();
  }
}

void Object::countPtsPerContour (std::vector<int> & points_per_contour)
{
  for (c_l_iterator i = contours.begin();i!=contours.end();i++)
  {
    points_per_contour.push_back(i->getNumRawPoints());
  }
}

void Object::fitSplines (Histogram & h)
{
  for (c_l_iterator i = contours.begin();i!=contours.end();i++)
  {
    //i->fitSplines(h,points_per_contour);
    i->fitSplines(h);
  }
}

void Object::computeHistogram (Histogram & h)
{
  for (c_l_iterator i = contours.begin();i!=contours.end();i++)
  {
    i->computeHistogram(h);
  }
}

void Object::printConfigFile (void)
{
  Controls & cs(Controls::instance());
  ///// print config file /////
  char filename[256],line[2048];
  FILE *F;
  sprintf(filename,"%s%s.config",cs.getOutputDir(),getName());
  // open file
  F = fopen(filename,"w");
  if (!F) { printf("Couldn't open output file %s\n",filename); exit(1); }
  // print file contents 
  sprintf(line,"PREFIX: %s\nSUFFIX: .pts\nOUTPUT_FILE_NAME: %s\n",getName() ,getName());
  fputs(line,F);
  if (cs.getCappingFlag())
  {
    sprintf(line,"SLICE_RANGE: %i %i\n",getMinSection()-1 ,getMaxSection()+1);
    sprintf(line,"%sMERGE_DISTANCE_SQUARE: 1E-24\n",line);
  }
  else
  {
    sprintf(line,"SLICE_RANGE: %i %i\n",getMinSection(),getMaxSection());
    sprintf(line,"%sMERGE_DISTANCE_SQUARE: 1E-24\n",line);
  }
  fputs(line,F);
  // close pts file
  fclose(F);
}

void Object::appendScriptFile (void)
{
  Controls & cs(Controls::instance());
  ///// append to script files /////
  char filename[256],line[2048];
  FILE *F;
  sprintf(filename,"%smesh.sh",cs.getOutputDir());
  // open file
  F = fopen(filename,"a");
  if (!F) { printf("Couldn't open output file %s\n",filename); exit(1); }
  // print file contents 
  sprintf(line,"echo ''\ncontour_tiler -f %s.config &> /dev/null\n",getName());
  sprintf(line,"%secho '%s meshed'\n",line,getName());
  fputs(line,F);
  // close pts file
  fclose(F);
  sprintf(filename,"%sconvert.sh",cs.getOutputDir());
  // open file
  F = fopen(filename,"a");
  if (!F) { printf("Couldn't open output file %s\n",filename); exit(1); }
  // print file contents 
  sprintf(line,"echo ''\npoly2mesh %s.poly > %s.mesh\n",getName(),getName());
  sprintf(line,"%secho '%s converted'\n",line,getName());
  fputs(line,F);
  // close pts file
  fclose(F);
}

void Object::printRawPtsFiles (void)
{
  Controls & cs(Controls::instance());
  char filename[256],line[2048];
  // for each contour in object
  for (c_l_iterator i =firstContourNonConst();i!=onePastLastContourNonConst();i++)
  {
    // create pts file
    sprintf(filename,"%sinput_%s%d.pts",cs.getOutputDir(),i->getName(),i->getSection());
    // open pts file
    FILE * F = fopen(filename,"a");
    if(!F){printf("Couldn't open output file %s\n",filename);exit(1);}
    // print number of contour points
    sprintf(line,"%i\n",i->getNumRawPoints());
    fputs(line,F);
    // for each interpolated point in contour
    for (c_p_l_iterator j=i->getFirstRawPoint();j!=i->getOnePastLastRawPoint();j++)
    {
      // print interpolated contour points
      sprintf(line,"%.15g %.15g %.15g\n",j->getX()*cs.getScale(),j->getY()*cs.getScale(),j->getZ()*cs.getScale());
      fputs(line,F);
    }
    // close pts file
    fclose(F);
  }
}

void Object::printPtsFiles (void)
{
  Controls & cs(Controls::instance());
  ///// print pts files of interpolated contour points /////
  char filename[256],line[2048];
  // for each contour in object
  for (c_l_iterator i =firstContourNonConst();i!=onePastLastContourNonConst();i++)
  {
    sprintf(filename,"%s%s%d.pts",cs.getOutputDir(),i->getName(),i->getSection());
    FILE * F = fopen(filename,"a");
    if (!F) { printf("Couldn't open output file %s\n",filename);exit(1); }
    // print number of contour points
    sprintf(line,"%i\n",i->getNumSplineSamples());
    fputs(line,F);
    // for each spline sample in contour
    for (c_p_iterator j=i->getFirstSplineSample();j!=i->getOnePastLastSplineSample();j++)
    {
      sprintf(line,"%.15g %.15g %.15g\n",j->getX()*cs.getScale(),j->getY()*cs.getScale(),j->getZ()*cs.getScale());
      fputs(line,F);
    }
    fclose(F);
  }
}

int Object::getNumRawPoints (void)
{
  int sum = 0;
  for (c_l_iterator i = contours.begin();i!=contours.end();i++)
  {
    sum += i->getNumRawPoints();
  }
  return sum;
}

