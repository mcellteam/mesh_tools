#include <iostream>
#include <string>
#include <stdlib.h>

#include "contour.h"
#include "controls.h"
#include "object.h"

Object::Object (char const * const str,const int & section)
:name(str),min_section(section),max_section(section),contours()
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

int Object::getRanges (std::vector<int> & ranges) const
{
  ranges.clear();
  int num_parts = 0;
  // initialize so 0 equals no output points for section
  std::vector<int> sections(max_section+1,0);
  // for each contour 
  for (c_c_l_iterator i=contours.begin();i!=contours.end();i++)
  {
    // record 1 if section has contour with output points
    if (i->getNumSamplePoints()>0) sections[i->getSection()]++; 
  }
  int part[2] = {-1,-1};
  // for each section
  for (int i=min_section;i<max_section+1;i++)
  {
    // if first section
    if (i==min_section)
    {
      // rise
      if (sections[i])
      {
        part[0] = i;
      }
      // no change
      else
      {
      }
    }
    else
    {
      // rise
      if (!sections[i-1] && sections[i])
      {
        part[0] = i;
      }
      // fall 
      else if (sections[i-1] && !sections[i])
      {
        part[1] = i-1;
      }
      // same
      else
      {
      }
    }
    // if complete range
    if (part[1]>=0)
    {
      ranges.push_back(part[0]);
      ranges.push_back(part[1]);
      num_parts++;
      // reset
      part[0] = -1;
      part[1] = -1;
    }
  }
  // if partial
  if (part[0]>=0)
  {
    ranges.push_back(part[0]);
    ranges.push_back(max_section);
    num_parts++;
    // reset
    part[0] = -1;
    part[1] = -1;
  }
  return num_parts;
}

/** Write output contour points to file.
 * \param[in] num_parts Number of separate parts of object.
 * \param[in] ranges Section numbers of object parts.
 */

void Object::writeOutputContours (const int & num_parts,
                                  const std::vector<int> & ranges)
{
  Controls & cs(Controls::instance());
  // for each object part
  for (int i=0;i<num_parts;i++)
  {
    char myname[256];
    // determine object part name
    if (num_parts==1)
    {
      sprintf(myname,"%s",name.c_str());
    }
    else
    {
      char temp[256];
      sprintf(temp,"%s%s",name.c_str(),cs.getMultiPartSuffix());
      sprintf(myname,temp,i);
    }
    const int min_sec = ranges[2*i];
    const int max_sec = ranges[2*i+1];
    assert((2*i+1)<static_cast<int>(ranges.size()));
    printConfigFile(myname,min_sec,max_sec);
    appendScriptFile(myname);
    printPtsFiles(myname,min_sec,max_sec);
    if (cs.getPrintDetailedInfo())
    {
      printRawPtsFiles(myname,min_sec,max_sec);
    }
    if (cs.getCappingFlag())
    {
      printCaps(myname,min_sec,max_sec);
    }
  }
}

std::string Object::getOutputContourSerStr(const int slice) {
  for (c_l_iterator c = contours.begin(); c != contours.end(); c++) {
    if (c->getSection() == slice) {
      return c->getSerString();
    }
  }

  return "";
}

/** Add new contour to object.
 * \param[in] mycontour New contour.
 */

void Object::addContour (Contour mycontour)
{
  contours.push_back(mycontour);
}

/** Remove duplicate points in each contour of object.
 */

void Object::removeDuplicates (void)
{
  // for each contour 
  for (c_l_iterator i=contours.begin();i!=contours.end();i++)
  {
    i->removeDuplicates();
  }
}

/** Print raw contour points to file
 */

void Object::printRawPoints (void) const
{
  char filename[256];
  // for each contour 
  for (c_c_l_iterator i=contours.begin();i!=contours.end();i++)
  {
    sprintf(filename,"%s%s_%i_raw.log",Controls::instance().getOutputDir(),
            i->getName(),i->getSection());
    FILE *F = fopen(filename,"w");
    fclose(F);
  }
  // for each contour 
  for (c_c_l_iterator i=contours.begin();i!=contours.end();i++)
  {
    sprintf(filename,"%s%s_%i_raw.log",Controls::instance().getOutputDir(),
            i->getName(),i->getSection());
    FILE *F = fopen(filename,"a");
    fclose(F);
    i->printRawPoints(filename);
  }
}

/** Erase contours with fewer points than user-specified threshold.
 */

void Object::purgeBadContours (void)
{
  // for each contour 
  c_l_iterator i=contours.begin();
  while (i!=contours.end())
  {
    int j = i->getNumRawPoints();
    if (j<Controls::instance().getPtPerContourThreshold())
    {
      cout << "Contour not processed. Too few points: " << j
           << " < " << Controls::instance().getPtPerContourThreshold()
           << " (threshold). <Contour " << i->getName()
           << " on section " << i->getSection() << ">\n";
      i = contours.erase(i); 
    }
    else
    {
      i++;
    }
  }
}

/** Determine if object is degenerate,
 *  i.e. if number of spline samples is less than three.
 * \return True if degenerate; false otherwise.
 */

bool Object::isDegenerate (void) const
{
  bool flag = false;
  // for each contour
  for (c_c_l_iterator i=contours.begin();i!=contours.end();i++)
  {
    if(i->isDegenerate())
    {
      flag=true;
    }
  }
  return flag;
}

/** Print pts files for this object corresponding to caps.
 * \param[in] myname Object name.
 * \param[in] min_sec Minimum section of object.
 * \param[in] max_sec Maximum section of object.
 */

void Object::printCaps (char const * myname,
                        const int & min_sec,
                        const int & max_sec) const
{
  Controls & cs(Controls::instance());
  // print min capping pts file
  char filename[256],line[2048];
  sprintf(filename,"%s%s%d.pts",cs.getOutputDir(),myname,min_sec-1);
  FILE * F = fopen(filename,"w");
  if (!F) { printf("Couldn't open output file %s\n",filename);exit(1);}
  assert(F);
  double z = (static_cast<double>(min_sec)-1.0)*
             cs.getSectionThickness()*
             cs.getOutputScaleFactor();
  sprintf(line,"1\n0.0 0.0 %g\n",z);
  fputs(line,F);
  fclose(F);
  // print max capping pts file
  sprintf(filename,"%s%s%d.pts",cs.getOutputDir(),myname,max_sec+1);
  F = fopen(filename,"w");
  if (!F) { printf("Couldn't open output file %s\n",filename);exit(1);}
  assert(F);
  z = (static_cast<double>(max_sec)+1.0)*
      cs.getSectionThickness()*cs.getOutputScaleFactor();
  sprintf(line,"1\n0.0 0.0 %g\n",z);
  fputs(line,F);
  fclose(F);
}

/** Fit splines to contour points in each contour in this object.
 * \param[in] h Instance of Histogram class for sample point deviations.
 * \param[in] si Instance of Histogram class for sample intervals
 *  after simulated annealing.
 * \param[in] si_before Instance of Histogram class for sample intervals
 *  before simulated annealing.
 */

void Object::processContour (Histogram & h,
                             Histogram & si,
                             Histogram & si_before)
{
  assert(si_before.getN()==si.getN());
  assert(si_before.getN()==h.getN());
  Controls & cs (Controls::instance()); 
  std::list<int> control_points;
  for (c_l_iterator i = contours.begin();i!=contours.end();i++)
  {
    if (Controls::instance().getPrintDetailedInfo())
    {
      i->clearControlLogFiles();
      i->clearSParamLogFiles();
      i->clearSampleLogFiles();
      i->clearUniformSampleLogFiles();
      i->clearArcLengthsLogFiles();
      i->clearRadiusCurvatureLogFiles();
    }
  }
  assert(si_before.getN()==si.getN());
  assert(si_before.getN()==h.getN());
  for (c_l_iterator i = contours.begin();i!=contours.end();i++)
  {
    assert(si_before.getN()==si.getN());
    assert(si_before.getN()==h.getN());
    // small number of raw points
    // but more than threshold
    // so keep
    if (cs.getReturnRawContourPoints() || i->fewButKeep())
    {
      assert(si_before.getN()==si.getN());
      assert(si_before.getN()==h.getN());
      assert(si_before.getN()==si.getN());
      assert(si_before.getN()==h.getN());
      i->setSamplesToRaw();
      for (c_p_iterator j=i->getFirstRawPoint();j!=i->getOnePastLastRawPoint();j++)
      {
        h.update(0);
      }
      i->updateSampleIntervals(si_before);
      i->updateSampleIntervals(si);
      assert(si_before.getN()==si.getN());
      assert(si_before.getN()==h.getN());
      if (cs.getPrintDetailedInfo())
      {
        i->printSimpleScripts();
      }
      assert(si_before.getN()==si.getN());
      assert(si_before.getN()==h.getN());
    }
    // small number of raw points
    // less than threshold
    // so omit
    else if (i->fewOmit()) continue;
    // adequate number of raw points
    else
    {
      assert(si_before.getN()==si.getN());
      assert(si_before.getN()==h.getN());
      std::vector<double> my_path_intervals;
      my_path_intervals.reserve(100);
      int cnt = cs.getMaxDeviationAdjustments();
      i->initializeControlPointVector(control_points);
      bool next_contour = false;
      int target_num_samples=0;
      while (true)
      {
        assert(si_before.getN()==si.getN());
        assert(si_before.getN()==h.getN());
        i->initializeControlPointMatrix(control_points);
        // fit and sample
        target_num_samples = i->processContour(my_path_intervals);
        if (target_num_samples<3)
        {
          i->setSamplesToRaw();
          i->calcSampleIntervals(my_path_intervals);
        }
        // too few control points for deviation control to work
        // so plan to break out of loop
        if (i->getNumSamplePoints()<cs.getPtPerContourThreshold())
        {
          next_contour=true;
          break;
        }
        if (target_num_samples<3)
        {
          std::cout << "Contour '" << i->getName()
                    << "' on section " << i->getSection()
                    << " passed unaltered since number of sample points ("
                    << target_num_samples << ") an insufficient  "
                    << " of points (4) for splining"
                    << " but number of raw points (" << i->getNumRawPoints()
                    << ") is sufficient and exceeds user threshold ("
                    << Controls::instance().getPtPerContourThreshold()
                    << ")." << std::endl;
          break;
        }
        // check deviations
        i->computeDeviations();
        assert(si_before.getN()==si.getN());
        assert(si_before.getN()==h.getN());
        // if linearly interpolating contour points
        // then deviation will be miniscule,
        // so skip deviation control
        if (cs.getReturnInterpolatedRawPoints()) break;;
        // if NOT interested in controlling deviations
        if (cs.getDeviationThreshold()==false) break;
        // if deviations are small
        if (i->deviationsLarge()==false) break;
        // if tried to reduce deviations too many times
        if (!cnt)
        {
          cout << "Warning: Deviation adjustment of contour (" << name
                << ") ended by threshold ("
                << cs.getMaxDeviationAdjustments() << ")\n";
          assert(si_before.getN()==si.getN());
          assert(si_before.getN()==h.getN());
          break;
        }
        // duplicate control point
        i->duplicateControlPoints(control_points);
        cnt--;
        assert(si_before.getN()==si.getN());
        assert(si_before.getN()==h.getN());
      }
      if (next_contour) continue;
      // update min and max deviation distances
      if (target_num_samples<3)
      {
        for (c_p_iterator j=i->getFirstRawPoint();j!=i->getOnePastLastRawPoint();j++)
        {
          h.update(0);
        }
        //assert(static_cast<int>(my_path_intervals.size())>0);
        //si_before.scan(my_path_intervals);
      }
      else
      {
        i->updateDeviations(h);
        //i->updateSampleIntervals(si_before);
      }
      // update interval statistics
      assert(static_cast<int>(my_path_intervals.size())>0);
      si_before.scan(my_path_intervals);
      i->updateSampleIntervals(si);
      assert(si_before.getN()==si.getN());
      assert(si_before.getN()==h.getN());
      if (cs.getPrintDetailedInfo())
      {
        i->writeFilesAfter();
        i->printScripts();
      }
      assert(si_before.getN()==si.getN());
      assert(si_before.getN()==h.getN());
    }
  }
  assert(si_before.getN()==si.getN());
  assert(si_before.getN()==h.getN());
}

/** Write to file the config information for tiling this object.
 * \param[in] myname Object name.
 * \param[in] min_sec Minimum section of object.
 * \param[in] max_sec Maximum section of object.
 */

void Object::printConfigFile (char  const * myname,
                              const int & min_sec,
                              const int & max_sec) const
{
  Controls & cs(Controls::instance());
  char filename[256],line[2048];
  sprintf(filename,"%s%s.config",cs.getOutputDir(),myname);
  FILE * F = fopen(filename,"w");
  if (!F) { printf("Couldn't open output file %s\n",filename); exit(1); }
  assert(F);
  sprintf(line,"PREFIX: %s\nSUFFIX: .pts\nOUTPUT_FILE_NAME: %s\n",myname,myname);
  fputs(line,F);
  if (cs.getCappingFlag())
  {
    sprintf(line,"SLICE_RANGE: %i %i\n",min_sec-1 ,max_sec+1);
    sprintf(line,"%sMERGE_DISTANCE_SQUARE: 1E-24\n",line);
  }
  else
  {
    sprintf(line,"SLICE_RANGE: %i %i\n",min_sec,max_sec);
    sprintf(line,"%sMERGE_DISTANCE_SQUARE: 1E-24\n",line);
  }
  fputs(line,F);
  fclose(F);
}

/** Append to bash script files commands to tile this object
 *  and convert to Hughes Hoppe mesh format.
 * \param[in] myname Object name.
 */

void Object::appendScriptFile (char const * myname) const
{
  Controls & cs(Controls::instance());
  char filename[256],line[2048];
  // command to tile 
  sprintf(filename,"%smesh.sh",cs.getOutputDir());
  FILE * F = fopen(filename,"a");
  if (!F) { printf("Couldn't open output file %s\n",filename); exit(1); }
  assert(F);
  sprintf(line,"echo ''\ncontour_tiler -f %s.config &> /dev/null\n",myname);
  sprintf(line,"%secho '%s meshed'\n",line,myname);
  fputs(line,F);
  fclose(F);
  // command to convert to mesh format
  sprintf(filename,"%sconvert.sh",cs.getOutputDir());
  F = fopen(filename,"a");
  if (!F) { printf("Couldn't open output file %s\n",filename); exit(1); }
  assert(F);
  sprintf(line,"echo ''\npoly2mesh %s.poly > %s.mesh\n",myname,myname);
  sprintf(line,"%secho '%s converted'\n",line,myname);
  fputs(line,F);
  fclose(F);
}

/** For each contour in object write input contour points to pts file.
 * \param[in] myname Object name.
 * \param[in] min_sec Minimum section of object.
 * \param[in] max_sec Maximum section of object.
 */

void Object::printRawPtsFiles (char const * myname,
                               const int & min_sec,
                               const int & max_sec) const
{
  // for each contour in object
  for (c_c_l_iterator i =contours.begin();i!=contours.end();i++)
  {
    if (i->getSection()<min_sec || i->getSection()>max_sec) continue;
    i->printRawPtsFiles(myname);
  }
}

/** For each contour in object write sampled contour points to pts file.
 * \param[in] myname Object name.
 * \param[in] min_sec Minimum section of object.
 * \param[in] max_sec Maximum section of object.
 */

void Object::printPtsFiles (char const * myname,
                            const int & min_sec,
                            const int & max_sec) const
{
  // for each contour in object
  for (c_c_l_iterator i =contours.begin();i!=contours.end();i++)
  {
    if (i->getSection()<min_sec || i->getSection()>max_sec) continue;
    i->printPtsFiles(myname);
  }
}

/** Get cumulative number of points in all contours in this object. 
 */

int Object::getNumRawPoints (void) const
{
  int sum = 0;
  for (c_c_l_iterator i = contours.begin();i!=contours.end();i++)
  {
    sum += i->getNumRawPoints();
  }
  return sum;
}

