// #####################################################
// #####################################################

int binarySearch(int* sortedArray, int first, int last, int key) {
    // function:
    //   Searches sortedArray[first]..sortedArray[last] for key.
    // returns: index of the matching element if it finds key,
    //         otherwise  -(index where it could be inserted)-1.
    // parameters:
    //   sortedArray in  array of sorted (ascending) values.
    //   first, last in  lower and upper subscript bounds
    //   key         in  value to search for.
    // returns:
    //   index of key, or -insertion_position -1 if key is not
    //                 in the array. This value can easily be
    //                 transformed into the position to insert it.

    int mid;

    while (first <= last) {
        mid = (first + last) / 2;  // compute mid point.
        if (key > sortedArray[mid])
            first = mid + 1;  // repeat search in top half.
        else if (key < sortedArray[mid])
            last = mid - 1; // repeat search in bottom half.
        else
            return mid;     // found it. return position /////
    }
    return -(first + 1);    // failed to find key
}

// #####################################################
// #####################################################

void copyControlFile (void) {

    char line[1024];
    char filein[128];
    char fileout[128];

    // open input data file
	sprintf(filein,"controls.cc");
    std::ifstream inFile(filein);
    if (inFile.fail()) // if stream cannot be opened
    { cout << "Can't open " << filein ; exit(1); }

	// open output data file
	sprintf(fileout,"%sparameter_values.dat",OUTPUT_DATA_DIR);
    std::ofstream outFile(fileout);
    if (outFile.fail()) // if stream cannot be opened
    { cout << "Can't open " << fileout ; exit(1); }

	// foreach line in input file
    while (inFile.getline(line,1024)) {
		outFile << line << endl;
	}

	// close files
    inFile.close();
    outFile.close();
}

// #####################################################
// #####################################################

int compare_double (const void* a, const void* b ) {

  return (int)( *(double*)a - *(double*)b );

}

// #####################################################
// #####################################################

int compare_int (const void* a, const void* b ) {

  return ( *(int*)a - *(int*)b );

}

// #####################################################
// #####################################################

time_t recordTime(ofstream& myfile, time_t before, char string[128])
{
	time_t after,diff;
	after = time (NULL);
	diff = after-before;
	myfile << string << " " << diff << " seconds\n";
	myfile.flush();
	return after;
}

// #####################################################
// #####################################################

int distinguishable(double a,double b)
{
  double diff=a-b;
  if (a<0) a=-a;
  if (b<0) b=-b;
  if (diff<0) diff=-diff;
  if (a>b) return (diff>a*DOUBLE_EPSILON);
  else return (diff>b*DOUBLE_EPSILON);
}

// #####################################################
// #####################################################

void biggest(double *x, int& big) {

    if ( x[0] < x[1] ) {
        // x[0] < x[1]
        if ( x[0] < x[2] ) {
            // x[0] < x[1]
            // x[0] < x[2]
            if (x[1] < x[2]) {
                // x[0] < x[1]
                // x[1] < x[2]
                // x[0] is smallest
                // x[2] is biggest
                big = 2;
            } else {
                // x[0] < x[1]
                // x[2] <= x[1]
                // x[0] is smallest
                // x[1] is biggest or tied with x[2]
                big = 1;
            }
        } else {
            // x[0] < x[1]
            // x[2] <= x[0]
            // x[2] is smallest or tied with x[0]
            // x[1] is biggest
            big = 1;
        }
    } else {
        // x[1] <= x[0]
        if ( x[0] < x[2] ) {
            // x[1] <= x[0]
            // x[0] < x[2]
            // x1 is smallest or tied with x[0]
            // x2 is biggest
            big = 2;
        } else {
        	big = 0;
        }
    }

}

// #####################################################
// #####################################################

void threeValueSort(double *x, double *biggest, double *smallest) {

    if ( x[0] < x[1] ) {
        // x[0] < x[1]
        if ( x[0] < x[2] ) {
            // x[0] < x[1]
            // x[0] < x[2]
            if (x[1] < x[2]) {
                // x[0] < x[1]
                // x[1] < x[2]
                // x[0] is smallest
                // x[2] is biggest
                *smallest = x[0];
                *biggest = x[2];
            } else {
                // x[0] < x[1]
                // x[2] <= x[1]
                // x[0] is smallest
                // x[1] is biggest or tied with x[2]
                *smallest = x[0];
                *biggest = x[1];
            }
        } else {
            // x[0] < x[1]
            // x[2] <= x[0]
            // x[2] is smallest or tied with x[0]
            // x[1] is biggest
            *smallest = x[2];
            *biggest = x[1];
        }
    } else {
        // x[1] <= x[0]
        if ( x[0] < x[2] ) {
            // x[1] <= x[0]
            // x[0] < x[2]
            // x1 is smallest or tied with x[0]
            // x2 is biggest
            *smallest = x[1];
            *biggest = x[2];
        } else {
            // x[1] <= x[0]
            // x[2] <= x[0]
            if (x[1] < x[2]) {
                // x[1] <= x[0]
                // x[2] <= x[0]
                // x[1] < x[2]
                // x[1] is smallest
                // x[0] is biggest
                *smallest = x[1];
                *biggest = x[0];
            } else {
                // x[1] <= x[0]
                // x[2] <= x[1]
                // x[2] is smallest or tied with x[1]
                // x[0] is biggest or tied with x[1]
                *smallest = x[2];
                *biggest = x[0];
            }
        }
    }

}

// #####################################################
// #####################################################

void checkLinePolygonIntersection(double* pvc[3],double pn[3],double lp[2][3],
									bool *line_flag, bool *poly_flag, bool *poly_edge_flag) {

	//pvc = polygon_vertex_coordinates
	//pn  = polygon_normal
	//lp  = line_points

	int i, j, a, b, tv[2];
    double num, den, pI[3], u , area[3], p, det[3];
	int big;

	*line_flag = false;
	*poly_flag = false;
	*poly_edge_flag = false;

	// do polygons intersect?
	// e.g. use plane-plane intersection and compute overlap
	// e.g. use edge-polygon intersection method (USED THIS HERE)

	// compute polygon area on each of three principal planes
	// xy
	area[0] = fabs(( (pvc[2][0]-pvc[1][0]) * pvc[0][1] +
					(pvc[0][0]-pvc[2][0]) * pvc[1][1] +
					(pvc[1][0]-pvc[0][0]) * pvc[2][1])/2.0);
	// yz
	area[1] = fabs(( (pvc[2][1]-pvc[1][1]) * pvc[0][2] +
					(pvc[0][1]-pvc[2][1]) * pvc[1][2] +
					(pvc[1][1]-pvc[0][1]) * pvc[2][2])/2.0);
	//zx
	area[2] = fabs(( (pvc[2][2]-pvc[1][2]) * pvc[0][0] +
					(pvc[0][2]-pvc[2][2]) * pvc[1][0] +
					(pvc[1][2]-pvc[0][2]) * pvc[2][0])/2.0);

	biggest(&area[0],big);
	if (big == 0) {
		tv[0] = 0;
		tv[1] = 1;
	} else if (big == 1) {
		tv[0] = 1;
		tv[1] = 2;
	} else {
		tv[0] = 0;
		tv[1] = 2;
	}


	// use line points, polygon normal and one polygon vertex
	// to compute where line intersects plane of polygon (if not parallel)

	// denominator of dot products
	den = pn[0]*(lp[1][0]-lp[0][0]) 
		+ pn[1]*(lp[1][1]-lp[0][1]) 
		+ pn[2]*(lp[1][2]-lp[0][2]);
	
	// if line and polygon plane are not parallel
	if  (den) {

		// numerator of dot products
		num = pn[0]*(pvc[0][0]-lp[0][0]) 
			+ pn[1]*(pvc[0][1]-lp[0][1]) 
			+ pn[2]*(pvc[0][2]-lp[0][2]);

		// point of intersection
		u = num/den;

		// if polygon cuts through line
		if (u > 0 && u < 1) {
			*line_flag = true;

			pI[0] = (1-u)*lp[0][0] + u*lp[1][0];
			pI[1] = (1-u)*lp[0][1] + u*lp[1][1];
			pI[2] = (1-u)*lp[0][2] + u*lp[1][2];

			////////// is point of intersection on other polygon? //////////
		
			// does point of intersection lie on polygon edge?
			*poly_edge_flag = false;

			if (detect_polygon_edge_intersection) {
				// if first coordinate of polygon vertex coordinates are distinguishable
				// else second coordinate of polygon vertex coordinates are distinguishable
				if (pvc[0][tv[0]] !=  pvc[1][tv[0]]) {
					p = pvc[0][tv[1]]+(pI[tv[0]]-pvc[0][tv[0]])/ 
						(pvc[1][tv[0]]-pvc[0][tv[0]])* 
						(pvc[1][tv[1]]-pvc[0][tv[1]]);
					if(!distinguishable(p,pI[tv[1]])){*poly_edge_flag = true;}
				} 
				else {
					p = pvc[0][tv[0]]+ (pI[tv[1]] -pvc[0][tv[1]])/ 
						(pvc[1][tv[1]] -pvc[0][tv[1]])* 
						(pvc[1][tv[0]] -pvc[0][tv[0]]);
					if(!distinguishable(p,pI[tv[0]])){*poly_edge_flag = true;}
				}
		
				if (pvc[1][tv[0]] != pvc[2][tv[0]]) {
					p = pvc[1][tv[1]]+ (pI[tv[0]] -pvc[1][tv[0]])/ 
						(pvc[2][tv[0]] -pvc[1][tv[0]])* 
						(pvc[2][tv[1]] -pvc[1][tv[1]]);
					if(!distinguishable(p,pI[tv[1]])){*poly_edge_flag = true;}
				} 
				else {
					p = pvc[1][tv[0]]+ (pI[tv[1]] -pvc[1][tv[1]])/ 
						(pvc[2][tv[1]] -pvc[1][tv[1]])* 
						(pvc[2][tv[0]] -pvc[1][tv[0]]);
					if(!distinguishable(p,pI[tv[0]])){*poly_edge_flag = true;}
				}
		
		
				if (pvc[2][tv[0]] != pvc[0][tv[0]]) {
					p = pvc[2][tv[1]]+ (pI[tv[0]] -pvc[2][tv[0]])/ 
						(pvc[0][tv[0]] -pvc[2][tv[0]])* 
						(pvc[0][tv[1]] -pvc[2][tv[1]]);
					if(!distinguishable(p,pI[tv[1]])){*poly_edge_flag = true;}
				} 
				else {
					p = pvc[2][tv[0]]+ (pI[tv[1]] -pvc[2][tv[1]])/ 
						(pvc[0][tv[1]] -pvc[2][tv[1]])* 
						(pvc[0][tv[0]] -pvc[2][tv[0]]);
					if(!distinguishable(p,pI[tv[0]])){*poly_edge_flag = true;}
				}
			}
		
			// if point of intersection is not on polygon edge
			if (!*poly_edge_flag) {

				// compute three determinants
				det[0] = (pvc[0][tv[0]] -pI[tv[0]]) *(pvc[1][tv[1]] -pI[tv[1]])
						-(pvc[1][tv[0]] -pI[tv[0]]) *(pvc[0][tv[1]] -pI[tv[1]]);
				det[1] = (pvc[1][tv[0]] -pI[tv[0]]) *(pvc[2][tv[1]] -pI[tv[1]])
						-(pvc[2][tv[0]] -pI[tv[0]]) *(pvc[1][tv[1]] -pI[tv[1]]);
	
				if ( ( (det[0] < 0) && (det[1] < 0) ) || ( (det[0] > 0) && (det[1] > 0) )){

					det[2] = (pvc[2][tv[0]] -pI[tv[0]]) *(pvc[0][tv[1]] -pI[tv[1]])
							-(pvc[0][tv[0]] -pI[tv[0]]) *(pvc[2][tv[1]] -pI[tv[1]]);
					if ( ( (det[0] < 0) && (det[1] < 0) && (det[2] < 0) ) || 
						( (det[0] > 0) && (det[1] > 0) && (det[2] > 0) ) ){
						// line intersects polygon plane inside polygon
					    *poly_flag = true;
					}
				}
			}
		}
	}
}

// #####################################################
// #####################################################

void addArrayElement(int*& ptr,int& num,int& max,int incr,int val) {

	int* temp;
	int i, temp_size;

	// before adding element to array
	// check if there is room
	if (num < max) {
		// add element
		ptr[num] = val;
		num++;
	} else {
		// first need to make more room in array
		// create temp array
		temp = new int[max];
		if (temp == NULL) { 
			cout << "Not enough memory for temp in addArrayElement\n";
			cout.flush();
			exit(1);
		}
		// store current data in temp array
		temp_size = max;
		for (i=0;i<temp_size;i++) {
			temp[i] = ptr[i];
		}
		// delete polygons array
		delete[] ptr;
		// compute new, bigger array size
		max = max + incr;
		// create bigger array;
		ptr = new int[max];
		if (ptr == NULL) { 
			cout << "Not enough memory for temp in addArrayElement\n";
			cout.flush();
			exit(1);
		}
		// copy temp data back to array
		for (i=0;i<temp_size;i++) {
			ptr[i] = temp[i];
		}
		// delete temp array
		delete[] temp;
		// add new element
		ptr[num] = val;
		num++;
	}
}

// #####################################################
// #####################################################

void addSortUniqueArrayElement(int*& ptr,int& num,int& max,int incr,int val) {

	if (num==0) {
		if (num < max ) {
			ptr[num] = val;
			num++;
		} else {
			// print error, because this case should never arise
			cout << "addSortUniqueArrayElement: ERROR: num and max are both 0.\n";
			exit(1);
		}
	} else {
        //find location of val in ptr
		int j;
		j = binarySearch(ptr,0,num-1,val);

		// if val was not found in ptr
		if (j<0) {
			// compute insert location of val
			j = -(j+1);

			// is j location in ptr available?
			if (j == num) {
				// insert val into ptr
				addArrayElement(ptr,num,max,incr,val);
			} else {
				// slide bar elements down
                if (num < max) {
                    // there is room to slide elements down
					int i;
                    for (i=num-1;i>=j;i--) {
                        ptr[i+1]=ptr[i];
                    }
                    // insert new element
                    ptr[j]=val;
                    num++;
                } else {
					// increase size of array
					addArrayElement(ptr,num,max,incr,ptr[num-1]);
                    // there is room to slide elements down
					int i;
                    for (i=num-3;i>=j;i--) {
                        ptr[i+1]=ptr[i];
                    }
                    // insert new element
                    ptr[j]=val;
                }
            }
		}
	}

}

// #####################################################
// #####################################################

int* uniqueArray(int* ptr,int& num, int max) {

	int* unique;
	unique = new int[max];
	if (unique == NULL) { 
		cout << "Not enough memory for unique in uniqueArray\n";
		cout.flush();
		exit(1);
	}
	int num_unique = 0;
	int i;

	if (num) {
		unique[0] = ptr[0];
		num_unique++;

		// for the 1st through penultimate elements of ptr (assumed sorted)
		for (i=0;i<num-1;i++) {
			// if next element is not the same as current
			if (ptr[i+1] != ptr[i]) {
				// add element to unique array
				unique[num_unique] = ptr[i+1];
				num_unique++;
			}
		}

		// update num
		num = num_unique;
	}

	// delete old array
	delete[] ptr;

	// return unique
	return unique;

}

// #####################################################
// #####################################################

void BoxClass::addPolygon(int& polygon_index) {

	int* temp;
	int i, temp_size;

	// before adding polygon to polygons array
	// check if there is room
	if (num_polygons < max_polygons) {
		// add polygon
		polygons[num_polygons] = polygon_index;
		num_polygons++;
	} else {
		// first need to make more room in polygons array
		// create temp array
		temp = new int[max_polygons];
		if (temp == NULL) { 
			cout << "Not enough memory for temp in addPolygon\n";
			cout.flush();
			exit(1);
		}
		// store current data in temp array
		temp_size = max_polygons;
		for (i=0;i<temp_size;i++) {
			temp[i] = polygons[i];
		}
		// delete polygons array
		delete[] polygons;
		// compute new, bigger polygons array size
		max_polygons = max_polygons+INCREMENT_OF_POLYGONS_IN_BOXES;
		// create bigger polygons array;
		polygons = new int[max_polygons];
		if (polygons == NULL) { 
			cout << "Not enough memory for polygon in addPolygons\n";
			cout.flush();
			exit(1);
		}
		// copy temp data back to polygons
		for (i=0;i<temp_size;i++) {
			polygons[i] = temp[i];
		}
		// delete temp array
		delete[] temp;
		// add new polygon
		polygons[num_polygons] = polygon_index;
		num_polygons++;
	}
}

// #####################################################
// #####################################################

BoxClass::~BoxClass(void) {

	delete[] polygons;

}

// #####################################################
// #####################################################

int PolygonClass::checkEdgeEdgeIntersection(	double* cpvc[3], int current_pairs[3][2], 
												double* opvc[3],  int current_polygon_index, 
												int other_polygon_index, bool share_edge_flag
											) {

	int i, j, k;
	double cv1[3], cv2[3], ov1[3], ov2[3], x[3], y[3];
	double num, den, u, q, qNum, qDen, uNum, uDen, cos_theta_square;
	double term1;
	bool intersect_flag, parallel_flag;

	// cpvc = current_polygon_vertex_coordinates
	// opvc = other_polygon_vertex_coordinates
	// cv   = current_vertex
	// ov   = other_vertex

	// initialize flag
	intersect_flag = false;

	// for each current polygon edge
	for (i=0;i<3;i++) {
		// for each other polygon edge
		for (j=0;j<3;j++) {
	
			// define two current and two other polygon vertices
			for (k=0;k<3;k++) {
				cv1[k] = cpvc[current_pairs[i][0]][k] ;
				cv2[k] = cpvc[current_pairs[i][1]][k] ;
				ov1[k] = opvc[current_pairs[j][0]][k] ;
				ov2[k] = opvc[current_pairs[j][1]][k] ;
			}
	
			// if the edges do not share a vertex
			if (
				(distinguishable(cv1[0],ov1[0]) &&
				 distinguishable(cv1[1],ov1[1]) &&
				 distinguishable(cv1[2],ov1[2]) ) 
				&&
				(distinguishable(cv1[0],ov2[0]) &&
				 distinguishable(cv1[1],ov2[1]) &&
				 distinguishable(cv1[2],ov2[2]) ) 
				&&
				(distinguishable(cv2[0],ov1[0]) &&
				 distinguishable(cv2[1],ov1[1]) &&
				 distinguishable(cv2[2],ov1[2]) ) 
				&&
				(distinguishable(cv2[0],ov2[0]) &&
				 distinguishable(cv2[1],ov2[1]) &&
				 distinguishable(cv2[2],ov2[2]) ) 
				) {

				// and the edges are not parallel
				parallel_flag = false;
				for (k=0;k<3;k++) {
					y[k] = ov2[k]-ov1[k];
				}

				term1 = (x[0]*y[0]+ x[1]*y[1]+ x[2]*y[2]);
				if ( !distinguishable(term1*term1,
				(x[0]*x[0]+ x[1]*x[1]+ x[2]*x[2])*(y[0]*y[0]+ y[1]*y[1]+ y[2]*y[2])) ) {parallel_flag = true;}
	
				if (!parallel_flag) {
					// compute scalars
					qDen = (cv2[0]-cv1[0])*(ov2[1]-ov1[1])-(cv2[1]-cv1[1])*(ov2[0]-ov1[0]);
					qNum = (cv2[0]-cv1[0])*(cv1[1]-ov1[1])-(cv2[1]-cv1[1])*(cv1[0]-ov1[0]);
					uNum = ov1[0]-cv1[0]+q*(ov2[0]-ov1[0]);
					uDen = cv2[0]-cv1[0];
	
					if(fabs(qDen)>DOUBLE_EPSILON && fabs(uDen)>DOUBLE_EPSILON) {
						q = qNum/qDen;
						u = uNum/uDen;
	
						if ( 
							((u > DOUBLE_EPSILON && u < (1-DOUBLE_EPSILON)) && 
							(q > DOUBLE_EPSILON && q < (1-DOUBLE_EPSILON)) ) ||
							( share_edge_flag && ((u > DOUBLE_EPSILON && u < (1-DOUBLE_EPSILON)) || 
												(q > DOUBLE_EPSILON && q < (1-DOUBLE_EPSILON))) )
							) {

						    return(1);
						}
					}
				}
			}
		}
	}

	return(0);

}

// #####################################################
// #####################################################

int PolygonClass::checkPolygonEdgePolygonIntersection(	double* cpvc[3],
														int cp[3][2], double on[3],
														double* opvc[3],  
														int current_polygon_index,
														int other_polygon_index) {

	// NOTE NEED ONLY CHECK EITHER (1) INTERSECTION OF CURRENT POLYGON EDGES
	// WITH OTHER POLYGON OR (2) VICE VERSA, SINCE THE ROLES OF CURRENT AND OTHER
	// POLYGON WILL GET REVERSED. CHOOSE OPTION 1.
	//
	// FURTHER NOTE THAT POLYGON INTERSECTION CHECK IS ASSYMETRIC
	// I.E. CHECK IS IF CURRENT POLYGON EDGE INTERSECTS OTHER POLYGON ONLY
	// NO CHECKING IS DONE IF CURRENT POLYGON IS INTERSECTED BY OTHER POLYGON EDGE

	//cpvc = current_polygon_vertex_coordinates
	//opvc = other_polygon_vertex_coordinates
	//tv   = target_vertices
	//cv1  = current_vertex 1
	//cv2  = current_vertex 2
	//on   = other_normal
	//cp   = current_pairs

	int i;
	double lp[2][3];
	bool line_flag, poly_flag, poly_edge_flag;

	line_flag = false;
	poly_flag = false;
	// for each current polygon edge
	for (i=0;i<3;i++) {
	
		if(!line_flag && !poly_flag) {

			lp[0][0] = cpvc[cp[i][0]][0];
			lp[0][1] = cpvc[cp[i][0]][1];
			lp[0][2] = cpvc[cp[i][0]][2];
			lp[1][0] = cpvc[cp[i][1]][0];
			lp[1][1] = cpvc[cp[i][1]][1];
			lp[1][2] = cpvc[cp[i][1]][2];

			checkLinePolygonIntersection(opvc,on,lp,&line_flag,&poly_flag,&poly_edge_flag);

		}
	}
	// do polygons intersect?
	if (line_flag && poly_flag) {return(1);}
	else {return(0);}
}

// #####################################################
// #####################################################

void PolygonClass::recordBoxes(int* list,int num_list) {

	int i;

	// record boxes in polygon struct
	if (boxes != NULL) {
		delete [] boxes;
	}
	boxes = new int[NUMBER_OF_BOXES];
	if (boxes == NULL) { 
		cout << "Not enough memory for boxes in PolygonClass::recordBoxes\n";
		cout.flush();
		exit(1);
	}
	num_boxes = 0;
	max_boxes = NUMBER_OF_BOXES;
	for (i=0;i<num_list;i++) {
		addArrayElement(boxes,num_boxes,max_boxes,INCREMENT_OF_BOXES,list[i]);
	}

}

// #####################################################
// #####################################################

PolygonClass::~PolygonClass(void) {

	delete[] boxes;
	delete[] intersections;

}

// #####################################################
// #####################################################

void VertexClass::computeMeanOutwardNormal(PolygonClass *polygons, 
											double normal[3]) {

	/////////////// ////////////////////////////
	// this code is being executed 
	// for each vertex structure in each object
	/////////////// ////////////////////////////

	int polygon_index, i;
	double normal_sum_x, normal_sum_y, normal_sum_z;

	// declare arrays for storing normal sums
	normal_sum_x =0;
	normal_sum_y =0;
	normal_sum_z =0;

	// for each adjacent polygon
	for (i=0;i<next_adjacent_polygon;i++) {

		// get adjacent polygon index
		polygon_index = adjacent_polygons_index[i];

		// get coordinates of polygon normal
		// and add to sum
		normal_sum_x += polygons[polygon_index].normal_components[0];
		normal_sum_y += polygons[polygon_index].normal_components[1];
		normal_sum_z += polygons[polygon_index].normal_components[2];
	}

	// average normals
	normal[0] = normal_sum_x/next_adjacent_polygon;
	normal[1] = normal_sum_y/next_adjacent_polygon;
	normal[2] = normal_sum_z/next_adjacent_polygon;

	// store normals
	outward_normal[0] = normal[0];
	outward_normal[1] = normal[1];
	outward_normal[2] = normal[2];

}

// #####################################################
// #####################################################

void VertexClass::computeEnergyForce(double* x, double* y, double* z, 
										int* edge_indices, int num_edge_indices,EdgeClass* edges) {

	int k, index;
	double deviation[3], scaled_spring_force;

	if (candidate) {

		////////// compute deviation vector //////////
		for (k=0;k<3;k++) {
			deviation[k] = closest[k]-new_coordinates[k];
		}

		////////// compute deviation distance //////////
		// if not nice, negate value of d
		deviation_distance_old = deviation_distance;
		deviation_distance = sqrt(deviation[0]*deviation[0]
								+deviation[1]*deviation[1]
								+deviation[2]*deviation[2]);
		if (!nice_data.nice) { deviation_distance = -deviation_distance;}

		// compute position error
		position_error = deviation_distance-TARGET_SEPARATION;
		// add energy contribution
		energy = SEPARATION_WEIGHT/2.0*position_error*position_error;
		// add force contribution
		// spring_force = spring constant * stretch
		// scaled_spring force = spring_force/deviation_distance
		scaled_spring_force = SEPARATION_WEIGHT*position_error/deviation_distance;

		// force cartesian component = spring_force * unit vector
		// where unit vector = (cartesian deviation_distance component) / deviation_distance
		// force cartesian component = scaled_spring_force * cartesian component deviation_distance
		force[0] = scaled_spring_force*deviation[0];
		force[1] = scaled_spring_force*deviation[1];
		force[2] = scaled_spring_force*deviation[2];
	
		////////// loop through all adjacent vertices of current vertex //////////
		// for each adjacent vertex of current vertex
		for (k = 0;k < next_adjacent_vertex;k++) {
	
			// energy contribution
			energy += EDGE_STRETCH_WEIGHT/2.0
						*(new_edge_lengths[k] -original_edge_lengths[k])
						*(new_edge_lengths[k] -original_edge_lengths[k]);
	
			// force contribution
			// spring_force = spring constant * stretch
			// scaled_spring force = spring_force/edge length
			scaled_spring_force = EDGE_STRETCH_WEIGHT
									*(new_edge_lengths[k] -original_edge_lengths[k])
									/new_edge_lengths[k];

			// force cartesian component = spring_force * unit vector
			// where unit vector = (adjacent vertex position - target vertex position) / edge length
			// force cartesian component = scaled_spring_force * cartesian component difference
			force[0] += scaled_spring_force *(x[k]-new_coordinates[0]);
			force[1] += scaled_spring_force *(y[k]-new_coordinates[1]);
			force[2] += scaled_spring_force *(z[k]-new_coordinates[2]);

		}

		////////// loop through all associated edges of current vertex //////////
		for (k=0;k<num_edge_indices;k++) {
	
			// energy contribution
			energy += edges[edge_indices[k]].energy/2;

			// identify which vertex is current
			if (edges[edge_indices[k]].other_vertex_indices[0] == vertex_index) { index = 0;}
			else { index = 1;}
	
			// force contribution
			force[0] += edges[edge_indices[k]].force[index][0];
			force[1] += edges[edge_indices[k]].force[index][1];
			force[2] += edges[edge_indices[k]].force[index][2];

		}

		// add intersection_force
		// NOTE: I DON'T HAVE A CORRESPONDING INTERSECTION ENERGY TERM
		force[0] += intersection_force[0];
		force[1] += intersection_force[1];
		force[2] += intersection_force[2];

	}
	
}

// #####################################################
// #####################################################

void VertexClass::setOdd(int num_odd_objects) {

	int i;

	// number of odd objects

	// if no odd
	if ( num_odd_objects == 0) {
		// vertex is nice
		nice_data.nice = 1;
	} else {
		// vertex is not nice
		nice_data.nice = 0;
	}
}

// #####################################################
// #####################################################

void VertexClass::computeNewCoords(double gain) {

	int k;
	double disp[3];

	// compute displacement vector
	for (k=0;k<3;k++) {
		disp[k] = gain*force[k];
	}

	// compute new vertex coordinates
	for (k=0;k<3;k++) {
		new_coordinates[k] = disp[k]+new_coordinates[k];
	}

	// write disp to file
	double net_disp;
	net_disp = sqrt(
				disp[0]*disp[0]+
				disp[1]*disp[1]+
				disp[2]*disp[2]
				);
}

// #####################################################
// #####################################################

void VertexClass::fileInit(int object_index) {

	char file[128];
	sprintf(file,"%s%s%d_%d.log",OUTPUT_DATA_DIR,VERTEX_LOG_FILE,object_index,vertex_index);
	Vertexfile.open(file);
    Vertexfile << "Object = " << object_index << ", "
			<< "Vertex = " << vertex_index << "\n";
    Vertexfile.width(15);
    Vertexfile << left << "Iteration";
    Vertexfile.width(15);
    Vertexfile << left << "candidate?";
    Vertexfile.width(15);
	Vertexfile << left << "Dev. Dist. (nm)";
    Vertexfile.width(15);
	Vertexfile << left << "Energy (J)";
    Vertexfile.width(15);
	Vertexfile << left << "nice?";
    Vertexfile.width(15);
	Vertexfile << left << "Force (nN)\n";
	Vertexfile.flush();
}

// #####################################################
// #####################################################

void VertexClass::updateLog(int iteration) {

	    Vertexfile.width(15);
		Vertexfile << left << iteration;
	    Vertexfile.width(15);
		Vertexfile << left << candidate;
	    Vertexfile.width(15);
		Vertexfile << left << deviation_distance;
	    Vertexfile.width(15);
		Vertexfile << left << energy;
	    Vertexfile.width(15);
		Vertexfile << left << nice_data.nice;
	    Vertexfile.width(15);
		Vertexfile << left << force[0] << " " << force[1] << " " << force[2] << "\n";
		Vertexfile.flush();
}

// #####################################################
// #####################################################

void VertexClass::fileOutit(void) {

	Vertexfile.close();

}

// #####################################################
// #####################################################

void VertexClass::addNeighbor(int& neighbor_index) {

	int* temp;
	int i, temp_size;

	// before adding neighbor to neighbors array
	// check if there is room
	if (next_neighbor_vertex < max_neighborhood_vertices_index) {
		// add neighbor
		neighborhood_vertices_index[next_neighbor_vertex] = neighbor_index;
		next_neighbor_vertex++;
	} else {
		// first need to make more room in neighbors array
		// create temp array
		temp = new int[max_neighborhood_vertices_index];
		if (temp == NULL) { 
			cout << "Not enough memory for temp in VertexClass::addNeighbor\n";
			cout.flush();
			exit(1);
		}
		// store current data in temp array
		temp_size = max_neighborhood_vertices_index;
		for (i=0;i<temp_size;i++) {
			temp[i] = neighborhood_vertices_index[i];
		}
		// delete neighbors array
		delete[] neighborhood_vertices_index;
		// compute new, bigger neighborss array size
		max_neighborhood_vertices_index = max_neighborhood_vertices_index+INCREMENT_OF_NEIGHBORS;
		// create bigger neighbors array;
		neighborhood_vertices_index = new int[max_neighborhood_vertices_index];
		if (neighborhood_vertices_index == NULL) { 
			cout << "Not enough memory for neighborhood_vertices_index in VertexClass::addNeighbor\n";
			cout.flush();
			exit(1);
		}
		// copy temp data back to neighbors
		for (i=0;i<temp_size;i++) {
			neighborhood_vertices_index[i] = temp[i];
		}
		// delete temp array
		delete[] temp;
		// add new neighbor
		neighborhood_vertices_index[next_neighbor_vertex] = neighbor_index;
		next_neighbor_vertex++;
	}
}

// #####################################################
// #####################################################

void VertexClass::addSortUniqueNeighbor(int& neighbor_index) {

	if (next_neighbor_vertex == 0) {
		if (next_neighbor_vertex < max_neighborhood_vertices_index) {
			// add neighbor
			neighborhood_vertices_index[next_neighbor_vertex] = neighbor_index;
			next_neighbor_vertex++;
		} else {
            // print error, because this case should never arise
            cout << "addSortUniqueNeighbor: ERROR: num and max are both 0.\n";
            exit(1);
		}
	} else {
		// find location of neighbor_index in vertices_index
		int j = 0;
        while (j<next_neighbor_vertex && (neighbor_index>neighborhood_vertices_index[j])) { j++; }

        // j = insert location of val into neighborhood_vertices_index
        // is j location in neighborhood_vertices_index available?
        if (j == next_neighbor_vertex) {
            // insert neighbor_index into neighborhood_vertices_index
            addNeighbor(neighbor_index);
        } else {
            // if neighbor_index equals neighborhood_vertices_index[j] do nothing
            // else insert neighbor_index into neighborhood_vertices_index
            if (neighbor_index != neighborhood_vertices_index[j]) {
                // slide bar elements down
                if (next_neighbor_vertex < max_neighborhood_vertices_index) {
                    // there is room to slide elements down
                    int i;
                    for (i=next_neighbor_vertex-1;i>=j;i--) {
                        neighborhood_vertices_index[i+1]=neighborhood_vertices_index[i];
                    }
                    // insert new element
                    neighborhood_vertices_index[j]=neighbor_index;
                    next_neighbor_vertex++;
                } else {
                    // increase size of array
                    addNeighbor(neighborhood_vertices_index[next_neighbor_vertex-1]);
                    // there is room to slide elements down
                    int i;
                    for (i=next_neighbor_vertex-3;i>=j;i--) {
                        neighborhood_vertices_index[i+1]=neighborhood_vertices_index[i];
                    }
                    // insert new element
                    neighborhood_vertices_index[j]=neighbor_index;
                }
			}
        }
    }
		
}

// #####################################################
// #####################################################

VertexClass::~VertexClass(void) {

	delete[] adjacent_polygons_index;
	delete[] adjacent_vertices_index;
	delete[] neighborhood_vertices_index;
	delete[] original_edge_lengths;
	delete[] new_edge_lengths;
	delete[] edges_indices;
}

// #####################################################
// #####################################################

void EdgeClass::initAngles(PolygonClass *polygons) {

	if (init_angle_to_flat) {
		original_angle_cosine = 1;
	} else {
		double n[2][3], dot_product, magn[2];

		// get polygon outward normals
		n[0][0] = polygons[polygon_indices[0]].normal_components[0];
		n[0][1] = polygons[polygon_indices[0]].normal_components[1];
		n[0][2] = polygons[polygon_indices[0]].normal_components[2];
		n[1][0] = polygons[polygon_indices[1]].normal_components[0];
		n[1][1] = polygons[polygon_indices[1]].normal_components[1];
		n[1][2] = polygons[polygon_indices[1]].normal_components[2];

		//// compute cosine of angle between normals ////
		// dot product
		dot_product = n[0][0]*n[1][0]+n[0][1]*n[1][1]+n[0][2]*n[1][2];
		// magnitude of normals
		magn[0] = sqrt(n[0][0]*n[0][0]+n[0][1]*n[0][1]+n[0][2]*n[0][2]);
		magn[1] = sqrt(n[1][0]*n[1][0]+n[1][1]*n[1][1]+n[1][2]*n[1][2]);

		// cosine of angle
		original_angle_cosine = dot_product/magn[0]/magn[1];
	}

}

// #####################################################
// #####################################################

void EdgeClass::computeAngleForceEnergy(PolygonClass *polygons) {

	double dot_product, mean[3], squared_magn[2], force_magn, normal_length[2],
			norm_cross[3], vec[3], force_dir[3], final_dot_product;
	double N1x, N1y, N1z, N2x, N2y, N2z;

	// get polygon outward normals
	N1x = polygons[polygon_indices[0]].normal_components[0];
	N1y = polygons[polygon_indices[0]].normal_components[1];
	N1z = polygons[polygon_indices[0]].normal_components[2];
	N2x = polygons[polygon_indices[1]].normal_components[0];
	N2y = polygons[polygon_indices[1]].normal_components[1];
	N2z = polygons[polygon_indices[1]].normal_components[2];
	normal_length[0] = sqrt(N1x*N1x+N1y*N1y+N1z*N1z);
	normal_length[1] = sqrt(N2x*N2x+N2y*N2y+N2z*N2z);

	//// force direction for other_vertex_indices[0]////
	// normals cross product = n[0] x n[1]
	norm_cross[0] = N1y*N2z-N1z*N2y;
	norm_cross[1] = N1z*N2x-N1x*N2z;
	norm_cross[2] = N1x*N2y-N1y*N2x;
	// vector: edge to other_vertex_indices[0]
	vec[0] = other_vertex_coords[0]-edge_vertex_coords[0];
	vec[1] = other_vertex_coords[1]-edge_vertex_coords[1];
	vec[2] = other_vertex_coords[2]-edge_vertex_coords[2];
	// force direction = normals cross product x vec
	force_dir[0] = norm_cross[1]*vec[2]-norm_cross[2]*vec[1];
	force_dir[1] = norm_cross[2]*vec[0]-norm_cross[0]*vec[2];
	force_dir[2] = norm_cross[0]*vec[1]-norm_cross[1]*vec[0];
	// dot product of force direction and normal
	final_dot_product = force_dir[0]*N1x+force_dir[1]*N1y+force_dir[2]*N1z;
	if (final_dot_product < 0 ) {flip = -1;}
	else {flip = 1;}

	//// compute cosine of angle between normals ////
	// dot product
	dot_product = N1x*N2x+N1y*N2y+N1z*N2z;
	// magnitude of normals
	squared_magn[0] = N1x*N1x+N1y*N1y+N1z*N1z;
	squared_magn[1] = N2x*N2x+N2y*N2y+N2z*N2z;
	// cosine of angle
	new_angle_cosine = sqrt(dot_product*dot_product/squared_magn[0]/squared_magn[1]);
	if (dot_product<0) { new_angle_cosine = -new_angle_cosine; }

	//// force and energy ////
	force_magn = ANGLE_STRETCH_WEIGHT*fabs(new_angle_cosine-original_angle_cosine);
	energy = ANGLE_STRETCH_WEIGHT*(new_angle_cosine-original_angle_cosine)
								*(new_angle_cosine-original_angle_cosine)/2;

	force[0][0] = force_magn/2*flip*N1x/normal_length[0];
	force[0][1] = force_magn/2*flip*N1y/normal_length[0];
	force[0][2] = force_magn/2*flip*N1z/normal_length[0];
	force[1][0] = force_magn/2*flip*N2x/normal_length[1];
	force[1][1] = force_magn/2*flip*N2y/normal_length[1];
	force[1][2] = force_magn/2*flip*N2z/normal_length[1];

}

// #####################################################
// #####################################################

void ObjectClass::fileInit(void) {

	char file[128];
	sprintf(file,"%s%s%d.log",OUTPUT_DATA_DIR,OBJECT_LOG_FILE,object_index);
	Objectfile.open(file);

    Objectfile.width(15);
    Objectfile << left << "Iteration";
    Objectfile.width(15);
    Objectfile << left << "#Nonnice Vert";
    Objectfile.width(15);
    Objectfile << left << "Pos. Error";
    Objectfile.width(15);
    Objectfile << left << "Force";
    Objectfile.width(15);
    Objectfile << left << "Energy" << endl;
}

// #####################################################
// #####################################################

void ObjectClass::fileOutit(void) {

	Objectfile.close();

}

// #####################################################
// #####################################################

void ObjectClass::InitVertices(void) {

	int i;

	// initialize vertex indices
	for (i=0;i<num_vertices;i++) {
		vertices[i].vertex_index = i;
//		if (vertex_print) {
		if (vertex_print && (i==2960)) {
			vertices[i].fileInit(object_index);
		}
	}

}

// #####################################################
// #####################################################

void ObjectClass::OutitVertices(void) {

	int i;

	// close vertex indices
	for (i=0;i<num_vertices;i++) {
		vertices[i].fileOutit();
	}

}

// ######################################
// ######################################

void ObjectClass::makeObject(int number_of_vertices, int number_of_polygons) {

    num_vertices = number_of_vertices;
    num_polygons = number_of_polygons;
    vertices = new VertexClass[number_of_vertices];
	if (vertices == NULL) { 
		cout << "Not enough memory for VertexClass\n";
		cout.flush();
		exit(1);
	}
	
}

// ######################################
// ######################################

void ObjectClass::addOriginalVertexCoordinates(int current_vertex,double *coordinates) {

    if (current_vertex >= num_vertices) {

        cout << "Error: Number of vertices exceeded in object# " << object_index
            << ". Number of vertices = " << num_vertices << "\n";
        cout << "Offending vertex# = " << current_vertex << ", coordinates = "
            << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << "\n";

    } else {
        vertices[current_vertex].original_coordinates[0] = coordinates[0];
        vertices[current_vertex].original_coordinates[1] = coordinates[1];
        vertices[current_vertex].original_coordinates[2] = coordinates[2];
    }

}

// ######################################
// ######################################

void ObjectClass::initializeNewCoordinatesAndDevDist(void) {

	int j, i;

    // for each vertex in parent object
    for (j = 0;j < num_vertices;j++) {
		for(i=0;i<3;i++) {
			vertices[j].new_coordinates[i] = vertices[j].original_coordinates[i];
		}
		vertices[j].deviation_distance = 0;
		vertices[j].deviation_distance_old = 0;
		vertices[j].new_edge_lengths = NULL;
		vertices[j].intersection_force[0] = 0;
		vertices[j].intersection_force[1] = 0;
		vertices[j].intersection_force[2] = 0;
	}
}

// #####################################################
// #####################################################

void ObjectClass::getPolygonCoordinates(int indices[3],double* vertex[3]) {

	//get coordinates of each vertex index
	vertex[0] = &vertices[indices[0]].new_coordinates[0];
	vertex[1] = &vertices[indices[1]].new_coordinates[0];
	vertex[2] = &vertices[indices[2]].new_coordinates[0];

}

// #####################################################
// #####################################################

void ObjectClass::findAdjacentPolygonsAndVertices(PolygonClass *polygons, 
												int polygon_count, int *polygons_array) {

    int cvi, m, k, next_polygon, start, end, i;
	int v[3];//, new_size;

	// cvi = current_vertex_index

	// restrict the set of polygons to search 
	// to those polygons on current object
	start = 0;
	for (i=0;i<object_index;i++) {
		start = start + polygons_array[i];
	}
	end = start + polygons_array[object_index];

    ////////// loop through all vertices structure elememts in parent object //////////
    // for each vertex in parent object
    for (cvi = 0;cvi < num_vertices;cvi++) {

        ////////// initialize variables //////////
        next_polygon = 0;
        vertices[cvi].next_adjacent_vertex = 0;
        vertices[cvi].adjacent_polygons_index = new int[NUMBER_OF_ADJACENT_POLYGONS];
		if (vertices[cvi].adjacent_polygons_index == NULL) { 
			cout << "Not enough memory for vertices[cvi].adjacent_polygons_index "
					<< "in ObjectClass::findAdjacentPolygonsAndVertices\n";
			cout.flush();
			exit(1);
		}
        vertices[cvi].next_adjacent_polygon = 0;
        vertices[cvi].max_adjacent_polygon = NUMBER_OF_ADJACENT_POLYGONS;
        vertices[cvi].adjacent_vertices_index = new int[NUMBER_OF_ADJACENT_VERTICES];
		if (vertices[cvi].adjacent_vertices_index == NULL) { 
			cout << "Not enough memory for vertices[cvi].adjacent_vertices_index "
					<< "in ObjectClass::findAdjacentPolygonsAndVertices\n";
			cout.flush();
			exit(1);
		}
        vertices[cvi].next_adjacent_vertex = 0;
        vertices[cvi].max_adjacent_vertex = NUMBER_OF_ADJACENT_VERTICES;

	}

	////////// loop through all polygons structure elements //////////

	// for each polygon structure in active set of master polygons array
	for (k=start;k<end;k++) {

		// get polygon vertex indices
		for (i=0;i<3;i++) {
			v[i] = polygons[k].vertices[i];
		}

		// add vertices as adjacencies to each other
//		addSortUniqueArrayElement(vertices[v[0]].adjacent_vertices_index, vertices[v[0]].next_adjacent_vertex,
//								vertices[v[0]].max_adjacent_vertex, INCREMENT_OF_ADJACENT_VERTICES, v[1]);
//		addSortUniqueArrayElement(vertices[v[0]].adjacent_vertices_index, vertices[v[0]].next_adjacent_vertex,
//								vertices[v[0]].max_adjacent_vertex, INCREMENT_OF_ADJACENT_VERTICES, v[2]);
//		addSortUniqueArrayElement(vertices[v[1]].adjacent_vertices_index, vertices[v[1]].next_adjacent_vertex,
//								vertices[v[1]].max_adjacent_vertex, INCREMENT_OF_ADJACENT_VERTICES, v[0]);
//		addSortUniqueArrayElement(vertices[v[1]].adjacent_vertices_index, vertices[v[1]].next_adjacent_vertex,
//								vertices[v[1]].max_adjacent_vertex, INCREMENT_OF_ADJACENT_VERTICES, v[2]);
//		addSortUniqueArrayElement(vertices[v[2]].adjacent_vertices_index, vertices[v[2]].next_adjacent_vertex,
//								vertices[v[2]].max_adjacent_vertex, INCREMENT_OF_ADJACENT_VERTICES, v[0]);
//		addSortUniqueArrayElement(vertices[v[2]].adjacent_vertices_index, vertices[v[2]].next_adjacent_vertex,
//								vertices[v[2]].max_adjacent_vertex, INCREMENT_OF_ADJACENT_VERTICES, v[1]);

		addArrayElement(vertices[v[0]].adjacent_vertices_index, vertices[v[0]].next_adjacent_vertex,
								vertices[v[0]].max_adjacent_vertex, INCREMENT_OF_ADJACENT_VERTICES, v[1]);
		addArrayElement(vertices[v[0]].adjacent_vertices_index, vertices[v[0]].next_adjacent_vertex,
								vertices[v[0]].max_adjacent_vertex, INCREMENT_OF_ADJACENT_VERTICES, v[2]);
		addArrayElement(vertices[v[1]].adjacent_vertices_index, vertices[v[1]].next_adjacent_vertex,
								vertices[v[1]].max_adjacent_vertex, INCREMENT_OF_ADJACENT_VERTICES, v[0]);
		addArrayElement(vertices[v[1]].adjacent_vertices_index, vertices[v[1]].next_adjacent_vertex,
								vertices[v[1]].max_adjacent_vertex, INCREMENT_OF_ADJACENT_VERTICES, v[2]);
		addArrayElement(vertices[v[2]].adjacent_vertices_index, vertices[v[2]].next_adjacent_vertex,
								vertices[v[2]].max_adjacent_vertex, INCREMENT_OF_ADJACENT_VERTICES, v[0]);
		addArrayElement(vertices[v[2]].adjacent_vertices_index, vertices[v[2]].next_adjacent_vertex,
								vertices[v[2]].max_adjacent_vertex, INCREMENT_OF_ADJACENT_VERTICES, v[1]);

		// add polygon index as adjacent polygon to vertices
		addArrayElement(vertices[v[0]].adjacent_polygons_index, vertices[v[0]].next_adjacent_polygon,
						vertices[v[0]].max_adjacent_polygon, INCREMENT_OF_ADJACENT_POLYGONS, k);
		addArrayElement(vertices[v[1]].adjacent_polygons_index, vertices[v[1]].next_adjacent_polygon,
						vertices[v[1]].max_adjacent_polygon, INCREMENT_OF_ADJACENT_POLYGONS, k);
		addArrayElement(vertices[v[2]].adjacent_polygons_index, vertices[v[2]].next_adjacent_polygon,
						vertices[v[2]].max_adjacent_polygon, INCREMENT_OF_ADJACENT_POLYGONS, k);

	}

	// sort arrays
	qsort(vertices[v[0]].adjacent_vertices_index,vertices[v[0]].next_adjacent_vertex,sizeof(int),compare_int);
	qsort(vertices[v[1]].adjacent_vertices_index,vertices[v[1]].next_adjacent_vertex,sizeof(int),compare_int);
	qsort(vertices[v[2]].adjacent_vertices_index,vertices[v[2]].next_adjacent_vertex,sizeof(int),compare_int);

//	Quicksort(vertices[v[0]].adjacent_vertices_index,0,vertices[v[0]].next_adjacent_vertex-1);
//	Quicksort(vertices[v[1]].adjacent_vertices_index,0,vertices[v[1]].next_adjacent_vertex-1);
//	Quicksort(vertices[v[2]].adjacent_vertices_index,0,vertices[v[2]].next_adjacent_vertex-1);

	// extract unique elements
	vertices[v[0]].adjacent_vertices_index = uniqueArray(vertices[v[0]].adjacent_vertices_index,
														vertices[v[0]].next_adjacent_vertex,
														vertices[v[0]].max_adjacent_polygon);
	vertices[v[1]].adjacent_vertices_index = uniqueArray(vertices[v[1]].adjacent_vertices_index,
														vertices[v[1]].next_adjacent_vertex,
														vertices[v[1]].max_adjacent_polygon);
	vertices[v[2]].adjacent_vertices_index = uniqueArray(vertices[v[2]].adjacent_vertices_index,
														vertices[v[2]].next_adjacent_vertex,
														vertices[v[2]].max_adjacent_polygon);


}

// #####################################################
// #####################################################

void ObjectClass::findNeighborhoodVertices(void) {

	int i, j, k, num_edges, radius;
	double mean_edge_length;
	int* list_last_last = NULL;
	int num_list_last_last;
	int max_list_last_last;
	int* list_last = NULL;
	int num_list_last;
	int max_list_last;
	int* list_new = NULL;
	int num_list_new;
	int max_list_new;
	int* temp1 = NULL;
	int num_temp1;
	int max_temp1;
	int* temp2 = NULL;
	int num_temp2;
	int max_temp2;
	int a, b;

	int NUMBER_LIST = 50;
	int INCREMENT_OF_LIST = 50;

	////////// compute mean edge length in object //////////
	mean_edge_length = 0;
	num_edges = 0;
	// for each vertex in object
    for (i=0;i<num_vertices;i++) {
		// for each adjacent edge
	    for (j=0;j<vertices[i].next_original_edge;j++) {
	        mean_edge_length += vertices[i].original_edge_lengths[j];
			num_edges++;
		}
	}
	// compute mean edge length
	mean_edge_length = mean_edge_length/num_edges; //nm

	cout << "mean_edge_length = " << mean_edge_length << endl;

	////////// compute neighborhood radius in edge lengths //////////
	radius = (int) ceil(NEIGHBORHOOD_RADIUS/mean_edge_length);

	////////// collect neighborhood vertices //////////
	// for each vertex in object
    for (i=0;i<num_vertices;i++) {

//		cout << "find neighbors: vertex " << i << " of " << num_vertices << endl;

		// initialize vectors
		vertices[i].neighborhood_vertices_index = new int[NUMBER_OF_NEIGHBORS];
		if (vertices[i].neighborhood_vertices_index == NULL) { 
			cout << "Not enough memory for neighborhood_vertices_index in ObjectClass::findNeighborhoodVertices\n";
			cout.flush();
			exit(1);
		}
		vertices[i].max_neighborhood_vertices_index = NUMBER_OF_NEIGHBORS;
		vertices[i].next_neighbor_vertex = 0;
		// add self to neighborhood
//		vertices[i].addSortUniqueNeighbor(i);
//		vertices[i].addNeighbor(i);
		addSortUniqueArrayElement(vertices[i].neighborhood_vertices_index, vertices[i].next_neighbor_vertex,
									vertices[i].max_neighborhood_vertices_index, INCREMENT_OF_NEIGHBORS,i);

		if (list_last_last != NULL) {
			delete [] list_last_last;
		}
		if (list_last != NULL) {
			delete [] list_last;
		}
		if (list_new != NULL) {
			delete [] list_new;
		}
		if (temp1 != NULL) {
			delete [] temp1;
		}
		if (temp2 != NULL) {
			delete [] temp2;
		}
		list_last_last = new int [NUMBER_LIST];
		if (list_last_last == NULL) { 
			cout << "Not enough memory for list_last_last (first) in ObjectClass::findNeighborhoodVertices\n";
			cout.flush();
			exit(1);
		}
		num_list_last_last = 0;
		max_list_last_last = NUMBER_LIST;
		list_last = new int [NUMBER_LIST];
		if (list_last == NULL) { 
			cout << "Not enough memory for list_last (first) in ObjectClass::findNeighborhoodVertices\n";
			cout.flush();
			exit(1);
		}
		num_list_last = 0;
		max_list_last = NUMBER_LIST;
		list_new = new int [NUMBER_LIST];
		if (list_new == NULL) { 
			cout << "Not enough memory for list_new (first) in ObjectClass::findNeighborhoodVertices\n";
			cout.flush();
			exit(1);
		}
		num_list_new = 0;
		max_list_new = NUMBER_LIST;
		temp1 = new int [NUMBER_LIST];
		if (temp1 == NULL) { 
			cout << "Not enough memory for temp1 (first) in ObjectClass::findNeighborhoodVertices\n";
			cout.flush();
			exit(1);
		}
		num_temp1 = 0;
		max_temp1 = NUMBER_LIST;
		temp2 = new int [NUMBER_LIST];
		if (temp2 == NULL) { 
			cout << "Not enough memory for temp2 (first) in ObjectClass::findNeighborhoodVertices\n";
			cout.flush();
			exit(1);
		}
		num_temp2 = 0;
		max_temp2 = NUMBER_LIST;
	
		// initialize list_last	
		addArrayElement(list_last,num_list_last,max_list_last,INCREMENT_OF_LIST,i);

		////// build neighborhood //////
		for (j=0;j<radius;j++) {

			//// build list_new ////
			// for each last vertex
			for (a=0;a<num_list_last;a++) {
				// for each adjacent vertex of last
			    for (k=0;k<vertices[list_last[a]].next_adjacent_vertex;k++) {
					// add adjacent vertex to new
					addSortUniqueArrayElement(list_new,num_list_new,max_list_new,INCREMENT_OF_LIST,
									vertices[list_last[a]].adjacent_vertices_index[k]);
				}
			}

			// extract elements from new that are not in last_last
			// and add to temp1
			if (j) {
				a = 0;
				b = 0;
				while ((a < num_list_last_last) && (b < num_list_new)) {
					if      (list_new[b] < list_last_last[a]) {
						addArrayElement(temp1,num_temp1,max_temp1,INCREMENT_OF_LIST,list_new[b] );
						b++;
					} else if (list_new[b] > list_last_last[a]) {a++;}
					else    {a++;b++; }
				}
				while (b < num_list_new) {
					addArrayElement(temp1,num_temp1,max_temp1,INCREMENT_OF_LIST,list_new[b] );
					b++;
				}
			} else {
				// add new to temp1
				for (a=0;a<num_list_new;a++) {
					addArrayElement(temp1,num_temp1,max_temp1,INCREMENT_OF_LIST,list_new[a] );
				}
	
			}

			// extract elements from new (i.e. temp1) that are not in last
			// and add to temp2
			if (j) {
				a = 0;
				b = 0;
				while ((a < num_list_last) && (b < num_temp1)) {
					if      (temp1[b] < list_last[a]) {
						addArrayElement(temp2,num_temp2,max_temp2,INCREMENT_OF_LIST,temp1[b] );
						b++;
					} else if (temp1[b] > list_last[a]) {a++;}
					else    {a++;b++; }
				}
				while (b < num_temp1) {
					addArrayElement(temp2,num_temp2,max_temp2,INCREMENT_OF_LIST,temp1[b] );
					b++;
				}
			} else {
				// add new to temp2
				for (a=0;a<num_list_new;a++) {
					addArrayElement(temp2,num_temp2,max_temp2,INCREMENT_OF_LIST,list_new[a] );
				}
	
			}

			// add new (i.e. temp2) to neighborhood
			for (a=0;a<num_temp2;a++) {
				// add adjacent vertex to new
//				vertices[i].addSortUniqueNeighbor(temp2[a]);
//				vertices[i].addNeighbor(temp2[a]);
				addSortUniqueArrayElement(vertices[i].neighborhood_vertices_index, 
											vertices[i].next_neighbor_vertex,
											vertices[i].max_neighborhood_vertices_index, 
											INCREMENT_OF_NEIGHBORS,temp2[a]);
			}

			//// set last_last to last ////
			delete[] list_last_last;
			list_last_last = new int[num_list_last];
			if (list_last_last == NULL) { 
				cout << "Not enough memory for list_last_last (second) in ObjectClass::findNeighborhoodVertices\n";
				cout.flush();
				exit(1);
			}
			max_list_last_last = num_list_last;
			num_list_last_last = 0;
			for (a=0;a<num_list_last;a++) {
				// add adjacent vertex to new
				list_last_last[a] = list_last[a];
				num_list_last_last++;
			}
			
			//// set last to new (i.e. temp2) ////
			delete[] list_last;
			list_last = new int[num_temp2];
			if (list_last == NULL) { 
				cout << "Not enough memory for list_last (second) in ObjectClass::findNeighborhoodVertices\n";
				cout.flush();
				exit(1);
			}
			max_list_last = num_temp2;
			num_list_last = 0;
			for (a=0;a<num_temp2;a++) {
				// add adjacent vertex to new
				list_last[a] = temp2[a];
				num_list_last++;
			}
			
			//// clear new ////
			delete [] list_new;
			delete [] temp1;
			delete [] temp2;
			list_new = new int [NUMBER_LIST];
			if (list_new == NULL) { 
				cout << "Not enough memory for list_new (second) in ObjectClass::findNeighborhoodVertices\n";
				cout.flush();
				exit(1);
			}
			max_list_new = NUMBER_LIST;
			num_list_new = 0;
			temp1 = new int [NUMBER_LIST];
			if (temp1 == NULL) { 
				cout << "Not enough memory for temp1 (second) in ObjectClass::findNeighborhoodVertices\n";
				cout.flush();
				exit(1);
			}
			max_temp1 = NUMBER_LIST;
			num_temp1 = 0;
			temp2 = new int [NUMBER_LIST];
			if (temp2 == NULL) { 
				cout << "Not enough memory for temp2 (second) in ObjectClass::findNeighborhoodVertices\n";
				cout.flush();
				exit(1);
			}
			max_temp2 = NUMBER_LIST;
			num_temp2 = 0;
		}

		cout << "vertices[i].next_neighbor_vertex = " << vertices[i].next_neighbor_vertex << endl;

	}

	delete [] list_last_last;
	delete [] list_last;
	delete [] list_new;
	delete [] temp1;
	delete [] temp2;

}

// ######################################
// ######################################

void ObjectClass::computeEdgeLengths(void) {

	int current_vertex_index, current_adjacent_vertex, adjacent_index;
	double xC, yC, zC, xA, yA, zA, diffX, diffY, diffZ, distance;

    ////////// loop through all vertices structure elememts in parent object //////////
    // for each vertex in parent object
    for (current_vertex_index = 0;current_vertex_index < num_vertices;current_vertex_index++) {

		// clear original edge lengths vector
		// initialize variable
		vertices[current_vertex_index].original_edge_lengths 
				= new double [vertices[current_vertex_index].next_adjacent_vertex];
		if (vertices[current_vertex_index].original_edge_lengths == NULL) { 
			cout << "Not enough memory for vertices[current_vertex_index].original_edge_lengths "
				<< "in ObjectClass::computeEdgeLengths\n";
			cout.flush();
			exit(1);
		}
		vertices[current_vertex_index].next_original_edge = 0;

    	////////// loop through all adjacent vertices of current vertex //////////
    	// for each adjacent vertex of current vertex
	    for (current_adjacent_vertex = 0;
			current_adjacent_vertex < vertices[current_vertex_index].next_adjacent_vertex;
			current_adjacent_vertex++) {

			// current adjacent vertex index
			adjacent_index = vertices[current_vertex_index].adjacent_vertices_index[current_adjacent_vertex];

			// get coordinates of vertices
			// current vertex
//			xC = vertices[current_vertex_index].original_coordinates[0];
//			yC = vertices[current_vertex_index].original_coordinates[1];
//			zC = vertices[current_vertex_index].original_coordinates[2];
			// adjacent
//			xA = vertices[adjacent_index].original_coordinates[0];
//			yA = vertices[adjacent_index].original_coordinates[1];
//			zA = vertices[adjacent_index].original_coordinates[2];

			// compute distance
//			diffX = xC-xA;
//			diffY = yC-yA;
//			diffZ = zC-zA;
			diffX = vertices[current_vertex_index].original_coordinates[0]
					-vertices[adjacent_index].original_coordinates[0];
			diffY = vertices[current_vertex_index].original_coordinates[1]
					-vertices[adjacent_index].original_coordinates[1];
			diffZ = vertices[current_vertex_index].original_coordinates[2]
					-vertices[adjacent_index].original_coordinates[2];
			distance = sqrt( pow(diffX,2) + pow(diffY,2) + pow(diffZ,2) );

			// add length to mesh
			vertices[current_vertex_index].original_edge_lengths[current_adjacent_vertex] = distance;
			vertices[current_vertex_index].next_original_edge++;

			if (print_flag) {
				cout << "computeEdgeLengths: "
					<< "object# = " << object_index 
					<< ", current vertex index = " << current_vertex_index 
					<< ", current adjacent vertex index  = " << adjacent_index 
					<< ", distance  = " << distance 
					<< ", C  = " << vertices[current_vertex_index].original_coordinates[0] << " " 
								<< vertices[current_vertex_index].original_coordinates[1] << " " 
								<< vertices[current_vertex_index].original_coordinates[2]
					<< ", A  = " << vertices[adjacent_index].original_coordinates[0] << " " 
								<< vertices[adjacent_index].original_coordinates[1] << " " 
								<< vertices[adjacent_index].original_coordinates[2]
					<< "\n";
			}

			// increment next original edge
//			vertices[current_vertex_index].next_original_edge++;
		}
	}
}

// ######################################
// ######################################

void ObjectClass::computeNewEdgeLengths(void) {

	int current_vertex_index, current_adjacent_vertex, adjacent_index;
	double xC, yC, zC, xA, yA, zA, diffX, diffY, diffZ, distance;

    ////////// loop through all vertices structure elememts in parent object //////////
    // for each vertex in parent object
    for (current_vertex_index = 0;current_vertex_index < num_vertices;current_vertex_index++) {

		// clear new edge lengths vector
		// initialize variable
		if (vertices[current_vertex_index].new_edge_lengths != NULL) {
			delete [] vertices[current_vertex_index].new_edge_lengths;
		}
		vertices[current_vertex_index].new_edge_lengths 
				= new double [vertices[current_vertex_index].next_adjacent_vertex];
		if (vertices[current_vertex_index].new_edge_lengths == NULL) { 
			cout << "Not enough memory for vertices[current_vertex_index].new_edge_lengths "
				<< "in ObjectClass::computeNewEdgeLengths\n";
			cout.flush();
			exit(1);
		}
		vertices[current_vertex_index].next_new_edge = 0;

    	////////// loop through all adjacent vertices of current vertex //////////
    	// for each adjacent vertex of current vertex
	    for (current_adjacent_vertex = 0;
			current_adjacent_vertex < vertices[current_vertex_index].next_adjacent_vertex;
			current_adjacent_vertex++) {

			// current adjacent vertex index
			adjacent_index = vertices[current_vertex_index]
							.adjacent_vertices_index[current_adjacent_vertex];

			// get coordinates of vertices
			// current vertex
			xC = vertices[current_vertex_index].new_coordinates[0];
			yC = vertices[current_vertex_index].new_coordinates[1];
			zC = vertices[current_vertex_index].new_coordinates[2];
			// adjacent
			xA = vertices[adjacent_index].new_coordinates[0];
			yA = vertices[adjacent_index].new_coordinates[1];
			zA = vertices[adjacent_index].new_coordinates[2];

			// compute distance
			diffX = xC-xA;
			diffY = yC-yA;
			diffZ = zC-zA;
			distance = sqrt( diffX*diffX + diffY*diffY + diffZ*diffZ );

			// add length to mesh
			vertices[current_vertex_index].new_edge_lengths[current_adjacent_vertex] = distance;
			vertices[current_vertex_index].next_new_edge++;

			if (print_flag) {
				cout << "computeNewEdgeLengths: "
					<< "object# = " << object_index 
					<< ", current vertex index = " << current_vertex_index 
					<< ", current adjacent vertex index  = " << adjacent_index 
					<< ", distance  = " << distance 
					<< ", C  = " << xC << " " << yC << " " << zC 
					<< ", A  = " << xA << " " << yA << " " << zA 
					<< "\n";
			}

			// increment next original edge
			vertices[current_vertex_index].next_new_edge++;
		}
	}
}

// #####################################################
// #####################################################

void ObjectClass::updateObjectLog(int iteration) {

    Objectfile.width(15);
    Objectfile << left << iteration;
    Objectfile.width(15);
    Objectfile << left << nonnice;
    Objectfile.width(15);
    Objectfile << left << position_error;
    Objectfile.width(15);
    Objectfile << left << force;
    Objectfile.width(15);
    Objectfile << left << energy << endl;
}

// #####################################################
// #####################################################

void ObjectClass::updateVertexLog(int iteration) {

	int i;

	// initialize vertex indices
	for (i=0;i<num_vertices;i++) {
		if (vertex_print && i==2960) {
//		if (vertex_print) {
			vertices[i].updateLog(iteration);
		}
	}

}

// #####################################################
// #####################################################

void ObjectClass::computeGlobalParams(void) {

	int i;
	double raw_force;

	// initialize object variables
	nonnice = 0;
	position_error = 0;
	force = 0;
	energy = 0;

	// for each vertex
	for (i=0;i<num_vertices;i++) {
		if (!vertices[i].nice_data.nice) {nonnice++;}
		if (vertices[i].candidate) {
			position_error += vertices[i].position_error;
			raw_force = sqrt( vertices[i].force[0]*vertices[i].force[0]+
							vertices[i].force[1]*vertices[i].force[1]+
							vertices[i].force[2]*vertices[i].force[2] );
			force += raw_force;
			energy += vertices[i].energy;
		}
	}	
}

// ######################################
// ######################################

void ObjectClass::initEdges(void) {

	num_edges = num_vertices+num_polygons-2+EDGE_COMPENSATION;

    edges = new EdgeClass[num_edges];
	if (edges == NULL) { 
		cout << "Not enough memory for EdgeClass\n";
		cout.flush();
		exit(1);
	}
	
}

// #####################################################
// #####################################################

void ObjectClass::buildEdges(PolygonClass *polygons, int polygon_count, int *polygons_array) {

	int next_edge, i, k, v1, v2, v3, v[3], start, end, j;
	int triplets[3][3], index, new_size;
	bool found_flag;

	// declare pairs
	triplets[0][0] = 0;
	triplets[0][1] = 1;
	triplets[0][2] = 2;
	triplets[1][0] = 1;
	triplets[1][1] = 2;
	triplets[1][2] = 0;
	triplets[2][0] = 2;
	triplets[2][1] = 0;
	triplets[2][2] = 1;

	next_edge = 0;

	// initialize valid
	for (i=0;i<num_edges;i++) {
		edges[i].valid = false;
    }

	// initialize vertex vector
	for (i=0;i<num_vertices;i++) {
		vertices[i].edges_indices = new int[NUMBER_OF_EDGES_VERTICES];
		if (vertices[i].edges_indices == NULL) { 
			cout << "Not enough memory for vertices[i].edges_indices in ObjectClass::buildEdges\n";
			cout.flush();
			exit(1);
		}
		vertices[i].num_edges_indices = 0;
		vertices[i].max_edges_indices = NUMBER_OF_EDGES_VERTICES;
    }


	////////// loop through all polygons structure elements //////////

	// restrict the set of polygons to search 
	// to those polygons on current object
	start = 0;
	for (i=0;i<object_index;i++) {
		start = start + polygons_array[i];
	}
	end = start + polygons_array[object_index];

	// for each polygon structure in active set of master polygons array
	for (k=start;k<end;k++) {

		// get polygon vertex indices
		for (i=0;i<3;i++) {
			v[i] = polygons[k].vertices[i];
		}

		// for each vertex pair
		for (i=0;i<3;i++) {

			v1 = v[triplets[i][0]];
			v2 = v[triplets[i][1]];
			v3 = v[triplets[i][2]];

			// check if pair is already in list
			found_flag = false;
			j = 0;


			while (!found_flag && (j < vertices[v1].num_edges_indices)) {
				if ((
					(edges[vertices[v1].edges_indices[j]].edge_vertex_indices[0] == v1) && 
					(edges[vertices[v1].edges_indices[j]].edge_vertex_indices[1] == v2)) || 
					((edges[vertices[v1].edges_indices[j]].edge_vertex_indices[0] == v2) && 
					(edges[vertices[v1].edges_indices[j]].edge_vertex_indices[1] == v1))  ) 
					{found_flag = true;}
				j++;
			}

			// if not found, then add pair to class and to list
			if (!found_flag) {

				// check size of array
				if (next_edge >= num_edges) {
					cout << "ERROR: increase size of EDGE_COMPENSATION\n";
					exit(1);
				}

				// add pair to class
				edges[next_edge].edge_vertex_indices[0] = v1;
				edges[next_edge].edge_vertex_indices[1] = v2;
				edges[next_edge].other_vertex_indices[0] = v3;
				edges[next_edge].object_index = object_index;
				edges[next_edge].polygon_indices[0] = k;
				edges[next_edge].valid = true;
				// add edge to vertices
				addArrayElement(vertices[v1].edges_indices, vertices[v1].num_edges_indices,
								vertices[v1].max_edges_indices, INCREMENT_OF_EDGES_VERTICES, next_edge);
				addArrayElement(vertices[v2].edges_indices, vertices[v2].num_edges_indices,
								vertices[v2].max_edges_indices, INCREMENT_OF_EDGES_VERTICES, next_edge);
				next_edge++;
			} else {
			// add second polygon index
				j--;
				edges[vertices[v1].edges_indices[j]].other_vertex_indices[1] = v3;
				edges[vertices[v1].edges_indices[j]].polygon_indices[1] = k;
			}
		}
    }
}

// #####################################################
// #####################################################

void ObjectClass::initEdgeAngles(PolygonClass *polygons) {

	int i;

	for (i=0;i<num_edges;i++) {
		if (edges[i].valid) {
			edges[i].initAngles(polygons);
		}
    }
}
	
// #####################################################
// #####################################################

void ObjectClass::computeEdgeAngleForceEnergy(PolygonClass *polygons) {

	int i;

	for (i=0;i<num_edges;i++) {

		if (edges[i].valid) {

			// initialize vertex coords of edge_vertex_indices[0]
			edges[i].edge_vertex_coords[0] = vertices[edges[i].edge_vertex_indices[0]].new_coordinates[0];
			edges[i].edge_vertex_coords[1] = vertices[edges[i].edge_vertex_indices[0]].new_coordinates[1];
			edges[i].edge_vertex_coords[2] = vertices[edges[i].edge_vertex_indices[0]].new_coordinates[2];

			// initialize vertex coords of other_vertex_indices[0]
			edges[i].other_vertex_coords[0] = vertices[edges[i].other_vertex_indices[0]].new_coordinates[0];
			edges[i].other_vertex_coords[1] = vertices[edges[i].other_vertex_indices[0]].new_coordinates[1];
			edges[i].other_vertex_coords[2] = vertices[edges[i].other_vertex_indices[0]].new_coordinates[2];

			// compute new angle, force, energy
			edges[i].computeAngleForceEnergy(polygons);

		}

    }
}

// #####################################################
// #####################################################

ObjectClass::~ObjectClass(void) {

	delete[] vertices;
	delete[] edges;

}

// #####################################################
// #####################################################

void ContainerClass::init(void) {

	// allocate memory for a specific number (full->object_count)
	// of array elements of type ObjectClass
	cout << "ContainerClass::init : object_count = " << object_count << "\n";
	objects = new ObjectClass[object_count];
	if (objects == NULL) { 
		cout << "Not enough memory for ObjectClass\n";
		cout.flush();
		exit(1);
	}

}

// #####################################################
// #####################################################

void ContainerClass::initObjects(int *vertices_array, int *polygons_array) {

	int k;

	// for each object in objects array allocate memory 
	// for an array called vertices of size scan.vertices[k]
	// with each array element being a class of type vertex.
	for (k=0;k<object_count;k++) {
		cout << "ContainerClass::initObjects : object = " << k
			<< ", vertices = " << vertices_array[k] 
			<< ", polygons = " << polygons_array[k] << "\n";
		objects[k].makeObject(vertices_array[k],polygons_array[k]);
		objects[k].object_index = k;
	}


}

// #####################################################
// #####################################################

void ContainerClass::initObjectAndVertexLogFiles(void) {

	int k;

	// initialize log files
	for (k=0;k<object_count;k++) {
		objects[k].fileInit();
		objects[k].InitVertices();
	}

}

// #####################################################
// #####################################################

void ContainerClass::OutitObjectAndVertexLogFiles(void) {

	int k;

	// close log files
	for (k=0;k<object_count;k++) {
		objects[k].fileOutit();
		if (vertex_print) {
			objects[k].OutitVertices();
		}
	}

}

// #####################################################
// #####################################################

void ContainerClass::initCoordsAndDist(void) {

	int k;

	// initialize new coordinates
	for (k=0;k<object_count;k++) {
		objects[k].initializeNewCoordinatesAndDevDist();
	}

	// check for frozen objects
//	if (FROZEN_OBJECTS) {
		// find frozen object
//		for (k=0;k<object_count;k++) {
//			// if found
//			if (!strcmp(objects[k].name,FROZEN_LIST)){
				// set frozen to true
//				objects[k].frozen = true;
//			} else {
				// set frozen to false
//				objects[k].frozen = false;
//			}
//		}
//	}

}

// #####################################################
// #####################################################

void ContainerClass::initFindAdj(PolygonClass* polygons,int polygon_count,
								int *polygons_array) {

	int k;

	// store adjacent polygons and nodes
	for (k=0;k<object_count;k++) {
		objects[k].findAdjacentPolygonsAndVertices(polygons,polygon_count,polygons_array);
	}

}

// #####################################################
// #####################################################

void ContainerClass::initFindNeighbor(void) {

	int k;

	// store adjacent polygons and nodes
	for (k=0;k<object_count;k++) {
		objects[k].findNeighborhoodVertices();
	}

}

// #####################################################
// #####################################################

void ContainerClass::computeEdges(void) {

	int k;

	// initialize new coordinates
	for (k=0;k<object_count;k++) {
		objects[k].computeEdgeLengths();
	}

}

// #####################################################
// #####################################################

void ContainerClass::computeNewEdges(void) {

	int k;

	// initialize new coordinates
	for (k=0;k<object_count;k++) {
		objects[k].computeNewEdgeLengths();
	}

}

// #####################################################
// #####################################################

void ContainerClass::updateObjectAndVertexLog(int iteration) {

	int k;

	for (k=0;k<object_count;k++) {
		objects[k].updateObjectLog(iteration);
		if (k==0) {
		objects[k].updateVertexLog(iteration);
		}
	}

}

// #####################################################
// #####################################################

void ContainerClass::fileInit(void) {

	char file[128];
	sprintf(file,"%s%s",OUTPUT_DATA_DIR,CONT_LOG_FILE);
	Containerfile.open(file);

    Containerfile.width(15);
    Containerfile << left << "Iteration";
    Containerfile.width(15);
    Containerfile << left << "#Nonnice Vert";
    Containerfile.width(15);
    Containerfile << left << "Pos. Error";
    Containerfile.width(15);
    Containerfile << left << "Force";
    Containerfile.width(15);
    Containerfile << left << "Energy" << endl;

}

// #####################################################
// #####################################################

void ContainerClass::fileOutit(void) {

	Containerfile.close();

}

// #####################################################
// #####################################################

void ContainerClass::updateFile(int iteration) {

    Containerfile.width(15);
    Containerfile << left << iteration;
    Containerfile.width(15);
    Containerfile << left << nonnice;
    Containerfile.width(15);
    Containerfile << left << position_error;
    Containerfile.width(15);
    Containerfile << left << force;
    Containerfile.width(15);
    Containerfile << left << energy << endl;

}

// #####################################################
// #####################################################

void ContainerClass::computeGlobalParams(void) {

	int i;

	// initialize container variables
	nonnice = 0;
	position_error = 0;
	force = 0;
	energy = 0;

	// for each object
	for (i=0;i<object_count;i++) {
		nonnice += objects[i].nonnice;
		position_error += objects[i].position_error;
		force += objects[i].force;
		energy += objects[i].energy;
	}

}

// #####################################################
// #####################################################

void ContainerClass::initObjectEdges(PolygonClass *polygons, int polygon_count, int *polygons_array) {

	int k;

	// initialize new coordinates
	for (k=0;k<object_count;k++) {
		objects[k].initEdges();
		objects[k].buildEdges(polygons,polygon_count,polygons_array);
	}

}

// #####################################################
// #####################################################

void ContainerClass::initObjectEdgeAngles(PolygonClass *polygons) {

	int k;

	// initialize new coordinates
	for (k=0;k<object_count;k++) {
		objects[k].initEdgeAngles(polygons);
	}

}

// #####################################################
// #####################################################

void ContainerClass::computeObjectEdgeAngleForceEnergy(PolygonClass* polygons) {

	int k;

	// store adjacent polygons and nodes
	for (k=0;k<object_count;k++) {
		objects[k].computeEdgeAngleForceEnergy(polygons);
	}

}

// #####################################################
// #####################################################

void ContainerClass::writeVertexClosestDistanceToFile(void) {

	char file[128];
	int i, j;

	double xmax = 6571;
	double xmin = 5171;
	double ymax = 6302;
	double ymin = 4902;
	double zmax = 6550;
	double zmin = 5150;

	// create output filename
	sprintf(file,"%sworld_object_separation.dat",OUTPUT_DATA_DIR);

	// open output file
	ofstream newfile (file,ios::out);

	if(newfile.is_open()){
	
		// for each object
		for(i=0;i<object_count;i++) {

			// for each vertex
			for(j=0;j<objects[i].num_vertices;j++) {

//				if (
//					(objects[i].vertices[j].new_coordinates[0] > xmin) &&
//					(objects[i].vertices[j].new_coordinates[0] < xmax) &&
//					(objects[i].vertices[j].new_coordinates[1] > ymin) &&
//					(objects[i].vertices[j].new_coordinates[1] < ymax) &&
//					(objects[i].vertices[j].new_coordinates[2] > zmin) &&
//					(objects[i].vertices[j].new_coordinates[2] < zmax) ) {

					// print deviation_distance
					newfile << objects[i].vertices[j].deviation_distance << endl;
//				}
			}

		}

		newfile.close();

	}
}

// #####################################################
// #####################################################

ContainerClass::~ContainerClass(void) {

	delete[] objects;

}

// #####################################################
// #####################################################

void ScanClass::scanDir(void) {

    std::string str;
    std::string::size_type found;

	DIR *pdir;								// pointer to a directory data structure
	struct dirent *pent;					// pointer to dirent structure

	cout << "input directory = " << INPUT_DATA_DIR << "\n";

    pdir = opendir(INPUT_DATA_DIR);
    if (!pdir) {
        printf ("opendir() failure; could not open %s terminating",INPUT_DATA_DIR);
        exit(1);
    }
    errno = 0;
	num_files = 0;
    while ((pent=readdir(pdir))){
		// copy char array to string
		str = pent->d_name;
		// if file of typ *.mesh
		found = str.find(".mesh",0);
        if (found != string::npos) {
			// save filename
			strcpy(files[num_files],str.c_str());
			// update index
			num_files++;
			if (!(num_files<NUMBER_OF_INPUT_FILES)) {
				printf("ERROR: increase NUMBER_OF_INPUT_FILES; current value = %d\n",NUMBER_OF_INPUT_FILES);
				exit(1);
			}
			// print file found to screen
	        cout << "file found: " << str << "\n";
		}
    }
    if (errno) {
        printf ("readdir() failure; terminating");
        exit(1);
    }
    closedir(pdir);
}

// #####################################################
// #####################################################

void ScanClass::scanFiles (ContainerClass& container,PolygonClass *polygons, bool scan_flag) {

	int count;
    std::string str,name;
    std::string::size_type pos1;

	// initialize cumulative counts
	vertex_count = 0;
	polygon_count = 0;

	// for each input file
	for (count=0;count<num_files;count++) {

		// copy char array to string
        str = files[count];

		if (print_flag) {
			if (scan_flag) { cout << "scanobject " << count << " ";}
			else { cout << "getdata " << count << "\n";}
		}
		if (count == NUMBER_OF_OBJECTS) {
			cout << "Error: allowed number of objects exceeded in input file.\n" <<
					"Increase NUMBER_OF_OBJECTS.\n\n";
			exit(1);
		}

        // record object name
		pos1 = str.find(".",0);
		if (!(pos1 ==std::string::npos)) {
			name = str.substr(0,pos1);
			if (!scan_flag) { strcpy(container.objects[count].name,name.c_str()); }
		}

		// scan file
		scanFile(container,polygons,files[count],scan_flag,count);
	}

	// store object count in class variable
	object_count = count;

}


// #####################################################
// #####################################################

void ScanClass::scanFile (ContainerClass& container,PolygonClass *polygons, char filename[32], bool scan_flag, int count) {

	char line[128],file[128];
    int vertex_num, polygon_num;
	int indices[3];
	double coor[3];
    std::string str,coor_str1,coor_str2,coor_str3;
    std::string::size_type pos1,pos2,pos3,pos4,pos5,len;

	sprintf(file,"%s%s",INPUT_DATA_DIR,filename);

    std::ifstream inFile(file);

	polygon_num = 0;
    vertex_num = 0;

    // open input data file
    if (inFile.fail()) // if stream cannot be opened
    {
        cout << "Can't open " << file ; // display error msg and
        exit(1); // terminate abnormally
    }

	// foreach line in file
    while (inFile.getline(line,1024)) {
        str = line;
	
        // if vertex
        if (!str.find("Vertex",0)) {

			// count vertex
			vertex_num++;
			vertex_count++;

			if (!scan_flag) {

				// split string into Vertex, index and three coordinates
				len = str.length();
	            pos1 = str.find(" ",0);
	            pos2 = str.find(" ",pos1+1);
	            pos3 = str.find(" ",pos2+1);
	            pos4 = str.find(" ",pos3+1);
	            pos5 = str.find(" ",pos4+1);
				// check if extra space between index and first coordinate
				if (pos3 == pos2+1) {
					coor_str1 = str.substr(pos3+1,pos4-pos3-1);
					coor_str2 = str.substr(pos4+1,pos5-pos4-1);
					coor_str3 = str.substr(pos5+1,len-pos5-1);
				} else {
					coor_str1 = str.substr(pos2+1,pos3-pos2-1);
					coor_str2 = str.substr(pos3+1,pos4-pos3-1);
					coor_str3 = str.substr(pos4+1,len-pos4-1);
				}

				// add vertices to object
				coor[0] = atof(coor_str1.c_str());
				coor[1] = atof(coor_str2.c_str());
				coor[2] = atof(coor_str3.c_str());
           		container.objects[count].addOriginalVertexCoordinates(vertex_num-1,coor);
	
				if (print_flag) {
					cout << "getdata: object# " << count << ", vertex# " << vertex_num-1 
						<< ", vertex coordinates = "   
						<< coor[0] << " " << coor[1] << " " << coor[2] << "\n";
				}
			}
		} else if (!str.find("Face",0)) {

			polygon_num++;
			polygon_count++;

			if (!scan_flag) {

				// split string into Face, index and three indices
				len = str.length();
	            pos1 = str.find(" ",0);
	            pos2 = str.find(" ",pos1+1);
	            pos3 = str.find(" ",pos2+1);
	            pos4 = str.find(" ",pos3+1);
				pos5 = str.find(" ",pos4+1);
				// check if extra space between index and first vertex index
				if (pos3 == pos2+1) {
					coor_str1 = str.substr(pos3+1,pos4-pos3-1);
					coor_str2 = str.substr(pos4+1,pos5-pos4-1);
					coor_str3 = str.substr(pos5+1,len-pos5-1);
				} else {
					coor_str1 = str.substr(pos2+1,pos3-pos2-1);
					coor_str2 = str.substr(pos3+1,pos4-pos3-1);
					coor_str3 = str.substr(pos4+1,len-pos4-1);
				}
	
				// add vertex indices to polygon
				// subtract 1 from each index to account for 0-based array indexing in c++
				indices[0] = atoi(coor_str1.c_str());
				indices[1] = atoi(coor_str2.c_str());
				indices[2] = atoi(coor_str3.c_str());
				polygons[polygon_count-1].object_index = count;
				polygons[polygon_count-1].vertices[0] = indices[0]-1;
				polygons[polygon_count-1].vertices[1] = indices[1]-1;
				polygons[polygon_count-1].vertices[2] = indices[2]-1;

				if (print_flag) {
					cout << "getdata: object# "
						<< polygons[polygon_count-1].object_index
						<< ", polygon# " << polygon_count-1 
						<< ", vertices indeces = " 
						<< polygons[polygon_count-1].vertices[0] << " "
						<< polygons[polygon_count-1].vertices[1] << " "
						<< polygons[polygon_count-1].vertices[2] << "\n";
				}
			}
		}
	}
	
	// close file
    inFile.close();

	if (print_flag) {
		if (scan_flag) { cout << "scanFile: #vertices = " << vertex_num << ", #polygons = " << polygon_num << "\n\n";}
		else { cout << "getdata: #vertices = " << vertex_num << ", #polygons = " << polygon_num << "\n\n";}
	}

	if (scan_flag) {
		// store number of vertices and polygons found in object
		vertices_array[count] = vertex_num;
		polygons_array[count] = polygon_num;
	}
}

// #####################################################
// #####################################################

void ScanClass::buildMeshAfter(ContainerClass& container ,PolygonClass *polygons) {

	int i, j, start, end, count;
	char file[128];

	// for each object
	for(i=0;i<object_count;i++) {

		// create output filename
		sprintf(file,"%s%s_MANIP.mesh",OUTPUT_DATA_DIR,container.objects[i].name);

		// open output file
		ofstream newfile (file,ios::out);

		if(newfile.is_open()){

			// for each vertex in object
			for(j=0;j<container.objects[i].num_vertices;j++) {

				// print index and final coordinates
				newfile << "Vertex "
						<< j+1 << " "
						<< container.objects[i].vertices[j].new_coordinates[0] << " "
						<< container.objects[i].vertices[j].new_coordinates[1] << " "
						<< container.objects[i].vertices[j].new_coordinates[2] << "\n";
			}

			// identify polygons associated with current object
			start = 0;
			for (j=0;j<i;j++) {
				start = start + polygons_array[j];
			}
			end = start + polygons_array[i];

			// for each polygon in object
			count = 1;
			for (j=start;j<end;j++) {
				newfile << "Face "
						<< count << " "
						<< polygons[j].vertices[0]+1 << " "
						<< polygons[j].vertices[1]+1 << " "
						<< polygons[j].vertices[2]+1 << "\n";
				count++;
			}
		}

		newfile.close();
	}
}

// #####################################################
// #####################################################

void ScanClass::fileInit(void) {

	char file[128];
	sprintf(file,"%s%s",OUTPUT_DATA_DIR,SCAN_LOG_FILE);
	Scanfile.open(file);

}

// #####################################################
// #####################################################

void ScanClass::fileOutit(void) {

	Scanfile.close();

}

// #####################################################
// #####################################################

void ScanClass::updateFile(void) {

	Scanfile << "Input data directory = " << INPUT_DATA_DIR << "\n"
			<< "Total number of input files = " << num_files << "\n"
			<< "Total number of (objects, vertices, polygons) = (" 
			<< object_count << ","
			<< vertex_count << ","
			<< polygon_count << ")\n\n";

    Scanfile.width(15);
	Scanfile << left << "Filename";	
    Scanfile.width(15);
	Scanfile << left << "object #";
    Scanfile.width(15);
	Scanfile << left << "# vertices";
    Scanfile.width(15);
	Scanfile << left << "#polygons" << endl;

	int i;
	for(i=0;i<num_files;i++) {

	    Scanfile.width(15);
		Scanfile << left << files[i];
	    Scanfile.width(15);
		Scanfile << left << i;
	    Scanfile.width(15);
		Scanfile << left << vertices_array[i];
	    Scanfile.width(15);
		Scanfile << left << polygons_array[i] << endl;
	}

}

// #####################################################
// #####################################################

void SpaceClass::init(void) {

	int x, y, z, loc, j;

	// subdivide space
	num_space[0] = (int) ceil( (world[1]-world[0])/SPACE_LENGTH );
	num_space[1] = (int) ceil( (world[3]-world[2])/SPACE_LENGTH );
	num_space[2] = (int) ceil( (world[5]-world[4])/SPACE_LENGTH );
	num_boxes = num_space[0]*num_space[1]*num_space[2];

	// allocate memory for boxes array with each 
	// array element of type box (boxes)
	if (boxes != NULL) {
//		for (j=0;j<num_boxes;j++) {
//			delete[] boxes[j].polygons;
//		}
		delete [] boxes;
	}
	boxes = new BoxClass[num_boxes];
	if (boxes == NULL) { 
		cout << "Not enough memory for BoxClass\n";
		cout.flush();
		exit(1);
	}

	// store box limits in boxes class
	// for each box
	for (z =0;z<num_space[2];z++) {
		for (y =0;y<num_space[1];y++) {
			for (x =0;x<num_space[0];x++) {
		
				loc = z*num_space[0]*num_space[1]+y*num_space[0]+x;
				boxes[loc].limits[0] = world[0]+x*SPACE_LENGTH;
				boxes[loc].limits[1] = world[0]+(x+1)*SPACE_LENGTH;
				boxes[loc].limits[2] = world[2]+y*SPACE_LENGTH;
				boxes[loc].limits[3] = world[2]+(y+1)*SPACE_LENGTH;
				boxes[loc].limits[4] = world[4]+z*SPACE_LENGTH;
				boxes[loc].limits[5] = world[4]+(z+1)*SPACE_LENGTH;
			}
		}
	}

	// initialize polygons arrays
	for (j=0;j<num_boxes;j++) {
		boxes[j].polygons = new int[NUMBER_OF_POLYGONS_IN_BOXES];
		if (boxes[j].polygons == NULL) { 
			cout << "Not enough memory for boxes[j].polygons in SpaceClass::init\n";
			cout.flush();
			exit(1);
		}
		boxes[j].num_polygons=0;
		boxes[j].max_polygons=NUMBER_OF_POLYGONS_IN_BOXES;
	}

}

// #####################################################
// #####################################################

SpaceClass::~SpaceClass(void) {

	delete[] boxes;
}

// #####################################################
// #####################################################

void SpaceClass::clear(void) {

	int j;
	// clear boxes polygon list
	for (j=0;j<num_boxes;j++) {
		delete[] boxes[j].polygons;
		boxes[j].num_polygons=0;
		boxes[j].max_polygons=0;
	}

}

// ######################################
// ######################################

void SpaceClass::boundWorld(ContainerClass& container, int object_count) {

    int i, j, k, count;
	double* x;
	double* y;
	double* z;

	// initialize vectors
	count = 0;
    // for each object
    for (i=0;i<object_count;i++) {
		count += container.objects[i].num_vertices;
    }

	x = new double[count];
	y = new double[count];
	z = new double[count];
	if (x==NULL || y==NULL || z==NULL) { 
		cout << "Not enough memory for x or y or z in SpaceClass::boundWorld\n";
		cout.flush();
		exit(1);
	}

	k = 0;
    ////////// loop through all objects //////////
    // for each object
    for (i=0;i<object_count;i++) {

        ////////// loop through all vertices structure elememts in parent object //////////
        // for each vertex in parent object
        for (j=0;j < container.objects[i].num_vertices;j++) {

            ////////// extract coordinates //////////
			x[k] = container.objects[i].vertices[j].new_coordinates[0];
			y[k] = container.objects[i].vertices[j].new_coordinates[1];
			z[k] = container.objects[i].vertices[j].new_coordinates[2];

			k++;
		}
    }

	/////////// extract min and max //////////
	// '*1.01' added so volume extends past outer most vertices
	// helps in collecting nice vertices
	//

	qsort(x,count,sizeof(double),compare_double);
	qsort(y,count,sizeof(double),compare_double);
	qsort(z,count,sizeof(double),compare_double);
//	Quicksort_double(x,0,count-1);
//	Quicksort_double(y,0,count-1);
//	Quicksort_double(z,0,count-1);

	if (x[0] < 0) {
		world[0] = x[0] * 1.01;
	} else {
		world[0] = x[0] * 0.99;
	}

	if (x[count-1] < 0) {
		world[1] = x[count-1] * 0.99;
	} else {
		world[1] = x[count-1] * 1.01;
	}

	if (y[0] < 0) {
		world[2] = y[0] * 1.01;
	} else {
		world[2] = y[0] * 0.99;
	}

	if (y[count-1] < 0) {
		world[3] = y[count-1] * 0.99;
	} else {
		world[3] = y[count-1] * 1.01;
	}

	if (z[0] < 0) {
		world[4] = z[0] * 1.01;
	} else {
		world[4] = z[0] * 0.99;
	}

	if (z[count-1] < 0) {
		world[5] = z[count-1] * 0.99;
	} else {
		world[5] = z[count-1] * 1.01;
	}


	if (print_flag) {
		cout << "\nworld bounds = [" 
			<< world[0] << " "
			<< world[1] << " "
			<< world[2] << " "
			<< world[3] << " "
			<< world[4] << " "
			<< world[5] << "]\n";
	}

	delete[] x;
	delete[] y;
	delete[] z;

}

// #####################################################
// #####################################################

//void SpaceClass::computeBoxesToCheck(double pvc[3][3], int* &boxes_to_check,int& num) {
void SpaceClass::computeBoxesToCheck(double* pvc[3], int* &boxes_to_check,int& num) {

	// pvc = polygon_vertex_coordinates
	
	int x,y,z, box_range[6], loc;
	double poly_limits[6], pvc_x[3], pvc_y[3], pvc_z[3];

	// identify polygon box limits
	pvc_x[0] = pvc[0][0];
	pvc_x[1] = pvc[1][0];
	pvc_x[2] = pvc[2][0];
	pvc_y[0] = pvc[0][1];
	pvc_y[1] = pvc[1][1];
	pvc_y[2] = pvc[2][1];
	pvc_z[0] = pvc[0][2];
	pvc_z[1] = pvc[1][2];
	pvc_z[2] = pvc[2][2];
	threeValueSort(&pvc_x[0],&poly_limits[1], &poly_limits[0]);
	threeValueSort(&pvc_y[0],&poly_limits[3], &poly_limits[2]);
	threeValueSort(&pvc_z[0],&poly_limits[5], &poly_limits[4]);

	// compute box index range that contains polygon box
	// note this range is zero lower-bounded (lowest range is zeroth box)
	// total range is 0..num_space[i]-1
	box_range[0] = (int) floor( (poly_limits[0]-world[0])/SPACE_LENGTH );   // -x
	box_range[1] = (int) floor( (poly_limits[1]-world[0])/SPACE_LENGTH );	//  x
	box_range[2] = (int) floor( (poly_limits[2]-world[2])/SPACE_LENGTH );	// -y
	box_range[3] = (int) floor( (poly_limits[3]-world[2])/SPACE_LENGTH );	//  y
	box_range[4] = (int) floor( (poly_limits[4]-world[4])/SPACE_LENGTH );   // -z
	box_range[5] = (int) floor( (poly_limits[5]-world[4])/SPACE_LENGTH );	//  z

	int count = 0;
	num = (box_range[1]-box_range[0]+1)*(box_range[3]-box_range[2]+1)*(box_range[5]-box_range[4]+1);
	boxes_to_check = new int[num];
	if (boxes_to_check == NULL) { 
		cout << "Not enough memory for boxes_to_check\n";
		cout.flush();
		exit(1);
	}

	for (z = box_range[4];z<box_range[5]+1;z++){
		for (y = box_range[2];y<box_range[3]+1;y++){
			for (x = box_range[0];x<box_range[1]+1;x++){
				loc = z*num_space[0]*num_space[1]+y*num_space[0]+x;
				boxes_to_check[count] = loc;
				count++;
			}
		}
	}
}

// #####################################################
// #####################################################

void SpaceClass::checkVertexOverlap(double pvc[3][3],int box_index,bool *vertex_overlap) {

	int k;

	// pvc = polygon_vertex_coordinates

	// check if the coordinates of any vertex is between 
	// box_limits for all three coordinate axes
	*vertex_overlap = false;

	for (k=0;k<3;k++) {
		if (
			(
			(pvc[k][0] > boxes[box_index].limits[0]) && 
			(pvc[k][0] < boxes[box_index].limits[1])  
			) && 
			(
			(pvc[k][1] > boxes[box_index].limits[2]) && 
			(pvc[k][1] < boxes[box_index].limits[3])  
			) && 
			(
			(pvc[k][2] > boxes[box_index].limits[4]) && 
			(pvc[k][2] < boxes[box_index].limits[5])  
			)  
			) {
			*vertex_overlap = true;
		}
	}
	if(print_flag) {
	cout << "checkVertexOverlap: vertex_overlap = " << *vertex_overlap << "\n";
	}
}

// #####################################################
// #####################################################
void SpaceClass::checkPolygonPlaneCrossesBox( int box_index ,double pvc[3][3] ,double *poly_normal
											,bool *pos_dot ,bool *neg_dot) {
	// pvc = polygon_vertex_coordinates

	double dot, A, B, C, limits[6];
	int i, j, k;

	*pos_dot = false;
	*neg_dot = false;

	A = pvc[0][0]*poly_normal[0];
	B = pvc[0][1]*poly_normal[1]; 
	C = pvc[0][2]*poly_normal[2];

	for (i=0;i<6;i++) {
		limits[i] = boxes[box_index].limits[i];
	}

	for (i=0;i<2;i++) {
		for (j=3;j>1;j--) {
			for (k=5;k>3;k--) {

				// compute dot product of polygon normal and vertex vectors
				dot = limits[i]*poly_normal[0]-A+limits[j]*poly_normal[1]-B+limits[k]*poly_normal[2]-C;

				// if any two dot products have opposite sign
				// then polygon plane crosses box
				if      (dot > 0 ) { *pos_dot = true;}
				else if (dot < 0 ) { *neg_dot = true;}
			}
		}
	}

}
// #####################################################
// #####################################################

void SpaceClass::checkPolygonEdgeIntersectBoxFace( int box_index ,double pvc[3][3]
											,bool *poly_edge,bool *box_face) {

	// pvc = polygon_vertex_coordinates

	int box_normal[3], m, n, i, val;
	double num, den, u, box_vertex[6][3], pI[3], 
			poly_x3[3], poly_y3[3], poly_z3[3], 
			poly_x4[3], poly_y4[3], poly_z4[3],
			box_limit[6];

	*poly_edge = false;
	*box_face = false;

	poly_x4[0] = pvc[1][0];
	poly_x4[1] = pvc[2][0];
	poly_x4[2] = pvc[0][0];

	poly_y4[0] = pvc[1][1];
	poly_y4[1] = pvc[2][1];
	poly_y4[2] = pvc[0][1];

	poly_z4[0] = pvc[1][2];
	poly_z4[1] = pvc[2][2];
	poly_z4[2] = pvc[0][2];

	box_limit[0] = boxes[box_index].limits[0];
	box_limit[1] = boxes[box_index].limits[1];
	box_limit[2] = boxes[box_index].limits[2];
	box_limit[3] = boxes[box_index].limits[3];
	box_limit[4] = boxes[box_index].limits[4];
	box_limit[5] = boxes[box_index].limits[5];

	box_vertex[0][0] = box_limit[0]; 
	box_vertex[0][1] = box_limit[2]; 
	box_vertex[0][2] = box_limit[4]; 

	box_vertex[1][0] = box_limit[1];
	box_vertex[1][1] = box_limit[3];
	box_vertex[1][2] = box_limit[5];

	box_vertex[2][0] = box_vertex[0][0]; 
	box_vertex[2][1] = box_vertex[0][1];
	box_vertex[2][2] = box_vertex[0][2];

	box_vertex[3][0] = box_vertex[1][0]; 
	box_vertex[3][1] = box_vertex[1][1];
	box_vertex[3][2] = box_vertex[1][2];

	box_vertex[4][0] = box_vertex[0][0]; 
	box_vertex[4][1] = box_vertex[0][1];
	box_vertex[4][2] = box_vertex[0][2];

	box_vertex[5][0] = box_vertex[1][0]; 
	box_vertex[5][1] = box_vertex[1][1];
	box_vertex[5][2] = box_vertex[1][2];

	//for each polygon edge
	for (m=0;m<3;m++) {

		// for each box face
		for (n=0;n<6;n++) {

			if (!*poly_edge  || !*box_face) {

			//use two polygon vertices, box face normal, 
			// and a box face vertex
			// to compute where line of polygon edge intersects
			// plane of box face (if not parallel)
			box_normal[0] = 0;
			box_normal[1] = 0;
			box_normal[2] = 0;
			val = -1;
			for (i=0;i<n;i++) {
				val = val*-1;
			}
			box_normal[n/2] = val;
	
			// dot products
			den = box_normal[0]*(poly_x4[m]-pvc[m][0])
				+ box_normal[1]*(poly_y4[m]-pvc[m][1])
				+ box_normal[2]*(poly_z4[m]-pvc[m][2]);

			// compute point of intersection
			if (den) {
				// point of intersection
				num = box_normal[0]*(box_vertex[n][0]-pvc[m][0])
					+ box_normal[1]*(box_vertex[n][1]-pvc[m][1])
					+ box_normal[2]*(box_vertex[n][2]-pvc[m][2]);
				u = num/den;
				pI[0] = pvc[m][0] + u*(poly_x4[m] -pvc[m][0]);
				pI[1] = pvc[m][1] + u*(poly_y4[m] -pvc[m][1]);
				pI[2] = pvc[m][2] + u*(poly_z4[m] -pvc[m][2]);
	
				// is intersection point between two polygon vertices 
				// and the two coordinate axes ranges?
				// if yes, then polygon is in box

				if ((u > 0) && (u < 1)) {	
					// point of intersection is between polygon vertices
					*poly_edge = true;
				}

				// is point of intersection between box face limits?
				if (n==0 || n==1) { // -x, +x
					if (
						(box_limit[2] < pI[1]) &&
						(box_limit[3] > pI[1]) &&
						(box_limit[4] < pI[2]) &&
						(box_limit[5] > pI[2]) 
						) { 
						*box_face = true;
					}
				} else if (n==2 || n==3) { // -y, y
					if (
						(box_limit[0] < pI[0]) &&
						(box_limit[1] > pI[0]) &&
						(box_limit[4] < pI[2]) &&
						(box_limit[5] > pI[2]) 
						) { 
						*box_face = true;
					}
				} else { // -z, z
					if (
						(box_limit[0] < pI[0]) &&
						(box_limit[1] > pI[0]) &&
						(box_limit[2] < pI[1]) &&
						(box_limit[3] > pI[1]) 
						) { 
						*box_face = true;
					}
				}
			}
			}		
		}	
	}
	if(print_flag) {
	cout << "checkPolygonEdgeIntersectBoxFace: poly_edge = " << *poly_edge
									<< ", box_face = " << *box_face << "\n";
	}
}

// #####################################################
// #####################################################

void SpaceClass::recordPolygon(int* list, int num_list,int polygon_index) {

	int* temp;
	int i, temp_size;

	// record polygon in boxes struct
	for (i=0;i<num_list;i++) {
		boxes[list[i]].addPolygon(polygon_index);
	}

}

// #####################################################
// #####################################################

void SpaceClass::checkBoxOverlap(double* pvc[3],int box_index,bool *box_overlap) {


    int z, y, x, k;
	double xvec[3], yvec[3], zvec[3];
	double biggestx, biggesty, biggestz, smallestx, smallesty, smallestz;
	double poly_limits[6], box_limits[6];

	xvec[0] = pvc[0][0];
	xvec[1] = pvc[1][0];
	xvec[2] = pvc[2][0];
	yvec[0] = pvc[0][1];
	yvec[1] = pvc[1][1];
	yvec[2] = pvc[2][1];
	zvec[0] = pvc[0][2];
	zvec[1] = pvc[1][2];
	zvec[2] = pvc[2][2];

    // sort object_index_list
	threeValueSort(&xvec[0], &biggestx, &smallestx);
	threeValueSort(&yvec[0], &biggesty, &smallesty);
	threeValueSort(&zvec[0], &biggestz, &smallestz);

    // identify polygon box limits
    poly_limits[0] = smallestx;
    poly_limits[1] = biggestx;
    poly_limits[2] = smallesty;
    poly_limits[3] = biggesty;
    poly_limits[4] = smallestz;
    poly_limits[5] = biggestz;

    // get box limits
	box_limits[0] = boxes[box_index].limits[0];
	box_limits[1] = boxes[box_index].limits[1];
	box_limits[2] = boxes[box_index].limits[2];
	box_limits[3] = boxes[box_index].limits[3];
	box_limits[4] = boxes[box_index].limits[4];
	box_limits[5] = boxes[box_index].limits[5];


    ////////// do boxes overlap //////////
    // yes if box max or min is between poly_limits for all three coordinate axes
    *box_overlap = false;

    if (
        (
        (box_limits[0] <= poly_limits[0]) && (box_limits[1] > poly_limits[0]) ||
        (box_limits[0] > poly_limits[0]) && (box_limits[0] < poly_limits[1])
        ) &&
        (
        (box_limits[2] <= poly_limits[2]) && (box_limits[3] > poly_limits[2]) ||
        (box_limits[2] > poly_limits[2]) && (box_limits[2] < poly_limits[3]) 
        ) &&
        (
        (box_limits[4] <= poly_limits[4]) && (box_limits[5] > poly_limits[4]) ||
        (box_limits[4] > poly_limits[4]) && (box_limits[4] < poly_limits[5]) 
        )
        ) {
        *box_overlap = true;
    }

}

// #####################################################
// #####################################################

void SpaceClass::fileInit(void) {

	char file[128];
	sprintf(file,"%s%s",OUTPUT_DATA_DIR,SPACE_LOG_FILE);
	Spacefile.open(file);

    Spacefile.width(15);
	Spacefile << left << "Iteration";
    Spacefile.width(15);
	Spacefile << left << "#boxes(x,y,z)";
    Spacefile.width(15);
	Spacefile << left << "total #boxes";	
    Spacefile.width(30);
	Spacefile << left << "world bounds (x,y,z)[min max]" << endl;


}

// #####################################################
// #####################################################

void SpaceClass::fileOutit(void) {

	Spacefile.close();

}

// #####################################################
// #####################################################

void SpaceClass::updateFile(int iteration) {

	char str[15];
	sprintf(str,"%d,%d,%d",num_space[0],num_space[1],num_space[2]);

    Spacefile.width(15);
	Spacefile << left << iteration;
    Spacefile.width(15);
	Spacefile << left << str;
    Spacefile.width(15);
	Spacefile << left << num_boxes;
	Spacefile << "["
			<< world[0] << " " 
			<< world[1] << " "
			<< world[2] << " "
			<< world[3] << " "
			<< world[4] << " "
			<< world[5] << "]\n";

}

// #####################################################
// #####################################################

void ManipClass::fileInit(void) {

	char file[128];
	sprintf(file,"%s%s",OUTPUT_DATA_DIR,MANIP_LOG_FILE);
	Manipfile.open(file);

}

// #####################################################
// #####################################################

void ManipClass::fileOutit(void) {

	Manipfile.close();

}

// #####################################################
// #####################################################

void ManipClass::computePolygonNormals(ContainerClass& container ,PolygonClass *polygons) {

	int j, i, object_index, vertices[3];
	double uX, uY, uZ, vX, vY, vZ;
	double* vertex[3];

    // for each polygon structure in master polygons array
	for (j=0;j<polygon_count;j++) {

	    // get vertex_index of each vertex of polygon
		for (i=0;i<3;i++) {
			vertices[i] = polygons[j].vertices[i];
		}

		//get object index of polygon
		object_index = polygons[j].object_index;

		//get coordinates of each vertex index
		container.objects[object_index].getPolygonCoordinates(vertices,vertex);

		// compute vectors 01 and 12
		uX = vertex[1][0]-vertex[0][0];
		uY = vertex[1][1]-vertex[0][1];
		uZ = vertex[1][2]-vertex[0][2];
		vX = vertex[2][0]-vertex[0][0];
		vY = vertex[2][1]-vertex[0][1];
		vZ = vertex[2][2]-vertex[0][2];

		// add cross product to normal in polygons
		polygons[j].normal_components[0] = uY*vZ-uZ*vY;
		polygons[j].normal_components[1] = uZ*vX-uX*vZ;
		polygons[j].normal_components[2] = uX*vY-uY*vX;

	}
}

// #####################################################
// #####################################################

void ManipClass::findNice(ContainerClass& container ,PolygonClass *polygons, SpaceClass& space) {

	int i, j;
	double normal[3], ray[3], nice_vol[6];
	int* polygon_flags;
	int num_polygon_flags;
	int max_polygon_flags;
	int num_odd_objects;
	int* crossed_polygons;
	int num_crossed_polygons;
	int max_crossed_polygons;
	int* unique_polygons;
	int num_unique_polygons;
	int max_unique_polygons;

	bool FLAG;

//	time_t after,before, diff1=0,diff2=0,diff3=0,diff4=0,diff5=0,diff6=0;


    ////////// loop through all objects //////////
    // for each object
    for (i=0;i<object_count;i++) {

        ////////// loop through all vertices structure elememts in parent object //////////
        // for each vertex in parent object

        for (j=0;j < container.objects[i].num_vertices;j++) {

			if ( (container.objects[i].vertices[j].deviation_distance < NICE_THRESHOLD) &&
				(container.objects[i].vertices[j].deviation_distance_old < NICE_THRESHOLD) ) {

//				if (i==0 && j==2960) { FLAG = true;}
//				else {FLAG = false;} ///////////////////////////////////////////////////////////////////
				FLAG = false;

				// initialize variables
				// initialize polygon_flags	
				polygon_flags = new int[NUMBER_OF_CROSSED_POLYGONS];
				if (polygon_flags == NULL) { 
					cout << "Not enough memory for polygon_flags\n";
					cout.flush();
					exit(1);
				}
				num_polygon_flags = 0;
				max_polygon_flags = NUMBER_OF_CROSSED_POLYGONS;
	
				// initialize crossed_polygons	
				crossed_polygons = new int[NUMBER_OF_CROSSED_POLYGONS];
				if (crossed_polygons == NULL) { 
					cout << "Not enough memory for crossed_polygons\n";
					cout.flush();
					exit(1);
				}
				num_crossed_polygons = 0;
				max_crossed_polygons = NUMBER_OF_CROSSED_POLYGONS;
	
				// initialize unique_polygons	
				unique_polygons = new int[NUMBER_OF_UNIQUE_POLYGONS];
				if (unique_polygons == NULL) { 
					cout << "Not enough memory for unique_polygons\n";
					cout.flush();
					exit(1);
				}
				num_unique_polygons = 0;
				max_unique_polygons = NUMBER_OF_UNIQUE_POLYGONS;
	
				num_odd_objects = 0;
	
				// pick an arbitrary outward vector
				// this function was already written, so used it
//				before = time (NULL);
				container.objects[i].vertices[j].computeMeanOutwardNormal(polygons,normal); // returns normal
//				after = time (NULL);
//				diff1 += after-before;
				findClosestAxis(normal,ray); // returns ray
//				before = time (NULL);
//				diff2 += before-after;
				collectNicePolygons(container,polygons,space,ray,i,unique_polygons,
									num_unique_polygons,max_unique_polygons,j,FLAG);
//				after = time (NULL);
//				diff3 += after-before;
				findIntersectedPolygons(container,polygons,space,i,j,ray,unique_polygons,
										num_unique_polygons,crossed_polygons,
										num_crossed_polygons,max_crossed_polygons,polygon_flags,
											num_polygon_flags,max_polygon_flags,FLAG);
//				before = time (NULL);
//				diff4 += before-after;
				findOddMeshes(polygons,crossed_polygons,num_crossed_polygons,
										polygon_flags,num_odd_objects,FLAG);
//				after = time (NULL);
//				diff5 += after-before;
		
				container.objects[i].vertices[j].setOdd(num_odd_objects);
//				before = time (NULL);
//				diff6 += before-after;
	
				// free memory
				delete [] unique_polygons;
				delete [] crossed_polygons;
				delete [] polygon_flags;

			}
		}
	}
//	cout << "computeMeanOutwardNormal = " << diff1 << endl;
//	cout << "findClosestAxis = " << diff2 << endl;
//	cout << "collectNicePolygons = " << diff3 << endl;
//	cout << "findIntersectedPolygons = " << diff4 << endl;
//	cout << "findOddMeshes = " << diff5 << endl;
//	cout << "setOdd = " << diff6 << endl;

//	cout << "object 0, vertex 2960, has nice value = "
//		<< container.objects[0].vertices[2960].nice_data.nice << endl;

}

// #####################################################
// #####################################################

void ManipClass::findClosestAxis(double normal[3], double ray[3]) {

	int i, max_index;
	double dot[6], max;

	int axis[6][3] = {
	    			{-1, 0, 0},
					{ 1, 0, 0},
					{ 0,-1, 0},
					{ 0, 1, 0},
					{ 0, 0,-1},
					{ 0, 0, 1}
					};

	// compute dot product of normal and six coordinate axes
	for (i=0;i<6;i++) {
		dot[i] = normal[0]*axis[i][0]+normal[1]*axis[i][1]+normal[2]*axis[i][2];
	}

	// find largest dot product
	max_index = 0;
	max = dot[0];
	for (i=1;i<6;i++) {
		if (dot[i] > max) {
			max = dot[i];
			max_index = i;
		}
	}

	// load ray
	for (i=0;i<3;i++) {
		ray[i] = axis[max_index][i];
	}
}

// #####################################################
// #####################################################

void ManipClass::collectNicePolygons(ContainerClass& container, PolygonClass *polygons, SpaceClass &space, 
									double ray[3], int object_index,int*& unique_polygons, 
								int& num_unique_polygons,int& max_unique_polygons,int vertex_index, bool FLAG) {

	int box_range[6], target_axis, x, y, z, loc_z, loc_zy,num, i, j;
	double coordinates[3];
	int* boxes_to_check;
	int num_boxes_to_check;
	int poi; // polygon_object_index
	int limit;

	// get vertex coordinate associated with ray
	coordinates[0] = container.objects[object_index].vertices[vertex_index].new_coordinates[0];
	coordinates[1] = container.objects[object_index].vertices[vertex_index].new_coordinates[1];
	coordinates[2] = container.objects[object_index].vertices[vertex_index].new_coordinates[2];

	// convert ray direction to array index
	if (ray[0]) {
		if (ray[0]==1) { target_axis = 1;}
		else 			{ target_axis = 0;}
	} else if (ray[1]) {
		if (ray[1]==1) { target_axis = 3;}
		else 			{ target_axis = 2;}
	} else {
		if (ray[2]==1) { target_axis = 5;}
		else 			{ target_axis = 4;}
	}

	if (FLAG) {
		cout << "\nManipClass::collectNicePolygons: "
			<< "current vertex coordinates = "
			<< coordinates[0] << " "
			<< coordinates[1] << " "
			<< coordinates[2]
			<< ", ray = "
			<< ray[0] << " "
			<< ray[1] << " "
			<< ray[2]
			<< ", target_axis = " << target_axis << endl;
	}

	// compute box index range that contains ray
	// note this range is zero lower-bounded (lowest range is zeroth box)
	// total range is 0..num_space[i]-1
	if (target_axis==0) {
		// -x
		box_range[0] = 0; 
		box_range[1] = (int) floor( (coordinates[0]-space.world[0])/SPACE_LENGTH );
		box_range[2] = (int) floor( (coordinates[1]-space.world[2])/SPACE_LENGTH );
		box_range[3] = box_range[2];
		box_range[4] = (int) floor( (coordinates[2]-space.world[4])/SPACE_LENGTH );
		box_range[5] = box_range[4];
	} else if (target_axis==1) {
		// +x
		box_range[0] = (int) floor( (coordinates[0]-space.world[0])/SPACE_LENGTH );
		box_range[1] = (int) floor( (space.world[1]-space.world[0])/SPACE_LENGTH );
		box_range[2] = (int) floor( (coordinates[1]-space.world[2])/SPACE_LENGTH );
		box_range[3] = box_range[2];
		box_range[4] = (int) floor( (coordinates[2]-space.world[4])/SPACE_LENGTH );
		box_range[5] = box_range[4];
	} else if (target_axis==2) {
		// -y
		box_range[0] = (int) floor( (coordinates[0]-space.world[0])/SPACE_LENGTH );
		box_range[1] = box_range[0];
		box_range[2] = 0;
		box_range[3] = (int) floor( (coordinates[1]-space.world[2])/SPACE_LENGTH );
		box_range[4] = (int) floor( (coordinates[2]-space.world[4])/SPACE_LENGTH );
		box_range[5] = box_range[4];
	} else if (target_axis==3) {
		// +y
		box_range[0] = (int) floor( (coordinates[0]-space.world[0])/SPACE_LENGTH );
		box_range[1] = box_range[0];
		box_range[2] = (int) floor( (coordinates[1]-space.world[2])/SPACE_LENGTH );
		box_range[3] = (int) floor( (space.world[3]-space.world[2])/SPACE_LENGTH );
		box_range[4] = (int) floor( (coordinates[2]-space.world[4])/SPACE_LENGTH );
		box_range[5] = box_range[4];
	} else if (target_axis==4) {
		// -z
		box_range[0] = (int) floor( (coordinates[0]-space.world[0])/SPACE_LENGTH );
		box_range[1] = box_range[0];
		box_range[2] = (int) floor( (coordinates[1]-space.world[2])/SPACE_LENGTH );
		box_range[3] = box_range[2];
		box_range[4] = 0;
		box_range[5] = (int) floor( (coordinates[2]-space.world[4])/SPACE_LENGTH );
	} else {
		// +z
		box_range[0] = (int) floor( (coordinates[0]-space.world[0])/SPACE_LENGTH );
		box_range[1] = box_range[0];
		box_range[2] = (int) floor( (coordinates[1]-space.world[2])/SPACE_LENGTH );
		box_range[3] = box_range[2];
		box_range[4] = (int) floor( (coordinates[2]-space.world[4])/SPACE_LENGTH );
		box_range[5] = (int) floor( (space.world[5]-space.world[4])/SPACE_LENGTH );
	}


	////////// collect boxes to check //////////
	// compute # boxes to check
	num = (box_range[1]-box_range[0]+1)*(box_range[3]-box_range[2]+1)*(box_range[5]-box_range[4]+1);
	boxes_to_check = new int[num];
	if (boxes_to_check == NULL) { 
		cout << "Not enough memory for boxes_to_check in ManipClass::collectNicePolygons\n";
		cout.flush();
		exit(1);
	}
	num_boxes_to_check = 0;
	
	for (z = box_range[4];z<box_range[5]+1;z++){
		loc_z = z*space.num_space[0]*space.num_space[1];
		for (y = box_range[2];y<box_range[3]+1;y++){
			loc_zy = loc_z+y*space.num_space[0];
			for (x = box_range[0];x<box_range[1]+1;x++){
				boxes_to_check[num_boxes_to_check] = loc_zy+x;
				num_boxes_to_check++;
			}
		}
	}

	if (FLAG) {
		cout << "ManipClass::collectNicePolygons: "
			<< "box_range = "
			<< box_range[0] << " "
			<< box_range[1] << " "
			<< box_range[2] << " "
			<< box_range[3] << " "
			<< box_range[4] << " "
			<< box_range[5]
			<< ", num_boxes_to_check = " << num_boxes_to_check << endl;
	}
	////////// gather polygons in boxes //////////

	// for each box
	for (j=0;j<num_boxes_to_check;j++) {

		// for each polygon in box
		limit = space.boxes[boxes_to_check[j]].num_polygons;
		for (i=0;i<limit;i++) {

			// add polygon if not from same object as current vertex
			poi = space.boxes[boxes_to_check[j]].polygons[i];
			if (polygons[poi].object_index != object_index) {
//				addSortUniqueArrayElement(unique_polygons,num_unique_polygons,max_unique_polygons,
//								INCREMENT_OF_UNIQUE_POLYGONS,space.boxes[boxes_to_check[j]].polygons[i]);
				addArrayElement(unique_polygons,num_unique_polygons,max_unique_polygons,
								INCREMENT_OF_UNIQUE_POLYGONS,space.boxes[boxes_to_check[j]].polygons[i]);
			}
		}

	}

	// sort array
//	Quicksort(unique_polygons,0,num_unique_polygons-1);
	qsort(unique_polygons,num_unique_polygons,sizeof(int),compare_int);

	// extract unique array
	unique_polygons = uniqueArray(unique_polygons,num_unique_polygons,max_unique_polygons);

	if (FLAG) {
		cout << "ManipClass::collectNicePolygons: "
			<< "unique_polygons = ";
		for (i=0;i<num_unique_polygons;i++) {
			cout << unique_polygons[i] << endl;
		}
	}

	delete[] boxes_to_check;

}

// #####################################################
// #####################################################

void ManipClass::findIntersectedPolygons(ContainerClass& container ,PolygonClass *polygons, SpaceClass& space,
											int object_index, int vertex_index, double ray[3], 
											int*& unique_polygons, int& num_unique_polygons, 
											int*& crossed_polygons, int& num_crossed_polygons, 
											int& max_crossed_polygons, 
											int*& polygon_flags, int& num_polygon_flags,
											int& max_polygon_flags,bool FLAG) {


	int i,j,k,index[3], polygon_object_index;
	double lp[2][3], coordinates[3];
	double* pvc[3];
	bool line_flag, poly_flag, poly_edge_flag;

	// get vertex coordinate associated with ray
	coordinates[0] = container.objects[object_index].vertices[vertex_index].new_coordinates[0];
	coordinates[1] = container.objects[object_index].vertices[vertex_index].new_coordinates[1];
	coordinates[2] = container.objects[object_index].vertices[vertex_index].new_coordinates[2];

	// extend ray to world bounds
	if (ray[0]) {
		if (ray[0]==1) { ray[0] = space.world[1]-coordinates[0];}
		else 			{ ray[0] = space.world[0]-coordinates[0];}
	} else if (ray[1]) {
		if (ray[1]==1) { ray[1] = space.world[3]-coordinates[1];}
		else 			{ ray[1] = space.world[2]-coordinates[1];}
	} else {
		if (ray[2]==1) { ray[2] = space.world[5]-coordinates[2];}
		else 			{ ray[2] = space.world[4]-coordinates[2];}
	}

	for (i=0;i<3;i++) {
	// tweak ray slightly
	  	ray[i] += RAY_EPSILON;
		// points on ray
		lp[0][i] = coordinates[i];
		lp[1][i] = lp[0][i]+ray[i];
	}

	if (FLAG) {
		cout << "\nManipClass::ManipClass::findIntersectedPolygons: "
			<< "ray = "
			<< ray[0] << " "
			<< ray[1] << " "
			<< ray[2]
			<< ", points on ray = "
			<< lp[0][0] << " "
			<< lp[0][1] << " "
			<< lp[0][2] << ", "
			<< lp[1][0] << " "
			<< lp[1][1] << " "
			<< lp[1][2] << "\n";
	}


	// for each unique polygon
	for (i=0;i<num_unique_polygons;i++) {

		// get polygon object index
		polygon_object_index = polygons[unique_polygons[i]].object_index;

		// get polygon vertex indices
		for (k=0;k<3;k++) {
			index[k] = polygons[unique_polygons[i]].vertices[k];
		}

		// get polygon vertex coordinates
		for (k=0;k<3;k++) {
			pvc[k] = &container.objects[polygon_object_index].vertices[index[k]].new_coordinates[0];
		}

		checkLinePolygonIntersection(pvc, polygons[unique_polygons[i]].normal_components,
										lp,&line_flag,&poly_flag,&poly_edge_flag);

		// does point intersect polygon
		if (poly_flag) {

			// add polygon_index to crossed array
			addArrayElement(crossed_polygons,num_crossed_polygons,max_crossed_polygons,
							INCREMENT_OF_CROSSED_POLYGONS,unique_polygons[i]);

			if (detect_polygon_edge_intersection) {
				// if intersection point falls on edge, signal with flag
				if ( poly_edge_flag ) {
					// add polygon flag
					addArrayElement(polygon_flags,num_polygon_flags,max_polygon_flags,
									INCREMENT_OF_CROSSED_POLYGONS,1);
				} else {
					// add polygon flag
					addArrayElement(polygon_flags,num_polygon_flags,max_polygon_flags,
									INCREMENT_OF_CROSSED_POLYGONS,0);
				}
			}
		}
	}

	if (FLAG) {
		cout << "\nManipClass::ManipClass::findIntersectedPolygons: "
			<< "crossed_polygon index = ";
		for (i=0;i<num_crossed_polygons;i++) {
			cout << crossed_polygons[i]
				<< ", object_index = " << polygons[crossed_polygons[i]].object_index
				<< ", polygon vertex indices = "
				<< polygons[crossed_polygons[i]].vertices[0] << " "
				<< polygons[crossed_polygons[i]].vertices[1] << " "
				<< polygons[crossed_polygons[i]].vertices[2] << "\n";
		}
	}


}

// #####################################################
// #####################################################

void ManipClass::findOddMeshes(PolygonClass *polygons, int* crossed_polygons, int num_crossed_polygons, 
								int* polygon_flags, int& num_odd_objects,bool FLAG) {

	int i, j, k, count, parity,
		polygon_object_index, edge_crossed_polygon_object_index, edge_num;
	int* object_index_list;
	int num_object_index_list;
	int max_object_index_list;
	int sum;

	object_index_list = new int [NUMBER_OF_CROSSED_POLYGONS];
	if (object_index_list == NULL) { 
		cout << "Not enough memory for object_index_list in ManipClass::findOddMeshes\n";
		cout.flush();
		exit(1);
	}
	num_object_index_list = 0;
	max_object_index_list = NUMBER_OF_CROSSED_POLYGONS;

    // for each crossed polygon
    for (i=0;i<num_crossed_polygons;i++) {

        // get polygon object index
        polygon_object_index = polygons[crossed_polygons[i]].object_index;

		if (detect_polygon_edge_intersection) {
			// if ray intersected polygon edge
			if (polygon_flags[i]){
				sum = 0;
				//look for another polygon with edge crossed
	    		for (j=i+1;j<num_crossed_polygons;j++) {
	        		// get polygon object index
			        edge_crossed_polygon_object_index = polygons[crossed_polygons[j]].object_index;
					// if ray intersected polygon edge
					if (polygon_flags[j] && (polygon_object_index==edge_crossed_polygon_object_index)){
						//increment count
						sum++;
					}
				}
	
				// compute parity of sum
				parity = sum % 2;
				// if even, add instance to object list
				if (!parity) {
					addArrayElement(object_index_list,num_object_index_list,max_object_index_list,
									INCREMENT_OF_CROSSED_POLYGONS,polygon_object_index);
				}
			}
		} else {	
			// add instance to object list
			addArrayElement(object_index_list,num_object_index_list,max_object_index_list,
							INCREMENT_OF_CROSSED_POLYGONS,polygon_object_index);
		}
	}

	if (FLAG) {
		cout << "\nManipClass::findOddMeshes: "
			<< "object_index_list = ";
		for (i=0;i<num_object_index_list;i++) {
			cout << object_index_list[i] << endl;
		}
	}

	// sort object_index_list
	qsort(object_index_list,num_object_index_list,sizeof(int),compare_int);
//	Quicksort(object_index_list,0,num_object_index_list-1);

	if (FLAG) {
		cout << "\nManipClass::findOddMeshes: "
			<< "object_index_list = ";
		for (i=0;i<num_object_index_list;i++) {
			cout << object_index_list[i] << endl;
		}
	}

	i=0;
	while (i<num_object_index_list) {
		if (i+1 != num_object_index_list) {
			if (object_index_list[i] == object_index_list[i+1]) {
	    	// skip pairs of objects
				i++;i++;
			} else {
				num_odd_objects++;
				i++;
			}
		} else {
		// add remaining objects to odd object list
			num_odd_objects++;
			i++;
		}
	}

	if (FLAG) {
		cout << "\nManipClass::findOddMeshes: "
			<< "num_odd_objects = " << num_odd_objects << endl;
	}

	delete[] object_index_list;

}

// #####################################################
// #####################################################

int ManipClass::checkPolygonPolygonIntersections(ContainerClass &container ,PolygonClass *polygons,
													int cpi,int opi) {

	int i, j, k, object_index, other_object_index, opvi[3], cpvi[3], current_pairs[3][2], num_unique;
	double cn[3], on[3], X[3], Y[3], colinear, cos_theta_square;
	double* opvc[3];
	double* cpvc[3];
	double term1;
	bool coplanar_flag, colinear_flag, parallel_flag, share_edge_flag, identical_flag;

	// cpvi = current_polygon_vertex_indices
	// cpvc = current_polygon_vertex_coordinates
	// opvc = other_polygon_vertex_coordinates
	// opvi = other_polygon_vertex_indices
	// cpi = current_polygon_index
	// opi = other_polygon_index
	// cn   = current_normal
	// on   = other_normal

	// get current polygon object index
	object_index = polygons[cpi].object_index;

	// get current polygon normal
	cn[0] = polygons[cpi].normal_components[0];
	cn[1] = polygons[cpi].normal_components[1];
	cn[2] = polygons[cpi].normal_components[2];

	// get current polygon vertex indices
	cpvi[0] = polygons[cpi].vertices[0];
	cpvi[1] = polygons[cpi].vertices[1];
	cpvi[2] = polygons[cpi].vertices[2];

	// get current polygon vertex coordinates
	container.objects[object_index].getPolygonCoordinates(polygons[cpi].vertices, cpvc);

	// describe three pairs of current polygon vertices
    current_pairs[0][0] = 0;
    current_pairs[0][1] = 1;
    current_pairs[1][0] = 1;
    current_pairs[1][1] = 2;
    current_pairs[2][0] = 2;
    current_pairs[2][1] = 0;

	// get other_polygon object index
	other_object_index = polygons[opi].object_index;

	// get other_polygon normal
	on[0] = polygons[opi].normal_components[0];
	on[1] = polygons[opi].normal_components[1];
	on[2] = polygons[opi].normal_components[2];

	// get other_polygon vertex indices
	opvi[0] = polygons[opi].vertices[0] ;
	opvi[1] = polygons[opi].vertices[1] ;
	opvi[2] = polygons[opi].vertices[2] ;

	// get other_polygon vertex coordinates
	container.objects[other_object_index].getPolygonCoordinates(polygons[opi].vertices, opvc);

	// initialize flags
	coplanar_flag = false;
	colinear_flag = false;
	parallel_flag = false;
	share_edge_flag = false;
	identical_flag = false;

	// are current polygon and other polygon parallel
	// i.e. is angle between normals equal to zero?
	// i.e. is the square of the cosine of the angle equal to 1?
	term1 = cn[0]*on[0]+ cn[1]*on[1]+ cn[2]*on[2];

	if ( !distinguishable(term1*term1,
	(cn[0]*cn[0]+cn[1]*cn[1]+cn[2]*cn[2])*(on[0]*on[0]+on[1]*on[1]+on[2]*on[2])) ) {parallel_flag = true;}

	// dot product of other polygon normal and 
	// line connecting point on other polygon to point on current polygon
	// try to choose vertices not shared by both polygons

	if ( (opvi[0] != cpvi[0]) && (opvi[0] != cpvi[1]) && (opvi[0] != cpvi[2]) ) { 
		X[0] = opvc[0][0]; 
		X[1] = opvc[0][1]; 
		X[2] = opvc[0][2];
	} else if ( (opvi[1] != cpvi[0]) && (opvi[1] != cpvi[1]) && (opvi[2] != cpvi[2]) ) { 
		X[0] = opvc[1][0]; 
		X[1] = opvc[1][1]; 
		X[2] = opvc[1][2];
	} else { 
		X[0] = opvc[2][0]; 
		X[1] = opvc[2][1]; 
		X[2] = opvc[2][2];
	}

	if ( (cpvi[0] != opvi[0]) && (cpvi[0] != opvi[1]) && (cpvi[0] != opvi[2]) ) {
		Y[0] = cpvc[0][0]; 
		Y[1] = cpvc[0][1]; 
		Y[2] = cpvc[0][2]; 
	} else if ( (cpvi[1] != opvi[0]) && (cpvi[1] != opvi[1]) && (cpvi[1] != opvi[2]) ) {
		Y[0] = cpvc[1][0]; 
		Y[1] = cpvc[1][1]; 
		Y[2] = cpvc[1][2]; 
	} else { 
		Y[0] = cpvc[2][0]; 
		Y[1] = cpvc[2][1]; 
		Y[2] = cpvc[2][2]; 
	}

	colinear = on[0]*(X[0]-Y[0]) + on[1]*(X[1]-Y[1]) + on[2]*(X[2]-Y[2]);

	// if polygons colinear, then normal and line are orthogonal
	// and dot product will be ~zero

	if (fabs(colinear)<DOUBLE_EPSILON) {colinear_flag = true;}
	if (parallel_flag && colinear_flag) {coplanar_flag = true;}
		
	// how many vertices are shared between current and other polygon
	num_unique = 0;
	// for each current polygon vertex
	for (j=0;j<3;j++) {
		if (
			distinguishable(cpvc[j][0], opvc[0][0]) ||
			distinguishable(cpvc[j][1], opvc[0][1]) ||
			distinguishable(cpvc[j][2], opvc[0][2])
			) {num_unique++;}
		if (
			distinguishable(cpvc[j][0], opvc[1][0]) ||
			distinguishable(cpvc[j][1], opvc[1][1]) ||
			distinguishable(cpvc[j][2], opvc[1][2])
			) {num_unique++;}
		if (
			distinguishable(cpvc[j][0], opvc[2][0]) ||
			distinguishable(cpvc[j][1], opvc[2][1]) ||
			distinguishable(cpvc[j][2], opvc[2][2])
			) {num_unique++;}
	}
	if (9-num_unique == 2) {share_edge_flag = true;}
	else if (9-num_unique == 3) {identical_flag = true;}

	////////// begin decision tree //////////

	if (coplanar_flag) {
		if (identical_flag) {
			// polygons are identical
			return(1);
		} else {
			//do polygon edges intersect?	
			// if yes, intersect
			// if no, do not intersect
			if (polygons[cpi].checkEdgeEdgeIntersection(cpvc, current_pairs, opvc,
													cpi, opi, share_edge_flag)) {return(1);}
			else {return(0);}
		}
	} else {
		if (!share_edge_flag) {
			// do polygons intersect?
			if (polygons[cpi].checkPolygonEdgePolygonIntersection( cpvc, current_pairs, on,  
																opvc, cpi, opi)) {return(1);}
			else {return(0);}
		} else { return(0);}
	}
}

// #####################################################
// #####################################################

void ManipClass::assignPolygonsToBoxes(ContainerClass& container ,PolygonClass *polygons, 
											SpaceClass& space) {

	int j;
	int* list;
	int num_list;
	int max_list;

	////////// identify in which boxes each polygon exists ////////

    // for each polygon structure in master polygons array
	for (j=0;j<polygon_count;j++) {

		list = new int[NUMBER_OF_LIST_ENTRIES];
		if (list == NULL) { 
			cout << "Not enough memory for list in ManipClass::assignPolygonsToBoxes\n";
			cout.flush();
			exit(1);
		}
		num_list = 0;
		max_list = NUMBER_OF_LIST_ENTRIES;

		// identify boxes containing part of polygon
		identifyBoxes(container,polygons,space,j,list,num_list,max_list);

		if (num_list == 0) {
			cout << "ERROR: NO BOXES FOR POLYGON INDEX " << j << "\n";
		}

		// record boxes in polygon struct
		polygons[j].recordBoxes(list,num_list);

		// record polygon in boxes struct
		space.recordPolygon(list,num_list,j);

		// free memory
		delete[] list;
	}
}

// #####################################################
// #####################################################

void ManipClass::identifyBoxes(ContainerClass& container ,PolygonClass *polygons,SpaceClass &space
								,int polygon_index,int*& list,int& num_list,int& max_list) {

	int i, num_boxes, j, object_index, polygon_vertex_index, k;
	double poly_normal[3];
	double* pvc[3];
	bool space_overlap, vertex_overlap, pos_dot, neg_dot, poly_edge, box_face, poly_face, box_edge;
	bool box_overlap;

	// pvc = polygon_vertex_coordinates

	// get polygon object index
	object_index = polygons[polygon_index].object_index;

	// get polygon normal
	for (j=0;j<3;j++) {
		poly_normal[j] = polygons[polygon_index].normal_components[j];
	}

	container.objects[object_index].getPolygonCoordinates(polygons[polygon_index].vertices,pvc);

	////////// compute boxes to check//////////
	int* boxes_to_check;
	space.computeBoxesToCheck(pvc,boxes_to_check,num_boxes);

	////////// for each box //////////
	for (k=0;k<num_boxes;k++) {
		i = boxes_to_check[k];
		space.checkBoxOverlap(pvc,i,&box_overlap);
		if (box_overlap) {
			addArrayElement(list,num_list,max_list,INCREMENT_OF_LIST_ENTRIES,i);
		}
	}
	delete[] boxes_to_check;
}

// #####################################################
// #####################################################

void ManipClass::setDeviationDistance(ContainerClass& container ,PolygonClass *polygons, 
										SpaceClass &space) {

	int i, j;

    ////////// loop through all objects //////////
    // for each object
    for (i=0;i<object_count;i++) {

        ////////// loop through all vertices structure elememts in parent object //////////
        // for each vertex in parent object
        for (j=0;j < container.objects[i].num_vertices;j++) {

			////////// find closest vertex to current vertex //////////
			////////// on object other than current object //////////
			////////// returns object index and vertex index for closest vertex //////////
			findClosestVertex(container, polygons, space, i, j);
		}
	}

}

// #####################################################
// #####################################################

void ManipClass::findClosestVertex(ContainerClass& container ,PolygonClass *polygons, SpaceClass &space, 
											int coi, int cvi) {

	int cbi[3], offset, box_range[6], x, y, z, closest_polygon_index , loc_z, loc_zy, loc,i, j, k, m, p;
	double current_coordinates[3], closest_coordinates[3], squareD, dotproduct, temp;
	bool closest_flag, candidates_flag, safe_flag; 

	int* closest_polygons;
	int num_closest_polygons;
	int count;
	int nnv;
	int poly_vertex, neighbor_vertex;
	bool niceness = container.objects[coi].vertices[cvi].nice_data.nice;
	int* pItem;

	// cvi = current_vertex_index
	// coi = current_object_index
	// cbi = current_box_index


	// get current vertex number of neighbor vertices
	nnv = container.objects[coi].vertices[cvi].next_neighbor_vertex;

	// get current vertex coordinate
	for (i=0;i<3;i++){
		current_coordinates[i] = container.objects[coi].vertices[cvi].new_coordinates[i];
	}

	// compute box index that contains current_vertex
	// note this index range is zero lower-bounded (lowest range is zeroth box)
	// total range is 0..num_space[i]-1
	for (i=0;i<3;i++){
		cbi[i] = (int) floor( (current_coordinates[i]-space.world[2*i])/SPACE_LENGTH );
	}

	// initialize candidate flag
	container.objects[coi].vertices[cvi].candidate = false;

	// define box offset from current_box to include in polygon check
	closest_flag = false;
	squareD = 0;
	for (offset=0;offset<(NUM_ADJACENT_BOXES+1);offset++) {
		// if closest vertex hasn't already been found
		if (!closest_flag) {
			box_range[0] = cbi[0]-offset;
			box_range[1] = cbi[0]+offset;
			box_range[2] = cbi[1]-offset;
			box_range[3] = cbi[1]+offset;
			box_range[4] = cbi[2]-offset;
			box_range[5] = cbi[2]+offset;
	
			if (box_range[0]<0) {box_range[0] = 0;}
			if (box_range[2]<0) {box_range[2] = 0;}
			if (box_range[4]<0) {box_range[4] = 0;}
			if (box_range[1]>(space.num_space[0]-1)) {box_range[1] = space.num_space[0]-1;}
			if (box_range[3]>(space.num_space[1]-1)) {box_range[3] = space.num_space[1]-1;}
			if (box_range[5]>(space.num_space[2]-1)) {box_range[5] = space.num_space[2]-1;}

			num_closest_polygons = 0;

			for (z = box_range[4];z<box_range[5]+1;z++){
				loc_z = z*space.num_space[0]*space.num_space[1];
				for (y = box_range[2];y<box_range[3]+1;y++){
					loc_zy = loc_z+y*space.num_space[0];
					for (x = box_range[0];x<box_range[1]+1;x++){
						num_closest_polygons += space.boxes[loc_zy+x].num_polygons;
					}
				}
			}

			closest_polygons = new int[num_closest_polygons];
			if (closest_polygons == NULL) { 
				cout << "Not enough memory for closest_polygons in ManipClass::findClosestVertex\n";
				cout.flush();
				exit(1);
			}

//			if (coi == 0 && cvi==2960) {
//				cout << "num_closest_polygons = " << num_closest_polygons << endl;
//			}

			////////// for each box //////////
			candidates_flag = false;
			count = 0;
			for (z = box_range[4];z<box_range[5]+1;z++){
				loc_z = z*space.num_space[0]*space.num_space[1];
				for (y = box_range[2];y<box_range[3]+1;y++){
					loc_zy = loc_z+y*space.num_space[0];
					for (x = box_range[0];x<box_range[1]+1;x++){
						loc = loc_zy+x;

						// for each polygon in box
						for (k=0;k<space.boxes[loc].num_polygons;k++) {
							safe_flag = true;
							// check if object index of polygon is different 
							// than object index of current vertex
							if (polygons[space.boxes[loc].polygons[k]].object_index == coi) {
								//// check if polygon vertices are not in neighborhood of current vertex ////

								// for each polygon vertex
								for (p=0;p<3;p++) {
									if (safe_flag) {
										poly_vertex = polygons[space.boxes[loc].polygons[k]].vertices[p];

										// look for polygon vertex index in sorted neighborhood list
//										m = 0;
//										search_flag = true;
//										while (m<nnv && search_flag) {
//											neighbor_vertex = container.objects[coi].vertices[cvi].neighborhood_vertices_index[m];
//											if (poly_vertex == neighbor_vertex) {safe_flag = false; search_flag = false;}
//											else if (poly_vertex < neighbor_vertex) {search_flag = false;}
//											m++;
//										}
										pItem = (int*) bsearch(&poly_vertex,
																container.objects[coi].vertices[cvi].neighborhood_vertices_index,
																nnv,sizeof(int),compare_int);
										if (pItem != NULL) {safe_flag = false;}
									}
								}
							}
							if (safe_flag) { 
								closest_polygons[count] = space.boxes[loc].polygons[k];
								count++;
								// signal that candidate polygons were found
								candidates_flag = true;
							}
						}
					}
				}
			}

			if (candidates_flag){
				// for each candidate polygon
				for (j=0;j<count;j++) {

					// find closest point on polygon to current vertex
					computeClosest(container,polygons,space,closest_polygons[j], 
									&current_coordinates[0],&closest_coordinates[0]);

//					if (coi == 0 && cvi==2960) {
//							cout << closest_coordinates[0] << " "
//								<< closest_coordinates[1] << " "
//								<< closest_coordinates[2]
//								<< " 1 0 0 1\n";
//					}
					// is closest point located correctly with respect to vertex outward normal vector?
					// i.e. if not nice, then dot product of deviation vector 
					// and outward_normal must be negative
					dotproduct = 0;
//					if (!container.objects[coi].vertices[cvi].nice_data.nice) {
					if (!niceness) {
						// not nice
						for (i=0;i<3;i++) {
//							dotproduct += closest_coordinates[i]*container.objects[coi].vertices[cvi].outward_normal[i];
							dotproduct += (closest_coordinates[i]-current_coordinates[i])
											*(container.objects[coi].vertices[cvi].outward_normal[i]-current_coordinates[i]);
						}
//						if (coi == 0 && cvi==2960) {
//							cout << "outward normal = "
//								<< container.objects[coi].vertices[cvi].outward_normal[0] << " "
//								<< container.objects[coi].vertices[cvi].outward_normal[1] << " "
//								<< container.objects[coi].vertices[cvi].outward_normal[2] << "\n";
//						}

					}
					if ( (!niceness && (dotproduct < 0)) || !dotproduct ){

						// compute square of deviation distance
						temp = 0;
						for (i=0;i<3;i++) {
							temp +=	(closest_coordinates[i]-current_coordinates[i])
									*(closest_coordinates[i]-current_coordinates[i]);
						}

						// is square of deviation distance closer than square of SEARCH_RADIUS
						if (temp<(SEARCH_RADIUS*SEARCH_RADIUS)) {

							// if yes, is square of deviation distance closer than
							// previously saved squares of deviation distance
							//
							if (temp<squareD || !squareD) {
	
								// if yes, 
								// save 
								for (i=0;i<3;i++) {
									container.objects[coi].vertices[cvi].closest[i] = closest_coordinates[i];
								}
								// and distance
								squareD = temp;
								// set closest_flag = true
								closest_flag = true;
								container.objects[coi].vertices[cvi].candidate = true;
								if (coi==1) {
									container.objects[coi].vertices[cvi].candidate = false;
								}
	
							}
						}
					}
				}
			}

			delete[] closest_polygons;
		}
	}
}
// #####################################################
// #####################################################

void ManipClass::computeClosest(ContainerClass& container ,PolygonClass *polygons, SpaceClass &space,
				int polygon_index, double *current_vertex_coordinates, double *closest_coordinates) {

	double normal[3], den, num, u, diff[3], candidates[7][3], squareD, temp; 
	double* P[3];
	double Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz, uNum, uDen, AdotA, AdotB, AdotC, BdotC, BdotB;
	double Ix, Iy, Iz, IxNum, IxDen, IyNum, IyDen, IzNum, IzDen, intersect[3], lp[2][3];
	int indices[3], polygon_object_index, i, j, current_pairs[3][2], count = 0, store;
    bool line_flag, poly_flag, poly_edge_flag;

	// get polygon vertex indices
	for (i=0;i<3;i++) {
		indices[i] = polygons[polygon_index].vertices[i];
	}

	// get polygon object index
	polygon_object_index = polygons[polygon_index].object_index;	

	// get polygon normal
	for (i=0;i<3;i++) {
		normal[i] = polygons[polygon_index].normal_components[i];
	}

	// get polygon vertex coordinates
	for (i=0;i<3;i++) {
		P[i] = &container.objects[polygon_object_index].vertices[indices[i]].new_coordinates[0];
	}

	// compute vector connecting arbitrary polygon vertex and current vertex
	for (i=0;i<3;i++) {
		diff[i] = P[0][i]-current_vertex_coordinates[i];
	}

	// compute indicators
	num = normal[0]*diff[0]   + normal[1]*diff[1]   + normal[2]*diff[2];
	den = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];

	////////// check current vertex/polygon plane coincidence //////////
	if (!(num < DOUBLE_EPSILON) ) {
		//current vertex does not lie on polygon plane
		//compute point on polygon plane that is closest to current vertex
		// i.e. intersection of plane with polygon normal through current vertex
		u = num/den;

		for (i=0;i<3;i++) {
			intersect[i] = current_vertex_coordinates[i]+u*normal[i];
		}

	    line_flag = false;
	    poly_flag = false;

		lp[0][0] = current_vertex_coordinates[0];
		lp[0][1] = current_vertex_coordinates[1];
		lp[0][2] = current_vertex_coordinates[2];
		lp[1][0] = intersect[0];
		lp[1][1] = intersect[1];
		lp[1][2] = intersect[2];

		checkLinePolygonIntersection(P,normal,lp,&line_flag,&poly_flag,&poly_edge_flag);

		// add polygon perpendicular intersection point
		if (poly_flag) {
			for (i=0;i<3;i++) {
				candidates[count][i] = current_vertex_coordinates[i]+u*normal[i];
			}
			count++;
		}

	}

	////////// add points of minimum distance between current vertex and each polygon edge//////////

	// describe three pairs of polygon vertices
    current_pairs[0][0] = 0;
    current_pairs[0][1] = 1;
    current_pairs[1][0] = 1;
    current_pairs[1][1] = 2;
    current_pairs[2][0] = 2;
    current_pairs[2][1] = 0;

	// for each polygon edge
	for (i=0;i<3;i++) {
		Ax = P[current_pairs[i][0]][0];
		Ay = P[current_pairs[i][0]][1];
		Az = P[current_pairs[i][0]][2];
		Bx = P[current_pairs[i][1]][0];
		By = P[current_pairs[i][1]][1];
		Bz = P[current_pairs[i][1]][2];
		Cx = current_vertex_coordinates[0];
		Cy = current_vertex_coordinates[1];
		Cz = current_vertex_coordinates[2];

		AdotA = Ax*Ax+Ay*Ay+Az*Az;
		AdotB = Ax*Bx+Ay*By+Az*Bz;
		AdotC = Ax*Cx+Ay*Cy+Az*Cz;
		BdotC = Bx*Cx+By*Cy+Bz*Cz;
		BdotB = Bx*Bx+By*By+Bz*Bz;

		uDen = AdotA+BdotB-2*AdotB;
		if(uDen) {
			uNum = AdotA-AdotB-AdotC+BdotC;
			u = uNum/uDen;
			if (u>0 && u<1) {
				// no need to check for u ==0 and u ==1,
				// since current vertex/polygon plane
				// coincidence was checked above.
				//
				// closest point on polygon edge line to current vertex
				// occurs on polygon edge between polygon vertices

				Ix = Ax+u*(Bx-Ax);
				Iy = Ay+u*(By-Ay);
				Iz = Az+u*(Bz-Az);
				candidates[count][0] = Ix;
				candidates[count][1] = Iy;
				candidates[count][2] = Iz;
				count++;
			}
		}
		
	}


	////////// add each polygon vertex//////////

	// for each polygon vertex
	for (i=0;i<3;i++) {
		candidates[count][0] = P[i][0];
		candidates[count][1] = P[i][1];
		candidates[count][2] = P[i][2];
		count++;
	}

	////////// find closest point //////////

	// for first candidate point
	// calculate the square of the distance to current vertex
	squareD = 0;
	for (i=0;i<3;i++) {
		squareD += 
					(current_vertex_coordinates[i]-candidates[0][i])
					*(current_vertex_coordinates[i]-candidates[0][i]);
	}
	store = 0;
	// for each other candidate
	// calculate distance to current vertex
	// if distance is less than stored, replace stored
//	for (i=1;i<count-1;i++) {
	for (i=1;i<count;i++) {
		temp = 0;
		for (j=0;j<3;j++) {
			temp += (current_vertex_coordinates[j]-candidates[i][j])
					*(current_vertex_coordinates[j]-candidates[i][j]);
		}
		if (temp<squareD) {
			store = i;
			squareD = temp;
		}
	}

	for (i=0;i<3;i++) {
		closest_coordinates[i] = candidates[store][i];
	}

}

// #####################################################
// #####################################################

void ManipClass::computeVertexEnergyForce(ContainerClass& container) {

	int i, j, k, adjacent_index;
	int** edge_indices;
	int* num_edge_indices;
	int* max_edge_indices;
	double *x,*y,*z;

    ////////// loop through all objects //////////
    // for each object
    for (i=0;i<object_count;i++) {

		// initialize variables
		edge_indices = new int*[container.objects[i].num_vertices]; 
		if (edge_indices == NULL) { 
			cout << "Not enough memory for edge_indices in ManipClass::computeVertexEnergyForce\n";
			cout.flush();
			exit(1);
		}
		num_edge_indices = new int[container.objects[i].num_vertices];
		if (num_edge_indices == NULL) { 
			cout << "Not enough memory for num_edge_indices in ManipClass::computeVertexEnergyForce\n";
			cout.flush();
			exit(1);
		}
		max_edge_indices = new int[container.objects[i].num_vertices];
		if (max_edge_indices == NULL) { 
			cout << "Not enough memory for max_edge_indices in ManipClass::computeVertexEnergyForce\n";
			cout.flush();
			exit(1);
		}
        for (j=0;j < container.objects[i].num_vertices;j++) {
			edge_indices[j] = new int [NUMBER_OF_EDGES_VERTICES];
			if (edge_indices[j] == NULL) { 
				cout << "Not enough memory for edge_indices[j] in ManipClass::computeVertexEnergyForce\n";
				cout.flush();
				exit(1);
			}
			num_edge_indices[j] = 0;
			max_edge_indices[j] = NUMBER_OF_EDGES_VERTICES;
		}

        ////////// loop through all edges in object //////////
        for (j=0;j < container.objects[i].num_edges;j++) {

			if (container.objects[i].edges[j].valid) {
				addArrayElement(edge_indices[container.objects[i].edges[j].other_vertex_indices[0]],
								num_edge_indices[container.objects[i].edges[j].other_vertex_indices[0]],
								max_edge_indices[container.objects[i].edges[j].other_vertex_indices[0]], 
								INCREMENT_OF_EDGES_VERTICES,j);
				addArrayElement(edge_indices[container.objects[i].edges[j].other_vertex_indices[1]],
								num_edge_indices[container.objects[i].edges[j].other_vertex_indices[1]],
								max_edge_indices[container.objects[i].edges[j].other_vertex_indices[1]], 
								INCREMENT_OF_EDGES_VERTICES,j);
			}

		}

        ////////// loop through all vertices structure elememts in parent object //////////
        // for each vertex in parent object
        for (j=0;j < container.objects[i].num_vertices;j++) {

			if (container.objects[i].vertices[j].candidate) {
				x = new double[container.objects[i].vertices[j].next_adjacent_vertex];
				y = new double[container.objects[i].vertices[j].next_adjacent_vertex];
				z = new double[container.objects[i].vertices[j].next_adjacent_vertex];
			if (x == NULL || y==NULL || z==NULL) { 
				cout << "Not enough memory for x or y or z in ManipClass::computeVertexEnergyForce\n";
				cout.flush();
				exit(1);
			}

				////////// loop through all adjacent vertices of current vertex //////////
				// for each adjacent vertex of current vertex
				for (k = 0;k < container.objects[i].vertices[j].next_adjacent_vertex;k++) {
	
					// current adjacent vertex index
					adjacent_index = container.objects[i].vertices[j].adjacent_vertices_index[k];
	
					// get coordinates of adjacent vertex
					x[k] = container.objects[i].vertices[adjacent_index].new_coordinates[0];
					y[k] = container.objects[i].vertices[adjacent_index].new_coordinates[1];
					z[k] = container.objects[i].vertices[adjacent_index].new_coordinates[2];
	
				}

				////////// compute energy and force //////////

				container.objects[i].vertices[j].computeEnergyForce(x,y,z,edge_indices[j],num_edge_indices[j],
																	container.objects[i].edges);

				// delete x, y, z pointers
				delete[] x;
				delete[] y;
				delete[] z;

			}
		}

		// compute object parameters
		container.objects[i].computeGlobalParams();

		// delete edge_indices pointer
        for (j=0;j < container.objects[i].num_vertices;j++) {
			delete[] edge_indices[j];
		}
		delete[] edge_indices;
		delete[] num_edge_indices;
		delete[] max_edge_indices;

	}

}

// #####################################################
// #####################################################

void ManipClass::computeNewVertexCoords(ContainerClass& container, double gain) {

	int i, j, k;
	double disp[3];

    ////////// loop through all objects //////////
    // for each object
    for (i=0;i<object_count;i++) {

        ////////// loop through all vertices structure elememts in parent object //////////
        // for each vertex in parent object
        for (j=0;j < container.objects[i].num_vertices;j++) {

			if (container.objects[i].vertices[j].candidate) {
 				container.objects[i].vertices[j].computeNewCoords(gain); 
			}

		}
	}
}

// #####################################################
// #####################################################

void ManipClass::computePolygonIntersectionForce(ContainerClass& container,
												PolygonClass *polygons, SpaceClass& space) {

	int i, v1, v2, v3, object_index;
	double nX, nY, nZ, normal_length, fX, fY, fZ;

	// set intersection bool for each polygon
	setPolygonIntersection(container,polygons,space);

	///// set force /////
	// for each polygon
    for (i=0;i<polygon_count;i++) {

		// if intersecting
		if (polygons[i].intersection) {

			// get vertex indices
			v1 = polygons[i].vertices[0];
			v2 = polygons[i].vertices[1];
			v3 = polygons[i].vertices[2];

			// get polygon normal
			nX = polygons[i].normal_components[0];
			nY = polygons[i].normal_components[1];
			nZ = polygons[i].normal_components[2];

			// normal length
			normal_length = sqrt(nX*nX+nY*nY+nZ*nZ);

			// get object index
			object_index = polygons[i].object_index;

			// compute force (opposite direction to normal)
			fX= -nX/normal_length*INTERSECTION_WEIGHT;
			fY= -nY/normal_length*INTERSECTION_WEIGHT;
			fZ= -nZ/normal_length*INTERSECTION_WEIGHT;

			// for each vertex
			// set intersection force
			container.objects[object_index].vertices[v1].intersection_force[0] = fX;
			container.objects[object_index].vertices[v1].intersection_force[1] = fY;
			container.objects[object_index].vertices[v1].intersection_force[2] = fZ;
			container.objects[object_index].vertices[v2].intersection_force[0] = fX;
			container.objects[object_index].vertices[v2].intersection_force[1] = fY;
			container.objects[object_index].vertices[v2].intersection_force[2] = fZ;
			container.objects[object_index].vertices[v3].intersection_force[0] = fX;
			container.objects[object_index].vertices[v3].intersection_force[1] = fY;
			container.objects[object_index].vertices[v3].intersection_force[2] = fZ;
			
		}
	}
}

// #####################################################
// #####################################################

void ManipClass::setPolygonIntersection(ContainerClass& container, PolygonClass *polygons, SpaceClass& space) {

	int i, j, k;
	int* combinations[2];
	int num_combinations;
	int* temp[2];
	int num_temp;
	int next_pair;

	// initialize polygon intersection status
    for (i=0;i<polygon_count;i++) {
		polygons[i].intersection = false;
	}

    // for each subspace
    for (i=0;i<space.num_boxes;i++) {

		// compute number of all possible polygon pairs
		num_combinations = 0;
        for (j=0;j<space.boxes[i].num_polygons;j++) {
			num_combinations += j;
		}

		// allocate memory for combination array
		combinations[0] = new int[num_combinations];
		combinations[1] = new int[num_combinations];

		// load combinations arrays
		num_combinations = 0;
        for (j=0;j<space.boxes[i].num_polygons-1;j++) {
        	for (k=j+1;k<space.boxes[i].num_polygons;k++) {
				combinations[0][num_combinations] = space.boxes[i].polygons[j];
				combinations[1][num_combinations] = space.boxes[i].polygons[k];
				num_combinations++;
			}
		}
		

        // while there are polygon combinations to check 
		next_pair = 0;
		while (next_pair<num_combinations) {

			////// check for intersection of next two polygons //////
			// if polygons intersect
			if (checkPolygonPolygonIntersections(container,polygons,
												combinations[0][next_pair],
												combinations[1][next_pair])) {
				// record intersection
				polygons[combinations[0][next_pair]].intersection = true;
				polygons[combinations[1][next_pair]].intersection = true;
				
				///// remove all instances of each polygon from combinations /////
				// allocate memory for temp array
				temp[0] = new int[num_combinations];
				temp[1] = new int[num_combinations];
				num_temp = 0;

				// for all remaining combinations
				for (k=next_pair+1;k<num_combinations;k++) {
					// if combination does not contain either of the intersection polygons
					if ((combinations[0][k] != combinations[0][next_pair]) &&
						(combinations[1][k] != combinations[0][next_pair]) &&
						(combinations[0][k] != combinations[1][next_pair]) &&
						(combinations[1][k] != combinations[1][next_pair]) ) {
						// add combination to temp
						temp[0][num_temp] = combinations[0][k];
						temp[1][num_temp] = combinations[1][k];
						num_temp++;
					}
				}

				// delete combinations array
				delete[] combinations[0];
				delete[] combinations[1];

				// allocate memory for combination array
				combinations[0] = new int[num_temp];
				combinations[1] = new int[num_temp];
				num_combinations = num_temp;
				next_pair = 0;

				// copy temp to combination
				for (k=0;k<num_temp;k++) {
					combinations[0][k] = temp[0][k];
					combinations[1][k] = temp[1][k];
				}

				// delete temp
				delete[] temp[0];
				delete[] temp[1];

			} else {
				// polygons do not intersect
				next_pair++;
			}

		}

		delete[] combinations[0];
		delete[] combinations[1];

	}
}

// #####################################################
// #####################################################
