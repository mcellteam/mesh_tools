class void_list
{
public:
  void_list *previous;
  void_list *next;
  void *data;
};
/*
class ParametersOld
{
public:
  double plot_rad_int;	// radius of curvature sampling interval for plotting sampling function
  double dmin;			// dmin = minimum spline sampling distance
  double dmax;			// dmax = maximum spline sampling distance
  double tau;			// decrease tau to increase steepness of sampling function
  double rI;			// inflection point of sampling function = tau*rI
						// decrease rI to sample finer
  double max_rad;		// calculated radius of curvature of spline will saturate at this value
  bool diag;			// set true to print diagnostic files
  int num;				// num is the # samples of the splines between contour points
						// sequential triplets of sampled points are used to
						// compute the radius of curvature , e.g. num/2 between contour points	
};
*/
class Parameters
{
public:
  double plot_rad_int;	// radius of curvature sampling interval for plotting sampling function
  double dmin;			// minimum spline sampling distance
  double dmax;			// maximum spline sampling distance
  double max_rad;		// calculated radius of curvature of spline will saturate at this value
  bool diag;			// set true to print diagnostic files
  int num;				// num is the # samples of the splines between contour points
						// sequential triplets of sampled points are used to
						// compute the radius of curvature , e.g. num/2 between contour points
  double T;				// sample period (= time to traverse dmax)
  double amax;			// max radial acceleration
};

class Object
{
public:
  char name[128];
  int min_section,max_section;
  void_list *contours;
  Object(char *str,int sec);
  ~Object(void);
};

Object::Object(char *str,int sec)
{
	strcpy(name,str);
	contours = NULL;
	min_section = sec;
	max_section = sec;
}

Object::~Object(void){
	void_list *p,*q;
	p=contours;
	while (p!=NULL) {
		q=p->next;
		delete p;
		p=q;
	}
}

class Point
{
public:
  double x,y,z;
  Point(char *str,int section,double thickness);
  Point(double xval,double yval,double zval);
};

Point::Point(double xval, double yval, double zval)
{
	x = xval;
	y = yval;
	z = zval;
}

Point::Point(char *str, int section, double thickness)
{
	char val[80];
	char *eptr;
	int i;

	// set z coordinate
	z = section*thickness;

	// get past 'points'
  	while (strchr(" points=\"\t",*str)!=NULL) {str++;}

	// grab x coordinate
	i=0;
	while (strchr("0123456789+-eE.",*str)!=NULL)
	{
		val[i++] = *str++;
	}
	val[i]=0;
	x = strtod(val,&eptr);
	if (val==eptr)
	{
		x=y=z=0;
		printf("Error in reading x coordinate\n");
	printf("str =%s\n",str);
		return;
	}

	// grab y coordinate
	while (strchr(" \t,",*str)!=NULL) { str++; }
	i=0;
	while (strchr("0123456789+-eE.",*str)!=NULL)
	{
		val[i++] = *str++;
	}
	val[i]=0;
	y = strtod(val,&eptr);
	if (val==eptr)
	{
		x=y=z=0;
		printf("Error in reading y coordinate\n");
		return;
	}
}

class Contour
{
public:
  char name[128];
  int section,num_interp_points;
  void_list *raw_points,*rawend;
  void_list *interp_points,*interpend;
  double *deviations;
  Contour(char* str, int sec);
  ~Contour(void);
  void removeDuplicates(void);
  void clearSpline(void);
  void linearlyInterp(double,double);
  void addPreviousRaw(void);
  void addPreviousInterp(void);
};

Contour::~Contour(void){
	void_list *p,*q;
	p=raw_points;
	while (p!=NULL) {
		q=p->next;
		delete (Point*)p->data;
		delete p;
		p=q;
	}
	p=interp_points;
	while (p!=NULL) {
		q=p->next;
		delete (Point*)p->data;
		delete p;
		p=q;
	}
	delete[] deviations;
}

Contour::Contour(char *str, int sec)
{
	char val[80];
	int i;
	// grab name
	i=0;
	while (strchr("\"",*str)==NULL){val[i++]=*str++;}
	val[i]=0;
	strcpy(name,val);
	section = sec;
	raw_points = NULL;
	rawend=NULL;
	interp_points = NULL;
}

void Contour::clearSpline(void)
{
	void_list *p,*q;
	p=interp_points;
	while(p!=NULL) {
		delete (Point*)p->data;
		q=p->next;
		delete p;
		p=q;
	}
	interp_points = NULL;
	interpend = NULL;
}

void Contour::addPreviousRaw(void){
	void_list *q,*prev;
	prev=NULL;
	rawend=NULL;
	// for each point
	for (q=raw_points;q!=NULL;q=q->next) {
		q->previous = prev;
		prev = q;
		rawend = q;
	}
}

void Contour::addPreviousInterp(void){
	void_list *q,*prev;
	prev=NULL;
	interpend=NULL;
	// for each point
	for (q=interp_points;q!=NULL;q=q->next) {
		q->previous = prev;
		prev = q;
		interpend = q;
	}
}

void_list * removeLink(void_list* L) {
	void_list *q;
    // and remove face from candidate face list
    if (L->previous!=NULL) {
        if (L->next!=NULL) {
            // if both previous and next exist
            (L->previous)->next = L->next;
            (L->next)->previous = L->previous;
        } else {
            // if previous exists and next does not
            (L->previous)->next = NULL;
        }
    } else {
        if (L->next!=NULL) {
            // if previous does not exist and next does
            (L->next)->previous = NULL;
        } // else { // if neither previous nor next exists }
    }
	// update pointer
	q=L->next;
	delete L;
	return q;
}

void_list * deletePoint(void_list *q,void_list *p,void_list *&ptr){
	Point *pt1,*pt2;
	pt1=(Point*)q->data;
	pt2=(Point*)p->data;
	// if points are identical
	if ((pt1->x==pt2->x)&&(pt1->y==pt2->y)){
		// delete point
		delete pt1;
		// adjust list pointer
		if (q==ptr) {ptr = q->next; }
		// remove current point from list
		q=removeLink(q);
	}
	return p;
}

void Contour::removeDuplicates(void)
{
	void_list *q,*ptr;
	Point *pt1,*pt2;
	ptr = raw_points;
	q=raw_points;
	while (q->next!=NULL) {q=deletePoint(q,q->next,ptr);}
    // adjust pointer
    raw_points = ptr;

	// compare first and last link in list
	q=deletePoint(raw_points,rawend,ptr);
    // adjust pointer
    raw_points = ptr;
}

void_list * interpPoints(void_list *q,void_list *p,void_list *&ptr,
							double maxdev,double scale,int flag){
	void_list *pp;
	double dist,distx,disty,x,y,count;
	int num;
	Point *v,*pt1,*pt2;
	pt1=(Point*)q->data;
	pt2=(Point*)p->data;
	// compute distance between points
	distx = pt2->x-pt1->x;
	disty = pt2->y-pt1->y;
	dist = sqrt(distx*distx+disty*disty);
	num = (int)(dist/(maxdev/scale)/3);
	if (num) {
		// linearly interpolate num evenly spaced points
		count = num+1;
		ptr = q;
		while(num){
			// insert point
			x = pt1->x+(count-(double)num)/count*distx;
			y = pt1->y+(count-(double)num)/count*disty;
			v = new Point(x,y,pt1->z);
			pp = new void_list();
			pp->next = p;
			pp->previous = ptr;
			pp->data = (void*)v;
			ptr->next = pp;
			p->previous = pp;
			// decrement num
			num--;
			ptr = pp;
		}
		if(flag){ptr->next=NULL;}
	}
	return p;
}

void Contour::linearlyInterp(double maxdev,double scale)
{
	void_list *q,*ptr=NULL;
	q=raw_points;
	while (q->next!=NULL) {
		q=interpPoints(q,q->next,ptr,maxdev,scale,0);
	}
	// compare first and last link in list
	q=interpPoints(rawend,raw_points,ptr,maxdev,scale,1);
	// adjust pointer
	rawend=ptr;
}

class Histogram
{
public:
  double min,max,mean,stddev,sum;
  int count[16],num;
  Histogram(void);
  void update(void_list*,int);
};

Histogram::Histogram(void)
{
  int i;
  for(i=0;i<16;i++){
	count[i]=0;
  }
  min=1e30;
  max=0;
  mean=stddev=sum=0;
  num=0;
}

void Histogram::update(void_list* q,int count) {
	///// update deviation distance statistics /////
	int i;
	double d;
	Contour* c=(Contour*)q->data;
	for (i=0;i<count;i++) {
		d=c->deviations[i];
		// update min and max deviation distance
		if(d<min) {min=d;}
		if(d>max) {max=d;}
		num++;
		sum+=d;
	}
}

class SplinePoint
{
public:
  double t,x,y,r,intfac;
};

class Weights
{
public:
	double *bx,*by;
};

