#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define EPS_C 5e-7
#define EPS_C 1e-15

class void_list
{
public:
  struct void_list *next;
  void *data;
};

int distinguishable(double a,double b,double eps)
{

  double c;

  c=a-b;

  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;

  if (b<a) eps*=a;
  else eps*=b;
  return (c>eps);
}

struct void_list* void_list_sort_by(struct void_list *vl,int (*leq)(void*,void*))
{
  struct void_list *stack[64];
  int stack_n[64];
  struct void_list *left,*right,*merge,*tail;
  int si = 0;
  
  while (vl != NULL)
  {
    if (vl->next == NULL)
    {
      stack[si] = vl;
      stack_n[si] = 1;
      vl = NULL;
      si++;
    }
    else if ( (*leq)(vl->data , vl->next->data) )
    {
      stack[si] = vl;
      stack_n[si] = 2;
      vl = vl->next->next;
      stack[si]->next->next = NULL;
      si++;
    }
    else
    {
      stack[si] = vl->next;
      stack_n[si] = 2;
      left = vl;
      vl = vl->next->next;
      stack[si]->next = left;
      left->next = NULL;
      si++;
    }
    while (si > 1 && stack_n[si-1]*2 >= stack_n[si-2])
    {
      stack_n[si-2] += stack_n[si-1];

      left = stack[si-2];
      right = stack[si-1];
      if ((*leq)(left->data , right->data)) { merge = left; left = left->next; }
      else { merge = right; right = right->next; }
      merge->next = NULL;
      tail = merge;

      while (1)
      {
        if (left==NULL)
        {
          tail->next = right; tail = right;
          break;
        }
        if (right==NULL)
        {
          tail->next = left; tail = left;
          break;
        }

        if ((*leq)(left->data , right->data))
        { 
          tail->next = left; tail = left; left = left->next;
        }
        else
        { 
          tail->next = right; tail = right; right = right->next; 
        }
      }
      
      stack[si-2] = merge;
      si--;   
    }
  }
  
  while (si > 1)  /* Exact duplicate of code in loop--keep it this way! */
  {
    stack_n[si-2] += stack_n[si-1];

    left = stack[si-2];
    right = stack[si-1];
    if ((*leq)(left->data , right->data)) { merge = left; left = left->next; }
    else { merge = right; right = right->next; }
    merge->next = NULL;
    tail = merge;

    while (1)
    {
      if (left==NULL)
      {
        tail->next = right; tail = right;
        break;
      }
      if (right==NULL)
      {
        tail->next = left; tail = left;
        break;
      }

      if ((*leq)(left->data , right->data))
      { 
        tail->next = left; tail = left; left = left->next;
      }
      else
      { 
        tail->next = right; tail = right; right = right->next; 
      }
    }
    
    stack[si-2] = merge;
    si--;   
  }
  
  return stack[0];
}


class Vertex
{
public:
  double x,y,z;
  int my_idx;
  void_list *owners;
  Vertex(char *triplet);
  bool near(Vertex &v);
  void steal(Vertex *v);
  int cmp(Vertex *v);
};

Vertex::Vertex(char *triplet)
{
  char val[80];
  char *eptr;
  int i;
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
  owners = NULL;
}

bool Vertex::near(Vertex &v)
{
  return (!distinguishable(x,v.x,EPS_C) &&
          !distinguishable(y,v.y,EPS_C) &&
          !distinguishable(z,v.z,EPS_C) );
}

int Vertex::cmp(Vertex *v)
{
  if (x < v->x && distinguishable(x,v->x,EPS_C)) return -1;
  else if (x > v->x && distinguishable(x,v->x,EPS_C)) return 1;
  else
  {
    if (y < v->y && distinguishable(y,v->y,EPS_C)) return -1;
    else if (y > v->y && distinguishable(y,v->y,EPS_C)) return 1;
    else
    {
      if (z < v->z && distinguishable(z,v->z,EPS_C)) return -1;
      else if (z > v->z && distinguishable(z,v->z,EPS_C)) return 1;
      return 0;
    }
  }
}


int leq(void *a,void *b)
{
  Vertex *v1 = (Vertex*)a;
  Vertex *v2 = (Vertex*)b;
  return v1->cmp(v2) < 1;
}

class Triangle
{
public:
  Vertex *v1,*v2,*v3;
};

void Vertex::steal(Vertex *v)
{
  Triangle *t;
  void_list *ll;
  for (ll=v->owners ; ll!=NULL ; ll=ll->next)
  {
    t = (Triangle*)ll->data;
    if (t->v1==v) t->v1=this;
    if (t->v2==v) t->v2=this;
    if (t->v3==v) t->v3=this;
  }
}

int compare (const void* a, const void* b ) {

	double i;
	void_list *j,*k;

	j = *(void_list**)a;
	k = *(void_list**)b;

	i = ((Vertex*)j->data)->x - ((Vertex*)k->data)->x;
	if (i<0) {	return 1;}
	else if (i>0) {	return -1;}
	else {
		i = ((Vertex*)j->data)->y - ((Vertex*)k->data)->y;
		if (i<0) {	return 1;}
		else if (i>0) {	return -1;}
		else {
			i = ((Vertex*)j->data)->z - ((Vertex*)k->data)->z;
			if (i<0) {	return 1;}
			else if (i>0) {	return -1;}
			else { return (0);}
		}
	}
}

int main(int argc,char *argv[])
{
  char *infile;
  char line[2048];
  char *str;
  void_list *vlh,*vl;
  Vertex *v;
  void_list *tlh,*tl;
  Triangle *t;
  void_list *q,*qq;
  FILE *F;
  if (argc != 2)
  {
    printf("Syntax: poly2mesh input_file\n");
    return 1;
  }
  
  infile = argv[1];
  F = fopen(infile,"r");
  if (!F)
  {
    printf("Couldn't open input file %s\n",infile);
    return 1;
  }
  
  vlh = NULL;
  tlh = NULL;
  int n = 0;
  int num_vert=0;  // JPK
  for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F))
  {
    if (!strcmp(str,"3\n"))
    {
      tl = new void_list();
      tl->next = tlh;
      tlh = tl;
      t = new Triangle();
      tl->data = (void*)t;
      
      str=fgets(line,2048,F);
      if (str==NULL) { printf("Error!\n"); return 1; }
      vl = new void_list();
      vl->next = vlh;
      v = new Vertex(str);
      vl->data = (void*)v;
      q = new void_list();
      q->next = v->owners;
      v->owners = q;
      q->data = (void*)t;
      t->v1 = v;
      vlh=vl;
      
      str=fgets(line,2048,F);
      if (str==NULL) { printf("Error!\n"); return 1; }
      vl = new void_list();
      vl->next = vlh;
      v = new Vertex(str);
      vl->data = (void*)v;
      q = new void_list();
      q->next = v->owners;
      v->owners = q;
      q->data = (void*)t;
      t->v2 = v;
      vlh = vl;
      
      str=fgets(line,2048,F);
      if (str==NULL) { printf("Error!\n"); return 1; }
      vl = new void_list();
      vl->next = vlh;
      v = new Vertex(str);
      vl->data = (void*)v;
      q = new void_list();
      q->next = v->owners;
      v->owners = q;
      q->data = (void*)t;
      t->v3 = v;
      vlh=vl;

		num_vert+=3; //JPK
    }
  }
  fflush(stdout);

	// create array of pointers to vertex list
	//printf("# vertices %i\n",num_vert);
	void_list** vert_array;
	vert_array = new void_list*[num_vert];
	num_vert=0;
	for (q=vlh;q!=NULL;q=q->next) {
		vert_array[num_vert++]=q;
	}


	qsort(vert_array,num_vert,sizeof(void_list*),compare);

	///// turn array into linked list /////
	void_list *sort_h,*sort;
	sort_h = NULL;
	int j;
	for (j=0;j<num_vert;j++) {
      	sort = new void_list();
		sort->next = sort_h;
		sort->data = (void*)((void_list*)vert_array[j]->data);
		sort_h = sort;
	}
 
  for (q=sort_h;q!=NULL;q=q->next)
  {
    if (q->next==NULL) break;
    while(((Vertex*)q->data)->cmp((Vertex*)q->next->data)==0)
    {
      ((Vertex*)(q->data))->steal((Vertex*)q->next->data);
      q->next = q->next->next;
      if (q->next==NULL) break;
    }
  }
  
  int i;
  for (q=sort_h,i=1;q!=NULL;q=q->next,i++)
  {
    v = (Vertex*)q->data;
    v->my_idx = i;
    printf("Vertex %d %.15g %.15g %.15g\n",v->my_idx,v->x,v->y,v->z);
  }
  
  
  for (tl=tlh,i=1;tl!=NULL;tl=tl->next,i++)
  {
    t = (Triangle*)tl->data;
    printf("Face %d %d %d %d\n",i,t->v1->my_idx,t->v2->my_idx,t->v3->my_idx);
  }
  fclose(F);
  return 0;
}
