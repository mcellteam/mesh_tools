#include "face.h"
#include "group.h"
#include "meshstitch.h"

Group::Group (void)
:count(0),c(NULL)
{
  count=0;
  c=NULL;
};

Group::~Group (void)
{
  delete[] c;
};

void Group::extraVertices (std::list<Vertex> & S)
{
  // for each group
  for (int i=0;i<count;i++)
  {
    // for each vertex in group
    for (iter_list_v q=c[i].verts.begin();q!=c[i].verts.end();q++)
    {
      bool found=false;
      // for each shared vertex
      iter_list_v p = S.begin();
      while(p!=S.end() && !found)
      {
        // if group vertex is shared
        if (q->index==p->index)
        {
          found=true;
        }
        p++;
      }
      // if group vertex is not shared
      if (!found)
      {
        // then it is extra vertex
        ExtraVertex EV;
        EV.v = &(*q);
        c[i].extra.push_back(EV);
      }
    }
  }
}


void Group::groupExtraVertices (void)
{
  // group extra vertices in each contour into sites
  ExtraVertex *ev1=NULL,*ev2=NULL,*ev3=NULL,*eva=NULL,*evb=NULL;
  ///// handle sites with two or more extra vertices
  //for each contour group
  for (int i=0;i<count;i++)
  {
    bool flag=true;
    int next_group=0;
    // while changes are made
    while(flag)
    {
      flag=false;
      // for each face in contour group
      for (iter_list_f q=c[i].faces.begin();q!=c[i].faces.end();q++)
      {
        int g1=0,g2=0,g3=0,ga=0,gb=0;
        // for each extra vertex in contour group
        for (iter_vec_ev p=c[i].extra.begin();p!=c[i].extra.end();p++)
        {
          // are any vertices extra?
          if ( q->v1 == p->v->index)
          {
            ev1=&(*p);
            g1=ev1->g;
          }
          else if ( q->v2 == p->v->index)
          {
            ev2=&(*p);
            g2=ev2->g;
          }
          else if ( q->v3 == p->v->index)
          {
            ev3=&(*p);
            g3=ev3->g;
          }
        }
        // were at least two vertices extra?
        if ((g1&&g2)||(g1&&g3)||(g2&&g3))
        {
          if (g1)
          {
            ga=g1;
            eva=ev1;
          }
          if (g2)
          {
            if (!ga)
            {
              ga=g2;
              eva=ev2;
            }
            else
            {
              gb=g2;
              evb=ev2;
            }
          }
          if (g3)
          {
            gb=g3;
            evb=ev3;
          }
          // check to see either extra vertex is in a group
          if (ga<0 && gb<0)
          {
            // if both ungrouped, then add both to same group
            eva->g=next_group;
            evb->g=next_group;
            next_group++;
            flag=true;
          }
          else if (ga==gb)
          {
            // if same and greater than zero, then do nothing
          }
          else if (ga<0 && gb>0)
          {
            eva->g=evb->g;
            flag=true;
          }
          else if (ga>0 && gb<0)
          {
            evb->g=eva->g;
            flag=true;
          }
          else
          {
            int largest,smallest;
            flag=true;
            // if different and greater than zero
            // then identify smallest group #
            if (ga<gb)
            {smallest=ga;largest=gb;}
            else {smallest=gb;largest=ga;}
            // then set any extra vertex with largest group # to smallest group #
            for (iter_vec_ev p=c[i].extra.begin();p!=c[i].extra.end();p++)
            {
              if (p->g==largest)
              {
                p->g=smallest;
              }
            }
          }
        }
      }
    }
    ///// handle sites with only a single extra vertex
    flag=true;
    // while changes are made
    while(flag)
    {
      flag=false;
      // for each face in contour group
      for (iter_list_f q=c[i].faces.begin();q!=c[i].faces.end();q++)
      {
        int g1=0,g2=0,g3=0;
        // for each extra vertex in contour group
        for (iter_vec_ev p=c[i].extra.begin();p!=c[i].extra.end();p++)
        {
          // are any vertices extra?
          if ( q->v1 == p->v->index)
          {
            ev1=&(*p);
            g1=ev1->g;
          }
          else if ( q->v2 == p->v->index)
          {
            ev2=&(*p);
            g2=ev2->g;
          }
          else if ( q->v3 == p->v->index)
          {
            ev3=&(*p);
            g3=ev3->g;
          }
        }
        // were any vertices extra?
        if (g1||g2||g3)
        {
          // check to see if extra vertex is in a group
          if (g1<0)
          {
            // if both ungrouped, then add both to same group
            ev1->g=next_group;
            next_group++;
            flag=true;
          } else if (g2<0)
          {
            // if both ungrouped, then add both to same group
            ev2->g=next_group;
            next_group++;
            flag=true;
          } else if (g3<0)
          {
            // if both ungrouped, then add both to same group
            ev3->g=next_group;
            next_group++;
            flag=true;
          }
        }
      }
    }
  }

  ///// count number of sites per contour
  // for each contour group
  for (int i=0;i<count;i++)
  {
    int size=0;	
    int num=0;
    // for each vertex in contour
    for (iter_list_v q=c[i].verts.begin();q!=c[i].verts.end();q++)
    {size++;}
    // create int array of size = number of vertices per contour+1
    int *array = new int[size+1];
    for (int j=0;j<size+1;j++)
    {
      array[j]=0;
    }
    // for each extra vertex in contour group
    for (iter_vec_ev p=c[i].extra.begin();p!=c[i].extra.end();p++)
    {
      // record group number in array
      array[p->g]=1;
    }
    for (int j=0;j<size+1;j++)
    {
      if (array[j])
      {num++;}
    }
    // record number of sites in contour
    c[i].setNum(num);
    // cleanup
    delete[] array;
  }
}

void Group::initializeSites (void)
{
  // for each contour
  for (int i=0;i<count;i++)
  {
    // initialize sites
    c[i].initSites();
  }
}

void Group::orderSites (void)
{
  Vertex *target = NULL;
  // for each contour
  for (int i=0;i<count;i++)
  {
    // for each site in contour
    for (int j=0;j<c[i].getNum();j++)
    {
      ///// identify site ev_0 and ev_n-1
      // for each extra vertex
      for (iter_vec_ev p=c[i].extra.begin();p!=c[i].extra.end();p++)
      {
        ExtraVertex *ev=&(*p);
        bool found=false;
        // if extra vertex group number matches site group
        if (ev->g==j)
        {
          // if extra vertex laterals, la or lb, match site laterals l1 or l2
          // then assign extra vertex to site ev in appropriate place if ev==NULL
          if (((Vertex*)ev->la)->index==((Vertex*)c[i].s[j].l1)->index)
          {
            if (c[i].s[j].ev[0]==NULL)
            {
              c[i].s[j].ev[0]=ev;
              found=true;
              target=(Vertex*)ev->lb;
            }
          }
          if (!found && ((Vertex*)ev->la)->index==((Vertex*)c[i].s[j].l2)->index)
          {
            if (c[i].s[j].ev[c[i].s[j].n-1]==NULL)
            {
              c[i].s[j].ev[c[i].s[j].n-1]=ev;	
              found=true;
            }
          }
          if (!found && ((Vertex*)ev->lb)->index==((Vertex*)c[i].s[j].l1)->index)
          {
            if (c[i].s[j].ev[0]==NULL)
            {
              c[i].s[j].ev[0]=ev;
              found=true;
              target=(Vertex*)ev->la;
            }
          }
          if (!found && ((Vertex*)ev->lb)->index==((Vertex*)c[i].s[j].l2)->index)
          {
            if (c[i].s[j].ev[c[i].s[j].n-1]==NULL)
            {
              c[i].s[j].ev[c[i].s[j].n-1]=ev;	
              found=true;
            }
          }
        }
      }

      ///// identify site ev_1 to ev_n-2
      // for site extra vertex slots ev_1 to ev_n-2
      for (int k=1;k<c[i].s[j].n-1;k++)
      {
        // for each extra vertex
        for (iter_vec_ev p=c[i].extra.begin();p!=c[i].extra.end();p++)
        {
          ExtraVertex *ev=&(*p);
          // if ev->v matches target 
          if (((Vertex*)ev->v)->index==target->index)
          {
            // record extra vertex in site
            c[i].s[j].ev[k]=ev;
            // set new target
            // if ev->la matches k-1
            if (((Vertex*)ev->la)->index==((Vertex*)((ExtraVertex*)c[i].s[j].ev[k-1])->v)->index)
            {
              // then ev->lb is new target
              target=(Vertex*)ev->lb;
            }
            // else ev->lb is new target
            else {target=(Vertex*)ev->la;}
          }
        }
      }
    }
  }
}

void Group::lateralVertices (std::list<Face> & S,int z_value,Vertex **vert_array)
{
  // group extra vertices in each contour into sites
  groupExtraVertices();
  // initialize sites
  initializeSites();
  ///// load sites
  // for each contour group
  for (int i=0;i<count;i++)
  {
    // for each extra vertex in contour group
    for (iter_vec_ev p=c[i].extra.begin();p!=c[i].extra.end();p++)
    {
      ExtraVertex *ev=&(*p);
      ///// identify lateral vertices
      // for each group face
      for (iter_list_f q=c[i].faces.begin();q!=c[i].faces.end();q++)
      {
        Vertex *v1=vert_array[q->v1];
        Vertex *v2=vert_array[q->v2];
        Vertex *v3=vert_array[q->v3];
        // if face contains extra vertex
        if (q->vertexInFace(p->v))
        {
          // if z_value matches and index is not same as extra vertex
          // remember need to account for 2 pairs of duplicate lateral vertices
          if ((v1->z==z_value)&&(q->v1!=((Vertex*)ev->v)->index))
          {
            // this vertex is lateral
            if (ev->la==NULL)
            {
              ev->la=v1;
            }
            else if (ev->lb==NULL&&(ev->la)->index!=q->v1)
            {
              ev->lb=v1;
            }
          }
          if ((v2->z==z_value)&&(q->v2!=(ev->v)->index))
          {
            // this vertex is lateral
            if (ev->la==NULL)
            {
              ev->la=v2;
            }
            else if (ev->lb==NULL&&(ev->la)->index!=q->v2)
            {
              ev->lb=v2;
            }
          }
          if ((v3->z==z_value)&&(q->v3!=(ev->v)->index))
          { 
            // this vertex is lateral
            if (ev->la==NULL)
            {
              ev->la=v3;
            }
            else if (ev->lb==NULL&&(ev->la)->index!=q->v3)
            {
              ev->lb=v3;
            }
          }
        }
      }

      // for each shared vertex
      bool found = true;
      iter_list_f q=S.begin();
      while(q!=S.end() && found)
      {
        // if lateral is shared
        if (((Vertex*)(ev->la))->index==q->index)
        {found=false;}
        q++;
      }
      // if lateral is shared
      if (!found)
      {
        // then save extra vertex lateral as site lateral 1 or 2
        if (c[i].s[ev->g].l1==NULL)
        {c[i].s[ev->g].l1=(Vertex*)ev->la;}
        else {c[i].s[ev->g].l2=(Vertex*)ev->la;}
      }
      // for each shared vertex
      found = true;
      q=S.begin();
      while(q!=S.end() && found)
      {
        // if lateral is shared
        if (((Vertex*)(ev->lb))->index==q->index)
        {found=false;}
        q++;
      }
      // if lateral is shared
      if (!found)
      {
        // then save extra vertex lateral as site lateral 1 or 2
        if (c[i].s[ev->g].l1==NULL)
        {c[i].s[ev->g].l1=(Vertex*)ev->lb;}
        else {c[i].s[ev->g].l2=(Vertex*)ev->lb;}
      }
    }
  }

  // determine the sequence of extra vertices in sites
  orderSites();
}

void Group::addVertex (std::list<Vertex> & V,Vertex *EV,int MV)
{
  // add extra vertex
  char buffer[128];
  sprintf (buffer,"Vertex %d %g %g %g",MV,EV->x,EV->y,EV->z);
  V.push_back(Vertex(buffer));
  if (0) {printf("Vertex %i %.15g %.15g %.15g\n",MV,EV->x,EV->y,EV->z);}
}

void Group::addFaces (std::list<Face> & F,int v1,int v2, int th,int orient,int &MF)
{
  // add faces
  MF++;
  char buffer[128];
  if (orient)
  {
    sprintf (buffer, "Face %d %d %d %d",MF,v1,v2,th);
    if (0) {printf("Face %i %i %i %i\n",MF,v1,v2,th);}
  }
  else
  {
    sprintf (buffer, "Face %d %d %d %d",MF,v2,v1,th);
    if (0) {printf("Face %i %i %i %i\n",MF,v2,v1,th);}
  }
  F.push_back(Face(buffer));
}

int Group::faceContainsLaterals (Vertex *v1,Vertex *v2,Face *q,int val,Vertex **vert_array)
{
  return ((v1->index==q->v1||v1->index==q->v2||v1->index==q->v3)&&
          (v2->index==q->v1||v2->index==q->v2||v2->index==q->v3)
          &&( vert_array[q->v1]->z!=val && vert_array[q->v2]->z!=val && vert_array[q->v3]->z!=val)
         );
}

void Group::convertLateral (std::list<Vertex> & V,double epsilon)
{
  // convert lateral vertices to other object
  // for each contour
  for (int i=0;i<count;i++)
  {
    // for each site in contour
    for (int j=0;j<c[i].getNum();j++)
    {
      Vertex *v1=c[i].s[j].l1;
      Vertex *v2=c[i].s[j].l2;
      // for each vertex of other object
      bool f1=false,f2=false;
      iter_list_v p=V.begin();
      while (p!=V.end() && (!f1 || !f2))
      {
        Vertex *o=&(*p);
        // for both site lateral vertices
        if ( !distinguishable(o->x,v1->x,epsilon) &&
             !distinguishable(o->y,v1->y,epsilon) &&
             !distinguishable(o->z,v1->z,epsilon)
           )
        {c[i].s[j].l1 = o;f1=true;}
        if ( !distinguishable(o->x,v2->x,epsilon) &&
             !distinguishable(o->y,v2->y,epsilon) &&
             !distinguishable(o->z,v2->z,epsilon)
           )
        {c[i].s[j].l2= o;f2=true;}
        p++;
      }
    }
  }
}

void Group::gatherThirdDeleteFace (std::list<Face> & F,int val,Vertex **vert_array)
{
  // for each contour
  for (int i=0;i<count;i++)
  {
    // for each site in contour
    for (int j=0;j<c[i].getNum();j++)
    {
      Vertex *v1=c[i].s[j].l1;
      Vertex *v2=c[i].s[j].l2;
      // for each face in other object
      bool found = false;
      iter_list_f p=F.begin();
      while(p!=F.end() && !found)
      {
        Face *f=&(*p);
        // if face contains both lateral vertices and z != val
        if (faceContainsLaterals(v1,v2,f,val,vert_array))
        {
          found = true;
          // identify third vertex
          if ((v1->index!=f->v1)&&(v2->index!=f->v1))
          {c[i].s[j].th=vert_array[f->v1];}
          else if ((v1->index!=f->v2)&&(v2->index!=f->v2))
          {c[i].s[j].th=vert_array[f->v2];}
          else {c[i].s[j].th=vert_array[f->v3];}
          // identify orientation
          if (((v1->index==f->v1)&&(v2->index==f->v2))||((v1->index==f->v3)&&
                                                         (v2->index==f->v1))||((v1->index==f->v2)&&(v2->index==f->v3)))
          {c[i].s[j].orient= 1;}
          else {c[i].s[j].orient=0;}
          // delete face link
          delete f;
          F.erase(p);
        }
        else
        {
          p++;
        }
      }
    }
  }
}

void Group::addFacesVertices (std::list<Face> & F,std::list<Vertex> & V,int max_vertex,int max_face)
{
  // for each contour
  for (int i=0;i<count;i++)
  {
    // for each site in contour
    for (int j=0;j<c[i].getNum();j++)
    {
      Vertex *th=c[i].s[j].th;
      // add face with l1, site ev0, and third
      int orient=c[i].s[j].orient;
      Vertex *v1=c[i].s[j].l1;
      Vertex *v2=c[i].s[j].ev[0]->v;
      max_vertex++;
      addFaces(F,v1->index,max_vertex,th->index,orient,max_face);
      addVertex(V,v2,max_vertex);
      iter_list_v v = V.end();
      c[i].s[j].ev[0]->v=&(*(--v));
      // add face with site ev_n-1, l2, and third
      if (c[i].s[j].n>1)
      {
        v1=c[i].s[j].ev[c[i].s[j].n-1]->v;
        v2=c[i].s[j].l2;
        max_vertex++;
        addFaces(F,max_vertex,v2->index,th->index,orient,max_face);
        addVertex(V,v1,max_vertex);
        v = V.end();
        c[i].s[j].ev[c[i].s[j].n-1]->v=&(*(--v));
      }
      else
      {
        v1=c[i].s[j].ev[c[i].s[j].n-1]->v;
        v2=c[i].s[j].l2;
        addFaces(F,max_vertex,v2->index,th->index,orient,max_face);
      }
      // add more faces
      if (c[i].s[j].n>1)
      {
        for (int k=0;k<c[i].s[j].n-2;k++)
        {
          v1=c[i].s[j].ev[k]->v;
          v2=c[i].s[j].ev[k+1]->v;
          max_vertex++;
          addFaces(F,v1->index,max_vertex,th->index,orient,max_face);
          addVertex(V,v2,max_vertex);
          v = V.end();
          c[i].s[j].ev[k+1]->v=&(*(--v));
        }
        int k=c[i].s[j].n-2;
        v1=c[i].s[j].ev[k]->v;
        v2=c[i].s[j].ev[k+1]->v;
        addFaces(F,v1->index,max_vertex,th->index,orient,max_face);
      }
    }
  }
}

void Group::getContours (std::list<Face> & evfsp,
                         Vertex **vert_array,
                         int max_vert,
                         int z_value,
                         std::list<Vertex> & S)
{
  // Collect vertices and faces for each contour that has shared vertices
  // between the two mesh files
  bool flag=true;
  //create array of ints
  int array[max_vert+1];
  for (int i=0;i<max_vert+1;i++)
  {
    array[i]=0;
  }
  int next_group=1;

  while (flag)
  {
    flag=false;
    //for each face
    for (iter_list_f p=evfsp.begin();p!=evfsp.end();p++)
    {
      int va=0,vb=0;
      // are any vertices at z_value?
      if (vert_array[p->v1]->z==z_value)
      {
        va=p->v1;
      }
      if (vert_array[p->v2]->z==z_value)
      {
        if (!va) {va=p->v2;}
        else{vb=p->v2;}
      }
      if (vert_array[p->v3]->z==z_value)
      {
        vb=p->v3;
      }
      if (va && vb)
      {
        // check to see if any are in a group
        int ga=array[va];
        int gb=array[vb];
        if (!ga && !gb)
        {
          // if both zero, then add both to same group
          array[va]=array[vb]=next_group;
          next_group++;
          flag=true;
        }
        else if (ga==gb)
        {
          // if same and nonzero, then do nothing
        }
        else if (!ga && gb)
        {
          array[va]=array[vb];
          flag=true;
        }
        else if (ga && !gb)
        {
          array[vb]=array[va];
          flag=true;
        }
        else
        {
          flag=true;
          // if different and nonzero
          // then identify smallest group #
          int smallest,largest;
          if (ga<gb)
          {smallest=ga;largest=gb;}
          else {smallest=gb;largest=ga;}
          // then set any vertex with other group #s to smallest group #
          for (int i=0;i<max_vert+1;i++)
          {
            if (array[i]==largest) array[i]=smallest;
          }
        }
      }
    }
  }

  ///// check each group for shared vertices
  // create array of size next_group
  int shared[next_group];
  for (int i=0;i<next_group;i++)
  {
    shared[i]=0;
  }

  int num=1;
  // for each array
  for (int i=0;i<max_vert+1;i++)
  {
    // if group nonzero
    if (array[i])
    {
      flag = true;
      iter_list_v q=S.begin();
      while (q!=S.end() && flag)
      {
        // if vertex is shared
        if (q->index==i)
        {
          // then set shared[group] to 1
          if (!shared[array[i]])
          {
            shared[array[i]]=num++;
          }
          flag = false;
        }
        q++;
      }
    }
  }
  num--;

  // create contours
  count = num;
  c = new Contour[num];
  if (c[0].getNum()!=0)
  {printf("THIS SHOULD NOT HAVE PRINTED!!!!");}
  if (c[0].extra.size()!=0)
  {printf("THIS SHOULD NOT HAVE PRINTED!!!!");}

  // load vertices into contours
  // for each vertex 
  for (int i=0;i<max_vert+1;i++)
  {
    if (shared[array[i]])
    {
      c[shared[array[i]]-1].verts.push_back(*vert_array[i]);
    }
  }

  // load face into contours
  //for each face
  for (iter_list_f p=evfsp.begin();p!=evfsp.end();p++)
  {
    flag = false;
    // for each vertex group
    int i=0;
    while(i<num && !flag)
    {
      // for each vertex in group
      iter_list_v q=c[i].verts.begin();
      while( q!=c[i].verts.end() && !flag)
      {
        // if vertex is part of face
        if (p->vertexInFace(&(*q)))
        {
          flag=true;
          // then add face to group
          c[i].faces.push_back(*p);
        }
        q++;
      }
      i++;
    }
  }
}

void Group::printContourFacesVertices (char *str)
{
  printf("%s\n",str);
  for (int i=0;i<count;i++)
  {
    printf("\nContour %i\n",i);
    for (iter_list_v q=c[i].verts.begin();q!=c[i].verts.end();q++)
    {
      printf("vertex %i %.15g %.15g %.15g\n",q->index,q->x,q->y,q->z);
    }
    for (iter_list_f q=c[i].faces.begin();q!=c[i].faces.end();q++)
    {
      printf("face %i %i %i %i\n",q->index,q->v1,q->v2,q->v3);
    }
  }
}

void Group::printOtherContourData (char *str)
{
  printf("%s\n",str);
  for (int i=0;i<count;i++)
  {
    printf("\nContour group %i\n",i);
    printf("\nnumber of sites %i\n",c[i].getNum());
    for (iter_vec_ev q=c[i].extra.begin();q!=c[i].extra.end();q++)
    {
      ExtraVertex *ev=&(*q);
      printf("extra vertex group %i\n",ev->g);
      printf("extra vertex %i %.15g %.15g %.15g\n",
             ((Vertex*)ev->v)->index,
             ((Vertex*)ev->v)->x,
             ((Vertex*)ev->v)->y,
             ((Vertex*)ev->v)->z);
      if (q->la!=NULL)
      {
        printf("extra vertex la %i %.15g %.15g %.15g\n",
               ((Vertex*)ev->la)->index,
               ((Vertex*)ev->la)->x,
               ((Vertex*)ev->la)->y,
               ((Vertex*)ev->la)->z);
      }
      if (q->lb!=NULL)
      {
        printf("extra vertex lb %i %.15g %.15g %.15g\n",
               ((Vertex*)ev->lb)->index,
               ((Vertex*)ev->lb)->x,
               ((Vertex*)ev->lb)->y,
               ((Vertex*)ev->lb)->z);
      }
    }
    for (int j=0;j<c[i].getNum();j++)
    {
      printf("Site %i lateral1 %i %.15g %.15g %.15g\n",j,
             ((Vertex*)c[i].s[j].l1)->index,
             ((Vertex*)c[i].s[j].l1)->x,
             ((Vertex*)c[i].s[j].l1)->y,
             ((Vertex*)c[i].s[j].l1)->z);
      printf("Site %i lateral2 %i %.15g %.15g %.15g\n",j,
             ((Vertex*)c[i].s[j].l2)->index,
             ((Vertex*)c[i].s[j].l2)->x,
             ((Vertex*)c[i].s[j].l2)->y,
             ((Vertex*)c[i].s[j].l2)->z);
      printf("Site %i third %i %.15g %.15g %.15g\n",j,
             ((Vertex*)c[i].s[j].th)->index,
             ((Vertex*)c[i].s[j].th)->x,
             ((Vertex*)c[i].s[j].th)->y,
             ((Vertex*)c[i].s[j].th)->z);
      printf("Site %i orientation %i\n",j,c[i].s[j].orient);
      for (int k=0;k<c[i].s[j].n;k++)
      {
        ExtraVertex *ev=c[i].s[j].ev[k];
        printf("Site %i ev %i vertex %i %.15g %.15g %.15g\n",j,k,
               ((Vertex*)ev->v)->index,
               ((Vertex*)ev->v)->x,
               ((Vertex*)ev->v)->y,
               ((Vertex*)ev->v)->z);
      }
    }
  }
}

