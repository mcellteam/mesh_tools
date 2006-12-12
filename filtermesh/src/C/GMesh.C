// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include <cctype>
#include "GMesh.h"
#include "Polygon.h"
#include "A3dStream.h"
#include "Array.h"

using std::istream;
using std::ostream;
using std::flush;
using std::ws;

const int BUFLEN=2000;

GMesh::GMesh() : ios(0) { }

GMesh::~GMesh() { }

void GMesh::copy(const GMesh& m)
{
	Mesh::copy(m);
	gmodFlag(m.gflag(ALL),1);
	ForMeshVertex(m,v) {
		setPoint(idvertex(m.vertexid(v)),m.point(v));
	} EndFor;
	ForMeshVertex(m,v) {
		modFlag(idvertex(m.vertexid(v)),m.flag(v,ALL),1);
	} EndFor;
	ForMeshFace(m,f) {
		modFlag(idface(m.faceid(f)),m.flag(f,ALL),1);
	} EndFor;
	ForMeshEdge(m,e) {
		modFlag(edge(idvertex(m.vertexid(m.vertex1(e))),
			     idvertex(m.vertexid(m.vertex2(e)))),
			m.flag(e,ALL),1);
	} EndFor;
}

void GMesh::copyVertices(const GMesh& m)
{
	Mesh::copyVertices(m);
	ForMeshVertex(m,v) {
		setPoint(idvertex(m.vertexid(v)),m.point(v));
	} EndFor;
	ForMeshVertex(m,v) {
		modFlag(idvertex(m.vertexid(v)),m.flag(v,ALL),1);
	} EndFor;
}

void GMesh::setPoint(Vertex v, const Point& pp)
{
	v->point=pp;
	if (ios) *ios << hform("MVertex %d  %15lg %15lg %15lg\n",
			       vertexid(v),pp[0],pp[1],pp[2]);
}

void GMesh::polygon(Face f, Polygon& poly) const
{
	poly.init(0);
	ForFaceVertex(*this,f,v) {
		poly+=point(v);
	} EndFor;
}

void GMesh::points(Face f, Point pa[3]) const
{
	Vertex va[3];
	vertices(f,va);
	for (int i=0;i<3;i++) pa[i]=point(va[i]);
}

double GMesh::length2(Edge e) const
{
	return dist2(point(vertex1(e)),point(vertex2(e)));
}

double GMesh::length(Edge e) const
{
	return mysqrt(length2(e));
}

// Info

MeshInfoString::MeshInfoString(const char* ps) : s(0) { setString(ps); }

MeshInfoString::~MeshInfoString() { delete[] s; }

const char* MeshInfoString::getString() const { return s; }

void MeshInfoString::setString(const char* ps)
{
	delete[] s,s=newString(ps);
}

// I/O

void GMesh::read(istream& is)
{
	char str[BUFLEN];
	for (;is;) {
		is >> ws;
		is.get(str,sizeof(str));
		if (!is) break;
		assertw1(is.peek()=='\n');
		is.ignore(9999,'\n');
		readline(str);
	}
	if (sdebug>=1) OK();
}

void GMesh::readline(char* str)
{
	if (str[0]=='#') return;
	char* sinfo=strchr(str,'{');
	if (sinfo) {
		*sinfo++=0;
		char* s=strchr(sinfo,'}');
		if (assertw1(s)) sinfo=0; else *s=0;
	}
	if (!strncmp(str,"Vertex ",7)) {
		Point p; int vi;
		int nf=sscanf(str,"Vertex %d %lg %lg %lg",&vi,&p[0],&p[1],&p[2]);
		assertx(nf==4);
		Vertex v=createVertexI(vi); setPoint(v,p);
		if (sinfo) {
			setInfo(v,new MeshInfoString(sinfo));
			if (strstr(sinfo,"cusp")) modFlag(v,CUSPV,1);
		}
	} else if (!strncmp(str,"Face ",5)) {
		static Array<Vertex> va; va.init(0);
		char* s=str+4;
		int fi=-1;
		for (;;) {
			while (*s && isspace(*s)) s++;
			if (!*s) break;
			char* beg=s;
			while (*s && isdigit(*s)) s++;
			if (assertw1(!*s || isspace(*s))) continue;
			int j=atoi(beg);
			if (fi<0) { fi=j; continue; }
			va+=idvertex(j);
		}
		if (assertw1(va.num()>=3)) return;
		Face f=fi?createFaceI(fi,va,va.num()):createFace(va,va.num());
		if (assertw1(f)) return;
		if (sinfo) setInfo(f,new MeshInfoString(sinfo));
	} else if (!strncmp(str,"MVertex ",8)) {
		Point p; int vi;
		int nf=sscanf(str,"MVertex %d %lg %lg %lg",&vi,&p[0],&p[1],&p[2]);
		assertx(nf==4);
		setPoint(idvertex(vi),p);
		if (sinfo) setInfo(idvertex(vi),new MeshInfoString(sinfo));
	} else if (!strncmp(str,"Edge ",5)) {
		int vi1, vi2;
		assertx(sscanf(str,"Edge %d %d",&vi1,&vi2)==2);
		Edge e=orderedEdge(idvertex(vi1),idvertex(vi2));
		if (sinfo) {
			setInfo(e,new MeshInfoString(sinfo));
			modFlag(e,SHARPE,strstr(sinfo,"sharp")?1:0);
		}
	} else if (!strncmp(str,"CVertex ",8)) {
		(void)createVertexI(atoi(str+8));
	} else if (!strncmp(str,"DVertex ",8)) {
		destroyVertex(idvertex(atoi(str+8)));
	} else if (!strncmp(str,"DFace ",6)) {
		destroyFace(idface(atoi(str+6)));
	} else if (!strncmp(str,"SubFV ",6)) {
		int fi, viold, vinew;
		assertx(sscanf(str,"SubFV %d %d %d",&fi,&viold,&vinew)==3);
		substituteFaceVertex(idface(fi),
				     idvertex(viold),idvertex(vinew));
	} else if (!strncmp(str,"Ecol ",5)) {
		int vi1, vi2;
		assertx(sscanf(str,"Ecol %d %d",&vi1,&vi2)==2);
		collapseEdge(orderedEdge(idvertex(vi1),idvertex(vi2)));
	} else if (!strncmp(str,"Eswa ",5)) {
		int vi1, vi2;
		assertx(sscanf(str,"Eswa %d %d",&vi1,&vi2)==2);
		assertx(swapEdge(orderedEdge(idvertex(vi1),idvertex(vi2))));
	} else if (!strncmp(str,"Espl ",5)) {
		int vi1, vi2, vi3;
		assertx(sscanf(str,"Espl %d %d %d",&vi1,&vi2,&vi3)==3);
		splitEdge(orderedEdge(idvertex(vi1),idvertex(vi2)),vi3);
	} else {
		if (Warning("GMesh::read: cannot parse line")) SHOWN(str);
	}
}

int GMesh::recognizeLine(const char* s)
{
	static const char* sar[]={"Vertex ","Face ","MVertex ","Edge ",
				  "CVertex ","DVertex ","DFace ","SubFV ",
				  "Ecol ","Eswa ","Espl ",0};
	static int lenar[20];
	if (!lenar[0]) {
		for (int i=0;sar[i];i++) lenar[i]=strlen(sar[i]);
	}
	for (int i=0;sar[i];i++)
		if (!strncmp(s,sar[i],lenar[i])) return 1;
	return 0;
}

void GMesh::write(ostream& os) const
{
	ForMeshOrderedVertex(*this,v) {
		const Point& p=point(v);
		os << hform("Vertex %d  %.15g %.15g %.15g", vertexid(v), p[0], p[1], p[2]);
		MeshInfo* in=info(v);
		const char* sinfo=in?in->getString():0;
		if (sinfo) os << " {" << sinfo << "}";
		os << "\n";
		assertx(os);
	} EndFor;
	ForMeshOrderedFace(*this,f) {
		os << "Face " << faceid(f) << " "; 
		ForFaceVertex(*this,f,v) {
			os << " " << vertexid(v);
		} EndFor;
		MeshInfo* in=info(f);
		const char* sinfo=in?in->getString():0;
		if (sinfo) os << " {" << sinfo << "}";
		os << "\n";
		assertx(os);
	} EndFor;
	ForMeshEdge(*this,e) {
		MeshInfo* in=info(e);
		const char* sinfo=in?in->getString():0;
		if (!sinfo) continue;
		os << "Edge " << vertexid(vertex1(e)) << " " <<
			vertexid(vertex2(e)) << " {" << sinfo << "}\n";
		assertx(os);
	} EndFor;
	os << flush;
}

void GMesh::write(WA3dStream& oa3d, const A3dVertexColor& col) const
{
	static A3dElem el;
	ForMeshOrderedFace(*this,f) {
		el.init(A3dElem::TPolygon);
		ForFaceVertex(*this,f,v) {
			el+=A3dVertex(point(v),Vector(0,0,0),col);
		} EndFor;
		oa3d.write(el);
	} EndFor;
	oa3d.flush();
}

ostream* GMesh::recordChanges(ostream* pos)
{
	ostream* oos=ios; ios=pos; return oos;
}

Vertex GMesh::createVertexI(int id)
{
	if (ios) *ios << "CVertex " << id << '\n';
	return Mesh::createVertexI(id);
}

// Override Mesh members
void GMesh::destroyVertex(Vertex v)
{
	if (ios) *ios << "DVertex " << vertexid(v) << '\n';
	Mesh::destroyVertex(v);
}

Face GMesh::createFaceI(int id, Vertex va[], int nv)
{
	Face f=Mesh::createFaceI(id,va,nv);
	if (!f) return 0;
	if (ios) {
		*ios << "Face " << id << ' ';
		for (int i=0;i<nv;i++) *ios << ' ' << vertexid(va[i]);
		*ios << '\n';
	}
	return f;
}

void GMesh::destroyFace(Face f)
{
	if (ios) *ios << "DFace " << faceid(f) << '\n';
	Mesh::destroyFace(f);
}

int GMesh::substituteFaceVertex(Face f, Vertex vold, Vertex vnew)
{
	if (!Mesh::substituteFaceVertex(f,vold,vnew)) return 0;
	if (ios) *ios << "SubFV " << faceid(f) << ' ' <<
		vertexid(vold) << ' ' << vertexid(vnew) << '\n';
	return 1;
}

void GMesh::collapseEdge(Edge e)
{
	if (sdebug>=1) valid(e);
	ostream* tos=ios; ios=0;
	Vertex v1=vertex1(e), v2=vertex2(e);
	int id1=vertexid(v1), id2=vertexid(v2);
	// Compute geometry for unified vertex (v1)
	int isb1=isBoundary(v1), isb2=isBoundary(v2), sumb=isb1+isb2;
	// GNUG 2.5.5 cant parse this
// 	Point p=(sumb==0 || sumb==2 ? interp(point(v1),point(v2)) :
// 		 isb1 ? point(v1) : point(v2));
	Point p;
	if (sumb==0 || sumb==2) p=interp(point(v1),point(v2));
	else p= isb1 ? point(v1) : point(v2);
	Set<Vertex> vsharp;
	ForEdgeVertex(*this,e,v) {
		ForVertexEdge(*this,v,ee) {
			if (ee!=e && flag(ee,SHARPE))
				vsharp.add(oppVertex(v,ee));
		} EndFor;
	} EndFor;
	Mesh::collapseEdge(e);	// v1 kept
	setPoint(v1,p);		// (ios==0)
	ForSet(vsharp,Vertex,v) {
		Edge e=edge(v1,v);
		if (!isBoundary(e)) modFlag(e,SHARPE,1);
	} EndFor;
	if (tos) ios=tos, *ios << "Ecol " << id1 << ' ' << id2 << '\n';
}

Vertex GMesh::splitEdge(Edge e, int id)
{
	if (sdebug>=1) valid(e);
	ostream* tos=ios; ios=0;
	Vertex v1=vertex1(e), v2=vertex2(e);
	Vertex vo1=sideVertex1(e), vo2=sideVertex2(e);
	int fl11, fl12, fl21=0, fl22=0, fle=flag(e,ALL);
	fl11=flag(edge(vo1,v1),ALL),fl12=flag(edge(vo1,v2),ALL);
	if (vo2) fl21=flag(edge(vo2,v1),ALL),fl22=flag(edge(vo2,v2),ALL);
	Vertex v=Mesh::splitEdge(e,id);
	modFlag(edge(vo1,v1),fl11,1),modFlag(edge(vo1,v2),fl12,1);
	if (vo2) modFlag(edge(vo2,v1),fl21,1),modFlag(edge(vo2,v2),fl22,1);
	modFlag(edge(v1,v),fle,1),modFlag(edge(v2,v),fle,1);
	setPoint(v,interp(point(v1),point(v2))); // (ios==0)
	if (tos) ios=tos, *ios << "Espl " << vertexid(v1) << ' ' <<
		vertexid(v2) << ' ' << vertexid(v) << '\n';
	return v;
}

Edge GMesh::swapEdge(Edge e)
{
	if (sdebug>=1) valid(e);
	ostream* tos=ios; ios=0;
	Vertex v1=vertex1(e), v2=vertex2(e);
	Vertex vo1=sideVertex1(e), vo2=sideVertex2(e); assertx(vo1 && vo2);
	int fl11, fl12, fl21, fl22;
	fl11=flag(edge(vo1,v1),ALL),fl12=flag(edge(vo1,v2),ALL);
	fl21=flag(edge(vo2,v1),ALL),fl22=flag(edge(vo2,v2),ALL);
	Edge ne=Mesh::swapEdge(e);
	modFlag(edge(vo1,v1),fl11,1),modFlag(edge(vo1,v2),fl12,1);
	modFlag(edge(vo2,v1),fl21,1),modFlag(edge(vo2,v2),fl22,1);
	if (tos) ios=tos, *ios << "Eswa " << vertexid(v1) << ' ' <<
		vertexid(v2) << '\n';
	return ne;
}
