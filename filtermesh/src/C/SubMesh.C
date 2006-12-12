// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe

#include "Hh.h"
#include "SubMesh.h"
#include "Homogeneous.h"
#include "MeshOp.h"		// for EdgeDihedralAngleCos()
#include "HashStruct.h"

// On 04/18/94, eliminated capacity to get higher order splines on creases.
// This did not seem to make sense in conjunction with extraordinary crease
// vertices.
// Also, this greatly reducing bookkeeping (msharpvv disappears)

// mofif and mfindex work by assuming that two faces--whose vertices are in
// the same order and have the same ids--will be subdivided in the same order
// regardless of the mesh they belong to.
// This is true because of Mesh::vertices and because of ForMeshOrderedFace.

struct vertexinfo : MeshInfo {
	vertexinfo(int pnume, int pnumsharpe)
		: nume(pnume), numsharpe(pnumsharpe) { }
	int nume;
	int numsharpe;
};

static int debug=GetenvValue("DEBUG");

//*** helper

// translate Vertex from one Mesh to another Mesh
static Vertex trvmm(Vertex v, const Mesh& mf, const Mesh& mt)
{
	return mt.idvertex(mf.vertexid(v));
}

// translate Face from one Mesh to another Mesh
static Face trfmm(Face f, const Mesh& mf, const Mesh& mt)
{
	return mt.idface(mf.faceid(f));
}

#if 0
// translate Edge from one Mesh to another Mesh
static Edge tremm(Edge e, const Mesh& mf, const Mesh& mt)
{
	return mt.edge(trvmm(mf.vertex1(e),mf,mt),
		       trvmm(mf.vertex2(e),mf,mt));
}
#endif

//*** Combvh

ALLOCATEPOOL(Combvh);

int Combvh::iscombination() const
{
	return fabs(c.sum()+h[3]-1.)<1e-6;
}

//*** Mvcvh

Mvcvh::~Mvcvh()
{
	clear();
}

void Mvcvh::shallowclear()
{
	Map<Vertex,Combvh*>::clear();
}

void Mvcvh::clear()
{
	ForMapValue(*this,Vertex,Combvh*,comb) { delete comb; } EndFor;
	shallowclear();
}

int Mvcvh::isconvolution() const
{
	ForMapValue(*this,Vertex,Combvh*,comb) {
		if (!comb->h.iszero()) return 0;
	} EndFor;
	return 1;
}

// co=ci*this
void Mvcvh::compose(const Combvh& ci, Combvh& co) const
{
	// assertx(co.c.empty() && co.h.iszero());
	co.h=ci.h;
	ForCombination(ci.c,Vertex,v,val) {
		const Combvh* comb=retrieve(v);
		if (!comb) {
			// missing entry -> assume identity map
			co.c[v]+=val;
		} else {
			co.h+=comb->h*val;
			if (comb->c.num()) {
				ForCombination(comb->c,Vertex,v2,val2) {
					co.c[v2]+=val2*val;
				} EndFor;
			}
		}
	} EndFor;
}

// this=mconv*this
void Mvcvh::compose(const Mvcvh& mconv)
{
	Mvcvh nthis;
	ForMapKeyValue(mconv,Vertex,v,Combvh*,pcomb) {
		const Combvh& comb=*pcomb;
		assertx(comb.h.iszero());
		if (!comb.c.num()) {
			// identity assumed, keep unchanged
		} else if (comb.c.num()==1) {
			Warning("Waste of space");
			assertw1(comb.c[v]==1);
		} else {
			Combvh& ncomb=*new Combvh;
			compose(comb,ncomb);
			nthis.enter(v,&ncomb);
			if (debug) assertx(ncomb.iscombination());
		}
	} EndFor;
	ForMapKeyValue(nthis,Vertex,v,Combvh*,comb) {
		Combvh* c=replace(v,comb);
		if (c) delete c; else enter(v,comb);
	} EndFor;
	nthis.shallowclear();	// prevent delete of combinations
}

//*** SubMesh

// weight of center vertex in position vertex mask, Loop scheme
// It is the weight that appears in the subdivision matrix,
// not the weight used when splitting+averaging.
inline double subdiv_a(int n)
{
	return n==6?.625 : 3./8.+square((3+2*cos(2*PI/n))/8);
}
// 3,0.4375  4,0.515625  5,0.579534  6,0.625  7,0.656826  1000,0.765619

SubMesh::SubMesh(GMesh& mesh)
: omesh(mesh), allvvar(1)
{
	init();
	m.copy(omesh);		// vertex ids will match, flags copied
	ForMeshFace(m,f) {
		mforigf.enter(f,trfmm(f,m,omesh));
		mfindex.enter(f,0);
	} EndFor;
	ForMeshFace(omesh,f) {
		Array<Face>& array=*new Array<Face>;
		mofif.enter(f,&array);
		array+=trfmm(f,omesh,m);
	} EndFor;
	ForMeshVertex(m,v) {
		int nsharpe=0;
		ForVertexEdge(m,v,e) { if (sharp(e)) nsharpe++; } EndFor;
		m.setInfo(v,new vertexinfo(m.degree(v),nsharpe));
	} EndFor;
	ForMeshVertex(m,v) {
		Combvh& comb=*new Combvh;
		cmvcvh.enter(v,&comb);
		if (m.flag(v,VARIABLEV)) {
			comb.c[trvmm(v,m,omesh)]=1;
		} else {
			comb.h=omesh.point(trvmm(v,m,omesh));
			allvvar=0;
		}
	} EndFor;
}

SubMesh::~SubMesh()
{
	clear();
}

void SubMesh::clear()
{
	m.clear();
	cmvcvh.clear();
	mforigf.clear();
	mfindex.clear();
	ForMapValue(mofif,Face,Array<Face>*,a) { delete a; } EndFor;
	mofif.clear();
}

//*** subdivide

void SubMesh::subdivide(double cosang)
{
	if (allvvar) {
		Mvcvh mconv;
		subdivideaux(cosang,&mconv);
		convolveself(mconv);
	} else {
		subdivideaux(cosang,0);
	}
}

void SubMesh::subdividen(int nsubdiv, int limit, double cosang)
{
	// whichever is faster
	if (allvvar) {
		Mvcvh mconv;
		ForIndex(i,nsubdiv) {
			Mvcvh mconv1;
			subdivideaux(cosang,&mconv1);
			mconv.compose(mconv1);
		} EndFor;
		if (limit) {
			Mvcvh mconv1;
			createconv(mconv1,&SubMesh::limitmask);
			mconv.compose(mconv1);
		}
		convolveself(mconv);
	} else {
		ForIndex(i,nsubdiv) {
			subdivideaux(cosang,0);
		} EndFor;
		if (limit) {
			Mvcvh mconv;
			createconv(mconv,&SubMesh::limitmask);
			convolveself(mconv);
		}
	}
}

void SubMesh::subdivideaux(double cosang, Mvcvh* pmconv)
{
	Mvcvh* mconv=pmconv?pmconv:new Mvcvh;
	if (cosang<.99999) {
		// udpate geometry so that EdgeDihedralCos makes sense
		updatevertexpositions();
		selrefine(*mconv,cosang);
	} else {
		refine(*mconv);
	}
	if (!pmconv) convolveself(*mconv),delete mconv,mconv=0;
	Mvcvh mconv2;
	createconv(mconv2,&SubMesh::averagingmask);
	if (!pmconv) convolveself(mconv2);
	else pmconv->compose(mconv2);
}

//*** compute convolutions

void SubMesh::refine(Mvcvh& mconv)
{
        int i;
	// Save current mesh objects for later iteration
	Stack<Vertex> stackv; Stack<Face> stackf; Stack<Edge> stacke;
	ForMeshVertex(m,v) { stackv.push(v); } EndFor;
	ForMeshOrderedFace(m,f) { stackf.push(f); } EndFor;
	ForMeshEdge(m,e) { stacke.push(e); } EndFor;
	Map<Edge,Vertex> menewv;
	// Create new vertices and make them midpoints of old edges
	ForStack(stacke,Edge,e) {
		Vertex v=m.createVertex();
		m.setInfo(v,new vertexinfo(m.isBoundary(e)?4:6,sharp(e)?2:0));
		menewv.enter(e,v);
		Combvh& comb=*new Combvh;
		mconv.enter(v,&comb);
		comb.c[m.vertex1(e)]=.5;
		comb.c[m.vertex2(e)]=.5;
	} EndFor;
	// Create new triangulation in situ with old one.
	int oldv=m.setDebug(-1);
	ForStack(stackf,Face,f) {
		Vertex va[3], vs[3];
		m.vertices(f,va);
//		for (int i=0;i<3;i++)
		for (i=0;i<3;i++)
			vs[i]=menewv.get(m.edge(va[i],va[(i+1)%3]));
		Face forig=mforigf.get(f);
		Array<Face>& array=*mofif.get(forig);
		for (i=0;i<3;i++) {
			Face fn=assertv(m.createFace(va[i],vs[i],vs[(i+2)%3]));
			mforigf.enter(fn,forig);
			mfindex.enter(fn,array.num());
			array+=fn;
		}
		Face fn=assertv(m.createFace(vs,3));
		mforigf.enter(fn,forig);
		mfindex.enter(fn,array.num());
		array+=fn;
	} EndFor;
	m.setDebug(oldv);
	// Update sharp edges
	ForStack(stacke,Edge,e) {
		if (!m.flag(e,GMesh::SHARPE)) continue;
		Vertex vnew=menewv.get(e);
		m.modFlag(m.edge(vnew,m.vertex1(e)),GMesh::SHARPE,1);
		m.modFlag(m.edge(vnew,m.vertex2(e)),GMesh::SHARPE,1);
	} EndFor;
	// Remove old triangulation
	ForStack(stackf,Face,f) {
		m.destroyFace(f);
		assertx(mforigf.remove(f));
		mfindex.remove(f); // can be index 0
	} EndFor;
}

//*** selrefine

static GMesh* gmesh;		// unfortunate

struct vvnewv {
	vvnewv(Vertex pv1, Vertex pv2, Vertex pvnew, int psharp)
		: vnew(pvnew), sharp(psharp)
	{ if (pv1<pv2) v1=pv1,v2=pv2; else v1=pv2,v2=pv1; }
	Vertex v1, v2, vnew;
	int sharp;
};

static Univ hashvvnewv(const vvnewv* s)
{
	return Conv<int>::e(gmesh->vertexid(s->v1)+gmesh->vertexid(s->v2)*761);
}

static int cmpvvnewv(const vvnewv* s1, const vvnewv* s2)
{
	return s1->v1!=s2->v1 || s1->v2!=s2->v2;
}

void SubMesh::selrefine(Mvcvh& mconv, double cosang)
{
        int i;
	// mforigf, mfindex, and mofif are not supported with this scheme!
	gmesh=&m;
	Stack<Face> stackf;
	ForMeshFace(m,f) { stackf.push(f); } EndFor;
	// Determine which edges will be subdivided
	Set<Edge> subde;	// edges to subdivide
	ForMeshEdge(m,e) {
		if (sharp(e) || EdgeDihedralAngleCos(m,e)<=cosang)
			subde.enter(e);
	} EndFor;
	Queue<Face> queuef;
	ForStack(stackf,Face,f) { queuef.enqueue(f); } EndFor;
	while (!queuef.empty()) {
		Face f=queuef.dequeue();
		int nnew=0;
		ForFaceEdge(m,f,e) {
			if (subde.contains(e)) nnew++;
		} EndFor;
		if (nnew!=2) continue; // ok, no propagating changes
		ForFaceEdge(m,f,e) {
			if (!subde.add(e)) continue;
			Face f2=m.oppFace(f,e);
			if (f2) queuef.enqueue(f2);
		} EndFor;
	}
	// Introduce new vertices at midpoints
	HashStruct<vvnewv> mvvnewv(hashvvnewv,cmpvvnewv);
	ForSet(subde,Edge,e) {
		Vertex v=m.createVertex();
		m.setInfo(v,new vertexinfo(m.isBoundary(e)?4:6,sharp(e)?2:0));
		mvvnewv.enter(new vvnewv(m.vertex1(e),m.vertex2(e),v,
					 m.flag(e,GMesh::SHARPE)));
		Combvh& comb=*new Combvh;
		mconv.enter(v,&comb);
		comb.c[m.vertex1(e)]=.5;
		comb.c[m.vertex2(e)]=.5;
	} EndFor;
	// Subdivide faces, destroys validity of flags(e)
	ForStack(stackf,Face,f) {
		Vertex va[3], vs[3]; int nnew=0, i0=-1;
		m.vertices(f,va);
//		for (int i=0;i<3;i++) {
		for (i=0;i<3;i++) {
			Edge e=m.edge(va[i],va[(i+1)%3]);
			vvnewv si(m.vertex1(e),m.vertex2(e),0,0);
			vvnewv* s=mvvnewv.retrieve(&si);
			vs[i]=s?s->vnew:0;
			if (vs[i]) nnew++,i0=i;
		}
		if (!nnew) continue;
		assertx(nnew==1 || nnew==3);
		m.destroyFace(f);
		if (nnew==1) {
			int i1=(i0+1)%3, i2=(i0+2)%3; Face fn;
			fn=assertv(m.createFace(va[i0],vs[i0],va[i2]));
			fn=assertv(m.createFace(vs[i0],va[i1],va[i2]));
		} else {
			for (i=0;i<3;i++)
				assertx(m.createFace(va[i],vs[i],vs[(i+2)%3]));
		}
	} EndFor;
	ForHashStruct(mvvnewv,vvnewv,s) {
		if (s->sharp) {
			m.modFlag(m.edge(s->vnew,s->v1),GMesh::SHARPE,1);
			m.modFlag(m.edge(s->vnew,s->v2),GMesh::SHARPE,1);
		}
		delete s;
	} EndFor;
	mvvnewv.clear();	// avoid warning message
}

void SubMesh::createconv(Mvcvh& mconv, FVMASK fsubdivision)
{
	assertx(mconv.empty());
	Combvh* comb=0;
	ForMeshVertex(m,v) {
		if (!comb) comb=new Combvh;
		(this->*fsubdivision)(v,*comb);
		assertx(comb->h.iszero());
		if (comb->c.empty()) continue;
		if (debug) assertx(comb->iscombination());
		mconv.enter(v,comb);
		comb=0;
	} EndFor;
	delete comb;
}

//*** averaging masks and limit masks

// Central weight for the subdivision masks at a cone, as function of n
const int TABCONEN=10;		// valid up to and including this
static const double tableconesubd[TABCONEN+1]=
{ 0, 0, 0, .625, .75, .827254, .875, .905872, .926777, .941511, .952254 };
static const double tableconelimit[TABCONEN+1]=
{ 0, 0, 0, .5, .6, .684624, .75, .799356, .836637, .865074, .887058 };

void SubMesh::averagingmask(Vertex v, Combvh& comb) const
{
	int ne=nume(v), nesharp=numsharpe(v);
	double wa=-1;		// central weight
	if (m.flag(v,GMesh::CUSPV)) { // cusp vertex
		if (nesharp>=2) return;	// becomes corner vertex
		int nuse=ne;
		if (assertw1(nuse<=TABCONEN)) nuse=TABCONEN;
		wa=tableconesubd[nuse]*2-1;
		if (weighta) wa=weighta*2-1;
	} else if (nesharp<=1) { // interior or dart vertex
		// normal case, quartic bspline surface
		// was wa=subdiv_a(ne); after refinement, wa=subdiv_a(ne)*2-1
		// s222 : mask was n/3 --- 1 (n times); after refinement, .25
		wa=s222?.25 : subdiv_a(ne)*2-1;
		if (weighta && ne!=6) wa=weighta*2-1;
	} else if (nesharp==2) { // on crease
		creaseaveragingmask(v,comb);
		return;
	} else if (nesharp>=3) { // corner vertex held constant
		return;
	} else { assertnever(""); }
	double wc=(1-wa)/ne;
	if (wa==1) return;
	if (wa) comb.c[v]=wa;
	if (wc) ForVertexVertex(m,v,vv) { comb.c[vv]=wc; } EndFor;
}

// Subdivision algorithm:
// - if dart-ord, dart-dart, dart-ecv, or dart-corner -> smooth mask
// - if ecv-ord or corner-ord -> special 3-5 mask
// - all other cases
//     (ord-ord, ecv-ecv, corner-corner, ecv-corner) -> crease mask
void SubMesh::creaseaveragingmask(Vertex v, Combvh& comb) const
{
	int adjdartvertices=0, adjcornervertices=0, adjecvertices=0, svi=0;
	static Array<Vertex> va; va.init(0);
	ForVertexEdge(m,v,e) {
		if (!sharp(e)) continue;
		Vertex vo=m.oppVertex(v,e);
		va+=vo;
		int nesharp=numsharpe(vo);
		if (nesharp==1) {
			adjdartvertices++;
		} else if (nesharp>=3) {
			adjcornervertices++,svi=va.num()-1;
		} else if (extraordinarycreasev(vo)) {
			adjecvertices++,svi=va.num()-1;
		}
	} EndFor;
	int adjspecial=adjdartvertices+adjcornervertices+adjecvertices;
	assertx(adjspecial<=2);
	if (adjdartvertices>=1) {
		// adjacent to single dart vertex, do smooth mask
		int ne=nume(v);
		double a=s222 ? .25 : subdiv_a(ne)*2-1;
		double wa=a, wc=(1-wa)/ne;
		comb.c[v]=wa;
		ForVertexVertex(m,v,vv) { comb.c[vv]=wc; } EndFor;
		return;
	}
	if (adjspecial!=1) {
		// 1,2,1 mask gives cubic spline
		comb.c[v]=.5; comb.c[va[0]]=.25; comb.c[va[1]]=.25;
		return;
	}
	if (adjecvertices==1 || adjcornervertices==1) {
		// adjacent to single extraordinary crease vertex or corner
		if (debug) SHOWL;
		comb.c[v]=.75; comb.c[va[1-svi]]=.25;
		return;
	}
	assertnever("");
}

int SubMesh::extraordinarycreasev(Vertex v) const
{
	assertx(numsharpe(v)==2); // optional
	if (m.isBoundary(v)) return nume(v)!=4;
	if (nume(v)!=6) return 1;
	int nside=0, sharpef=0;
	ForVertexCcwVertex(m,v,vv) {
		if (sharpef==1) nside++;
		if (sharp(m.edge(v,vv))) sharpef++;
	} EndFor;
	assertx(sharpef==2);
	return nside!=3;
}

void SubMesh::limitmask(Vertex v, Combvh& comb) const
{
	int ne=nume(v), nesharp=numsharpe(v);
	if (m.flag(v,GMesh::CUSPV)) { // cusp
		if (nesharp>=2) return;	// becomes corner vertex, held constant
		int nuse=ne;
		if (assertw1(nuse<=TABCONEN)) nuse=TABCONEN;
		double wa=tableconelimit[nuse];
		assertw1(!weighta);
		double wc=(1-wa)/ne;
		comb.c[v]=wa;
		ForVertexVertex(m,v,vv) { comb.c[vv]=wc; } EndFor;
		return;
	}
	if (nesharp>=3) return;	    // corner vertex held constant
	if (nesharp<=1) {	// interior or dart vertex
		if (weighta && ne!=6) { Warning("to do??"); return; }
		// was wa=3/(3+8*a1); wc=8*a1/(3+8*a1)/ne; wb=wc;
		// since there is no refinement here, should be the same
		// s222 : mask is n --- 1 (n times)
		double a=s222?.5 : subdiv_a(ne);
		double wa=3/(11-8*a), wc=(1-wa)/ne;
		if (wa) comb.c[v]=wa;
		if (wc) ForVertexVertex(m,v,vv) { comb.c[vv]=wc; } EndFor;
	} else if (nesharp==2) { // bspline curve
		// for cubic, 1-4-1 mask
		double wa=4./6., wb=1./6., wc=0;
		if (extraordinarycreasev(v))
			wa=3./5.,wb=1./5.;
		comb.c[v]=wa;
		ForVertexEdge(m,v,e) {
			double w=sharp(e)?wb:wc;
			if (w) comb.c[m.oppVertex(v,e)]=w;
		} EndFor;
	} else { assertnever(""); }
}

//*** misc

void SubMesh::convolveself(const Mvcvh& mconv)
{
	cmvcvh.compose(mconv);
}

const Combvh& SubMesh::combination(Vertex v) const
{
	return *cmvcvh.get(v);
}

void SubMesh::composecmvcvh(const Combvh& ci, Combvh& co) const
{
	cmvcvh.compose(ci,co);
}

void SubMesh::updatevertexposition(Vertex v)
{
	Combvh& comb=*cmvcvh.get(v);
	Homogeneous h(comb.h);
	ForCombination(comb.c,Vertex,vv,val) {
		h+=omesh.point(vv)*val;
	} EndFor;
	m.setPoint(v,toPoint(h));
}

void SubMesh::updatevertexpositions()
{
	ForMeshVertex(m,v) { updatevertexposition(v); } EndFor;
}

Face SubMesh::origface(Face f) const
{
	return mforigf.get(f);
}

void SubMesh::origfaceindex(Face fi, Face& of, int& pindex) const
{
	of=mforigf.get(fi);
	pindex=mfindex.get(fi);
}

Face SubMesh::getface(Face of, int index) const
{
	const Array<Face>& array=*mofif.get(of);
	return array[index];
}

//*** debug

void SubMesh::showmvcvh(const Mvcvh& mvcvh) const
{
	SHOWF("Mvcvh = {\n");
	ForMapKeyValue(mvcvh,Vertex,v,Combvh*,comb) {
		SHOWF(" vertex %d {\n",m.vertexid(v));
		SHOWF("  h[3]=%g\n",comb->h[3]);
		ForCombination(comb->c,Vertex,v,val) {
			SHOWF("  v%-4d  %g\n",m.vertexid(v),val);
		} EndFor;
		SHOWF(" }\n");
	} EndFor;
	SHOWF("}\n");
}

void SubMesh::showcmvcvh() const
{
	showmvcvh(cmvcvh);
}

//*** helper

int SubMesh::sharp(Edge e) const
{
	return m.isBoundary(e) || m.flag(e,GMesh::SHARPE);
}

int SubMesh::nume(Vertex v) const
{
	return (reinterpret_cast<vertexinfo*>(m.info(v)))->nume;
}

int SubMesh::numsharpe(Vertex v) const
{
	return (reinterpret_cast<vertexinfo*>(m.info(v)))->numsharpe;
}

Edge SubMesh::oppsharpe(Vertex v, Edge ee) const
{
	Edge er=0;
	ForVertexEdge(m,v,e) {
		if (!sharp(e)) continue;
		if (e==ee) { ee=0; } else { assertx(!er); er=e; }
	} EndFor;
	assertx(!ee && er);
	return er;
}

Vertex SubMesh::oppsharpv(Vertex v, Vertex v2) const
{
	return m.oppVertex(v,oppsharpe(v,m.edge(v,v2)));
}
