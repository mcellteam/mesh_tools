// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef SubMesh_h
#define SubMesh_h

#include "GMesh.h"
#include "Combination.h"
#include "Homogeneous.h"
#include "Array.h"

struct Combvh {
	Combination<Vertex> c;
	Homogeneous h;
	int iscombination() const;
	POOLALLOCATION(Combvh);
	Combvh() {} ~Combvh() {} // GNUG 2.5.8
};

class Mvcvh : public Map<Vertex,Combvh*> {
  public:
	~Mvcvh();		// automatically deletes combinations
	void shallowclear();	// no delete of combinations
	void clear();
	int isconvolution() const; // check combination is affine
	void compose(const Combvh& ci, Combvh& co) const; // co=ci*this
	// compose two maps to produce one, die unless mconv.isconvolution()
	void compose(const Mvcvh& mconv); // this=mconv*this
};

class SubMesh {
  public:
	SubMesh(GMesh& mesh);
	~SubMesh();
	void clear();
	enum { VARIABLEV=1<<28 }; // flag bit for GMesh vertices
	// mesh() may be modified if no more SubMesh operations will be done.
	GMesh& mesh() { return m; }
	const GMesh& mesh() const { return m; }
	GMesh& origmesh() { return omesh; }
	const GMesh& origmesh() const { return omesh; }
// subdivide (makes use of refine(),createconv(),convolveself(),...)
	void subdivide(double cosang=1.);
	void subdividen(int nsubdiv, int limit, double cosang=1.);
// Combinations
	// get a combination (expressing v of mesh() in terms of origmesh())
	const Combvh& combination(Vertex v) const;
	// Compose c1 with cmvcvh to get combination in terms of orig. verts.
	void composecmvcvh(const Combvh& ci, Combvh& co) const;
// update vertex positions on mesh() according to its mask
	void updatevertexposition(Vertex v);
	void updatevertexpositions();
// misc
	void maskparameters(int ps222, double pweighta)
	{ s222=ps222; weighta=pweighta; }
	void init() { s222=0; weighta=0; }
// omesh to and from mesh	
	Face origface(Face f) const;
	void origfaceindex(Face fi, Face& fo, int& pindex) const;
	Face getface(Face of, int index) const;
// split and compute splitting masks
	void refine(Mvcvh& mconv); // 4to1 split at edge medians
	// refine near creases, and refine edges with cosdihedral <cosang
	void selrefine(Mvcvh& mconv, double cosang);
// compute averaging masks
	typedef void (SubMesh::*FVMASK)(Vertex v, Combvh& comb) const;
	void createconv(Mvcvh& mconv, FVMASK f); // use a subdivision mask
// the masks
	void averagingmask(Vertex v, Combvh& comb) const;
	void limitmask(Vertex v, Combvh& comb) const;
// apply a convolution
	void convolveself(const Mvcvh& mconv); // cmvcvh=mconv*cmvcvh
// debug
	void showmvcvh(const Mvcvh& mvcvh) const;
	void showcmvcvh() const;
  private:
	GMesh& omesh;
	GMesh m;
	Mvcvh cmvcvh;		     // maps vertices of m to Combvh of omesh
	Map<Face,Face> mforigf;	     // face of m -> face of omesh
	Map<Face,int> mfindex;	     // face of m -> index within origf
	Map<Face,Array<Face>*> mofif; // face of omesh -> (int -> face of m)
	int allvvar;		     // 1 if all vertices are VARIABLEV
	//
	int s222;
	double weighta;
	//
	int sharp(Edge e) const;
	int nume(Vertex v) const;
	int numsharpe(Vertex v) const;
	Edge oppsharpe(Vertex v, Edge e) const;
	Vertex oppsharpv(Vertex v, Vertex v2) const;
	void subdivideaux(double cosang, Mvcvh* mconv);
	void creaseaveragingmask(Vertex v, Combvh& comb) const;
	int extraordinarycreasev(Vertex v) const;
	DISABLECOPY(SubMesh);
};

#endif
