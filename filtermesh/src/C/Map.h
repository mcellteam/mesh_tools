// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Hash_h
#define Hash_h

#if 0
{
	Map<Edge,Vertex> mev;
	ForMapKeyValue(mev,Edge,e,Vertex,v) { do(e,v); } EndFor;
	ForMapKey(mev,Edge,Vertex,e) { do(e); } EndFor;
	ForMapValue(mev,Edge,Vertex,v) { do(v); } EndFor;
}
#endif

#include "Pool.h"

// Base "map" class.  Maps pointer-sized values onto pointer-sized values.
class BMap {
  public:

	/* Create a new base map */
	BMap();

	/* Destroy a base map */
	~BMap();

	/* Clear this map */
	void clear();

	/* "Enter" a new value into the map */
	void enter(Univ k, Univ v); // k must be new

	/* Check if this map contains a given key */
	bool contains(Univ k) const { return find(k) != 0; }

	/*
	 * Retrieve a value from the map.  Success/failure is indicated
	 * by "present"
	 */
	Univ retrieve(Univ k, int& present) const;

	/*
	 * Retrieve a value from the map.  Success/failure is indicated
	 * by "present"
	 */
	Univ retrieve(Univ k, bool& present) const;

	/*
	 * Retrieve a value from the map.  Returns 0 if value is not found or
	 * if the value happens to be 0...
	 */
	Univ retrieve(Univ k) const;

	/* Get a value from the map.  Segfault if value not present. */
	Univ get(Univ k) const { return find(k)->v; }

	/*
	 * Remove a value from the map, returning the previous value or 0 if no
	 * such value existed.  (this is ambiguous)
	 */
	Univ remove(Univ k);

	/*
	 * Replace a value in the map, returning the previous value or 0 if no
	 * such value existed.  (this is ambiguous).  NOTE: If no such value
	 * existed before, no value is entered into the map!
	 */
	Univ replace(Univ k, Univ v);

	/* Get the number of entries in this map */
	int num() const { return inum; }

	/* Check if this map is empty */
	bool empty() const { return inum == 0; }

	/* Validate the state of this map. */
	void OK() const;

	/* An empty map in case you need one! */
	static const BMap EMPTY;

// not recommented:

	/*
	 * Find or create a 0-initialized item in the map and return the addr
	 * of its value.
	 */
	Univ* specialretrieveenteraddr(Univ k);

	/* Add an item only if the key does not yet exist in the map */
	int specialadd(Univ k, Univ v);	// ret was_there

	POOLALLOCATION(BMap);
#if !defined(__DECCXX)
  protected:
#endif
	struct Node {
		POOLALLOCATION(Node);
		Univ k;
		Univ v;
		Node* n;
	};
	// {b,isize,inum} could be made Array<Node*> but for efficiency
	Node** b;		// is zero if empty set!

	/* The (trivial) hash function for the map. */
	int hashk(Univ k) const { return Conv<unsigned>::d(k)%isize; }

	void removeaux(int buckn, Node* n, Node* last);

  private:
    friend class BMapIter;
	int isize;		// 0 if !b
	int inum;
	mutable int fbuckn;	// <= first non-zero index, 0 if !b ??
	void quickenter(Univ k, Univ v);
	Node* find(Univ k) const;
	void resize(int newsize);
	DISABLECOPY(BMap);
};

class Random;

class BMapIter {
  public:
	BMapIter(const BMap& hh);
	BMapIter(const BMap& hh, Random& r);
	~BMapIter() {}
	operator bool() const { return n != 0; }
	void next() { if ((n = n->n) == 0) advance(); }
	Univ key() const { return n->k; }
	Univ value() const { return n->v; }
  private:
	const BMap& h;
	int bn;
	const BMap::Node* n;
	void findrealfbuckn();
	void advance();
	void advancewrap();
	// shallow copy is safe
};

//----------------------------------------------------------------------------

// no inum++, no resize, no fbuckn
inline void BMap::quickenter(Univ k, Univ v)
{
	register int buckn=hashk(k);
	register Node* n=new Node;
	n->k=k; n->v=v; n->n=b[buckn]; b[buckn]=n;
}

inline void BMap::enter(Univ k, Univ v)
{
	inum++;
	if (!b) resize(5);
	else if (inum>isize*3) resize((isize-1)*5+3);
	quickenter(k,v);
	fbuckn=0;
}

inline BMap::Node* BMap::find(Univ k) const
{
	if (!b) return 0;
	for (Node* n=b[hashk(k)];n;n=n->n)
		if (n->k==k) return n;
	return 0;
}

inline Univ BMap::retrieve(Univ k, bool& present) const
{
	Node* n = find(k);
	present = (n != 0);
	return n ? n->v : 0;
}

inline Univ BMap::retrieve(Univ k, int& present) const
{
	bool success = false;
	Univ r = retrieve(k, success);
	present = success ? 1 : 0;
	return r;
}

inline Univ BMap::retrieve(Univ k) const
{
	Node* n=find(k);
	// return n?assertv(n->v):0;
	return n?n->v:0;
}

inline Univ BMap::replace(Univ k, Univ v)
{
	Node* n=find(k);
	Univ ov;
	return n?(ov=n->v,n->v=v,ov):0;
}

inline Univ* BMap::specialretrieveenteraddr(Univ k)
{
	Node* n = 0;
	if (!b) { enter(k,0); return &find(k)->v; }
	int buckn=hashk(k);
	for (n=b[buckn];n;n=n->n)
		if (n->k==k) return &n->v;
	n=new Node; n->k=k; n->v=0; n->n=b[buckn]; b[buckn]=n;
	fbuckn=0;
	if (++inum>isize*3) { resize((isize-1)*5+3); return &find(k)->v; }
	return &n->v;
}

inline int BMap::specialadd(Univ k, Univ v)
{
        Node* n;
	if (!b) { enter(k,v); return 1; }
	int buckn=hashk(k);
//	for (Node* n=b[buckn];n;n=n->n)
	for (n=b[buckn];n;n=n->n)
		if (n->k==k) return 0;
	n=new Node; n->k=k; n->v=v; n->n=b[buckn]; b[buckn]=n;
	fbuckn=0;
	if (++inum>isize*3) resize((isize-1)*5+3);
	return 1;
}

// Hand optimize this one
inline void BMapIter::advance()
{
	// assertx(!n);  too expensive
	// for (++bn;bn<h.isize;bn++) if (n=h.b[bn]) break;
	register int lbn=bn;
	register const int lisize=h.isize;
	register const BMap::Node** const lb = const_cast<BMap::Node const **>(h.b);
	register const BMap::Node* ln=0;
	for (++lbn;lbn<lisize;lbn++)
		if ((ln = lb[lbn]) != 0) break;
	bn=lbn; n=ln;
}

//----------------------------------------------------------------------------

template<class K, class V> class MapIter;

template<class K, class V>
class Map : public BMap {
  public:
	Map() { }
	~Map() { }
	inline void enter(K k, V v)
	{ BMap::enter(Conv<K>::e(k),Conv<V>::e(v)); }
	inline bool contains(K k) const
	{ return BMap::contains(Conv<K>::e(k)); }
	inline V retrieve(K k, int& present) const
	{ return Conv<V>::d(BMap::retrieve(Conv<K>::e(k),present)); }
	inline V retrieve(K k, bool& present) const
	{ return Conv<V>::d(BMap::retrieve(Conv<K>::e(k),present)); }
	inline V retrieve(K k) const
	{ return Conv<V>::d(BMap::retrieve(Conv<K>::e(k))); }
	inline V get(K k) const
	{ return Conv<V>::d(BMap::get(Conv<K>::e(k))); }
	inline V remove(K k)
	{ return Conv<V>::d(BMap::remove(Conv<K>::e(k))); }
	inline V replace(K k, V v)
	{ return Conv<V>::d(BMap::replace(Conv<K>::e(k),Conv<V>::e(v))); }
// 	typedef MapIter<K,V> Iter;
};

template<class K, class V>
class MapIter : public BMapIter {
  public:
	inline MapIter(const Map<K,V>& hh) : BMapIter(hh) { }
	inline MapIter(const Map<K,V>& hh, Random& r) : BMapIter(hh,r) { }
	inline ~MapIter() { }
	inline K key() const { return Conv<K>::d(BMapIter::key()); }
	inline V value() const { return Conv<V>::d(BMapIter::value()); }
};

#define ForMapKeyValue(S,T1,V1,T2,V2) \
{ for (MapIter<T1,T2> zz(S);zz;zz.next()) \
  { T1 V1=zz.key(); T2 V2=zz.value();
#define ForMapKey(S,T1,T2,V) { for (MapIter<T1,T2> zz(S);zz;zz.next()) \
			       { T1 V=zz.key();
#define ForMapValue(S,T1,T2,V) { for (MapIter<T1,T2> zz(S);zz;zz.next()) \
				 { T2 V=zz.value();

#endif
