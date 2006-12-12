// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef HashStruct_h
#define HashStruct_h

#if 0
{
	HashStruct hs<segment>;
	ForHashStruct(hs,segment,s) { do(s); } EndFor;
}
#endif

#include "Map.h"

// Maintains a set of structures.
// The structures are hashed according to a user-provided function hashf.
//  Note: hashf can return any value including zero.
// Also, equality of two structures is determined by user function cmpf.
//  Note: cmpf should return: does_not_match (as does strcmp)!
//  Note for hacker: when cmpf is invoked, its first arg is st (below)
//   Thus, st (structure template) can be fake and cmpf can be asymmetric.
// Note: BHashStruct should be empty prior to destruction.

class BHashStruct : public BMap {
  public:
	typedef Univ (*HASHF)(Univ);
	typedef int (*CMPF)(Univ, Univ);
	BHashStruct(HASHF phashf, CMPF pcmpf);
	~BHashStruct();
	void enter(Univ s);	// s!=0, s must be new
	int add(Univ s);	// ret: is_new
	Univ retrieve(Univ st) const; // ret sfound or 0
	int remove(Univ st);	      // ret wasfound
	Univ getone() const;	      // die if empty
	Univ removeone();	      // die if empty
  private:
	HASHF hashf;
	CMPF cmpf;
	Univ replace(Univ k, Univ v); // disable ?? permitted
};

//----------------------------------------------------------------------------

inline Univ BHashStruct::getone() const {
	BMapIter si(*this); return si.value();
}
inline Univ BHashStruct::removeone() {
	Univ e=getone(); remove(e); return e;
}

//----------------------------------------------------------------------------

template<class T> class HashStructIter;

// GNUG chokes on Conv<T*>, so I did it manually here
template<class T>
class HashStruct : public BHashStruct {
  public:
	typedef Univ (*HASHF)(const T*);
	typedef int (*CMPF)(const T*, const T*);
	HashStruct(HASHF phashf, CMPF pcmpf) :
	BHashStruct(reinterpret_cast<BHashStruct::HASHF>(phashf),reinterpret_cast<BHashStruct::CMPF>(pcmpf)) { }
	~HashStruct() { }
	inline void enter(T* e) { BHashStruct::enter(Conv<T *>::e(e)); }
	inline T* retrieve(const T* e) const
	{ return Conv<T *>::d(BHashStruct::retrieve(Conv<T const *>::e(e))); }
	inline int add(T* e) { return BHashStruct::add(Conv<T *>::e(e)); }
	inline int remove(const T* e)
	{ return BHashStruct::remove(Conv<T const *>::e(e)); }
	inline T* getone() const { return Conv<T *>::d(BHashStruct::getone()); }
	inline T* removeone() { return Conv<T *>::d(BHashStruct::removeone()); }
// 	typedef HashStructIter<T> Iter;
};

template<class T>
class HashStructIter : public BMapIter {
  public:
	HashStructIter(const HashStruct<T>& s) : BMapIter(s) { }
	HashStructIter(const HashStruct<T>& s, Random& r) : BMapIter(s,r) { }
	~HashStructIter() { }
	inline T* operator()() const { return Conv<T*>::d(BMapIter::value()); }
  private:
	void key() const;	// disable
	void value() const;	// disable
};

#define ForHashStruct(S,T,V) \
{ for (HashStructIter<T> zz(S);zz;zz.next()) { T* V=zz();

#endif
