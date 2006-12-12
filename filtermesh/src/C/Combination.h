// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Combination_h
#define Combination_h

#if 0
{
	Combination<Vertex> comb;
	ForCombination(comb,Vertex,v,val) { do(v,val); } EndFor;
}
#endif

#include "Map.h"
#include "Pool.h"

class BCombination {
  public:
	BCombination();
	~BCombination();
	void clear();
	double operator[](Univ e) const;
	double& operator[](Univ e); // use const version for lookup!
	int num() const;	   // should squeeze first
	int empty() const;	   // should squeeze first
	double sum() const;	// slow
	void squeeze() const;	// virtual const
	static const BCombination EMPTY;
	POOLALLOCATION(BCombination);
  private:
    friend class BCombinationIter;
	Map<Univ,double> m;
	DISABLECOPY(BCombination);
};

class BCombinationIter {
  public:
	BCombinationIter(const BCombination& c);
	~BCombinationIter() { }
	operator bool() const;
	void next();
	Univ elem() const;
	double value() const;
  private:
	MapIter<Univ,double> mi;
	DISABLECOPY(BCombinationIter);
};

//----------------------------------------------------------------------------

inline double BCombination::operator[](Univ e) const {
	// if not present double(0) returned, hack
	return m.retrieve(e);
}
inline double& BCombination::operator[](Univ e) {
	Univ* p=m.specialretrieveenteraddr(e);
	// if not there, adds 0, which happens to equal double(0), hack
	return * reinterpret_cast<double *>(p);
}
inline int BCombination::num() const { return m.num(); }
inline int BCombination::empty() const { return m.empty(); }


inline BCombinationIter::BCombinationIter(const BCombination& c)
: mi(c.m) { }
inline BCombinationIter::operator bool() const { return mi; }
inline void BCombinationIter::next() { mi.next(); }
inline Univ BCombinationIter::elem() const { return mi.key(); }
inline double BCombinationIter::value() const { return mi.value(); }

//----------------------------------------------------------------------------

template<class T> class CombinationIter;

template<class T>
class Combination : public BCombination {
  public:
	Combination() { }
	~Combination() { }
	inline double operator[](T e) const {
		return BCombination::operator[](Conv<T>::e(e));
	}
	inline double& operator[](T e) {
		return BCombination::operator[](Conv<T>::e(e));
	}
// 	typedef CombinationIter<T> Iter;
};

template<class T>
class CombinationIter : public BCombinationIter {
  public:
	inline CombinationIter(const Combination<T>& c)
		: BCombinationIter(c) { }
	inline ~CombinationIter() { }
	inline T elem() const { return Conv<T>::d(BCombinationIter::elem()); }
};

#define ForCombination(S,T1,V1,V2) \
{ for (CombinationIter<T1> zz(S);zz;zz.next()) \
  { T1 V1=zz.elem(); double V2=zz.value();
    
#endif
