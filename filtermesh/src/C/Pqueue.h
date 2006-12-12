// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Pqueue_h
#define Pqueue_h

#include "Map.h"
#include "Pool.h"

// Self-resizing priority queue
class BPqueue {
  public:
	BPqueue(int size=0);
	virtual ~BPqueue();
	virtual void clear();
	void enter(Univ e, double pri); // pri>=0 always
	int empty() const;
	int num() const;
	Univ min() const;	// ret e, else die
	double minpriority() const; // ret pri, else die
	virtual Univ removeminI();  // ret e, else die
	virtual void enterUnsorted(Univ e, double pri);
	virtual void sort();
  protected:
	struct Node {
		Univ e;
		double pri;
	};
	int isize;
	int inum;
	Node* ar;
	// {ar,isize,inum} could be made Array<Node> but for efficiency
	void nswitch(int n1, int n2);
	void resize(int newsize);
	virtual void vnswitch(int n1, int n2);
	virtual void vadjust(int n);
  private:
	void padjust(int n, int uptoo);
	DISABLECOPY(BPqueue);
};

// Hashed priority queue
class BHPqueue : public BPqueue {
  public:
	BHPqueue(int size=0);
	~BHPqueue();
	void clear();
	Univ removeminI();
	double retrieve(Univ e) const; // ret pri or <0
	double remove(Univ e);	      // ret pri or <0
	double update(Univ e, double pri); // ret prevpri or <0
	double enterupdate(Univ e, double pri); // ret prevpri or <0
	void enterUnsorted(Univ e, double pri);
	void sort();
  protected:
	void vnswitch(int n1, int n2);
	void vadjust(int n);
  private:
	Map<Univ,int> m;	// element -> index in array
	int find(Univ e) const;	// ret -1 if not there
	void padjust(int n, int uptoo);
	void nswitch(int n1, int n2);
	void adjust(int n);
};

//----------------------------------------------------------------------------

inline int BPqueue::empty() const { return !inum; }
inline int BPqueue::num() const { return inum; }
inline Univ BPqueue::min() const { assertx(inum); return ar[0].e; }
inline double BPqueue::minpriority() const { assertx(inum); return ar[0].pri; }

inline void BPqueue::enterUnsorted(Univ e, double pri)
{
	assertx(pri>=0);
	if (inum==isize) resize(isize*2);
	ar[inum].e=e;
	ar[inum].pri=pri;
	inum++;
}

inline void BPqueue::enter(Univ e, double pri)
{
	enterUnsorted(e,pri);
	vadjust(inum-1);	// virtual call
}

//----------------------------------------------------------------------------

template<class T>
class Pqueue : public BPqueue {
  public:
	Pqueue() { }
	~Pqueue() { }
	inline void enter(T e, double pri)
	{ BPqueue::enter(Conv<T>::e(e),pri); }
	inline T min() const { return Conv<T>::d(BPqueue::min()); }
	inline T removemin() { return Conv<T>::d(BPqueue::removeminI()); }
	inline void enterUnsorted(T e, double pri)
	{ BPqueue::enterUnsorted(Conv<T>::e(e),pri); }
};

template<class T>
class HPqueue : public BHPqueue {
  public:
	HPqueue() { }
	~HPqueue() { }
	inline void enter(T e, double pri)
	{ BHPqueue::enter(Conv<T>::e(e),pri); }
	inline T min() const { return Conv<T>::d(BHPqueue::min()); }
	inline T removemin() { return Conv<T>::d(BHPqueue::removeminI()); }
	inline void enterUnsorted(T e, double pri)
	{ BHPqueue::enterUnsorted(Conv<T>::e(e),pri); }
	inline double retrieve(T e) const
	{ return BHPqueue::retrieve(Conv<T>::e(e)); }
	inline double remove(T e) { return BHPqueue::remove(Conv<T>::e(e)); }
	inline double update(T e, double pri)
	{ return BHPqueue::update(Conv<T>::e(e),pri); }
	inline double enterupdate(T e, double pri)
	{ return BHPqueue::enterupdate(Conv<T>::e(e),pri); }
};

#endif
