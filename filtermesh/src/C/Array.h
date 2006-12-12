// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#ifndef Array_h
#define Array_h

#if 0
{
	static Array<Point> ar; ar.init(np);
	Array<Point> ar; ar+=Point(1,2,3);
}
#endif

// SArray<T> supports classes that have no public operator=
template<class T>
class SArray {
  public:
	inline SArray(int psize) : size(psize), a(size?new T[size]:0) { }
	inline ~SArray()			{ delete[] a; }
	// inline void ok(int i) const		{ assertx(i>=0 && i<size); }
	inline void ok(int i) const		{ }
	inline operator const T*() const	{ return a; }
	inline operator T*()			{ return a; }
	inline const T& operator[](int i) const	{ ok(i); return a[i]; }
	inline T& operator[](int i)		{ ok(i); return a[i]; }
  protected:
	int size;		// must come before a
	T* a;
  private:
	inline void operator=(SArray<T>&) { } // disable
	inline SArray(SArray<T>&) { }	      // disable
};

// for SArray<T>, T must have a public operator=
template<class T>
class Array : public SArray<T> {
  private:
	int n;
	inline void resize(int ns) {
		size=ns;
		// assertx(n<=size);
		T* na=size?new T[size]:0;
		// in next line, Array<T> need for GNUG 2.4.5 sometimes
		for (int i=0;i<(Array<T>::n);i++) na[i]=a[i];
		delete[] a,a=na;
	}
  public:
	inline Array(int ne=0)			: SArray<T>(ne), n(ne) { }
	inline ~Array()				{ }
	inline int num() const			{ return n; }
	inline void clear()			{ delete[] a,a=0; n=size=0; }
	inline void init(int ne) { // allocate ne, CLEAR if too small
		if (ne>size) { clear(); resize(ne); }
		n=ne;
	}
	inline void need(int ne) { // allocate at least ne, COPY if too small
		if (ne>size) resize(::max(int(n*1.5)+3,ne));
		n=ne;
	}
	inline int add(int ne) { // ret: previous num()
		int cn=n; need(n+ne); return cn;
	}
	inline Array<T>& operator+=(const T& e) {
		int i=add(1); a[i]=e; return *this; // add() may allocate a[]
	}
	inline void squeeze() 			{ if (n<size) resize(n); }
	// tighter assertions
	// inline void ok(int i) const		{ assertx(i>=0 && i<n); }
	inline void ok(int i) const		{ } 
	inline const T& operator[](int i) const	{ ok(i); return a[i]; }
	inline T& operator[](int i)		{ ok(i); return a[i]; }
#ifdef CXX
	// next for stupid cxx
	inline operator const T*() const
	{ return SArray<T>::operator const T*(); }
	inline operator T*()
	{ return SArray<T>::operator T*(); }
#endif
  private:
	inline void operator=(Array<T>&) { } // disable
	inline Array(Array<T>&) : SArray<T>(0) { } // disable
};

#endif
