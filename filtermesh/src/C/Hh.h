// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Hh_h
#define Hh_h

#include <cstdlib>
#include <cmath>
#include <unistd.h>		// UNIX system calls
#include <csignal>
#include <cstdio>
#include <cstring>
#include <netinet/in.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdint.h>

#if defined(sgi) && !defined(__GNUG__)
#define SGICC
#endif

#if defined(__DECCXX) && defined(__alpha)
typedef void (*HHSIG_PF)(int);
#elif defined(__GNUG__)
typedef sighandler_t HHSIG_PF;
// #define HHSIG_PF SignalHandler
#elif defined(__DECCXX)
#define HHSIG_PF sigvec_t
#elif defined(ATTCC)
#define HHSIG_PF SIG_PF
#elif defined(SGICC)
typedef void (*HHSIG_PF)(...);
#else
??
#endif

#ifdef __DECCXX
extern "C" void bcopy(const void* src, void* dst, int nbytes);
#endif
#if defined(sparc) || defined(__alpha) || defined(__FreeBSD__)
#define fsqrt sqrt
#elif defined(__GNUG__)
// extern "C" double fsqrt(double);
#endif

extern int hhassert(int warncode, const char* s, const char* file, int line);
template<class T> inline T assertvaux(T expr, const char* s,
				      const char* file, int line)
{ if (!expr) hhassert(-1,s,file,line); return expr; }
// Warning: ret: message_printed
#define Warning(s) hhassert(1, s, __FILE__, __LINE__)
// assertnever: ret: void
#define assertnever(str) \
(void)hhassert(-1, "assertnever(" str ")", __FILE__, __LINE__)
// assertx: ret: void
#define assertx(expr) \
((void)((expr)?0:hhassert(-1, "assertx(" #expr ")", __FILE__, __LINE__)))
// assertw: ret: expr_failed
#define assertw(expr) \
((expr)?0:(hhassert(+0, "assertw(" #expr ")", __FILE__, __LINE__),1))
// assertw1: ret: expr_failed
#define assertw1(expr) \
((expr)?0:(hhassert(+1, "assertw1(" #expr ")", __FILE__, __LINE__),1))
// assertv: ret: val
#define assertv(e) \
assertvaux(e, "assertv(" #e ")", __FILE__, __LINE__)

extern void FlushTimers();
extern void FlushStats();
extern void FlushWarnings();
inline void CleanUp() { FlushTimers(); FlushStats(); FlushWarnings(); }

typedef uint64_t Univ;

inline double i64tofloat(uint64_t const &i)
{ return * reinterpret_cast<double const *>(&i); }

inline uint64_t floattoi64(const double& f)
{ return * reinterpret_cast<uint64_t const *>(&f); }

template<class T>
class Conv {
  public:
	inline static Univ e(T e) { return static_cast<Univ>(e); }
	inline static T d(Univ e) { return static_cast<T>(e); }
};

template<class T>
class Conv<T *> {
  public:
	inline static Univ e(T *e) { return reinterpret_cast<Univ>(e); }
	inline static T *d(Univ e) { return reinterpret_cast<T *>(e); }
};

class Conv<double> {
  public:
	inline static Univ e(double e)
	    { return reinterpret_cast<Univ>(floattoi64(e)); }
	inline static double d(Univ e)
	{
	    uint64_t v = reinterpret_cast<uint64_t>(e);
	    return i64tofloat(v);
	}
};

extern void* hmalloc(void* oldp, int size);
extern void hfree(void* e);

#ifndef PI
const double PI=3.14159265358979323846;
#endif

template<class T>
inline void hqsort(T* ar, int n, int (*cmp)(const T*, const T*))
{ qsort(ar,n,sizeof(T),reinterpret_cast<int(*)(const void*,const void*)>(cmp)); }

// GNUG 2.5.0 sometimes requires :: before min() and max()
template<class T> inline void swap(T* e1, T* e2) { T e=*e1; *e1=*e2; *e2=e; }
template<class T> inline T sign(T e) { return e>0?1:e<0?-1:0; }
template<class T> inline T abs(T e) { return e<0?-e:e; }
template<class T> inline T square(T e) { return e*e; }
template<class T> inline T min(T a, T b) { return a<b?a:b; }
template<class T> inline T max(T a, T b) { return a>b?a:b; }

inline double torad(double deg) { return deg/180*PI; }
inline double todeg(double rad) { return rad/PI*180; }

// Useful for "NEST { int localvar; ... }" so that indentation is correct
#define NEST
#define bcase break; case
#define ocase case
#define bdefault break; default
#define ForIndex(i,upb) { int zz=upb; for (int i=0;i<zz;i++) {
#define DummyEndFor }} // for gemacs formatting
#define For {{ for
#define EndFor }}
#define DISABLECOPY(class) class& operator=(class&); class(class&)


// SHOWN(x*y);
// SHOW(point+vector);
// SHOWS("Could not open display");
// SHOWF("%s: Argument '%s' ambiguous, '%s' assumed\n",argv[0],arg,assumed);
// SHOWL;
// os << hform(" Endprogram: %dgons %dlines\n",ngons,nlines);

// Show an expression
#define SHOW(x) (std::cerr << #x << " = " << (x))
// Show an expression, and terminate with a newline
#define SHOWN(x) (std::cerr << #x << " = " << (x) << "\n")
// Show a string
#define SHOWS(s) (std::cerr << (s) << "\n")
// Formatted show, uses small buffer
extern void SHOWF(const char* format, ...);
// Output a string prefixed with "# " on cerr, and on cout if _is_a_file
extern void SHOWDF(const char* format, ...);
// Output a string prefixed with "# " on cout only if _is_a_file
extern void SHOWFF(const char* format, ...);
// Show current line number and file
#define SHOWL (std::cerr << "Now in " << __FILE__ << " at line " << \
	       __LINE__ << "\n")
// Format a string.  Uses a few rotating small buffers
// use { char s[bigv]; sprintf() or ostrstream() } if you need a big buffer
extern const char* hform(const char* format, ...);

// Allocate a duplicate of string, using new char[], use delete[]!
extern char* newString(const char* s);

extern int GetenvValue(const char* varname);

// Returns fcntl return value
extern int Setfdnodelay(int fd, int nodelay);

extern void SpawnOnSockets(int* fdi, int* fdo,
			   FILE** fi, FILE** fo, void (*spawn)());

// Return absolute time, in seconds
extern double GetPreciseTime();
// Return absolute time, in seconds
extern int GetTime();
// Return string of time, without a trailing '\n'
extern const char* CTime();
// Delay for some number of seconds
extern void Delay(double sec);

// These could be made inline but would require including some .h files
extern void FloatToStd(double const* in, char* ou);
extern void StdToFloat(char const* in, double* ou);
extern void IntToStd(unsigned int const* in, char* ou);
extern void StdToInt(char const* in, unsigned int* ou);
extern void ShortToStd(unsigned short const* in, char* ou);
extern void StdToShort(char const* in, unsigned short* ou);

// Prevent NaN's from appearing.
inline double myacos(double a)
{
	assertw(fabs(a)<1.001);
	if (a<-1) a=-1;
	if (a>1) a=1;
	return acos(a);
}

inline double mysqrt(double a)
{
	assertx(a>-1e-5);
	double b=a<0?0:a;
//	return fsqrt(b);
	return sqrtf(b);
}

// extern void unsetenv(const char* name);
// extern int setenv(const char* name, const char* value,
// 		  int overwrite); // 0=success

// Multi-line header, each line begins with "# " and ends with "\n"
extern const char* CreateHeader(int argc, char const** argv);

#endif
