// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Map.h"		// for hhwarning()

#include <sys/types.h>		// for socket struct, u_short in hton
#include <sys/socket.h>		// for socketpair
#include <fcntl.h>		// for fcntl
#include <cerrno>		// for errno with __GNUG__

#include <sys/time.h>		// gettimeofday(), timeval, setitimer()
#include <ctime>		// gettimeofday(), timeval, time(), ctime()
#include <new>			// set_new_handler()

#include <cstdarg>		// to use with hform()

#include <iostream>
using std::cout;
using std::cerr;

// Note: normally, these are defined to act on unsigned long (32bit),
//  but I assume that 'int' is 32 bit, so it's ok on alphas too
//extern "C" {			// <netinet/in.h>
//	unsigned int htonl(unsigned int),ntohl(unsigned int);
//	unsigned short htons(unsigned short),ntohs(unsigned short);
//}

extern "C" char** environ;

extern "C" {
#if defined(__DECCXX)
	int socketpair(int, int, int, int sv[2]);
	int gettimeofday(struct timeval*, struct timezone*);
#ifndef __alpha
	int putenv(const char* string);
#endif
#elif __GNUG__==2 && __GNUC_MINOR__>=5
	int socketpair(int, int, int, int sv[2]);
//	int putenv(const char* string);
	int setitimer(int, const struct itimerval*,struct itimerval*);
	int gettimeofday(struct timeval*, struct timezone*);
#elif defined(__alpha)
	int socketpair(int, int, int, int sv[2]);
#endif
}

static void new_failed()
{
	assertnever("new is out of memory");
}

static int init()
{
	std::set_new_handler(&new_failed);

// JW, 2006-10-26: These assumptions _should_ no longer hold.  Cross your
//                 fingers, folks!
//
//	// Assumptions present in my C++ code (Hugues Hoppe)
//	assertx(sizeof(int)>=4);
//	assertx(sizeof(void*)>=sizeof(int));
//	assertx(sizeof(void*)>=sizeof(double));
	return 0;
}

static int dummy=init();

void* hmalloc(void* oldp, int size)
{
	assertx(size>=0);
	void* r = oldp ? realloc(oldp, size) : malloc(size);
	return assertv(r);
}

void hfree(void* e)
{
	if (e) free(e);
}

int Setfdnodelay(int fd, int nodelay)
{
#ifdef sgi
	// on sgi, setting nodelay on terminal fd may cause window closure
	if (nodelay) assertx(!isatty(fd));
#endif
	return fcntl(fd,F_SETFL,nodelay?
#if defined(hpux) || defined(sgi)
		     O_NONBLOCK
#elif defined(FNDELAY)
		     FNDELAY
#elif defined(O_NDELAY)
		     O_NDELAY
#endif
		     :0)!=-1;
}

void SpawnOnSockets(int* fdi, int* fdo, FILE** fi, FILE** fo,
		    void (*spawn)())
{
	int s[2],child;
	assertx(!socketpair(AF_UNIX,SOCK_STREAM,0,s));
	assertx((child=fork())>=0);
	if (!child) {
		close(s[1]);
		dup2(s[0],0);
		dup2(s[0],1);
		if (s[0]<2)
			close(s[0]);
		(*spawn)();
		_exit(0);
	}
	close(s[0]);
	if (fi) *fi=fdopen(s[1],"r");
	if (fo) *fo=fdopen(s[1],"w");
	if (fdi) *fdi=s[1];
	if (fdo) *fdo=s[1];
}


static struct Warningclass {
	~Warningclass() { flush(); }
	void flush();
	Map<const char*,int> m; // string -> number of times printed
} Warnings;

void Warningclass::flush()
{
	if (m.empty()) return;
	SHOWDF("Summary of warnings:\n");
	ForMapKeyValue(m,const char*,s,int,n) {
		SHOWDF(" %5d '%s'\n",n,s);
	} EndFor;
	m.clear();
}

void FlushWarnings() { Warnings.flush(); }

int hhassert(int warncode, const char* ps, const char* file, int line)
{
	const char* s=(ps?ps:"?");
	if (warncode>=0) {
		if (Warnings.m.contains(s)) {
			Warnings.m.replace(s,Warnings.m.get(s)+1);
			if (warncode==1) return 0;
		} else {
			Warnings.m.enter(s,1);
		}
	}
	cerr << (warncode<0?"Fatal assertion error":"assertion warning") <<
		": " << s << " in line " << line <<
		" of file " << file << "\n";
	if (GetenvValue("ASSERTABORT") ||
	    warncode>=0 && GetenvValue("ASSERTWABORT") ||
	    warncode<0 && GetenvValue("ASSERTXABORT")) abort();
	if (warncode<0) {
		if (errno) perror("possible error");
		_exit(1);
	}
	return 1;
}

double GetPreciseTime()
{
	struct timeval ti;
	assertx(!gettimeofday(&ti,0));
	return double(ti.tv_sec)+double(ti.tv_usec)*1e-6;
}

int GetTime()
{
	return time(0);
}

const char* CTime()
{
	time_t ti=GetTime();
	char* s=ctime(&ti);
	s[24]=0;		// strip out trailing '\n'
	return s;
}


//------------------------------------------------------------------------
//*** conversions to/from standard byte order

static union {
	unsigned char little_endian;
	unsigned int  d;
} endianness;

static int init_endianness()
{
	endianness.d = 1;
        return 0;
}

static int _init_endianness = init_endianness();

void ReverseMemory(char *region, int len)
{
	char *from = region, *to = region + len;
	do {
		char c = *from;
		*from++ = *to;
		*--to = c;
	}
	while (to != region);
}

void FloatToStd(double const* in, char* ou)
{
	if (endianness.little_endian)
	{
		double d = *in;
		ReverseMemory(reinterpret_cast<char *>(&d), sizeof(double));
		memcpy(ou, &d, sizeof(double));
	}
	else
		memcpy(ou, in, sizeof(double));
}

void StdToFloat(char const* in, double* ou)
{
	if (endianness.little_endian)
	{
		double d;
		memcpy(&d, in, sizeof(double));
		ReverseMemory(reinterpret_cast<char *>(&d), sizeof(double));
		*ou = d;
	}
	else
		memcpy(ou, in, sizeof(double));
}

void IntToStd(unsigned int const* in, char* ou)
{
	* reinterpret_cast<unsigned int *>(ou) = htonl(* in);
}

void StdToInt(char const* in, unsigned int* ou)
{
	* ou = ntohl(* reinterpret_cast<unsigned int const *>(in));
}

void ShortToStd(unsigned short const* in, char* ou)
{
	* reinterpret_cast<unsigned short *>(ou) = htons(* in);
}

void StdToShort(char const* in, unsigned short* ou)
{
	* ou = ntohs(* reinterpret_cast<unsigned short const *>(in));
}

static void nullfunc(int)
{
	signal(SIGALRM,HHSIG_PF(nullfunc)); // for ATT unix
}

void Delay(double sec)
{
	struct itimerval t;
	signal(SIGALRM,HHSIG_PF(nullfunc));
	t.it_value.tv_sec=int(sec);
	t.it_value.tv_usec=int(sec*1e6)%1000000;
	t.it_interval.tv_sec=0;
	t.it_interval.tv_usec=10000; // 10 millisec, so pause() doesn't stick
	if (setitimer(ITIMER_REAL,&t,0)) perror("setitimer");
	pause();
	t.it_value.tv_sec=0;
	t.it_value.tv_usec=0;
	t.it_interval.tv_sec=0;
	t.it_interval.tv_usec=0;
	if (setitimer(ITIMER_REAL,&t,0)) perror("setitimer");
	// SIGALRM not reset to SIG_DFL just to be sure
}


const int bufsiz=512;

// Warning: I have found that routines which use va_list should not be
// compiled with -O under cxx.
// Could be because vprintf() routines are defined in terms of varargs.h
// and not stdarg.h

const char* hform(const char* format, ...) // not reentrant!!!!
{
	const int numbuf=10;
	static char buf[numbuf][bufsiz];
	static int bufn;
	va_list ap;
	va_start(ap,format);
	bufn=(bufn+1)%numbuf;
	vsnprintf(buf[bufn], bufsiz, format, ap);
	va_end(ap);
	return buf[bufn];
}

void SHOWF(const char* format, ...)
{
	char buf[bufsiz];
	va_list ap;
	va_start(ap,format);
	vsnprintf(buf,bufsiz,format,ap);
	va_end(ap);
	cerr << buf;
}

// old
//	isatty1	isatty2	cout	cerr
//	0	0	1	0
//	0	1	1	1
//	1	0	1	0
//	1	1	0	1
// new
//	isatty1	isatty2	cout	cerr
//	0	0	0	1
//	0	1	1	1
//	1	0	1	1
//	1	1	0	1
void SHOWDF(const char* format, ...)
{
	static int done=0, needcout, needcerr;
	if (!done++) {
		int is1=isatty(1), is2=isatty(2);
		// old
		needcout=!is1 || !is2;
		needcerr=is2;
		// new
		needcout=is1+is2==1;
		needcerr=1;
		// SHOWN(is1); SHOWN(is2); SHOWN(needcout); SHOWN(needcerr);
	}
	char buf[bufsiz];
	va_list ap;
	va_start(ap,format);
	vsnprintf(buf,bufsiz,format,ap);
	va_end(ap);
	//if (needcout) cout << "# " << buf;
	if (needcerr) cerr << "# " << buf;
}

void SHOWFF(const char* format, ...)
{
	static int done=0, wantcout;
	if (!done++) wantcout=!isatty(1);
	if (!wantcout) return;
	char buf[bufsiz];
	va_list ap;
	va_start(ap,format);
	vsnprintf(buf,bufsiz,format,ap);
	va_end(ap);
	//cout << "# " << buf;
	cerr << "# " << buf;
}

char* newString(const char* s)
{
	return strcpy(new char[strlen(s)+1],s);
}

int GetenvValue(const char* varname)
{
	const char* s=getenv(varname);
	if (!s || *s=='0' || *s=='n' || *s=='N' || *s=='f' || *s=='F')
		return 0;
	return max(atoi(s),1);
}



// static const char* _findenv(const char* name, int& offset)
// {
//         const char* c;
// 	int len=0;
// //	for (const char* c=name;*c && *c!='=';c++,len++);
// 	for (c=name;*c && *c!='=';c++,len++);
//         for (char** p=environ;*p;p++)
//                 if (!strncmp(*p,name,len))
//                         if (*(c=*p+len)=='=') {
//                                 offset=p-environ;
//                                 return ++c;
//                         }
//         return 0;
// }

#if 0
// Access to 'environ' is not reentrant
void unsetenv(const char* name)
{
	int offset;
	while (_findenv(name,offset)) // if set multiple times
                for (char** p=&environ[offset];;p++)
                        if (!(*p=*(p+1)))
                                break;
}
#endif

// return: 0=success, -1=out_of_mem
// int setenv(const char* name, const char* value, int overwrite)
// {
// 	assertx(overwrite);
// 	return putenv(newString(hform("%s=%s",name,value)));
// }

const char* CreateHeader(int argc, char const** argv)
{
	int len=0;
	const char* s="# ";
	for (int i=0;i<argc;i++) {
		len+=strlen(argv[i])+1;
		if (i && len>75) s=hform("%s\n#  ",s),len=strlen(argv[i])+1;
		s=hform("%s %s",s,argv[i]);
	}
	return hform("%s\n",s);
}
