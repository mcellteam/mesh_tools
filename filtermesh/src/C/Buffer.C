// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include <cerrno>
#include "Hh.h"
#include "Buffer.h"

#if __GNUG__==2 && __GNUC_MINOR__==5
extern "C" bcopy(const void*, void*, _G_size_t);
#endif

const int INITIALSIZE=4096;
const int READSIZE=1024;
const int WRITESIZE=8192;


//------------------------------------------------------------------------
//*** Buffer

Buffer::Buffer(int pfd, int pmaxsize)
{
	fd=pfd;
	maxsize=pmaxsize;
	size=0;
	ar=0;
	beg=0;
	n=0;
	steof=sterr=steot=0;
	blocked=0;
	detect_to=0;
	timev=0;
}

Buffer::~Buffer()
{
	delete[] ar;
}

void Buffer::shift()
{
	assertx(ar && beg);
	if (n) bcopy(&ar[beg],&ar[0],n);
	beg=0;
}

void Buffer::expand()
{
	assertw(!beg);		// not necessary, current implementation
	if (maxsize) assertx(size<maxsize);
	int nsize=!size?INITIALSIZE:size*2;
	if (maxsize && nsize>maxsize) nsize=maxsize;
	char* nar=new char[nsize];
	if (ar) bcopy(&ar[beg],&nar[beg],n);
	delete[] ar;
	ar=nar;
	size=nsize;
}

double Buffer::getcurtime()
{
	return timev?*timev:GetPreciseTime();
}

void Buffer::recordtime()
{
	if (detect_to)
		ti_lastop=getcurtime();
}

void Buffer::checkeot()
{
	if (!blocked || !detect_to)
		return;
	if (getcurtime()-ti_lastop>ti_to)
		steot=1;
}

// 0= use system calls
double* Buffer::usetime(double* timevar)
{
	double* old=timev;
	timev=timevar;
	return old;
}

// 0= no timeouts
double Buffer::settimeout(double timeinter)
{
	double old=ti_to;
	ti_to=timeinter;
	detect_to=ti_to!=0;
	if (blocked) recordtime();
	return old;
}

void Buffer::reset_timeout()
{
	steot=0;
	if (blocked) recordtime();
}


//------------------------------------------------------------------------
//*** RBuffer

// 0=no, 1=yes, -1=other
// if buffer is full, shift() || expand()
int RBuffer::fill()
{
	assertx(!eof() && !eot() && !err());
	if (beg+n==size && beg)
		shift();
	else if (n==size)
		expand();
	int ntoread=min(size-beg-n,READSIZE),nread;
	assertx(ntoread);
	for (;;) {
		nread=read(fd,&ar[beg+n],ntoread);
		if (nread<0) {
			if (errno==EINTR)
				continue; // for ATT UNIX (hpux)
			if (errno==EWOULDBLOCK || errno==EAGAIN) {
				if (!blocked) {
					blocked=1;
					recordtime();
				} else {
					checkeot();
				}
				return steot?-1:0;
			}
			// perror("RBuffer:read:");
			sterr=1;
			return -1;
		} else if (!nread) {
			steof=1;
			return -1;
		}
		break;
	}
	blocked=0;
	n+=nread;
	return 1;
}

// have read n bytes
void RBuffer::extract(int num)
{
	assertx(num && num<=n);
	n -= num;
	beg += num;
	if (!n) beg=0;
	// if (beg%4) shift();
}

int RBuffer::hasline() const
{
	int i=0;
	for (;;) {
		if (++i>=num()) return 0;
		if ((*this)[i]=='\n') return 1;
	}
}

int RBuffer::examineline(char* s, int len)
{
        int i = 0;
//	for (int i=0;i<len;i++) {
	for ( ; i < len; i++) {
		if (i >= n) return 0;
		s[i] = ar[beg + i];
		if (s[i]=='\n') break;
	}
	assertx(++i<len);
	s[i]=0;
	return 1;
}

int RBuffer::extractline(char* s, int len)
{
	if (! examineline(s,len)) return 0;
	extract(strlen(s));
	return 1;
}

char RBuffer::bgetchar(int bi) const
{
	return (*this)[bi];
}

int RBuffer::bgetint(int bi) const
{
	assertx(bi>=0 && bi+4<=n);
	char s[4];
	for (int i=0;i<4;i++) s[i]=ar[beg+bi+i];
	int r;
	StdToInt(s, reinterpret_cast<unsigned int *>(&r));
	return r;
}

short RBuffer::bgetshort(int bi) const
{
	assertx(bi>=0 && bi+2<=n);
	char s[2];
	for (int i=0;i<2;i++) s[i]=ar[beg+bi+i];
	short r;
	StdToShort(s, reinterpret_cast<unsigned short *>(&r));
	return r;
}

double RBuffer::bgetfloat(int bi) const
{
	assertx(bi>=0 && bi+4<=n);
	char s[4];
	for (int i=0;i<4;i++) s[i]=ar[beg+bi+i];
	double r;
	StdToFloat(s,&r);
	return r;
}

//------------------------------------------------------------------------
//*** WBuffer

// 0=part, 1=complete, -1=other
// no alignment problem to worry about since buffer is never word-accessed
// (at least not by user)
int WBuffer::flush(int nb)
{
	int nwritten;
	if (!nb) nb=n;
	for (;;) {
		assertx(nb<=n);
		if (!nb) return 1;
		nwritten=write(fd,&ar[beg],nb);
		if (nwritten<0) {
			if (errno==EINTR)
				continue;
			if (errno==EWOULDBLOCK || errno==EAGAIN) {
				if (!blocked) {
					blocked=1;
					recordtime();
				} else {
					checkeot();
				}
				return steot?-1:0;
			}
			// perror("WBuffer:write:");
			sterr=1;
			return -1;
		}
		assertx(nwritten>0 && nwritten<=nb);
		blocked=0;
		beg+=nwritten;
		n-=nwritten;
		if (!n) beg=0;
		nb-=nwritten;
	}
}

// if full, shift() || expand()
void WBuffer::bput(const void* buf, int nbytes)
{
	assertx(!eof() && !eot() && !err());
	if (beg+n+nbytes>size && beg)
		shift();
	while (n+nbytes>size)
		expand();
	bcopy(buf,&ar[beg+n],nbytes);
	n+=nbytes;
	if (n>=WRITESIZE)
		flush(WRITESIZE);
}

void WBuffer::bput(char c)
{
	bput(&c,1);
}

void WBuffer::bput(short i)
{
	unsigned short ti=i;
	char buffer[sizeof(unsigned short)];
	ShortToStd(&ti, buffer);
	bput(buffer, sizeof(unsigned short));
}

void WBuffer::bput(int i)
{
	char buffer[sizeof(unsigned int)];
	IntToStd(reinterpret_cast<unsigned int *>(&i), buffer);
	bput(buffer, sizeof(unsigned int));
}

void WBuffer::bput(double f)
{
	double tf=f;
	char buffer[sizeof(double)];
	FloatToStd(&tf, buffer);
	bput(buffer, sizeof(double));
}
