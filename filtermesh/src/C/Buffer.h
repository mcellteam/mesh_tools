// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Buffer_h
#define Buffer_h

class Buffer {
  public:
	Buffer(int fd, int maxsize=0); // default is infinite size
	~Buffer();
	int eof() const;	// end of file
	int err() const;	// error in system call
	int eot() const;	// end of time
	double* usetime(double* timevar); // 0= use system calls (default)
	double settimeout(double timeinter); // 0= no timeouts (default)
	void reset_timeout();		     // clear eot()
  protected:
	int fd;			// file descriptor associated
	int size;		// current size of ar[]
	char* ar;		// buffer
	int beg;		// index of first element in ar[]
	int n;			// number of elements in buffer (beg+n<=size)
	int steof;
	int sterr;
	int steot;
	int blocked;		// last operation would have blocked
	void checkeot();
	void recordtime();	// success, update ti_lastop
	void shift();		// shift data to beginning of buffer
	void expand();		// increase buffer size
  private:
	int maxsize;		// maximum size buffer can expand to (0=inf)
	int detect_to;		// if 0, next 3 are not used
	double ti_to;		// time-out interval
	double ti_lastop;	// time of last successful operation
	double* timev;		// pointer to time, 0= use system call
	double getcurtime();
	DISABLECOPY(Buffer);
};

class RBuffer : public Buffer {
  public:
	RBuffer(int fd, int maxsize=0) : Buffer(fd,maxsize) { }
	~RBuffer() { }
	int fill();		// ret: 0=no, 1=yes, -1=other
	void extract(int num);	// have read num bytes
	int num() const;
	char operator[](int bi) const;
	int hasline() const;
	// next 2: die if len not sufficient, include '\n', ret success
	int examineline(char* s, int len);
	int extractline(char* s, int len); // examineline() + extract()
	char  bgetchar(int bi) const; // same as operator[]
	int   bgetint(int bi) const;
	short bgetshort(int bi) const;
	double bgetfloat(int bi) const;
};

class WBuffer : public Buffer {
  public:
	WBuffer(int fd, int maxsize=0) : Buffer(fd,maxsize) { }
	~WBuffer() { assertw(!n); }
	int flush(int nb=0);	// (nb==0 is all) ret: 0=part, 1=all, -1=other
	void bput(const void* buf, int nbytes);
	void bput(char c);
	void bput(short i);
	void bput(int i);
	void bput(double f);
};

//----------------------------------------------------------------------------

inline int Buffer::eof() const { return steof; }

inline int Buffer::err() const { return sterr; }

inline int Buffer::eot() const { return steot; }


inline int RBuffer::num() const { return n; }

inline char RBuffer::operator[](int bi) const
{
	assertx(bi>=0 && bi<n);
	return ar[beg+bi];
}

#endif
