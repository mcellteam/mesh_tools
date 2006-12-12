
/* hh: I changed the #if structure */

typedef struct FCOMPLEX {double r,i;} fcomplex;
typedef struct IMMENSE {unsigned long l,r;} immense;
typedef struct GREAT {unsigned short l,c,r;} great;

#if defined(ANSI)
	void  adi(double **a, double **b, double **c, double **d, double **e,
		double **f, double **g, double **u, int jmax, int k,
		double alpha, double beta, double eps);
	void  amoeba(double **p, double *y, int ndim, double ftol,
		double (*funk)(double *), int *iter);
	void  anneal(double *x, double *y, int *iorder, int ncity);
	void  avevar(double *data, int n, double *ave, double *svar);
	void  balanc(double **a, int n);
	void  bcucof(double *y, double *y1, double *y2, double *y12, double d1,
		double d2, double **c);
	void  bcuint(double *y, double *y1, double *y2, double *y12, double x1l, 
		double x1u, double x2l, double x2u, double x1, double x2, 
		double *ansy, double *ansy1, double *ansy2);
	double bessi(int n, double x);
	double bessi0(double x);
	double bessi1(double x);
	double bessj(int n, double x);
	double bessj0(double x);
	double bessj1(double x);
	double bessk(int n, double x);
	double bessk0(double x);
	double bessk1(double x);
	double bessy(int n, double x);
	double bessy0(double x);
	double bessy1(double x);
	double beta(double z, double w);
	double betacf(double a, double b, double x);
	double betai(double a, double b, double x);
	double bico(int n, int k);
	void  bksub(int ne, int nb, int jf, int k1, int k2, double ***c);
	double bnldev(double pp, int n, int *idum);
	double brent(double ax, double bx, double cx, double (*f)(double), double tol,
		double *xmin);
	void  bsstep(double *y, double *dydx, int nv, double *xx, double htry,
		double eps, double *yscal, double *hdid, double *hnext, 
		void (*derivs)(double,double *,double *));
	void  caldat(long julian, int *mm, int *id, int *iyyy);
	double cel(double qqc, double pp, double aa, double bb);
	void  chder(double a, double b, double *c, double *cder, int n);
	double chebev(double a, double b, double *c, int m, double x);
	void  chebft(double a, double b, double *c, int n, double (*func)(double));
	void  chebpc(double *c, double *d, int n);
	void  chint(double a, double b, double *c, double *cint, int n);
	void  chsone(double *bins, double *ebins, int nbins, int knstrn, 
		double *df, double *chsq, double *prob);
	void  chstwo(double *bins1, double *bins2, int nbins, int knstrn, 
		double *df, double *chsq, double *prob);
	void  cntab1(int **nn, int n1, int nj, double *chisq, double *df, 
		double *prob, double *cramrv, double *ccc);
	void  cntab2(int **nn, int ni, int nj, double *h, double *hx, double *hy, 
		double *hygx, double *hxgy, double *uygx, double *uxgy,
		double *uxy);
	void  convlv(double *data, int n, double *respns, int m, int isign, 
		double *ans);
	void  correl(double *data1, double *data2, int n, double *ans);
	void  cosft(double *y, int n, int isign);
	void  covsrt(double **covar, int ma, int *lista, int mfit);
	void  crank(int n, double *w, double *s);
	double dbrent(double ax, double bx, double cx, double (*f)(double),
		double (*df)(double), double tol, double *xmin);
	void  ddpoly(double *c, int nc, double x, double *pd, int nd);
	void  des(immense inp, immense key, int *newkey, int isw, immense *out);
	void  ks(immense key, int n, great *kn);
	void  cyfun(unsigned long ir, great k, unsigned long *iout);
	double df1dim(double x);
	void  dfpmin(double *p, int n, double ftol, int *iter, double *fret, 
		double (*func)(double *), void (*dfunc)(double *,double *));
	void  difeq(int k, int k1, int k2, int jsf, int is1, int isf, 
		int *indexv, int ne, double **s, double **y);
	void  dlinmin(double *p, double *xi, int n, double *fret,
		double (*func)(double *), void (*dfunc)(double *,double *));
	void  eclass(int *nf, int n, int *lista, int *listb, int m);
	void  eclazz(int *nf, int n, int (*equiv)(int,int));
	void  eigsrt(double *d, double **v, int n);
	double el2(double x, double qqc, double aa, double bb);
	void  elmhes(double **a, int n);
	void  eulsum(double *sum, double term, int jterm, double *wksp);
	double evlmem(double fdt, double *cof, int m, double pm);
	double expdev(int *idum);
	double f1dim(double x);
	double factln(int n);
	double factrl(int n);
	void  fgauss(double x, double *a, double *y, double *dyda, int na);
	void  fit(double *x, double *y, int ndata, double *sig, int mwt, double *a, 
		double *b, double *siga, double *sigb, double *chi2, double *q);
	void  fixrts(double *d, int npoles);
	void  fleg(double x, double *pl, int nl);
	void  flmoon(int n, int nph, long *jd, double *frac);
	void  four1(double *data, int nn, int isign);
	void  fourn(double *data, int *nn, int ndim, int isign);
	void  fpoly(double x, double *p, int np);
	void  frprmn(double *p, int n, double ftol, int *iter, double *fret, 
		double (*func)(double *), void (*dfunc)(double *,double *));
	void  ftest(double *data1, int n1, double *data2, int n2, double *f, 
		double *prob);
	double gamdev(int ia, int *idum);
	double gammln(double xx);
	double gammp(double a, double x);
	double gammq(double a, double x);
	double gasdev(int *idum);
	void  gauleg(double x1, double x2, double *x, double *w, int n);
	void  gaussj(double **a, int n, double **b, int m);
	void  gcf(double *gammcf, double a, double x, double *gln);
	double golden(double ax, double bx, double cx, double (*f)(double), double tol, 
		double *xmin);
	void  gser(double *gamser, double a, double x, double *gln);
	void  hqr(double **a, int n, double *wr, double *wi);
	void  hunt(double *xx, int n, double x, int *jlo);
	void  indexx(int n, double *arrin, int *indx);
	int   irbit1(unsigned long int *iseed);
	int   irbit2(unsigned long int *iseed);
	void  jacobi(double **a, int n, double *d, double **v, int *nrot);
	long  julday(int mm, int id, int iyyy);
	void  kendl1(double *data1, double *data2, int n, double *tau, double *z,
		double *prob);
	void  kendl2(double **tab, int i, int j, double *tau, double *z,
		double *prob);
	void  ksone(double *data, int n, double (*func)(double), double *d,
		double *prob);
	void  kstwo(double *data1, int n1, double *data2, int n2, double *d,
		double *prob);
	void  laguer(fcomplex *a, int m, fcomplex *x, double eps, int polish);
	void  lfit(double *x, double *y, double *sig, int ndata, double *a, int ma, 
		int *lista, int mfit, double **covar, double *chisq,
		void (*funcs)(double,double *,int));
	void  linmin(double *p, double *xi, int n, double *fret, double (*func)(double));
	void  locate(double *xx, int n, double x, int *j);
	void  lubksb(double **a, int n, int *indx, double *b);
	void  ludcmp(double **a, int n, int *indx, double *d);
	void  mdian1(double *x, int n, double *xmed);
	void  mdian2(double *x, int n, double *xmed);
	void  medfit(double *x, double *y, int ndata, double *a, double *b,
		double *abdev);
	void  memcof(double *data, int n, int m, double *pm, double *cof);
	double midexp(double (*funk)(double), double aa, double bb, int n);
	double midinf(double (*funk)(double), double aa, double bb, int n);
	double midpnt(double (*func)(double), double a, double b, int n);
	double midsql(double (*funk)(double), double aa, double bb, int n);
	double midsqu(double (*funk)(double), double aa, double bb, int n);
	void  mmid(double *y, double *dydx, int nvar, double xs, double htot,
		int nstep, double *yout,
		void (*derivs)(double,double *,double *));
	void  mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, 
		double *fc, double (*func)(double));
	void  mnewt(int ntrial, double *x, int n, double tolx, double tolf);
	void  moment(double *data, int n, double *ave, double *adev, double *sdev, 
		double *svar, double *skew, double *curt);
	void  mprove(double **a, double **alud, int n, int *indx, double *b, 
		double *x);
	void  mrqcof(double *x, double *y, double *sig, int ndata, double *a, int ma, 
		int *lista, int mfit, double **alpha, double *beta, double
		*chisq, void (*funcs)(double,double *,double *,double *,int));
	void  mrqmin(double *x, double *y, double *sig, int ndata, double *a,
		int ma, int *lista, int mfit, double **covar, double **alpha, 
		double *chisq, void (*funcs)(double,double *,double *,double *,
		int),double *alamda);
	void  odeint(double *ystart, int nvar, double x1, double x2, double eps,
		double h1, double hmin, int *nok, int *nbad,
		void (*derivs)(double,double *,double *),
		void  (*rkqc)(double *,double *,int,double *,double,double,double
		*,double *,double *,void (*)(double,double *,double *)));
	void  pcshft(double a, double b, double *d, int n);
	void  pearsn(double *x, double *y, int n, double *r, double *prob, double *z);
	void  piksr2(int n, double *arr, double *brr);
	void  piksrt(int n, double *arr);
	void  pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
		double ***c, double **s);
	double plgndr(int l, int m, double x);
	double poidev(double xm, int *idum);
	void  polcoe(double *x, double *y, int n, double *cof);
	void  polcof(double *xa, double *ya, int n, double *cof);
	void  poldiv(double *u, int n, double *v, int nv, double *q, double *r);
	void  polin2(double *x1a, double *x2a, double **ya, int m, int n, double x1, 
		double x2, double *y, double *dy);
	void  polint(double *xa, double *ya, int n, double x, double *y, double *dy);
	void  powell(double *p, double **xi, int n, double ftol, int *iter,
		double *fret, double (*func)(double *));
	void  predic(double *data, int ndata, double *d, int npoles, 
		double *future, int nfut);
	double probks(double alam);
	void  pzextr(int iest, double xest, double *yest, double *yz, double *dy,
		int nv, int nuse);
	void  qcksrt(int n, double *arr);
	double qgaus(double (*func)(double), double a, double b);
	double qromb(double (*func)(double), double a, double b);
	double qromo(double (*func)(double), double a, double b,
		double (*choose)(double (*)(double),double,double,int));
	void  qroot(double *p, int n, double *b, double *c, double eps);
	double qsimp(double (*func)(double), double a, double b);
	double qtrap(double (*func)(double), double a, double b);
	double quad3d(double (*func)(double,double,double), double x1, double x2);
	double ran0(int *idum);
	double ran1(int *idum);
	double ran2(long *idum);
	double ran3(int *idum);
	double ran4(int *idum);
	void  rank(int n, int *indx, int *irank);
	void  ratint(double *xa, double *ya, int n, double x, double *y, double *dy);
	void  realft(double *data, int n, int isign);
	void  red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
		int ic1, int jc1, int jcf, int kc, double ***c, double **s);
	void  rk4(double *y, double *dydx, int n, double x, double h, double *yout,
		void (*derivs)(double,double *,double *));
	void  rkdumb(double *vstart, int nvar, double x1, double x2, int nstep, 
		void (*derivs)(double,double *,double *));
	void  rkqc(double *y, double *dydx, int n, double *x, double htry, 
		double eps, double *yscal, double *hdid, double *hnext, 
		void (*derivs)(double,double *,double *));
	double rofunc(double b);
	double rtbis(double (*func)(double), double x1, double x2, double xacc);
	double rtflsp(double (*func)(double), double x1, double x2, double xacc);
	double rtnewt(void (*funcd)(double,double *,double *), double x1, double x2,
		double xacc);
	double rtsafe(void (*funcd)(double,double *,double *), double x1, double x2,
		double xacc);
	double rtsec(double (*func)(double), double x1, double x2, double xacc);
	void  rzextr(int iest, double xest, double *yest, double *yz, double *dy,
		int nv, int nuse);
	void  scrsho(double (*fx)(double));
	void  shell(int n, double *arr);
	void  shoot(int nvar, double *v, double *delv, int n2, double x1, double x2,
		double eps, double h1, double hmin, double *f, double *dv);
	void  shootf(int nvar, double *v1, double *v2, double *delv1, double *delv2,
		int n1, int n2, double x1, double x2, double xf, double eps, 
		double h1, double hmin, double *f, double *dv1, double *dv2);
	void  simp1(double **a, int mm, int *ll, int nll, int iabf, int *kp,
		double *bmax);
	void  simp2(double **a, int n, int *l2, int nl2, int *ip, int kp,
		double *q1);
	void  simp3(double **a, int i1, int k1, int ip, int kp);
	void  simplx(double **a, int m, int n, int m1, int m2, int m3, 
		int *icase, int *izrov, int *iposv);
	void  sinft(double *y, int n);
	void  smooft(double *y, int n, double pts);
	void  sncndn(double uu, double emmc, double *sn, double *cn, double *dn);
	void  solvde(int itmax, double conv, double slowc, double *scalv,
		int *indexv, int ne, int nb, int m, double **y, double ***c,
		double **s);
	void  sor(double **a, double **b, double **c, double **d, double **e, 
		double **f, double **u, int jmax, double rjac);
	void  sort(int n, double *ra);
	void  sort2(int n, double *ra, double *rb);
	void  sort3(int n, double *ra, double *rb, double *rc);
	void  sparse(double *b, int n, double *x, double *rsq);
	void  spctrm(FILE *fp, double *p, int m, int k, int ovrlap);
	void  spear(double *data1, double *data2, int n, double *d, double *zd,
		double *probd, double *rs, double *probrs);
	void  splie2(double *x1a, double *x2a, double **ya, int m, int n,
		double **y2a);
	void  splin2(double *x1a, double *x2a, double **ya, double **y2a, int m,
		int n, double x1, double x2, double *y);
	void  spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
	void  splint(double *xa, double *ya, double *y2a, int n, double x, double *y);
	void  svbksb(double **u, double *w, double **v, int m, int n, double *b,
		double *x);
	void  svdcmp(double **a, int m, int n, double *w, double **v);
	void  svdfit(double *x, double *y, double *sig, int ndata, double *a, 
		int ma, double **u, double **v, double *w, double *chisq,
		void (*funcs)(double,double *,int));
	void  svdvar(double **v, int ma, double *w, double **cvm);
	void  toeplz(double *r, double *x, double *y, int n);
	void  tptest(double *data1, double *data2, int n, double *t, double *prob);
	void  tqli(double *d, double *e, int n, double **z);
	double trapzd(double (*func)(double), double a, double b, int n);
	void  tred2(double **a, int n, double *d, double *e);
	void  tridag(double *a, double *b, double *c, double *r, double *u, int n);
	void  ttest(double *data1, int n1, double *data2, int n2, double *t,
		double *prob);
	void  tutest(double *data1, int n1, double *data2, int n2, double *t,
		double *prob);
	void  twofft(double *data1, double *data2, double *fft1, double *fft2,
		int n);
	void  vander(double *x, double *w, double *q, int n);
	int   zbrac(double (*func)(double), double *x1, double *x2);
	void  zbrak(double (*fx)(double), double x1, double x2, int n, double *xb1,
		double *xb2, int *nb);
	double zbrent(double (*func)(double), double x1, double x2, double tol);
	void  zroots(fcomplex *a, int m, fcomplex *roots, int polish);
#elif defined(LINT_ARGS)
	void adi(double **, double **, double **, double **, 
		double **, double **, double **, double **, 
		int, int, double, double, double);
	void amoeba(double **, double *, int, double, double (*)(), int *);
	void anneal(double *, double *, int *, int);
	void avevar(double *, int, double *, double *);
	void balanc(double **, int);
	void bcucof(double *, double *, double *, double *, double, 
		double, double **);
	void bcuint(double *, double *, double *, double *, double, 
		double, double, double, double, double, double *, 
		double *, double *);
	double bessi(int, double);
	double bessi0(double);
	double bessi1(double);
	double bessj(int, double);
	double bessj0(double);
	double bessj1(double);
	double bessk(int, double);
	double bessk0(double);
	double bessk1(double);
	double bessy(int, double);
	double bessy0(double);
	double bessy1(double);
	double beta(double, double);
	double betacf(double, double, double);
	double betai(double, double, double);
	double bico(int, int);
	void bksub(int, int, int, int, int, double ***);
	double bnldev(double, int, int *);
	double brent(double, double, double, double (*)(), double, double *);
	void bsstep(double *, double *, int, double *, double,
		double, double *, double *, double *, void (*)());
	void caldat(long, int *, int *, int *);
	double cel(double, double, double, double);
	void chder(double, double, double *, double *, int);
	double chebev(double, double, double *, int, double);
	void chebft(double, double, double *, int, double (*)());
	void chebpc(double *, double *, int);
	void chint(double, double, double *, double *, int);
	void chsone(double *, double *, int, int, double *, 
		double *, double *);
	void chstwo(double *, double *, int, int, double *, 
		double *, double *);
	void cntab1(int **, int, int, double *, double *, 
		double *, double *, double *);
	void cntab2(int **, int, int, double *, double *, 
		double *, double *, double *, double *, double *, double *);
	void convlv(double *, int, double *, int, int, double *);
	void correl(double *, double *, int, double *);
	void cosft(double *, int, int);
	void covsrt(double **, int, int *, int);
	void crank(int, double *, double *);
	double dbrent(double, double, double, double (*)(), double (*)(),
		double, double *);
	void ddpoly(double *, int, double, double *, int);
	void des(immense, immense, int *, int, immense *);
	void ks(immense, int, great *);
	void cyfun(unsigned long, great, unsigned long *);
	double df1dim(double);
	void dfpmin(double *, int, double, int *, double *, double (*)(),
		void (*)());
	void difeq(int, int, int, int, int, int, int *, int, 
		double **, double **);
	void dlinmin(double *, double *, int, double *, double (*)(), void (*)());
	void eclass(int *, int, int *, int *, int);
	void eclazz(int *, int, int (*)());
	void eigsrt(double *, double **, int);
	double el2(double, double, double, double);
	void elmhes(double **, int);
	void eulsum(double *, double, int, double *);
	double evlmem(double, double *, int, double);
	double expdev(int *);
	double f1dim(double);
	double factln(int);
	double factrl(int);
	void fgauss(double, double *, double *, double *, int);
	void fit(double *, double *, int, double *, int, double *, 
		double *, double *, double *, double *, double *);
	void fixrts(double *, int);
	void fleg(double, double *, int);
	void flmoon(int, int, long *, double *);
	void four1(double *, int, int);
	void fourn(double *, int *, int, int);
	void fpoly(double, double *, int);
	void frprmn(double *, int, double, int *, double *, double (*)(),
		void (*)());
	void ftest(double *, int, double *, int, double *, double *);
	double gamdev(int, int *);
	double gammln(double);
	double gammp(double, double);
	double gammq(double, double);
	double gasdev(int *);
	void gauleg(double, double, double *, double *, int);
	void gaussj(double **, int, double **, int);
	void gcf(double *, double, double, double *);
	double golden(double, double, double, double (*)(), double, double *);
	void gser(double *, double, double, double *);
	void hqr(double **, int, double *, double *);
	void hunt(double *, int, double, int *);
	void indexx(int, double *, int *);
	int irbit1(unsigned long int *);
	int irbit2(unsigned long int *);
	void jacobi(double **, int, double *, double **, int *);
	long julday(int, int, int);
	void kendl1(double *, double *, int, double *, double *, double *);
	void kendl2(double **, int, int, double *, double *, double *);
	void ksone(double *, int, double (*)(), double *, double *);
	void kstwo(double *, int, double *, int, double *, double *);
	void laguer(fcomplex *, int, fcomplex *, double, int);
	void lfit(double *, double *, double *, int, double *, int, 
		int *, int, double **, double *, void (*)());
	void linmin(double *, double *, int, double *, double (*)());
	void locate(double *, int, double, int *);
	void lubksb(double **, int, int *, double *);
	void ludcmp(double **, int, int *, double *);
	void mdian1(double *, int, double *);
	void mdian2(double *, int, double *);
	void medfit(double *, double *, int, double *, double *, double *);
	void memcof(double *, int, int, double *, double *);
	double midexp(double (*)(), double, double, int);
	double midinf(double (*)(), double, double, int);
	double midpnt(double (*)(), double, double, int);
	double midsql(double (*)(), double, double, int);
	double midsqu(double (*)(), double, double, int);
	void mmid(double *, double *, int, double, double, int, double *,
		void (*)());
	void mnbrak(double *, double *, double *, double *, 
		double *, double *, double (*)());
	void mnewt(int, double *, int, double, double);
	void moment(double *, int, double *, double *, double *, 
		double *, double *, double *);
	void mprove(double **, double **, int, int *, double *, double *);
	void mrqcof(double *, double *, double *, int, double *, int, 
		int *, int, double **, double *, double *, void (*)());
	void mrqmin(double *, double *, double *, int, double *, int, 
		int *, int, double **, double **, double *, void (*)(),
		double *);
	void odeint(double *, int, double, double, double, double, 
		double, int *, int *, void (*)(), void (*)());
	void pcshft(double, double, double *, int);
	void pearsn(double *, double *, int, double *, double *, double *);
	void piksr2(int, double *, double *);
	void piksrt(int, double *);
	void pinvs(int, int, int, int, int, int, double ***, double **);
	double plgndr(int, int, double);
	double poidev(double, int *);
	void polcoe(double *, double *, int, double *);
	void polcof(double *, double *, int, double *);
	void poldiv(double *, int, double *, int, double *, double *);
	void polin2(double *, double *, double **, int, int, double, 
		double, double *, double *);
	void polint(double *, double *, int, double, double *, double *);
	void powell(double *, double **, int, double, int *, double *, double (*)());
	void predic(double *, int, double *, int, double *, int);
	double probks(double);
	void pzextr(int, double, double *, double *, double *, int, int);
	void qcksrt(int, double *);
	double qgaus(double (*)(double), double, double);
	double qromb(double (*)(double), double, double);
	double qromo(double (*)(double), double, double, double (*)());
	void qroot(double *, int, double *, double *, double);
	double qsimp(double (*)(double), double, double);
	double qtrap(double (*)(double), double, double);
	double quad3d(double (*)(double), double, double);
	double ran0(int *);
	double ran1(int *);
	double ran2(long *);
	double ran3(int *);
	double ran4(int *);
	void rank(int, int *, int *);
	void ratint(double *, double *, int, double, double *, double *);
	void realft(double *, int, int);
	void red(int, int, int, int, int, int, int, int, int, int, 
		int, double ***, double **);
	void rk4(double *, double *, int, double, double, double *, void (*)());
	void rkdumb(double *, int, double, double, int, void (*)());
	void rkqc(double *, double *, int, double *, double, double, 
		double *, double *, double *, void (*)());
	double rofunc(double);
	double rtbis(double (*)(), double, double, double);
	double rtflsp(double (*)(), double, double, double);
	double rtnewt(void (*)(), double, double, double);
	double rtsafe(void (*)(), double, double, double);
	double rtsec(double (*)(), double, double, double);
	void rzextr(int, double, double *, double *, double *, int, int);
	void scrsho(double (*)());
	void shell(int, double *);
	void shoot(int, double *, double *, int, double, double, double, 
		double, double, double *, double *);
	void shootf(int, double *, double *, double *, double *, int, int, 
		double, double, double, double, double, double, double *, 
		double *, double *);
	void simp1(double **, int, int *, int, int, int *, double *);
	void simp2(double **, int, int *, int, int *, int, double *);
	void simp3(double **, int, int, int, int);
	void simplx(double **, int, int, int, int, int, 
		int *, int *, int *);
	void sinft(double *, int);
	void smooft(double *, int, double);
	void sncndn(double, double, double *, double *, double *);
	void solvde(int, double, double, double *, int *, int, int, 
		int, double **, double ***, double **);
	void sor(double **, double **, double **, double **, double **, 
		double **, double **, int, double);
	void sort(int, double *);
	void sort2(int, double *, double *);
	void sort3(int, double *, double *, double *);
	void sparse(double *, int, double *, double *);
	void spctrm(FILE *, double *, int, int, int);
	void spear(double *, double *, int, double *, 
		double *, double *, double *, double *);
	void splie2(double *, double *, double **, int, int, double **);
	void splin2(double *, double *, double **, double **, int, int, 
		double, double, double *);
	void spline(double *, double *, int, double, double, double *);
	void splint(double *, double *, double *, int, double, double *);
	void svbksb(double **, double *, double **, int, int,
		double *, double *);
	void svdcmp(double **, int, int, double *, double **);
	void svdfit(double *, double *, double *, int, double *, int, 
		double **, double **, double *, double *, void (*)());
	void svdvar(double **, int, double *, double **);
	void toeplz(double *, double *, double *, int);
	void tptest(double *, double *, int, double *, double *);
	void tqli(double *, double *, int, double **);
	double trapzd(double (*)(double), double, double, int);
	void tred2(double **, int, double *, double *);
	void tridag(double *, double *, double *, double *, double *, int);
	void ttest(double *, int, double *, int, double *, double *);
	void tutest(double *, int, double *, int, double *, double *);
	void twofft(double *, double *, double *, double *, int);
	void vander(double *, double *, double *, int);
	int zbrac(double (*)(), double *, double *);
	void zbrak(double (*)(), double, double, int, double *, double *, int *);
	double zbrent(double (*)(), double, double, double);
	void zroots(fcomplex *, int, fcomplex *, int);
#else
	void adi();
	void amoeba();
	void anneal();
	void avevar();
	void balanc();
	void bcucof();
	void bcuint();
	double bessi();
	double bessi0();
	double bessi1();
	double bessj();
	double bessj0();
	double bessj1();
	double bessk();
	double bessk0();
	double bessk1();
	double bessy();
	double bessy0();
	double bessy1();
	double beta();
	double betacf();
	double betai();
	double bico();
	void bksub();
	double bnldev();
	double brent();
	void bsstep();
	void caldat();
	double cel();
	void chder();
	double chebev();
	void chebft();
	void chebpc();
	void chint();
	void chsone();
	void chstwo();
	void cntab1();
	void cntab2();
	void convlv();
	void correl();
	void cosft();
	void covsrt();
	void crank();
	double dbrent();
	void ddpoly();
	void des();
	void ks();
	void cyfun();
	double df1dim();
	void dfpmin();
	void dlinmin();
	void difeq();
	void eclass();
	void eclazz();
	void eigsrt();
	double el2();
	void elmhes();
	void eulsum();
	double evlmem();
	double expdev();
	double f1dim();
	double factln();
	double factrl();
	void fgauss();
	void fit();
	void fixrts();
	void fleg();
	void flmoon();
	void four1();
	void fourn();
	void fpoly();
	void frprmn();
	void ftest();
	double gamdev();
	double gammln();
	double gammp();
	double gammq();
	double gasdev();
	void gauleg();
	void gaussj();
	void gcf();
	double golden();
	void gser();
	void hqr();
	void hunt();
	void indexx();
	int irbit1();
	int irbit2();
	void jacobi();
	long julday();
	void kendl1();
	void kendl2();
	void ksone();
	void kstwo();
	void laguer();
	void lfit();
	void linmin();
	void locate();
	void lubksb();
	void ludcmp();
	void mdian1();
	void mdian2();
	void medfit();
	void memcof();
	double midexp();
	double midinf();
	double midpnt();
	double midsql();
	double midsqu();
	void mmid();
	void mnbrak();
	void mnewt();
	void moment();
	void mprove();
	void mrqcof();
	void mrqmin();
	void odeint();
	void pcshft();
	void pearsn();
	void piksr2();
	void piksrt();
	void pinvs();
	double plgndr();
	double poidev();
	void polcoe();
	void polcof();
	void poldiv();
	void polin2();
	void polint();
	void powell();
	void predic();
	double probks();
	void pzextr();
	void qcksrt();
	double qgaus();
	double qromb();
	double qromo();
	void qroot();
	double qsimp();
	double qtrap();
	double quad3d();
	double ran0();
	double ran1();
	double ran2();
	double ran3();
	double ran4();
	void rank();
	void ratint();
	void realft();
	void red();
	void rk4();
	void rkdumb();
	void rkqc();
	double rofunc();
	double rtbis();
	double rtflsp();
	double rtnewt();
	double rtsafe();
	double rtsec();
	void rzextr();
	void scrsho();
	void shell();
	void shoot();
	void shootf();
	void simp1();
	void simp2();
	void simp3();
	void simplx();
	void sinft();
	void smooft();
	void sncndn();
	void solvde();
	void sor();
	void sort();
	void sort2();
	void sort3();
	void sparse();
	void spctrm();
	void spear();
	void splie2();
	void splin2();
	void spline();
	void splint();
	void svbksb();
	void svdcmp();
	void svdfit();
	void svdvar();
	void toeplz();
	void tptest();
	void tqli();
	double trapzd();
	void tred2();
	void tridag();
	void ttest();
	void tutest();
	void twofft();
	void vander();
	int zbrac();
	void zbrak();
	double zbrent();
	void zroots();
#endif
