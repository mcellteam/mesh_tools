/*
  subroutine sqrdc(x,ldx,n,p,qraux,jpvt,work,job)
  integer ldx,n,p,job
  integer jpvt(1)
  real x(ldx,1),qraux(1),work(1)
  */
/*
  subroutine sqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
  integer ldx,n,k,job,info
  real x(ldx,1),qraux(1),y(1),qy(1),qty(1),b(1),rsd(1),xb(1)
  */

void sqrdc_();
void sqrsl_();

main()
{
	static double a[3][2]={ {1,0}, {0,1}, {.8,.2} };
	static double rh[3]={3,4,2};
	static double fa[2][3], qraux[2], work[2], qty[2], b[2];
	static int jpvt[2]={0,0}, pivot=1, modeb=100, info;
	int m=3, n=2, i,j;
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
			fa[j][i]=a[i][j];
	sqrdc_(fa,&m,&m,&n,qraux,jpvt,work,&pivot);
	sqrsl_(fa,&m,&m,&n,qraux,rh,0,qty,b,0,0,&modeb,&info);
	printf("info=%d  (should be 0)\n",info);
	printf("b[0]=%g  (should be 2.42857)\n",b[0]);
	printf("b[1]=%g  (should be 3.85714)\n",b[1]);
}
