int distinguishable(double a,double b,double eps)
{
  double c;
  c=a-b;
  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;
  if (b<a) eps*=a;
  else eps*=b;
  return (c>eps);
}

Weights* loadWeightArrays(void_list *p,Weights* w,int limit) {
	///// load arrays of weights /////
	// END OF RAW POINTS ARRAY WAS PASSED AS P 
	// AND LIST IS TRAVERSED IN 'PREVIOUS' DIRECTION
	int i=0,j,k;
	void_list *qq;
	// spline 0
	for (qq=p;qq!=NULL;qq=qq->previous) {
		if (qq->previous==NULL) {
			w->bx[i] = ((Point*)qq->data)->x;
			w->by[i] = ((Point*)qq->data)->y;
			i++;
			break;
		}
	}
	qq=p;
	w->bx[i] = ((Point*)qq->data)->x;
	w->by[i] = ((Point*)qq->data)->y;
	i++;
	qq=qq->previous;
	w->bx[i] = ((Point*)qq->data)->x;
	w->by[i] = ((Point*)qq->data)->y;
	i++;
	qq=qq->previous;
	w->bx[i] = ((Point*)qq->data)->x;
	w->by[i] = ((Point*)qq->data)->y;
	i++;
	// splines 1 through m-2
	for (j=1;j<(limit-2);j++) { // note '-2' because m = contour_array[i]-1
		k=0;
		for (qq=p;qq!=NULL;qq=qq->previous) {
			if( (k==(j-1)) || (k==j) || (k==(j+1)) || (k==(j+2)) ) {
				w->bx[i] = ((Point*)qq->data)->x;
				w->by[i] = ((Point*)qq->data)->y;
				i++;
			}
			k++;
			if (qq->previous==NULL) break;
		}
	}
	// spline m-1
	k=0;
	for (qq=p;qq!=NULL;qq=qq->previous) {
		if (k==(limit-3)) { break; }
		k++;
	}
	w->bx[i] = ((Point*)qq->data)->x;
	w->by[i] = ((Point*)qq->data)->y;
	i++;
	qq=qq->previous;
	w->bx[i] = ((Point*)qq->data)->x;
	w->by[i] = ((Point*)qq->data)->y;
	i++;
	qq=qq->previous;
	w->bx[i] = ((Point*)qq->data)->x;
	w->by[i] = ((Point*)qq->data)->y;
	i++;
	qq=p;
	w->bx[i] = ((Point*)qq->data)->x;
	w->by[i] = ((Point*)qq->data)->y;
	i++;
	// spline m
	k=0;
	for (qq=p;qq!=NULL;qq=qq->previous) {
		if (k==(limit-2)) { break; }
		k++;
	}
	w->bx[i] = ((Point*)qq->data)->x;
	w->by[i] = ((Point*)qq->data)->y;
	i++;
	qq=qq->previous;
	w->bx[i] = ((Point*)qq->data)->x;
	w->by[i] = ((Point*)qq->data)->y;
	i++;
	qq=p;
	w->bx[i] = ((Point*)qq->data)->x;
	w->by[i] = ((Point*)qq->data)->y;
	i++;
	qq=qq->previous;
	w->bx[i] = ((Point*)qq->data)->x;
	w->by[i] = ((Point*)qq->data)->y;

	return w;
}


SplinePoint* computeSplines(SplinePoint* sp,Weights *w,double dt,int limit,Parameters *pa) {
	int i,j,index;
	double inc,inc2,inc3,xdot,xdotdot,ydot,ydotdot,den,num_part,maxr=pa->max_rad;
	// for each point in contour
	for (i=0;i<limit;i++) {
		// for each parameter increment
		for(j=0;j<pa->num;j++) {
			index = i*pa->num+j;
			inc=j*dt;
			inc2 = inc*inc;
			inc3 = inc2*inc;
			// store time
			sp[index].t   = i+inc;
			sp[index].x   = 1.0/6.0*(w->bx[4*i+0]*(-inc3+3.0*inc2-3.0*inc+1.0) 
					+ w->bx[4*i+1]*(3.0*inc3-6.0*inc2+4.0) 
					+ w->bx[4*i+2]*(-3.0*inc3+3.0*inc2+3.0*inc+1.0) 
					+ w->bx[4*i+3]*inc3);
			sp[index].y   = 1.0/6.0*(w->by[4*i+0]*(-inc3+3.0*inc2-3.0*inc+1.0) 
					+ w->by[4*i+1]*(3.0*inc3-6.0*inc2+4.0) 
					+ w->by[4*i+2]*(-3.0*inc3+3.0*inc2+3.0*inc+1.0) 
					+ w->by[4*i+3]*inc3);
			xdot    = 1.0/6.0*(w->bx[4*i+0]*(-3.0*inc2+6.0*inc-3.0)
					+ w->bx[4*i+1]*(9.0*inc2-12.0*inc)
					+ w->bx[4*i+2]*(-9.0*inc2+6.0*inc+3.0)
					+ w->bx[4*i+3]*3.0*inc2);
			ydot    = 1.0/6.0*(w->by[4*i+0]*(-3.0*inc2+6.0*inc-3.0)
					+ w->by[4*i+1]*(9.0*inc2-12.0*inc)
					+ w->by[4*i+2]*(-9.0*inc2+6.0*inc+3.0)
					+ w->by[4*i+3]*3.0*inc2);
			xdotdot = 1.0/6.0*(w->bx[4*i+0]*(-6.0*inc+6.0)
					+ w->bx[4*i+1]*(18.0*inc-12.0)
					+ w->bx[4*i+2]*(-18.0*inc+6.0)
					+ w->bx[4*i+3]*6.0*inc);
			ydotdot = 1.0/6.0*(w->by[4*i+0]*(-6.0*inc+6.0)
					+ w->by[4*i+1]*(18.0*inc-12.0)
					+ w->by[4*i+2]*(-18.0*inc+6.0)
					+ w->by[4*i+3]*6.0*inc);
			den     = fabs(xdot*ydotdot-ydot*xdotdot);
			num_part = sqrt(xdot*xdot+ydot*ydot);
			sp[index].intfac = num_part;
			if (den) { sp[index].r = num_part*num_part*num_part/den; }
			else { sp[index].r = maxr;}
			if (sp[index].r>maxr){sp[index].r=maxr;}
		}
	}
	return sp;
}
/*
void sampleSplinesOld(void_list* q,SplinePoint* sp,double dt,int limit, double thickness,
					char *outdir,int tag,Parameters *pa) {
	///// sample splines /////
	Contour *c;
	c=(Contour*)q->data;
	void_list *pp;
	Point *v;
	int i,j,k=0,count=0,num_sampled=0,myswitch,kvec[limit*pa->num/2];
	// dmin = minimum spline sampling distance
	// dmax = maximum spline sampling distance
	// decrease tau to increase steepness
	// inflection point of sampling function = tau*rI
	// decrease rI to sample finer
//	double dmin=.001,dmax=.050,tau=2,rI=-6,inc;
	double dl,l_accum=0,l=0,r_mean,ds,ds_mean=0,intfac_array[3],r_array[2],xval,yval,inc;
	double rvec[limit*pa->num/2];
	for (i=0;i<limit;i++) {
		for(j=0;j<pa->num;j++) {
			inc=(double)(i+j*dt);
			myswitch = j%2;
			if(inc && !myswitch) {
				intfac_array[0] = sp[(int)((inc-2.0*dt)/dt)].intfac;
				r_array[0] = sp[(int)((inc-2.0*dt)/dt)].r;
				intfac_array[1] = sp[(int)((inc-dt)/dt)].intfac;
				intfac_array[2] = sp[(int)(inc/dt)].intfac;
				r_array[1] = sp[(int)(inc/dt)].r;
				xval = sp[(int)(inc/dt)].x;
				yval = sp[(int)(inc/dt)].y;
				// increment along spline length
				dl = dt/6.0*(intfac_array[0]+4.0*intfac_array[1]+intfac_array[2]);
				// accumulated contour length traversal
				l += dl;
				// accumulated spline length traversal
				l_accum += dl;
				// mean radius of curvature
				r_mean = (r_array[0]+r_array[1])/2.0;
				// incremental sampling distance
				ds = (pa->dmax-pa->dmin)/(1.0+exp(-r_mean/pa->tau-pa->rI))+pa->dmin;
				count++;
				// mean sampling distance
				ds_mean = (ds_mean*(double)(count-1.0)+ds)/count;
				// sample data
				if (l_accum >= ds_mean) {
					// add interpolated point to contour
					pp = new void_list();
					pp->next = c->interp_points;
					v = new Point(xval,yval,c->section*thickness);
					pp->data = (void*)v;
					c->interp_points = pp;
					// clear variables
					l_accum = 0;
					ds_mean = 0;
					count = 0;
					num_sampled++;
				}
				if(pa->diag){rvec[k]=r_mean;kvec[k]=k;k++;}
			}
		}
	}

	// store number of sampled points
	c->num_interp_points = num_sampled;

	// add previous data
	c->addPreviousInterp();

	// diagnostics
	if(pa->diag) {
		char filename[256],line[2048];
		FILE *F;
		double r_test[201],d_r[201];
		r_test[0] = 0.1;
		d_r[0] = (pa->dmax-pa->dmin)/(1.0+exp(-r_test[0]/pa->tau-pa->rI))+pa->dmin;
		for(i=1;i<201;i++){
			r_test[i]=i*pa->plot_rad_int;
			d_r[i] = (pa->dmax-pa->dmin)/(1.0+exp(-r_test[i]/pa->tau-pa->rI))+pa->dmin;
		}
		///// print sampling function points /////
		sprintf(filename,"%ssampling_function_%i.dat",outdir,tag);
		F = fopen(filename,"w");
		if (!F) { printf("Couldn't open output file %s\n",filename); return; }
		// for each point
		for(i=0;i<201;i++){
			sprintf(line,"%.15g %.15g\n",r_test[i],d_r[i]);
			fputs(line,F);
		}
		fclose(F);
		///// print radius of curvature /////
		sprintf(filename,"%s%s_%i_%i.rad",outdir,c->name,c->section,tag);
		F = fopen(filename,"w");
		if (!F) { printf("Couldn't open output file %s\n",filename); return; }
		// for each point
		for(i=0;i<limit*pa->num/2-1;i++){
			sprintf(line,"%i %.15g\n",kvec[i],rvec[i]);
			fputs(line,F);
		}
		fclose(F);
	}
}
*/

void sampleSplines(void_list* q,SplinePoint* sp,double dt,int limit, double thickness,
					char *outdir,int tag,Parameters *pa) {
	///// sample splines /////
	Contour *c;
	c=(Contour*)q->data;
	void_list *pp;
	Point *v;
	int i,j,k=0,count=0,num_sampled=0,myswitch,kvec[limit*pa->num/2];
	// dmin = minimum spline sampling distance
	// dmax = maximum spline sampling distance
	// decrease tau to increase steepness
	// inflection point of sampling function = tau*rI
	// decrease rI to sample finer
//	double dmin=.001,dmax=.050,tau=2,rI=-6,inc;
	double dl,l_accum=0.0,l=0.0,r_mean,ds,ds_mean=0.0,intfac_array[3],r_array[2],xval,yval,inc;
	double vmax=pa->dmax/pa->T,v_mean,amax=pa->amax,delt,t_accum=0.0,T=pa->T;
	double rvec[limit*pa->num/2],vvec[limit*pa->num/2],tvec[limit*pa->num/2];
	for (i=0;i<limit;i++) {
		for(j=0;j<pa->num;j++) {
			inc=(double)(i+j*dt);
			myswitch = j%2;
			if(inc && !myswitch) {
				intfac_array[0] = sp[(int)((inc-2.0*dt)/dt)].intfac;
				r_array[0] = sp[(int)((inc-2.0*dt)/dt)].r;
				intfac_array[1] = sp[(int)((inc-dt)/dt)].intfac;
				intfac_array[2] = sp[(int)(inc/dt)].intfac;
				r_array[1] = sp[(int)(inc/dt)].r;
				xval = sp[(int)(inc/dt)].x;
				yval = sp[(int)(inc/dt)].y;
				// increment along spline length
				dl = 2*dt/3.0*(intfac_array[0]+4.0*intfac_array[1]+intfac_array[2]);
				// mean radius of curvature
				r_mean = (r_array[0]+r_array[1])/2.0;
				// mean velocity
				v_mean = sqrt(amax*r_mean);
				if (v_mean>vmax) {v_mean=vmax;}
				// time increment
				delt = dl/v_mean;
				// accumulated time
				t_accum += delt;
				// sample data
				if (t_accum >= T) {
					// add interpolated point to contour
					pp = new void_list();
					pp->next = c->interp_points;
					v = new Point(xval,yval,c->section*thickness);
					pp->data = (void*)v;
					c->interp_points = pp;
					// clear variables
					t_accum = 0.0;
					num_sampled++;
				}
				if(pa->diag){
					rvec[k]=r_mean;
					vvec[k]=v_mean;
					tvec[k]=delt;
					kvec[k]=k;k++;
				}
			}
		}
	}

	// store number of sampled points
	c->num_interp_points = num_sampled;

	// add previous data
	c->addPreviousInterp();

	// diagnostics
	if(pa->diag) {
		char filename[256],line[2048];
		FILE *F;
		///// print radius of curvature /////
		sprintf(filename,"%s%s_%i_%i.rad",outdir,c->name,c->section,tag);
		F = fopen(filename,"w");
		if (!F) { printf("Couldn't open output file %s\n",filename); return; }
		// for each point
		for(i=0;i<limit*pa->num/2-1;i++){
			sprintf(line,"%i %.15g\n",kvec[i],rvec[i]);
			fputs(line,F);
		}
		fclose(F);
		///// print velocity /////
		sprintf(filename,"%s%s_%i_%i.vel",outdir,c->name,c->section,tag);
		F = fopen(filename,"w");
		if (!F) { printf("Couldn't open output file %s\n",filename); return; }
		// for each point
		for(i=0;i<limit*pa->num/2-1;i++){
			sprintf(line,"%i %.15g\n",kvec[i],vvec[i]);
			fputs(line,F);
		}
		fclose(F);
		///// print incremental time /////
		sprintf(filename,"%s%s_%i_%i.tim",outdir,c->name,c->section,tag);
		F = fopen(filename,"w");
		if (!F) { printf("Couldn't open output file %s\n",filename); return; }
		// for each point
		for(i=0;i<limit*pa->num/2-1;i++){
			sprintf(line,"%i %.15g\n",kvec[i],tvec[i]);
			fputs(line,F);
		}
		fclose(F);
	}
}

void computeDeviation(void_list* q,SplinePoint* sp,int count,int num) {
	///// compute deviation distance between raw points and splines /////
	Point *P;
	Contour *c;
	void_list *p,*pp;
	c=(Contour*)q->data;
	c->deviations = new double[count];
	double diffx,diffy,dist0,distneg,distneg0,distpos,distpos0;
	int j=0,i=0,m,n;
	for (p=c->rawend;p!=NULL;p=p->previous) {
		P=(Point*)p->data;
		diffx = P->x-sp[j].x;
		diffy = P->y-sp[j].y;
		dist0 = sqrt(diffx*diffx+diffy*diffy);
		// check negative dir
		m=j;
		distneg=dist0;
		do {
			distneg0=distneg;
			if(!m){m=num*count-1;}
			else {m--;}
			diffx = P->x-sp[m].x;
			diffy = P->y-sp[m].y;
			distneg = sqrt(diffx*diffx+diffy*diffy);
		} while (distneg<distneg0);
		// check positive dir
		n=j;
		distpos=dist0;
		do {
			distpos0=distpos;
			if(n==num*count){n=0;}
			else {n++;}
			diffx = P->x-sp[n].x;
			diffy = P->y-sp[n].y;
			distpos = sqrt(diffx*diffx+diffy*diffy);
		} while (distpos<distpos0);

		if (dist0<distneg && dist0<distpos){
			c->deviations[i]=dist0;
		} else if (distneg0<dist0 && distneg0<distpos0) {
			c->deviations[i]=distneg0;
		} else if (distpos0<dist0 && distpos0<distneg0) {
			c->deviations[i]=distpos0;
		} else { 
			printf("\n\nweird error!\n");
			printf("Contour %s, section %d, #raw points %d, num %d\n",c->name,c->section,count,num);
			printf("rawx %.15g, rawy %.15g\n",P->x,P->y);
			printf("j %i, dist0 %.15g, dist0x %.15g, dist0y %.15g\n",j,dist0,sp[j].x,sp[j].y);
			printf("m %i, distneg0 %.15g, distneg %.15g," ,m,distneg0,distneg);
			printf(" distnegx %.15g, distnegy %.15g\n" ,sp[m].x,sp[m].y);
			printf("n %i, distpos0 %.15g, distpos %.15g," ,n,distpos0,distpos);
			printf(" distposx %.15g, distposy %.15g\n" ,sp[n].x,sp[n].y);
		}
		j+=num;
		i++;
	}
}

Weights* checkDeviation(void_list *q,void_list *p,Weights* w,int limit,double maxdev,double scale) {
	///// check deviation at each point in contour /////
	// END OF RAW POINTS ARRAY WAS PASSED AS P 
	// AND LIST IS TRAVERSED IN 'PREVIOUS' DIRECTION
	void_list *pp;
	int i;
	bool flag=false;
	// for each deviation
	for (i=0;i<limit;i++){
		// if deviation exceeds threshold
		if (((Contour*)q->data)->deviations[i]*scale>maxdev) {
			flag = true;
			// edit weight array
			if(!i) {
				w->bx[4*(limit-1)+3] = w->bx[4*(limit-1)+2];
				w->by[4*(limit-1)+3] = w->by[4*(limit-1)+2];
				w->bx[4*i+0] = w->bx[4*i+1];
				w->by[4*i+0] = w->by[4*i+1];
			} else {
				w->bx[4*(i-1)+3] = w->bx[4*(i-1)+2];
				w->by[4*(i-1)+3] = w->by[4*(i-1)+2];
				w->bx[4*i+0] = w->bx[4*i+1];
				w->by[4*i+0] = w->by[4*i+1];
			}
		}
	}
	if (flag) {	return w;}
	else {return NULL;}
}

void printDiagnostics(void_list* q,SplinePoint* sp,int num,char *outdir,int tag) {
	Contour *c;
	c=(Contour*)q->data;
	void_list *qq;
	int j;
	char filename[256],line[2048];
	FILE *F;
	///// print raw points /////
	sprintf(filename,"%s%s_%i_%i.raw",outdir,c->name,c->section,tag);
	F = fopen(filename,"w");
	if (!F) { printf("Couldn't open output file %s\n",filename); return; }
	// for each point
	for (qq=c->rawend;qq!=NULL;qq=qq->previous) {
		sprintf(line,"%.15g %.15g\n",((Point*)qq->data)->x,((Point*)qq->data)->y);
		fputs(line,F);
	}
	fclose(F);
	///// print spline points /////
	sprintf(filename,"%s%s_%i_%i.spline",outdir,c->name,c->section,tag);
	F = fopen(filename,"w");
	if (!F) { printf("Couldn't open output file %s\n",filename); return; }
	// for each spline point
	for(j=0;j<num;j++){
		sprintf(line,"%.15g %.15g\n",sp[j].x,sp[j].y);
		fputs(line,F);
	}
	fclose(F);
	///// print interpolated points /////
	sprintf(filename,"%s%s_%i_%i.interp",outdir,c->name,c->section,tag);
	F = fopen(filename,"w");
	if (!F) { printf("Couldn't open output file %s\n",filename); return; }
	// for each point
	for (qq=c->interpend;qq!=NULL;qq=qq->previous) {
		sprintf(line,"%.15g %.15g\n",((Point*)qq->data)->x,((Point*)qq->data)->y);
		fputs(line,F);
	}
	fclose(F);
	///// print grace script /////
	sprintf(filename,"%sbfile_%s_%i_%i",outdir,c->name,c->section,tag);
	F = fopen(filename,"w");
	if (!F) { printf("Couldn't open output file %s\n",filename); return; }
	sprintf(line,"#Obligatory descriptive comment\n");
	fputs(line,F);
	sprintf(line,"READ xy \"%s_%i_%i.raw\"\n",c->name,c->section,tag);
	fputs(line,F);
	sprintf(line,"READ xy \"%s_%i_%i.spline\"\n",c->name,c->section,tag);
	fputs(line,F);
	sprintf(line,"READ xy \"%s_%i_%i.interp\"\n",c->name,c->section,tag);
	fputs(line,F);
	sprintf(line,"legend on\ns0 legend ");
	sprintf(line,"%s\"Raw\"\ns1 legend \"Spline\"\ns2 legend \"Interpolated\"\n",line);
	fputs(line,F);
	sprintf(line,"s0 symbol 1\ns0 symbol color 2\ns0 symbol fill color 2\n");
	sprintf(line,"%ss0 line type 0\ns0 line color 2\ns0 errorbar color 2\n",line);
	fputs(line,F);
	sprintf(line,"s1 symbol color 1\ns1 symbol fill color 1\n");
	sprintf(line,"%ss1 line color 1\ns1 errorbar color 1\n",line);
	fputs(line,F);
	sprintf(line,"s2 symbol 2\ns2 symbol color 4\ns2 symbol fill color 4\n");
	sprintf(line,"%ss2 line type 0\ns2 line color 4\ns2 errorbar color 4\n",line);
	fputs(line,F);
	fclose(F);
	///// print shell script /////
	sprintf(filename,"%s%s_%i_%i.csh",outdir,c->name,c->section,tag);
	F = fopen(filename,"w");
	if (!F) { printf("Couldn't open output file %s\n",filename); return; }
	sprintf(line,"#!/bin/csh\n");
	fputs(line,F);
	sprintf(line,"xmgrace -noask -nosafe -batch bfile_%s_%i_%i\n",c->name,c->section,tag);
	fputs(line,F);
	fclose(F);
}

void_list * getContours(int argc,char *argv[],double thickness){
	int i,min_section,max_section;
	char *indir,infile[128],*str,line[2048],*name,filename[256],*eptr,*temp;
	FILE *F;
	void_list *c,*q,*ch;
	ch=NULL;
	Contour *cont;
	Point *v;
	bool contour_flag;
	indir = argv[1];
	min_section = (int) strtod(argv[3],&eptr);
	max_section = (int) strtod(argv[4],&eptr);
	// adjust indir
	strcpy(filename,indir);
	temp=strrchr(indir,'/');
	if(!temp) {strcat(filename,"/");}
	else if(*++temp) {strcat(filename,"/");}
	strcat(filename,argv[2]);
	// for each reconstruct input file
	printf("\n");
	for (i=min_section;i<max_section+1;i++) {
		// open file
		sprintf(infile,"%s.%d",filename,i);
		F = fopen(infile,"r");
		if (!F) { printf("Couldn't open input file %s\n",infile); return NULL;}
		else{ printf("Input file found: %s\n",infile); }
		contour_flag = false;
		// for every line in file
		for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) {
			// if start of contour
			if (strstr(str,"<Contour")!=NULL) {
				// find name
				name = strstr(str,"name=");
				name += 6; // advance pointer to start of name
				// create new contour
				cont = new Contour(name,i);
				c = new void_list();
				c->next = ch;
				c->data = (void*)cont;
				ch = c;
				printf("Contour found: %s\n",((Contour*)ch->data)->name);
				// set contour flag
				contour_flag = true;
			} 
			else if (strstr(str,"/>")!=NULL){contour_flag = false;}
			else if (contour_flag) {
				// add point to contour
				v = new Point(str,i,thickness);
				q = new void_list();
				q->next = ((Contour*)ch->data)->raw_points;
				q->data = (void*)v;
				((Contour*)ch->data)->raw_points = q;
			}
		}
		fclose(F);
	}
	printf("\n");
	return ch;
}

int countContours(void_list *ch){
	void_list *q;
	int i=0;
	// for each contour
	for (q=ch;q!=NULL;q=q->next) {i++;}
	return i;
}

int * createArray(int *c_array,void_list* ch,int num_contours){
	int i=0,num_points;
	void_list *q,*p,*qq;
	c_array=new int[num_contours];
	// for each contour
	for (q=ch;q!=NULL;q=q->next) {
		p = ((Contour*)q->data)->raw_points;
		// for each point
		num_points=0;
		for (qq=p;qq!=NULL;qq=qq->next){num_points++;}
		c_array[i]=num_points;
		i++;
	}
	return c_array;
}

void fitSplines(void_list *ch,Histogram *h,double t,int *c_array,double dev_thr,
				double scale,char *outdir,Parameters *pa,int num_contours){
	SplinePoint *sp;
	Contour *c;
	Weights *w1,*w2;
	w1 = new Weights;
	void_list *q,*p;
	double dt;
	int m=0,i=0;
	// for each contour
	for (q=ch;q!=NULL;q=q->next) {
		printf("Splining contour %i of %i\n",i,num_contours-1);
		c=(Contour*)q->data;
		if(c_array[m]<3){
			printf("contour has less than 3 points and ");
			printf("was skipped: contour %s, section %d,",c->name,c->section);
			printf(" num_points %d\n",c_array[m]);
		} else {
			///// load arrays of weights /////
			p = c->rawend;
			w1->bx = new double[4*c_array[m]];
			w1->by = new double[4*c_array[m]];
			w1 = loadWeightArrays(p,w1,c_array[m]);
			// W ARRAY IS CONSTRUCTED WITH FIRST WEIGHT CORRESPONDING TO FIRST RAW POINT IN INPUT FILE.
			///// compute splines for current contour /////
			dt = 1.0/(double)pa->num;
			sp = new SplinePoint[c_array[m]*pa->num];
			sp = computeSplines(sp,w1,dt,c_array[m],pa);
			// SP IS CONSTRUCTED WITH SAME ORIENTATION AS W.
			///// sample splines /////
			sampleSplines(q,sp,dt,c_array[m],t,outdir,i,pa);
			// INTERPOLATED POINTS HAVE SAME ORIENTATION AS RAW POINTS.
			// FIRST OFF THE LINKED LIST WAS LAST ADDED TO LIST.
			///// compute deviation distance between raw points and splines /////
			computeDeviation(q,sp,c_array[m],pa->num);
			///// check deviations /////
			if (dev_thr) {
				w2 = checkDeviation(q,p,w1,c_array[m],dev_thr,scale);
				if(w2!=NULL){
					// recompute splines
					delete[] sp;
					sp = new SplinePoint[c_array[m]*pa->num];
					sp = computeSplines(sp,w2,dt,c_array[m],pa);
					// clear inter_points in contour
					c->clearSpline();
					sampleSplines(q,sp,dt,c_array[m],t,outdir,i,pa);
					computeDeviation(q,sp,c_array[m],pa->num);
				}
			}
			///// update min and max deviation distances /////
			h->update(q,c_array[m]);
			if (pa->diag){printDiagnostics(q,sp,c_array[m]*pa->num,outdir,i);}
			// clean up	
			delete[] w1->bx;
			delete[] w1->by;
			delete[] sp;
		}
		m++;i++;
	}
	delete w1;
}

void computeHistogram(Histogram *h,void_list *ch,int *c_array){
	void_list *q;
	Contour *c;
	double foo=0;
	int i,m=0;
	h->mean = h->sum/h->num;
	// for each contour
	for (q=ch;q!=NULL;q=q->next) {
		c=(Contour*)q->data;
		if(c_array[m]<3){
			printf("contour has less than 3 points and was skipped:");
			printf(" contour %s, section %d, num_points %d\n",c->name,c->section,c_array[m]);
		} else {
			// for each raw point in contour
			for(i=0;i<c_array[m];i++){
				// compute stddev scratch work
				foo+=(c->deviations[i]-h->mean)*(c->deviations[i]-h->mean);
				// bin deviation
				h->count[(int)(c->deviations[i]/(h->max/14))]++;
			}
		}
		m++;
	}
	h->stddev=sqrt(foo/(h->num-1));
}

void_list* createObjects(void_list *ch){
	void_list *q,*qq,*pp,*objectsh,*target;
	objectsh = NULL;
	Contour *c;
	Object *o;
	// for each contour
	for (q=ch;q!=NULL;q=q->next) {
		c=(Contour*)q->data;
		// has object been created with same name?
		target=NULL;
		qq=objectsh;
		while (qq!=NULL && !target) {
			if(!strcmp(((Object*)qq->data)->name,c->name)){target=qq;}
			qq=qq->next;
		}
		// if not
		if (!target) {
			// create a new object and save pointer to new object
			pp = new void_list();
			pp->next = objectsh;
			o = new Object(c->name,c->section);
			pp->data = (void*)o;
			objectsh = pp;
			target = pp;
		}
		// add contour to pointer object
		o=(Object*)target->data;
		pp = new void_list();
		pp->next = o->contours;
		pp->data = (void*)c;
		o->contours = pp;
		// update pointer object min and max
		if (c->section<o->min_section){o->min_section=c->section;}
		if (c->section>o->max_section){o->max_section=c->section;}
	}
	return objectsh;
}

void cleanup(void_list *objectsh,Histogram *h,void_list *ch,int *c_array){
	void_list *p,*q;
	// delete objects
	p=objectsh;
	while (p!=NULL) {
		q=p->next;
		delete (Object*)p->data;
		delete p;
		p=q;
	}
	// delete histogram
	delete h;
	// delete contours
	p=ch;
	while (p!=NULL) {
		q=p->next;
		delete (Contour*)p->data;
		delete p;
		p=q;
	}
	// delete array
	delete[] c_array;
}

void clearPtsFiles(char *outdir,void_list *o){
	////////// clear any existing pts files //////////
	void_list *q,*qq;
	char filename[256],line[2048];
	FILE *F;
	Contour *c;

	// for each object
	for (q=o;q!=NULL;q=q->next) {
		// for each contour in object
		for (qq=((Object*)q->data)->contours;qq!=NULL;qq=qq->next) {
			c=(Contour*)qq->data;
			// create pts file
			sprintf(filename,"%s%s%d.pts",outdir,c->name,c->section);
			// open pts file
	        F = fopen(filename,"w");
	        if (!F) { printf("Couldn't open output file %s\n",filename); exit(0); }
			// close pts file
			fclose(F);
		}
	}

	////////// clear script files //////////
	sprintf(filename,"%smesh.csh",outdir);
	F = fopen(filename,"w");
	if (!F) { printf("Couldn't open output file %s\n",filename); exit(0); }
	sprintf(line,"#!/bin/csh\n\n");
	fputs(line,F);
	fclose(F);
	sprintf(filename,"%sconvert.csh",outdir);
	F = fopen(filename,"w");
	if (!F) { printf("Couldn't open output file %s\n",filename); exit(0); }
	sprintf(line,"#!/bin/csh\n\n");
	fputs(line,F);
	fclose(F);
}

void printConfigFile(char *outdir,Object* o,int capping_flag){
	///// print config file /////
	char filename[256],line[2048];
	FILE *F;
	sprintf(filename,"%s%s.config",outdir,o->name);
	// open file
	F = fopen(filename,"w");
	if (!F) { printf("Couldn't open output file %s\n",filename); exit(1); }
	// print file contents 
	sprintf(line,"PREFIX: %s\nSUFFIX: .pts\nOUTPUT_FILE_NAME: %s\n",o->name ,o->name);
	fputs(line,F);
	if(capping_flag) {
		sprintf(line,"SLICE_RANGE: %i %i\n",o->min_section-1 ,o->max_section+1);
		sprintf(line,"%sMERGE_DISTANCE_SQUARE: 1E-24\n",line);
	}else{
		sprintf(line,"SLICE_RANGE: %i %i\n",o->min_section ,o->max_section);
		sprintf(line,"%snMERGE_DISTANCE_SQUARE: 1E-24\n",line);
	}
	fputs(line,F);
	// close pts file
	fclose(F);
}

void appendScriptFile(char *outdir,Object* o){
	///// append to script files /////
	char filename[256],line[2048];
	FILE *F;
	sprintf(filename,"%smesh.csh",outdir);
	// open file
	F = fopen(filename,"a");
	if (!F) { printf("Couldn't open output file %s\n",filename); exit(1); }
	// print file contents 
	sprintf(line,"echo ''\ncontour_tiler -f %s.config >&! /dev/null\n",o->name);
	sprintf(line,"%secho '%s meshed'\n",line,o->name);
	fputs(line,F);
	// close pts file
	fclose(F);
	sprintf(filename,"%sconvert.csh",outdir);
	// open file
	F = fopen(filename,"a");
	if (!F) { printf("Couldn't open output file %s\n",filename); exit(1); }
	// print file contents 
	sprintf(line,"echo ''\npoly2mesh %s.poly >! %s.mesh\n",o->name ,o->name);
	sprintf(line,"%secho '%s converted'\n",line,o->name);
	fputs(line,F);
	// close pts file
	fclose(F);
}

void printPtsFiles(char *outdir,Object* o,double scale){
	///// print pts files of interpolated contour points /////
	void_list *qq,*pp;
	char filename[256],line[2048];
	FILE *F;
	Contour *c;
	Point *p;
	// for each contour in object
	for (qq=o->contours;qq!=NULL;qq=qq->next) {
		c=(Contour*)qq->data;
		// create pts file
		sprintf(filename,"%s%s%d.pts",outdir,c->name,c->section);
		// open pts file
        F = fopen(filename,"a");
		if(!F){printf("Couldn't open output file %s\n",filename);exit(1);}
		// print number of contour points
		sprintf(line,"%i\n",c->num_interp_points);
		fputs(line,F);
		// for each interpolated point in contour
		for (pp=c->interp_points;pp!=NULL;pp=pp->next) {
			p=(Point*)pp->data;
			// print interpolated contour points
			sprintf(line,"%.15g %.15g %.15g\n",p->x*scale,p->y*scale,p->z*scale);
			fputs(line,F);
		}
		// close pts file
		fclose(F);
	}
}

void printCaps(char *outdir,Object* o,double thickness,double scale){
	///// print min capping pts file /////
	char filename[256],line[2048];
	FILE *F;
	sprintf(filename,"%s%s%d.pts",outdir,o->name,o->min_section-1);
	// open file
	F = fopen(filename,"w");
	if (!F) { printf("Couldn't open output file %s\n",filename);exit(1);}
	// print file contents
	sprintf(line,"1\n0.0 0.0 %d\n",(int)((o->min_section-1)*thickness*scale));
	fputs(line,F);
	// close pts file
	fclose(F);
	///// print max capping pts file /////
	sprintf(filename,"%s%s%d.pts",outdir,o->name,o->max_section+1);
	// open file
 	F = fopen(filename,"w");
	if (!F) { printf("Couldn't open output file %s\n",filename);exit(1);}
	// print file contents 
	sprintf(line,"1\n0.0 0.0 %d\n",(int)((o->max_section+1)*thickness*scale));
	fputs(line,F);
	// close pts file
	fclose(F);
}

void createCallingScript(char *outdir,char *script){
	char filename[256],line[2048];
	FILE *F;
	sprintf(filename,"%s%s",outdir,script);
	F = fopen(filename,"w");
	if (!F) { printf("Couldn't open output file %s\n",filename);exit(1);}
	sprintf(line,"#!/bin/csh\n\n/bin/csh mesh.csh\n/bin/csh convert.csh\n");
	fputs(line,F);
	fclose(F);
}

void printStatistics(Histogram *h,double scale){
	////////// print deviation statistics //////////
	printf("\n\nSpline deviation statistics:\n\n");
	printf("  Avg deviation = %g +- %g\n",h->mean*scale,h->stddev*scale);
	printf("  Smallest deviation:  %10.5g   |  Largest deviation:  %10.5g\n\n",
			h->min*scale,h->max*scale);
	printf("  Deviation histogram:\n");
	printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
			h->max/14.0*0.0*scale,h->max/14.0*1.0*scale,h->count[0],
			h->max/14.0*8.0*scale,h->max/14.0*9.0*scale,h->count[9]);
	printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
			h->max/14.0*1.0*scale,h->max/14.0*2.0*scale,h->count[1],
			h->max/14.0*9.0*scale,h->max/14.0*10.0*scale,h->count[10]);
	printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
			h->max/14.0*2.0*scale,h->max/14.0*3.0*scale,h->count[2],
			h->max/14.0*10.0*scale,h->max/14.0*11.0*scale,h->count[11]);
	printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
			h->max/14.0*3.0*scale,h->max/14.0*4.0*scale,h->count[3],
			h->max/14.0*11.0*scale,h->max/14.0*12.0*scale,h->count[12]);
	printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
			h->max/14.0*4.0*scale,h->max/14.0*5.0*scale,h->count[4],
			h->max/14.0*12.0*scale,h->max/14.0*13.0*scale,h->count[13]);
	printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
			h->max/14.0*5.0*scale,h->max/14.0*6.0*scale,h->count[5],
			h->max/14.0*13.0*scale,h->max/14.0*14.0*scale,h->count[14]);
	printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
			h->max/14.0*6.0*scale,h->max/14.0*7.0*scale,h->count[6],
			h->max/14.0*14.0*scale,h->max/14.0*15.0*scale,h->count[15]);
	printf("  %6.5g - %-6.5g       : %9d    | %6.5g -                : %9d\n",
			h->max/14.0*7.0*scale,h->max/14.0*8.0*scale,h->count[7],
			h->max/14.0*15.0*scale,0.0);
	printf("\n");
}
