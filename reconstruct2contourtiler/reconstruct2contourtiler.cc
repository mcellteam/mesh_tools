#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "classes.cc"
#include "functions.cc"

int main(int argc,char *argv[]){

	if (argc != 10)
	{
		printf("\nSyntax: reconstruct2contourtiler input_directory ");
		printf("filename_prefix min_section max_section section_thickess ");
		printf("scale output_directory capping_flag deviation_threshold\n\n");
		printf("Description: Converts reconstruct contour format to contour_tiler input format.\n");
		printf("		All files in input directory are assumed to be ");
		printf("of the form filename_prefix.section#.\n");
		printf("		min_section to max_section is the section range to be converted.\n");
		printf("		Section_thickness should be in same scale as x,y contour points.\n");
		printf("		x,y and z coordinates of sampled splines will be multipled by scale in output.\n");
		printf("		capping_flag=1 to attempt end capping.capping_flag=0 to leave ends open.\n");
		printf("		deviation_threshold is the maximum allowed deviation of the spline from raw\n");
		printf("		contour points in scaled units. Set ");
		printf("deviation_threshold to 0 to disable thresholding.\n\n");
		return 1;
	}
		fprintf(stderr,"\n\nMore sophisticated capping must be ");
		fprintf(stderr,"performed as a postprocess on the pts output files.\n\n");


	////////// declare variables /////////
	char outdir[128],output_script[32] = "mesh_and_convert.csh",*eptr,*temp;
	void_list *q,*ch,*objectsh;
	Histogram *h;
	Parameters *pa;
	Object *o;	
	Contour *c;
	h = new Histogram();
	pa = new Parameters;
	int i,capping_flag,*contour_array,num_contours;
	double thickness,scale,deviation_threshold;

	// spline parameters
	pa->plot_rad_int=.1;// radius of curvature sampling interval for plotting sampling function
	pa->num=100;		// num is the # samples of the splines between contour points
						// sequential triplets of sampled points are used to
						// compute the radius of curvature , e.g. num/2 between contour points
	pa->diag=false;		// set true to print diagnostic files
	pa->max_rad=1E10;	// calculated radius of curvature of spline will saturate at this value
	pa->dmin=.001;		// dmin = minimum spline sampling distance
	pa->dmax=.050;		// dmax = maximum spline sampling distance
	pa->T=1.0;			// sample period (= time to traverse dmax)
	pa->amax=1E-2;		// max radial acceleration
//	pa->tau=2;			// decrease tau to increase steepness of sampling function
//	pa->rI=-6;			// inflection point of sampling function = tau*rI
						// decrease rI to sample finer

	////////// get data /////////
	thickness = strtod(argv[5],&eptr);
	scale = strtod(argv[6],&eptr);
	strcpy(outdir,argv[7]);
	capping_flag = (int)strtod(argv[8],&eptr);
	deviation_threshold = strtod(argv[9],&eptr);

	// adjust outdir
	temp=strrchr(argv[7],'/');
	if(!temp) {strcat(outdir,"/");}
	else if(*++temp) {strcat(outdir,"/");}

	// create contours
	ch = getContours(argc,argv,thickness);

	////////// add previous data //////////
	for (q=ch;q!=NULL;q=q->next){((Contour*)q->data)->addPreviousRaw();}

	///// check contours for duplicate points /////
	for (q=ch;q!=NULL;q=q->next){((Contour*)q->data)->removeDuplicates();}

	// count contours
	num_contours=countContours(ch);

	///// linearly interpolate raw contour points /////
	if (deviation_threshold) {
		i=0;
		for (q=ch;q!=NULL;q=q->next) {
			c=(Contour*)q->data;
			printf("Contour %s, section %d\n",c->name,c->section);
			c->linearlyInterp(deviation_threshold,scale);
		}
	}

	////////// add previous data //////////
	for (q=ch;q!=NULL;q=q->next){
		((Contour*)q->data)->addPreviousRaw();
	}

	// NOTE: RAW POINTS IN CONTOURS ARE STORED IN REVERSE ORDER
	// FIRST OFF THE LINKED LIST WAS LAST ADDED TO LIST.


	////////// create array of number of points in each contour //////////
	contour_array = createArray(contour_array,ch,num_contours);

	////////// fit splines //////////
	fitSplines(ch,h,thickness,contour_array,deviation_threshold,scale,outdir,pa,num_contours);

	///// compute histogram /////
	computeHistogram(h,ch,contour_array);

	////////// create objects //////////
	objectsh=createObjects(ch);

	////////// clear any existing pts and script files //////////
	clearPtsFiles(outdir,objectsh);

	////////// print each object //////////
	// for each object
	for (q=objectsh;q!=NULL;q=q->next) {
		o=(Object*)q->data;
		///// print config file /////
		printConfigFile(outdir,o,capping_flag);
		/////  append to script files /////
		appendScriptFile(outdir,o);
		///// print pts files of interpolated contour points /////
		printPtsFiles(outdir,o,scale);
		if(capping_flag){printCaps(outdir,o,thickness,scale);}

	}

	////////// create calling script //////////
	createCallingScript(outdir,output_script);

	////////// print deviation statistics //////////
	printStatistics(h,scale);

	cleanup(objectsh,h,ch,contour_array);

	return 0;
}
