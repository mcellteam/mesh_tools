 
/*****************************************************************************/
/*****************************************************************************/
/**                                                                         **/
/** This SHASTRA software is not in the Public Domain. It is distributed on **/
/** a person to person basis, solely for educational use and permission is  **/
/** NOT granted for its transfer to anyone or for its use in any commercial **/
/** product.  There is NO warranty on the available software and neither    **/
/** Purdue University nor the Applied Algebra and Geometry group directed   **/
/** by C.  Bajaj accept responsibility for the consequences of its use.     **/
/**                                                                         **/
/*****************************************************************************/
/*****************************************************************************/

/*---------------------------------------------------------------
  proj1.c  --  Main programs. 

major functions:
1. parse command line input (argv)
---------------------------------------------------------------- */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "common.h"
#include "parse_config.h"
#include "contour.h"
#include "myutil.h"
#include "contour_read.h"
#include "main.h"

int MATH_FUNC = 0,DEBUG=0,SILENT=0;
int CONTOURS_ONLY=0;
int FILTER_TYPE=0;
int FROM_CONTOURS=1;
int SECURITY_CHECK=1;
int SAVE_CONTOURS=0;
int DO_APPROX_ONLY=0;
int STORE_1ST_PASS=0;
int CHANGE_Z=0;
int SEPERATE_CONVEX_FILE=0;
int ALL_CORRESPOND=0;
int CHRISTIANSEN_METHOD=0;
int COUNT_TRIANGLE_ONLY=0;
int NUMERICALLY_STABLE=1;
int NO_OUTPUT_FILE=0;
int MESH_FILE=0; /* prepare file for mesh generation */

int SAME_LEVEL=1; /* if 1: only tile between the same level */
int DIFF_LEVEL_FILES=0; /* if 1: put diff. level into diff. files */ 
int Current_level=0;
int CALC_VOLUME_ONLY=0;
int NON_VALID=0;
int Beg_slice=-1000;

double TURBERANCE=0.001;
float Scale_factor=1.0;
double Tolerance=1.0;
int WriteIndividualFileInPoly=0;
double XSHIFT;
double Read_file_time=0.0;
float FIRST_Z=0;
Interpo_struct Istru;
Name_struct Nstru;
double MERGE_DISTANCE_SQUARE=1e-4;
float Scale_z=1.0;

extern int get_no_of_poly_triangles();
static char *Config_file=NULL;

////////////////
int STATIC_CHQ;

extern void my_clear_approximate();
extern void my_clear_branch();
extern void my_clear_branch_util();
extern void my_clear_contour();
extern void my_clear_contour_read();
extern void my_clear_correspond();
extern void my_clear_cross();
extern void my_clear_CTSliceH2DInit();
extern void my_clear_decom_util();
extern void my_clear_decompose();
extern void my_clear_filter();
extern void my_clear_gen_data();
extern void my_clear_group();
extern void my_clear_image2cont2D();
extern void my_clear_mymarch();
extern void my_clear_parse_config();
extern void my_clear_qsort();
extern void my_clear_read_slice();
extern void my_clear_scale();
extern void my_clear_slice();
extern void my_clear_tile_hash();
extern void my_clear_tile_util();
extern void my_clear_tiling();

////////////////


void clear_all_statics()
{

	MATH_FUNC = 0;
	DEBUG=0,SILENT=0;
	CONTOURS_ONLY=0;
	FILTER_TYPE=0;
	FROM_CONTOURS=1;
	SECURITY_CHECK=1;
	SAVE_CONTOURS=0;
	DO_APPROX_ONLY=0;
	STORE_1ST_PASS=0;
	CHANGE_Z=0;
	SEPERATE_CONVEX_FILE=0;
	ALL_CORRESPOND=0;
	CHRISTIANSEN_METHOD=0;
	COUNT_TRIANGLE_ONLY=0;
	NUMERICALLY_STABLE=1;
	NO_OUTPUT_FILE=0;
	MESH_FILE=0; /* prepare file for mesh generation */

	SAME_LEVEL=1; /* if 1: only tile between the same level */
	DIFF_LEVEL_FILES=0; /* if 1: put diff. level into diff. files */ 
	Current_level=0;
	CALC_VOLUME_ONLY=0;
	NON_VALID=0;
	Beg_slice=-1000;


	TURBERANCE=0.001;
	Scale_factor=1.0;
	Tolerance=1.0;
	WriteIndividualFileInPoly=0;
	XSHIFT;
	Read_file_time=0.0;
	FIRST_Z=0;
	MERGE_DISTANCE_SQUARE=1e-4;
	Scale_z=1.0;


	STATIC_CHQ =-10000;

	my_clear_approximate();
	my_clear_branch();
	my_clear_branch_util();
	my_clear_contour();
	my_clear_contour_read();
	my_clear_correspond();
	my_clear_cross();
	my_clear_CTSliceH2DInit();
	my_clear_decom_util();
	my_clear_decompose();
	my_clear_filter();
	my_clear_gen_data();
	my_clear_group();
	my_clear_image2cont2D();
	my_clear_mymarch();
	my_clear_parse_config();
	my_clear_qsort();
	my_clear_read_slice();
	my_clear_scale();
	my_clear_slice();
	my_clear_tile_hash();
	my_clear_tile_util();
	my_clear_tiling();

}

#ifndef i860  /* for Intel Paragon machine compatibility*/
	int mynode(){return 0;}
	int numnodes(){return 1;}
	void do_process(){return;}
#endif



void print_usage(char *p)
{
    if(mynode())return;

	printf("Usage: %s -f config_file [-V]\n",p);
	printf("config_file: the configuration file, see document\n");
	/*
	printf("-d ? : debug level\n");
	printf("-g : Write individual group into .poly file\n");
	printf("-c : contour file only, wrote into .pts\n");
	printf("-s : save cotours\n");
	printf("-v : separate convex file\n");
	printf("-p : store the 1st pass of 1st tiling in name_1st.poly\n");
	printf("-a ? : scale & shift, so all contours correspond to each other,\n");
	printf("       float between -0.5 to 0.5, for shift factor\n");
	printf("-C : Christiansen_method, only 1-1 and two slices\n");
	printf("-m : count marching cubes triangles only, input must be images\n");
	printf("-n : numerically unstable, more sharp triangles, shorter time\n");
	printf("-t : no output file will be produced\n");
	printf("-l : silence\n");
	*/
	printf("-V : calculate volume area only, no tiling is performed\n");
	exit(1);
}

int parse_argv(int argc, char **argv)
{ 
	/* Is it a group host */
    char c,*p,**argv0;
    int status=0;

    argv0=argv;
    if(argc == 1) {print_usage(argv0[0]); return 0; }
    argv++;         /* skip one command */
    while(*argv)
    {
        p=*argv;
        if(**argv != '-')
        {
            switch(status)
            {
                case 2: /* debug */
					sscanf(*argv,"%lf",&XSHIFT);
					if(XSHIFT>0.5 || XSHIFT<-0.5)
					{
						fprintf(stderr,"Warning, wrong XSHIFT, set to 0.0\n");
						XSHIFT=0.0;
					}
					else fprintf(stderr," XSHIFT =%lf \n",XSHIFT);
	                break;
        
				case 5: /* debug */
                    DEBUG =atoi(*argv);
                    break;
                case 7: 
				    Config_file=*argv;
                    break;
				default:
					if(!mynode())fprintf(stderr,"Unknown command: %s\n",p);
					break;
            }
            status =0;
        }
        else
        {
            if(status!=0) /* previous option needs a value */
            {
                if(!mynode())
					fprintf(stderr,"%s needs a value\n",*(argv-1));
			    print_usage(argv0[0]);
				return 0;
            }

            (*argv)++;      /* skip the '-' */
            c = **argv;     /* command arg */
            (*argv)++;      /* skip  */
            
			if((**argv)!=0)
            {
                fprintf(stderr,"Unknown option %s\n",p);
                print_usage(argv0[0]);
				return 0;
            }

            switch(c)
            {
                case 'd': status=5; break; /* debugger */
                case 'g': status=0; WriteIndividualFileInPoly=1; break;
                case 'a': status=2; ALL_CORRESPOND=1;
					  break; 

                case 'f': status=7; break; /* configuration file */
                case 'c': CONTOURS_ONLY=1; break; /* contours only */
                case 'C': CHRISTIANSEN_METHOD=1; break; /* contours only */
                case 's': SAVE_CONTOURS=1; break; /*  */
                case 'p': STORE_1ST_PASS=1; break; /*  */
                case 'v': SEPERATE_CONVEX_FILE=1; break; /*  */
                case 'V': CALC_VOLUME_ONLY=1; break; /*  */
                case 'm': COUNT_TRIANGLE_ONLY=1; break; /*  */
                case 'n': NUMERICALLY_STABLE=0; break; /*  */
                case 't': NO_OUTPUT_FILE=1; break; /*  */
                case 'l': SILENT=1; break; /*  */
				
				default: 
				fprintf(stderr,"Unknown command: -%c\n",c);
				break;
            }
        }
        argv++;         /* skip one command */
    }

    return 1;
}

void handle_fatal_error()
{
    fprintf(stderr,"fatal error\n");
    exit(1);
}

int main(int argc,char ** argv)
{
	#ifdef i860
		double btime;
	#endif
	/*    int i,err,tri_num; */
    if(argc==1)	
		print_usage(*argv);

    parse_argv(argc, argv);

	STATIC_CHQ = -100; //later, set it to -10000;
    
    if(Config_file==NULL) 
	{
	    printf("Error or no config file\n");
	    exit(1);
    }
    
	read_database(Config_file,&Istru,&Nstru);
    if(!strcmp(Nstru.output_file_name,"no_writing"))
		NO_OUTPUT_FILE=1;

    if(FROM_CONTOURS && COUNT_TRIANGLE_ONLY)
	{
		fprintf(stderr,"ERROR! -m option needs images, no FROM_CONTOURS\n");
		exit(1);
    }

    Istru.thres+=0.001123f; /* randomize it to avoid .... */

    if(Istru.detz>0)	CHANGE_Z=1;

    if(MATH_FUNC>4)		MATH_FUNC=0;
    
	if(MATH_FUNC)		fprintf(stderr,"The input data is generated as :");

    if(MATH_FUNC==1)	fprintf(stderr,"sphere\n");

    else if(MATH_FUNC==2)	fprintf(stderr,"cylinder\n");
    else if(MATH_FUNC==3)	fprintf(stderr,"x^2 + y^2 - z^2\n");
    else if(MATH_FUNC==4)	fprintf(stderr,"special testing functions\n");

    if(Istru.epsilon!=0.0)	Tolerance=Istru.epsilon;
    if(Tolerance>1.0)
	{
		if(!mynode())
			fprintf(stderr,"Warning! change epsilon from %lf to 1.0\n",Tolerance);
		Tolerance=1.0;
    }
    else if(Tolerance<0.0)
	{
		if(!mynode())
			fprintf(stderr,"Warning! change epsilon from %lf to 0.5\n",Tolerance);
		Tolerance=0.5;
    }
	#ifdef i860
		btime=dclock();
	#endif

	if(numnodes()>1)
		do_process();
    else 
	{
		if(DO_APPROX_ONLY && !FROM_CONTOURS)
		{
			if(!mynode())
				fprintf(stderr,"DO_APPROX_ONLY should have FROM_CONTOURS\n");
			exit(1);
		}

		if(DO_APPROX_ONLY)
		{
			do_approx_only(Nstru.output_file_name,Istru.beg, 
							Istru.end - Istru.beg, Istru.detz,Nstru.prefix, Nstru.suffix);
			exit(0);
		}

		if(FROM_CONTOURS)
		{
			Beg_slice = Istru.beg;
			tile_all_from_contours(Nstru.output_file_name,Istru.beg, 
									Istru.end - Istru.beg, Istru.detz,Nstru.prefix,Nstru.suffix);
		#ifdef i860
			btime=dclock()-btime;
			fprintf(stderr,"Total time %lf sec.\n",btime);
		#endif
			printf("no of poly triangles is %d\n", get_no_of_poly_triangles());
//			clear_all_statics();
			exit(0);
		}

		if(!FROM_CONTOURS)
		{
			find_contours_from_images();
			/* The output.pts contours becomes the input data */
			free(Nstru.suffix);
			free(Nstru.prefix);
			malloc_copy(&(Nstru.suffix),".pts");
			malloc_copy(&(Nstru.prefix),Nstru.output_file_name);
			/* the Zunit is controlled by contour, not by user */
			Istru.detz=0.0; 
			if(CONTOURS_ONLY && !COUNT_TRIANGLE_ONLY)
			{
				printf("It generated %d %s?.pts contour files\n",Istru.end - Istru.beg +1,
						Nstru.output_file_name);
				exit(0);
			}
		}

		if(COUNT_TRIANGLE_ONLY)		exit(1);

		tile_all_from_contours(Nstru.output_file_name,Istru.beg, 
								Istru.end - Istru.beg, Istru.detz, Nstru.prefix,Nstru.suffix);

		if(!FROM_CONTOURS)
		{
			printf("It generated %d %s?.pts contour files\n",Istru.end - Istru.beg +1,
					Nstru.output_file_name);
		}
    }

	#ifdef i860
		btime=dclock()-btime;
		fprintf(stderr,"Total time %lf sec.\n",btime);
	#endif
	return 0;
}
	


