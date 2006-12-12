/*  Transform a wrl (PRO/E) to a raw file */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>


/* global variables */
extern int no_of_tris , orig_verts;

extern char *ifname, *ofname;

void Put_Position (char * str, FILE *fp_out)
{
	int i;

    for ( i=0; i<(int)strlen(str); i++)  {
         if (str[i]==',')  str[i] = ' ';
    }

    orig_verts++;
    fputs(str, fp_out);
}


void Put_Triangle (char *str, FILE *fp_out)
{
    char  tmp[40];
    int i, kk=0;

    for ( i=0; i<(int)strlen(str); i++)  {
         if (str[i]==',')  str[i] = ' ';
    }

    for ( i=0; i<(int)strlen(str); i++)  {
        if (str[i]=='-') {
           tmp[kk]='\n'; 
           tmp[kk+1]='\0'; 
           no_of_tris++;
           fputs(tmp, fp_out);
           kk=0; i++;
        } else {
            tmp[kk] = str[i];
            kk++;
        }
    }
}

void wrl2raw()
{
   char   *substr;
   char   str[120];
   FILE  *fp_in, *fp_out;

   int    inptstate=0, outptstate=0;
   int    inidxstate=0, outidxstate=0;
   no_of_tris = orig_verts =0;

   fp_in = fopen(ifname,"r");
   if (fp_in==NULL) {
      fprintf(stderr, "Cannot open the file: %s\n",ifname);
      exit(1);
   }

   fp_out = fopen(ofname,"w");
   if (fp_out==NULL) {
      fprintf(stderr, "Cannot open the file: %s\n",ofname);
      exit(1);
   }

   sprintf(str,"%6d %6d\n",orig_verts,no_of_tris);
   fputs(str, fp_out);

   while ( fgets(str, 120, fp_in) != NULL) {
       if (inptstate==0 && outptstate==0) {
          substr = strstr(str, "point");
          if (substr!=NULL) {
              inptstate = 1;
              continue;
          }
       }

       if (inptstate==1) {
          substr = strchr (str, ']');
          if (substr!=NULL) {
              inptstate = 0;
              outptstate = 1;
          } else 
            /* fputs(str, fp_out);  */
            Put_Position (str, fp_out);
          continue;
       }
 
       if (inidxstate==0) {
          substr = strstr(str, "coordIndex");
          if (substr!=NULL) {
              inidxstate = 1;
              continue;
          }
       }

       if (inidxstate==1) {
          substr = strchr (str, ']');
          if (substr!=NULL) {
              inidxstate = 0;
          } else
            // fputs(str, fp_out);
            Put_Triangle (str, fp_out);
          continue;
       }
   }    

   rewind(fp_out);

   sprintf(str,"%6d %6d\n",orig_verts,no_of_tris);
   fputs(str, fp_out);

/* fprintf(fp_out,"%d %d\n ",orig_verts,no_of_tris); */

   fclose (fp_in);
   fclose (fp_out);
}
