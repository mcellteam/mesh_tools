/*  Transform a iv (like engine) to a raw file */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


/* global variables */
int  idx=0, polygon[100];

extern char *ifname, *ofname;
extern int no_of_tris , orig_verts;

void Process_Position (char * str, FILE *fp_out)
{
    int    i, hasdigit;
    char   *substr1, *substr2;


    /* deal with the line containing "point [" or "]" */
    substr1 = strstr(str, "point");
    substr2 = strstr(str, "]");
    if (substr1!=NULL || substr2!=NULL) {
        hasdigit=0;
	for (i=0; i<(int)strlen(str); i++) {
	    if (str[i]>='0' && str[i]<='9') 
                  hasdigit=1;
            else if (str[i] != '.' && str[i] != '\n')
                str[i] = ' ';
        }
        if (hasdigit==1) {
           orig_verts++;
           fputs(str, fp_out);
        }
        return;
    }

    for (i=0; i<(int)strlen(str); i++)  {
         if (str[i]==',')  str[i] = ' ';
    }
    orig_verts++;
    fputs(str, fp_out);
    return;
}


int IsRegular(int v0, int v1, int v2) 
{
    if (v0!=v1 && v0!=v2 && v1!=v2)
        return 1;
    else
        return 0;
}


void Add_Triangles(FILE *fp_out)
{

     if (idx<3) 
         return;

     if (idx==3) {
         if (IsRegular(polygon[0], polygon[1], polygon[2])) {
            fprintf(fp_out, "\t%d  %d  %d \n", polygon[0], polygon[1], polygon[2]);
            no_of_tris++;
         }
         idx = 0;
         return;
     }

     // Do 2D triangularation
     for (int i=1; i<idx-1; i++) {
         if (IsRegular(polygon[0], polygon[i], polygon[i+1])) {
            fprintf(fp_out, "\t%d  %d  %d \n", polygon[0], polygon[i], polygon[i+1]);
            no_of_tris++;
         }
     }

/*
     for (int i=0; i<idx-2; i++) {
         if ((i%2)==0)
            fprintf(fp_out, "\t%d  %d  %d \n", polygon[i], polygon[i+1], polygon[i+2]);
         else
            fprintf(fp_out, "\t%d  %d  %d \n", polygon[i], polygon[i+2], polygon[i+1]);
         no_of_tris++;
     }
*/
     

     idx=0;
     return;
}


void Get_Triangles(char *str, FILE *fp_out)
{
    int k=0, tmp;
    int length = strlen(str);

    while (k<length) {
        if (str[k]==' ' || str[k]=='\t') {
           k++;
           continue;
        }

        if (str[k]=='\n' || str[k]=='\0') 
            return;

        if (str[k]=='-') {
           Add_Triangles(fp_out);
           k = k+2;
           continue;
        }
        
        tmp=str[k] - '0';
        k++;
        while (str[k]!=' ' && str[k]!='\t' && str[k]!='-') {
            tmp = 10*tmp + str[k] - '0';
            k++;
        }
        polygon[idx] = tmp;
        idx++;
    }
    return;
}


void Process_Triangle (char *str, FILE *fp_out)
{

    int   i, hasdigit;
    char  *substr1, *substr2;


    /* deal with the line containing "point [" or "]" */
    substr1 = strstr(str, "coordIndex");
    substr2 = strstr(str, "]");
    if (substr1!=NULL || substr2!=NULL) {
        hasdigit=0;
        for (i=0; i<(int)strlen(str); i++) {
            if (str[i]>='0' && str[i]<='9')
                  hasdigit=1;
            else if (str[i] != '-' && str[i] != '\n' && str[i] != '\0' )
                str[i] = ' ';
        }
        if (hasdigit==1) {
           Get_Triangles (str, fp_out);
        }
        return;
    }


    for (i=0; i<(int)strlen(str); i++)  {
         if (str[i]==',')  str[i] = ' ';
    }
    Get_Triangles (str, fp_out);
    return;
}



void vrml2raw()
{
   char   *substr;
   char   str[120];
   FILE  *fp_in, *fp_out;

   int    inptstate=0, outptstate=0;
   int    inidxstate=0;
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
              Process_Position(str, fp_out);
              inptstate = 1;
              continue;
          }
       }

       if (inptstate==1) {
          substr = strchr (str, ']');
          if (substr!=NULL) {
              Process_Position(str, fp_out);
              inptstate = 0;
              outptstate = 1;
          } else 
            Process_Position (str, fp_out);
          continue;
       }
 
       if (inidxstate==0) {
          substr = strstr(str, "coordIndex");
          if (substr!=NULL) {
              Process_Triangle(str, fp_out);
              inidxstate = 1;
              continue;
          }
       }

       if (inidxstate==1) {
          substr = strchr (str, ']');
          if (substr!=NULL) {
              Process_Triangle(str, fp_out);
              inidxstate = 0;
          } else
              Process_Triangle(str, fp_out);
          continue;
       }
   }    

   rewind(fp_out);

   sprintf(str,"%6d %6d\n",orig_verts,no_of_tris);
   fputs(str, fp_out);

   fclose (fp_in);
   fclose (fp_out);
}
