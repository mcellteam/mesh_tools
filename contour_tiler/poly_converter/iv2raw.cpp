
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

extern char *ifname, *ofname;
extern int no_of_tris , orig_verts;

int comp_comma(char  *str)
{
    int tmp=0;
    int nd=0;

    for (int i=0; i<(int)strlen(str); i++) 
         if (str[i]==',')  tmp++;

    if (tmp==3) tmp=4;

    if (tmp==0) {
       for (int i=0; i<(int)strlen(str); i++) 
           if ('0' <= str[i] && str[i] <= '9' ) 
              nd++;
       if (nd>=3) tmp=1;
    }

    return tmp;
}

void iv2raw()
{
   char   str[120];
   FILE  *fp_in, *fp_out;
   int    nc;

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
	   nc = comp_comma(str);
	   if (nc!=1 && nc!=4) continue;

       for (int i=0; i<(int)strlen(str)-1; i++) {
		  if (str[i]=='-' && nc==4) 
			 { str[i] = '\0';  break; }
          if ( str[i] == ',')  str[i] = ' ';
       }

       if (nc==1) {
		  orig_verts++;
          fputs(str, fp_out);
       } else {
		  no_of_tris++;
          fputs(str, fp_out);
		  fputc('\n', fp_out);
	   }
   }    
   rewind(fp_out);

   sprintf(str,"%6d %6d\n",orig_verts,no_of_tris);
   fputs(str, fp_out);

   fclose (fp_in);
   fclose (fp_out);
}
