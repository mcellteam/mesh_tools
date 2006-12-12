/*
%%%%%%%%%%%%%%%%%%%%%%%%%
The SHASTRA software is not in the Public Domain. It is distributed 
on a person to person basis, solely for educational use and 
permission is NOT granted for its transfer to anyone or 
for its use in any commercial product. 
There is NO warranty on the available software and neither Purdue 
University nor the Applied Algebra and Geometry group directed by C. Bajaj
accept responsibility for the consequences of its use.
%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*
 *  volume.c  --  ct volume interface
 *
 *  programmer: Dan Schikore
 *        date: 6/25/93
 */

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

//#include <utils\utils.h> KLC :-> removed as all the calls to the fzopen() and fzclose()
//were replaced by fopen() and fclose () resply. The original calls were to read from
//compressed files.. ".z" or ".gz" formats..

//#include <interp\interp.h> //KLC

#include "ct/ct.h"
#include "ct/sunras.h"
#include "ct/slice.h"
#include "ct/slc.h"

#define DEBUGx

#define MAX(x, y) ((x)>(y)?(x):(y))
#define MIN(x, y) ((x)<(y)?(x):(y))

/*****************************/
static void flipl(unsigned char *p)
{
   unsigned char t;

   t = p[0];  p[0]=p[3];  p[3] = t;
   t = p[1];  p[1]=p[2];  p[2] = t;
}

void CTVolumeFree(CTVolume vol)
{
   free(vol->prefix);
   free(vol->suffix);
   free(vol);
}

/****************************************************
 * CTVolumeCreate - create a volume
 *         type = one of CTRAS, CTCTF, CTPGM, CTVOL, CTSLC
 ****************************************************
 */
CTVolume CTVolumeCreate(CTVolType type, char *prefix, char *suffix, int xdata,
                        int ydata, int zdata, double xunit, double yunit,
                        double zunit, int first, int last, double noise,
                        CTByteOrder byteorder, CTVoxelBits voxelbits)
{
   CTVolume vol;

   vol = (CTVolume)malloc(sizeof(CTVolumeS));
   vol->voltype = type;
   vol->prefix = (char *)malloc(strlen(prefix)+1);
   strcpy(vol->prefix, prefix);
   vol->suffix = (char *)malloc(strlen(suffix)+1);
   strcpy(vol->suffix, suffix);
   vol->first = first;
   vol->last = last;
   vol->numslice = (last-first+1);
   vol->xdata = xdata;
   vol->ydata = ydata;
   vol->zdata = zdata;
   vol->xunit = xunit;
   vol->yunit = yunit;
   vol->zunit = zunit;
   vol->noise = noise;
   vol->outtype = type;
   vol->byteorder = byteorder;
   vol->voxelbits = voxelbits;
   vol->outfp = NULL;
   return(vol);
}

/*
 * CTVolumeInitSLC - initialize an SLC volume
 */
CTVolume CTVolumeInitSLC(char *file, double noise)
{
   FILE *fp;
   int magic, xsize, ysize, zsize, bits;
   float xunit, yunit, zunit;
//   C_UnitType units;
//   C_DataOrigin origin;
//   C_DataModification datamod;

   //KLC..
   //fzopen() is used to open compressed files..... ".z or .gz" declared in utils\fileutil.h
   //if ((fp=fzopen(file, "r")) == NULL)
   //   return(NULL);

   if ((fp=fopen(file, "r")) == NULL)
      return(NULL);

   fread(&magic, sizeof(int), 1, fp);

   if (magic != SLC_MAGIC)
      return(NULL);

   fread(&xsize, sizeof(int), 1, fp);
   fread(&ysize, sizeof(int), 1, fp);
   fread(&zsize, sizeof(int), 1, fp);
   fread(&bits, sizeof(int), 1, fp);
   fread(&xunit, sizeof(float), 1, fp);
   fread(&yunit, sizeof(float), 1, fp);
   fread(&zunit, sizeof(float), 1, fp);

	fclose(fp);
   //fzclose(fp); //KLC

   return(CTVolumeCreate(CTSLC, file, "", xsize, ysize, zsize, xunit, yunit,
                       zunit, 1, zsize, noise, CT_MSBF,
                       bits==8?CT_8BIT:CT_16BIT));
}

/*
 * CTVolumeInitVOL - initialize a VOL volume
 */
CTVolume CTVolumeInitVOL(char *file, int first, int last, double zsep,
                         double noise,
                         CTByteOrder byteorder, CTVoxelBits voxelbits)
{
   return(CTVolumeCreate(CTVOL, file, "", 256, 256, last-first+1, 1.0, 1.0,
                         zsep, first, last, noise, byteorder, voxelbits));
}

/*
 * CTVolumeInitRAS - initialize a RAS volume
 */
CTVolume CTVolumeInitRAS(char *prefix, char *suffix, int first, int last,
                         double zsep, double noise)
{
   int flipit=0, xsize, ysize;
   struct rasterfile sunheader;
   FILE *fp;

   if ((fp=CTOpenSliceNum(prefix, first, suffix)) == NULL)
      return(NULL);

   fread(&sunheader, sizeof(struct rasterfile), 1, fp);
   
   if (sunheader.ras_magic != RAS_MAGIC) 
   {
      flipl( (unsigned char *) &sunheader.ras_magic);
   
	  if (sunheader.ras_magic == RAS_MAGIC)
         flipit = 1;
      else
         flipl( (unsigned char *) &sunheader.ras_magic);
   }

   if (sunheader.ras_magic != RAS_MAGIC)
      return(NULL);

   if (flipit) 
   {
      flipl((unsigned char *) &sunheader.ras_width);
      flipl((unsigned char *) &sunheader.ras_height);
   }

   xsize = sunheader.ras_width;
   ysize = sunheader.ras_height;

   fclose(fp);
   //fzclose(fp); //KLC

   return(CTVolumeCreate(CTRAS, prefix, suffix, xsize, ysize, last-first+1,
          1.0, 1.0, zsep, first, last, noise,flipit?CT_MSBF:CT_LSBF, CT_8BIT));
}

/*
 * CTVolumeInitCTF - initialize a CTF volume
 */
CTVolume CTVolumeInitCTF(char *prefix, char *suffix, int first, int last,
                         double zsep, double noise, CTByteOrder byteorder,
                         CTVoxelBits voxelbits)
{
   return(CTVolumeCreate(CTCTF, prefix, suffix, 512, 512, last-first+1,
          1.0, 1.0, zsep, first, last, noise, byteorder, voxelbits));
}

/*
 * CTVolumeInitPGM - initialize a PGM volume
 */
CTVolume CTVolumeInitPGM(char *prefix, char *suffix, int first, int last,
                         double zsep, double noise)
{
   int xsize, ysize, c;
   char buf[200];
   FILE *fp;

   if ((fp=CTOpenSliceNum(prefix, first, suffix)) == NULL)
      return(NULL);

   if (getc(fp) != 'P' || ((c=getc(fp)) != '2' && (c != '5'))) 
   {
      /* magic number failed */
      return(NULL);
   }

   if (c=='2' || c=='5') 
   {
      if (fscanf(fp, "%d%d", &xsize, &ysize) != 2) 
	  {
         fgets(buf, 200, fp);
         if (fscanf(fp, "%d%d", &xsize, &ysize) != 2)
            return(NULL);
      }
   }

   fclose(fp);
   //fzclose(fp); //KLC

   if (c=='2')
      return(CTVolumeCreate(CTPGM, prefix, suffix, xsize, ysize, last-first+1,
             1.0, 1.0, zsep, first, last, noise, CT_MSBF, CT_16BIT));
   else
      return(CTVolumeCreate(CTPGM, prefix, suffix, xsize, ysize, last-first+1,
             1.0, 1.0, zsep, first, last, noise, CT_MSBF, CT_8BIT));
}

/*
 * CTVolumeInitH3D - initialize an H3D volume
 */
CTVolume CTVolumeInitH3D(char *file, double noise)
{
   double xunit, yunit, zunit;
   int xsize, ysize, zsize;
   FILE *fp;

   if ((fp=fopen(file, "r")) == NULL)
      return(NULL);

   //fzopen is used for compressed files.. ".z" or ".gz" files declared in utils\fileutil.h
   //if ((fp=fzopen(file, "r")) == NULL)
   //   return(NULL);

   fread(&xsize, sizeof(int), 1, fp);
   fread(&ysize, sizeof(int), 1, fp);
   fread(&zsize, sizeof(int), 1, fp);

   fread(&xunit, sizeof(double), 1, fp);
   fread(&yunit, sizeof(double), 1, fp);
   fread(&zunit, sizeof(double), 1, fp);

   return(CTVolumeCreate(CTH3D, file, "", xsize, ysize, zsize, xunit, yunit,
                       zunit, 1, zsize, noise, CT_MSBF, CT_FLOAT));
}

/*-------------------------------------------------------------------
 * CTVolumeWriteSLCHeader - write the SLC header for a volume
 *-------------------------------------------------------------------
 */
void CTVolumeWriteSLCHeader(CTVolume vol)
{
   int magic, bits;
   float xunit, yunit, zunit;
   C_UnitType unittype;
   C_DataOrigin origin;
   C_DataModification mod;
   C_CompressionType comp;

   magic = SLC_MAGIC;
   fwrite(&magic, sizeof(int), 1, vol->outfp);

   fwrite(&CTVolXSize(vol), sizeof(int), 1, vol->outfp);
   fwrite(&CTVolYSize(vol), sizeof(int), 1, vol->outfp);
   fwrite(&CTVolZSize(vol), sizeof(int), 1, vol->outfp);

   bits = (CTVolVoxelBits(vol)==CT_8BIT?8:16);
   fwrite(&bits, sizeof(int), 1, vol->outfp);

   xunit = (float)CTVolXSep(vol);
   fwrite(&xunit, sizeof(float), 1, vol->outfp);
   
   yunit = (float) CTVolYSep(vol);
   fwrite(&yunit, sizeof(float), 1, vol->outfp);
   
   zunit = (float) CTVolZSep(vol);
   fwrite(&zunit, sizeof(float), 1, vol->outfp);

   unittype = C_MILLIMETER_UNIT;
   fwrite(&unittype, sizeof(C_UnitType), 1, vol->outfp);

   origin = C_COMPUTED_TOMOGRAPHY_DATA;
   fwrite(&origin, sizeof(C_DataOrigin), 1, vol->outfp);

   mod = C_ORIGINAL_DATA;
   fwrite(&mod, sizeof(C_DataModification), 1, vol->outfp);

   comp = C_NO_COMPRESSION;
   fwrite(&comp, sizeof(C_CompressionType), 1, vol->outfp);
}

/*-------------------------------------------------------------------
 * CTVolumeWriteH3DHeader - write the H3D header for a volume
 *-------------------------------------------------------------------
 */
void CTVolumeWriteH3DHeader(CTVolume vol)
{
   fwrite(&CTVolXSize(vol), sizeof(int), 1, vol->outfp);
   fwrite(&CTVolYSize(vol), sizeof(int), 1, vol->outfp);
   fwrite(&CTVolZSize(vol), sizeof(int), 1, vol->outfp);

   fwrite(&CTVolXSep(vol), sizeof(double), 1, vol->outfp);
   fwrite(&CTVolYSep(vol), sizeof(double), 1, vol->outfp);
   fwrite(&CTVolZSep(vol), sizeof(double), 1, vol->outfp);
}

/*-------------------------------------------------------------------
 * CTVolumeOpen - open a volume for writing
 *-------------------------------------------------------------------
 */
int CTVolumeOpen(CTVolume vol, CTVolType type, char *prefix, char *suffix)
{
   char buf[250];

   if (vol->outfp != NULL)
      return(CT_ERROR);
   vol->outtype = type;
   
   switch (type) 
   {
      case CTSLC:
         sprintf(buf, "%s%s", prefix, suffix);
         //vol->outfp = fzopen(buf, "w"); //KLC
		 vol->outfp = fopen(buf, "w");
         CTVolumeWriteSLCHeader(vol);
         break;
   
	  case CTVOL:
         sprintf(buf, "%s%s", prefix, suffix);
         //vol->outfp = fzopen(buf, "w"); //KLC
		 vol->outfp = fopen(buf, "w");
         break;
      
	  case CTH3D:
         sprintf(buf, "%s%s", prefix, suffix);
         //vol->outfp = fzopen(buf, "w"); //KLC
		 vol->outfp = fopen(buf, "w");
         CTVolumeWriteH3DHeader(vol);
         break;
      
	  case CTPGM:
      case CTRAS:
      case CTCTF:
         break;
   }

   return(CT_OK);
}

/*-------------------------------------------------------------------
 * CTVolumeClose - close a volume
 *-------------------------------------------------------------------
 */
int CTVolumeClose(CTVolume vol)
{
   if (vol->outfp == NULL) 
   {
      /* some types don't open a volume file.. */
      if (vol->outtype==CTPGM || vol->outtype==CTRAS || vol->outtype==CTCTF)
         return(CT_OK);
      return(CT_ERROR);
   }
   //fzclose(vol->outfp); //KLC
   fclose(vol->outfp); 
   vol->outfp = NULL;
   return(CT_OK);
}

/*
 * CTVolumeSlice - slice a volume at arbitrary angle
 */
CTVolume CTVolumeSlice(CTVolume vol, CTSlice *slices, CTSlice *slice,
                       double p1[3], double p2[3], double p3[3], int dimx,
                       int dimy)
{
   CTVolume newvol = NULL;
   //KLC.... this function was not being called .... HOPEFULLY!!!
   //the tlerp() function call is not defined in the interp\interp.h or any associated files.
   /*double xv[3], yv[3], p[3], xf, yf;
   int i, x, y, xvox, yvox, zvox;
   double s, t, u;

   newvol=CTVolumeCreate(CTH3D, "", "", dimx, dimy, 1, 1.0, 1.0, 1.0,
                         1, 1, 1.0e50, CT_MSBF, CT_FLOAT);
   *slice = CTSliceCreate(dimx, dimy, 0, 0, newvol, 1);
   for (i=0; i<3; i++) {
      xv[i] = p2[i]-p1[i];
      yv[i] = p3[i]-p1[i];
   }
   for (x=0; x<dimx; x++)
      for (y=0; y<dimy; y++) {
         xf = x/(dimx-1.0);
         yf = y/(dimy-1.0);
         for (i=0; i<3; i++)
            p[i] = p1[i] + xf*xv[i] + yf*yv[i];
         xvox=(p[0]/CTVolXSep(vol));
         yvox=(p[1]/CTVolYSep(vol));
         zvox=(p[2]/CTVolZSep(vol));
         if (xvox < 0 || yvox < 0 || zvox < 0 || xvox > CTVolXSize(vol)-2 ||
             yvox > CTVolYSize(vol)-2 || zvox > CTVolZSize(vol)-2)
            CTSliceFloatData(*slice, x, y) = 0.0;
         else {
            s = (p[0] - xvox*CTVolXSep(vol));
            t = (p[1] - yvox*CTVolYSep(vol));
            u = (p[2] - zvox*CTVolZSep(vol));
            if (CTVolVoxelBits(vol) == CT_8BIT) {
               CTSliceFloatData(*slice, x, y) = tlerp(
                       CTSliceCharData(slices[zvox], xvox, yvox),
                       CTSliceCharData(slices[zvox], xvox+1, yvox),
                       CTSliceCharData(slices[zvox+1], xvox, yvox),
                       CTSliceCharData(slices[zvox+1], xvox, yvox),
                       CTSliceCharData(slices[zvox], xvox, yvox+1),
                       CTSliceCharData(slices[zvox], xvox+1, yvox+1),
                       CTSliceCharData(slices[zvox+1], xvox, yvox+1),
                       CTSliceCharData(slices[zvox+1], xvox+1, yvox+1),
                       s, t, u);
            }
            else if (CTVolVoxelBits(vol) == CT_16BIT) {
               CTSliceFloatData(*slice, x, y) = tlerp(
                       CTSliceShortData(slices[zvox], xvox, yvox),
                       CTSliceShortData(slices[zvox], xvox+1, yvox),
                       CTSliceShortData(slices[zvox+1], xvox, yvox),
                       CTSliceShortData(slices[zvox+1], xvox, yvox),
                       CTSliceShortData(slices[zvox], xvox, yvox+1),
                       CTSliceShortData(slices[zvox], xvox+1, yvox+1),
                       CTSliceShortData(slices[zvox+1], xvox, yvox+1),
                       CTSliceShortData(slices[zvox+1], xvox+1, yvox+1),
                       s, t, u);
            }
            else if (CTVolVoxelBits(vol) == CT_FLOAT) {
               CTSliceFloatData(*slice, x, y) = tlerp(
                       CTSliceFloatData(slices[zvox], xvox, yvox),
                       CTSliceFloatData(slices[zvox], xvox+1, yvox),
                       CTSliceFloatData(slices[zvox+1], xvox, yvox),
                       CTSliceFloatData(slices[zvox+1], xvox, yvox),
                       CTSliceFloatData(slices[zvox], xvox, yvox+1),
                       CTSliceFloatData(slices[zvox], xvox+1, yvox+1),
                       CTSliceFloatData(slices[zvox+1], xvox, yvox+1),
                       CTSliceFloatData(slices[zvox+1], xvox+1, yvox+1),
                       s, t, u);
            }
         }
      }*/
	printf("CTVolumeSlice not defined in volume.c \n\n\n");
    return(newvol);
}
