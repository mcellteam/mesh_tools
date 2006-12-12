#ifndef CT_H
#define CT_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CT_ERROR 0
#define CT_OK    1

typedef enum {
   CTCTF,            /* .dat file, 512 x 512   */
   CTRAS,            /* sun raster file, X x Y */
   CTPGM,            /* portable greymap X x Y */
   CTVOL,            /* ct volume file (raw)   */
   CTSLC,            /* sunysb SLC format      */
   CTH3D             /* hartley 3d transform   */
} CTVolType;

typedef enum {
   CT_MSBF,
   CT_LSBF
} CTByteOrder;

typedef enum {
   CT_8BIT,
   CT_16BIT,
   CT_FLOAT
} CTVoxelBits;


typedef struct {
   CTVolType voltype;         /* the type of volume data            */
   char *prefix, *suffix;     /* prefix and suffix to access volume */

   int xdata, ydata, zdata;   /* number of samples in each direction */
   double xunit, yunit, zunit;/* relative spacing in each direction  */
   int numslice, first, last;
   double noise;              /* noise threshold, values above are not */
                              /* considered in calculating max density */

   FILE *outfp;               /* a file pointe for writing volumes     */
   CTVolType outtype;         /* the output volume format              */

   CTByteOrder byteorder;
   CTVoxelBits voxelbits;
} CTVolumeS, *CTVolume;

typedef struct {
   CTVolume vol;
   int num;
   unsigned short *data;
   unsigned char  *cdata;
   float *fdata;
   int x1, y1, x2, y2,    /* corners of portion in memory   */
       width, height;     /* width,height of portion loaded */
   double mind, maxd;
   double ratio;          /* ratio of nonzero data to total */
} CTSliceS, *CTSlice;

CTSlice  CTSliceRead(CTVolume, int, int, int, int, int);
int      CTSliceWrite(char *, CTSlice);
void     CTSliceSetVals(CTSlice, int, int, int, int, int, int);
void     CTSliceSetVal(CTSlice, int, int, float);
CTSlice  CTSliceFree(CTSlice);
CTSlice  CTSliceCopy(CTSlice, int);
CTSlice  CTSliceCreate(int, int, float, float, CTVolume, int);
void     CTSliceCompMinMaxD(CTSlice);
CTSlice  CTSliceHartley2D(CTSlice);


CTVolume CTVolumeCreate(CTVolType, char *, char *, int, int, int, double,
           double, double, int, int, double, CTByteOrder, CTVoxelBits);
void     CTVolumeFree(CTVolume);
CTVolume CTVolumeInitSLC(char *, double);
CTVolume CTVolumeInitVOL(char *, int, int, double, double,
                         CTByteOrder, CTVoxelBits);
CTVolume CTVolumeInitRAS(char *, char *, int, int, double, double);
CTVolume CTVolumeInitCTF(char *, char *, int, int, double, double, CTByteOrder,
                         CTVoxelBits);
CTVolume CTVolumeInitPGM(char *, char *, int, int, double, double);
CTVolume CTVolumeInitH3D(char *, double);
int      CTVolumeOpen(CTVolume, CTVolType, char *, char *);
int      CTVolumeClose(CTVolume);
CTVolume CTVolumeSlice(CTVolume, CTSlice *, CTSlice *, double [3], double [3],
                       double [3], int, int);

void            CTSliceH2DInit(int, int);
CTSlice         CTSliceHartley2D(CTSlice);
CTVolume        CTVolumeHartley3D(CTVolume);


#define CT_VolType(vol)       ((vol)->voltype)
#define CTVolPrefix(vol)     ((vol)->prefix)
#define CTVolSuffix(vol)     ((vol)->suffix)
#define CTVolXSep(vol)       ((vol)->xunit)
#define CTVolYSep(vol)       ((vol)->yunit)
#define CTVolZSep(vol)       ((vol)->zunit)
#define CTVolNoise(vol)      ((vol)->noise)
#define CTVolMaxX(vol)       (((vol)->xdata-1)*(vol)->xunit)
#define CTVolMaxY(vol)       (((vol)->ydata-1)*(vol)->yunit)
#define CTVolMaxZ(vol)       (((vol)->zdata-1)*(vol)->zunit)
#define CTVolXSize(vol)      ((vol)->xdata)
#define CTVolYSize(vol)      ((vol)->ydata)
#define CTVolZSize(vol)      ((vol)->zdata)
#define CTVolFirstSlice(vol) ((vol)->first)
#define CTVolLastSlice(vol)  ((vol)->last)
#define CTVolNumSlices(vol)  ((vol)->last - (vol)->first + 1)
#define CTVolByteOrder(vol)  ((vol)->byteorder)
#define CTVolVoxelBits(vol)  ((vol)->voxelbits)

#define CTSliceWidth(slice)  ((slice)->width)
#define CTSliceHeight(slice) ((slice)->height)
#define CTSliceXSize(slice)  (CTVolXSize((slice)->vol))
#define CTSliceYSize(slice)  (CTVolYSize((slice)->vol))
#define CTSliceMinX(slice)   ((slice)->x1)
#define CTSliceMinY(slice)   ((slice)->y1)
#define CTSliceMaxX(slice)   ((slice)->x2)
#define CTSliceMaxY(slice)   ((slice)->y2)
#define CTSliceMaxD(slice)   ((slice)->maxd)
#define CTSliceMinD(slice)   ((slice)->mind)
#define CTSliceZVal(slice)   ((slice)->vol->zunit*(slice)->num)
#define CTSliceVol(slice)    ((slice)->vol)
#define CTSliceNum(slice)    ((slice)->num)
#define CTSliceRatio(slice)  ((slice)->ratio)

#define CTSliceChars(slice)  ((slice)->vol->voxelbits==CT_8BIT)
#define CTSliceShorts(slice) ((slice)->vol->voxelbits==CT_16BIT)
#define CTSliceFloats(slice) ((slice)->vol->voxelbits==CT_FLOAT)
#define CTSliceVoxelBits(slice) (CTVolVoxelBits((slice)->vol))

#define CTSliceSData(slice) ((slice)->data)
#define CTSliceCData(slice) ((slice)->cdata)
#define CTSliceFData(slice) ((slice)->fdata)

#define CTSliceShortData(slice, x, y) ((slice)->data[(x-CTSliceMinX(slice))* \
         CTSliceHeight(slice)+(y-CTSliceMinY(slice))])
#define CTSliceCharData(slice, x, y) ((slice)->cdata[(x-CTSliceMinX(slice))* \
         CTSliceHeight(slice)+(y-CTSliceMinY(slice))])
#define CTSliceFloatData(slice, x, y) ((slice)->fdata[(x-CTSliceMinX(slice))* \
         CTSliceHeight(slice)+(y-CTSliceMinY(slice))])
#define CTSliceData(slice, x, y) \
            ((slice)->vol->voxelbits==CT_16BIT?CTSliceShortData(slice, x, y): \
             (slice)->vol->voxelbits==CT_8BIT? CTSliceCharData(slice, x, y): \
             CTSliceFloatData(slice, x, y))

#ifdef __cplusplus
}
#endif

#endif
