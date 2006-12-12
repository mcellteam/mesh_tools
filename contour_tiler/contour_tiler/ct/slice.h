#ifndef SLICE_H
#define SLICE_H

#ifdef __cplusplus
extern "C" {
#endif

FILE *CTOpenFile(CTVolume, int);
FILE *CTOpenSliceNum(char *, int, char *);

#ifdef __cplusplus
}
#endif

#endif
