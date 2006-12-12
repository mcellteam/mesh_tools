#ifndef MYUTIL_H
#define MYUTIL_H

#include "ct/ct.h"

char *mymalloc(int size);

char *mycalloc(int size);

void malloc_copy(char **str1, char *str2);

void print_slice_structure(CTSlice p);

void print_volume_structure(CTVolume p);

CTVolume  InitVolume(int type, char *prefix, char *suffix, int first, int last,
					 double zunits);

#endif
