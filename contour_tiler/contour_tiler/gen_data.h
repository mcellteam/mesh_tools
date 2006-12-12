#ifndef GENDATA_H
#define GENDATA_H

CTVolume myInitVolume(int type, char *prefix, char *suffix, int beg, int end,
					  double detz);

CTVolume myInitVolume16bits(int type, char *prefix, char *suffix, int beg, int end,
							double detz);

CTSlice myCTSliceRead(CTVolume vol, int num, int x1, int y1, int x2, int y2);

void generate_ellipse(unsigned short *sbuf, int width, int height, int origx, int origy,
					  double a, double b);

CTSlice specialCTSliceRead(CTVolume vol, int num, int x1, int y1, int x2, int y2);


#endif
