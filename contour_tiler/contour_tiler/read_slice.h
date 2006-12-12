#ifndef READ_SLICE_H
#define READ_SLICE_H

void make_pos_table(CTSlice slice);

void pass_volume(CTVolume vol_in);

float valuefunction(int x, int y, int z);

void slicefunction(int z, int needs);

#endif
