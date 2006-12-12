#ifndef CONTOUR_READ_H
#define CONTOUR_READ_H

void close_files();

void file_open_error(char *s);

void open_files(char *name);

FILE* openAnUntiled(int level);

FILE* openATiled(int level);

void tile_all_from_contours(char *out_name, int beg_slice, int dimz, double detz, 
							char *prefix, char *suffix);

CTSlice init_volume_and_slice();

void init_find_contours_from_images();

void free_find_contours_from_images();

int  find_contours_from_images();

#endif
