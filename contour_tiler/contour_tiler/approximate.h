#ifndef APPROXIMATE_H
#define APPROXIMATE_H

void adaptive_approx(int start, int end, short *ind_ary,char *flag_ary);

int approximate_line(int num, short *short_group_ary, char *flag_ary);

int approximate_contours( double *toler_ary,VertType *lines_ary, int lines_ind,
					short **group_ary, short *group_ind_ary, int group_num,
					short **new_group_ary, short *new_group_ind_ary);

#endif
