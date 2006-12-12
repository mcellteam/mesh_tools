#ifndef CONTOUR_H
#define CONTOUR_H

static void alloc_Lines_ary_if_empty(int index);

static int put_lines_on_map(int index,int line_num, double lines[4][SPACE_VERT_DIM]);

int make_contour_counter_clockwise(signed char is_positive,VertType *vert_ary,
								   short *vert_ind_ary, short poly_num);

void save_a_group(FILE *ofp, VertType *lines_ary, short *group_tab, int group_size);

void save_contours(int topdown, char *s, int slice_num);

void save_all_groups(int topdown, char *s, int slice_num);

int print_group_connection_table(int group_num0, int group_num1);

void allocate_basic_memory();

void allocate_backup_memory();

void free_basic_memory();

void free_backup_memory();

void free_individual_backup_memory();

void group_reset(int index);

void read_contours_from_1_file(int slice_num, char *prefix, char *sufix,int ind);

int read_contours_from_1_file_sub(int slice_num, char *name, int ind);

int read_contours_from_2_files(int k, int ind0, int ind1,char *prefix, char *sufix);

void initialize_tile_table(TileType *tile_table, int num);

void allocate_best_tiling_table(int ind0, int ind1);

int do_major_work(char *out_name, int k, int ind0, int ind1);

void handle_correspondence(int ind0, int ind1, int flag);

void static_save_contours_pts_format(char *s, int slice_num);

int do_approx_only(char *outname, int base_slice_num, int dimz,double zunit, 
				   char *prefix, char *sufix);

void save_both_contours(FILE *ofp);

int tile_from_contours(char *out_name, int base_slice_num, int dimz,double zunit, 
					   char *prefix, char *sufix);

int tile_from_contours_phase1(int base_slice_num, int dimz,char *prefix,
							  char *sufix, int begnum, int endnum);

int closed_all_tiling_files();

int copy_input_to_backup();

void copy_backup_to_input(int level);

void copy_one_level_to_backup(int level);

int  save_backup_files(int num, int begnum, int endnum);

int make_file_consistant(int num);

int calc_volume(int base_slice_num, int dimz, double zunit, char *prefix,
				char *sufix, double *volume_ary);

#endif
