#ifndef PARSE_CONFIG_H
#define PARSE_CONFIG_H

char * read_parse_file(char *filename);

int parse_config(char *str);

int read_database(char *filename, Interpo_struct *pstru, Name_struct *nstru);

int read_database(char *filename, Interpo_struct *pstru, Name_struct *nstru);

int tcp_read_database(char *str, Interpo_struct *pstru, Name_struct *nstru);

int do_database(Interpo_struct *pstru, Name_struct *nstru);

#endif
