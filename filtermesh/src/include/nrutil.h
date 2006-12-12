/* hh: I added the #if structure */

#if defined(ANSI) || defined(LINT_ARGS)

double *vector(int, int);
double **matrix(int, int, int, int);
double **convert_matrix(double*, int, int, int, int);
double *dvector(int, int);
double **dmatrix(double**, int, int, int, int);
int *ivector(int, int);
int **imatrix(int**, int, int, int, int);
double **submatrix(double**, int, int, int, int, int, int);
void free_vector(double*, int, int);
void free_dvector(double*, int, int);
void free_ivector(int*, int, int);
void free_matrix(double**, int, int, int, int);
void free_dmatrix(double**, int, int, int, int);
void free_imatrix(int**, int, int, int, int);
void free_submatrix(double**, int, int, int, int);
void free_convert_matrix(double**, int, int, int, int);
void nrerror(const char*);

#else

double *vector();
double **matrix();
double **convert_matrix();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
double **submatrix();
void free_vector();
void free_dvector();
void free_ivector();
void free_matrix();
void free_dmatrix();
void free_imatrix();
void free_submatrix();
void free_convert_matrix();
void nrerror();

#endif
