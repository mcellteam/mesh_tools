/*
MCell (tm) Version 2.08 5/18/1998

Copyright (C) 1997,1998, The Salk Institute & Cornell University.
MCell was written jointly by T.M. Bartol Jr. & J.R. Stiles,
with design input for Monte Carlo algorithms from E.E. Salpeter.

  Acknowledgements:
    T.J. Sejnowski for development input and support
    (NSF Grant IBN-9603611), and M.M. Salpeter for fostering
    quantitative experimental applications.  Additional support
    from NIH Grant K08NS01776 (J.R. Stiles).

MCell is a scientific simulation software tool distributed freely to
registered beta-test research laboratories, and must be obtained as a
machine-specific executable program from the authors.  Copying and/or
editing of MCell without the authors' consent is expressly forbidden.
MCell is provided as is and the authors assume no responsibility for
simulation results obtained by users.  Any material published from
MCell simulations must acknowledge the authors and granting
institutions.
*/

/* Header file for 3D vector routines */

#ifndef VECTOR_MATH
#define VECTOR_MATH

#define MY_PI 3.14159265358979323846

struct vector2 {
	double x;
	double y;
};

struct vector3 {
	double x;
	double y;
	double z;
};

void normalize();
void init_matrix();
void scale_matrix();
void translate_matrix();
void rotate_matrix();
void tform_matrix();
void mult_matrix();
double vect_length();
double dot_prod();
void vectorize();
void cross_prod();

#endif
