#ifndef FUNCION_H
#define FUNCION_H

#include "constants.h"
//extern typedef Vec;

Comp determinant(const Mat&);
double frob_norm(const Mat&);
Mat get_array_response(double, std::size_t);
void set_rand_channel( Mat& );
void accept_random_step(Map_mat*, Map_mat*);
void set_random_step(Map_mat*);
double find_new_sum_rate(Map_mat*, Map_mat*, Map_mat*);
void write_idx( std::ostream* );
void simulated_annealing( double, std::ostream* );

#endif
