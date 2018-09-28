#ifndef FUNCION_H
#define FUNCION_H

#include "constants.h"
//extern typedef Vec;

Comp determinant(Mat);
double frob_norm(Mat&);
Mat get_array_response(double, std::size_t);
void set_rand_channel( Mat& );
void accept_random_step(Map_mat*, Map_mat*);
void set_random_step(Map_mat*);
double find_new_sum_rate(Map_mat*, Map_mat*, Map_mat*);

#endif
