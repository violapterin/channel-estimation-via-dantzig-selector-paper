#include <map>
#include <complex>
#include <Eigen/Dense>
#include "const.h"


const double Const::pi =3.14159265358979;
const std::complex<double> Const::ii =std::sqrt( static_cast<Comp>(-1) );
const double Const::array_t_phase =0.5;
const double Const::array_r_phase =0.6;
const int Const::nn_hh =24;
const int Const::nn_rr =12;
const int Const::nn_bb =4;
const int Const::max_idx_ssnnrr =10;
const double Const::step_ssnnrr =0.3;
