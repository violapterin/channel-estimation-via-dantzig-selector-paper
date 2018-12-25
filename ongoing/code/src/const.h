#ifndef CONST_H
#define CONST_H

#include <map>
#include <complex>
#include <Eigen/Dense>

typedef std::complex <double> Comp;
typedef Eigen::MatrixXcd Mat;
typedef Eigen::VectorXcd Vec;


class Const
{
public:
   static const double pi;
   static const std::complex<double> ii;
   static const double array_t_phase; // a phase factor in uniform linear array of Tx
   static const double array_r_phase; // a phase factor in uniform linear array of Rx
   static const int nn_hh; // # of antennae
   static const int nn_rr; // dim. of analog stage
   static const int nn_bb; // dim. of digital stage
   static const int max_idx_ssnnrr; // maximum index of SNR
   static const double step_ssnnrr; // maximum index of SNR

private:
};

#endif
