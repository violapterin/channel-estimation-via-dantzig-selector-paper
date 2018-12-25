#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <map>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
typedef std::complex <double> Comp;
typedef ublas::matrix <Comp> Mat;
//typedef std::map< std::size_t, Mat> Map_mat;

namespace cst
{
   extern const double pi;
   extern const Comp ii;
   extern const double a_t_phase;
   extern const double a_r_phase;
   extern const double spread;
   extern const std::size_t nn;
   extern const std::size_t nn_c;
   extern const std::size_t nn_s;
   extern const double max_iter;
   extern const double step_size;
}

#endif
