#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <map>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
typedef std::complex<double> Comp;
typedef ublas::matrix <Comp> Mat;
typedef std::map< std::size_t, Mat> Map_mat;

namespace cst
{
   extern const double pi;
   extern const Comp ii;
   extern const double a_t_phase;
   extern const double a_r_phase;
   extern const double spread;
   extern const std::size_t nn;
   extern const std::size_t kk;
   extern const std::size_t s_blk_kk;
   extern const std::size_t nn_blk_kk;
   extern const std::size_t uu;
   extern const std::size_t s_blk_uu;
   extern const std::size_t nn_blk_uu;
   extern const std::size_t ff;
   extern const std::size_t nn_c;
   extern const std::size_t nn_s;
   extern const double pp;
   extern const std::size_t t_record;
   extern const std::size_t t_change_temp;
   extern const double max_iter;
   extern const double step_size;
   extern const double delta;
   extern const double temp_init;
   extern const double arr_cooling_param[4];
}

#endif
