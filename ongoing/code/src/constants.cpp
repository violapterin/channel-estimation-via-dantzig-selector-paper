#include <map>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "constants.h"

namespace cst
{
   extern const double pi =3.14159265358979;
   extern const std::complex<double> ii =std::sqrt(-1);
   extern const double a_t_phase =0.5; // a phase factor in uniform linear array of Tx
   extern const double a_r_phase =0.6; // a phase factor in uniform linear array of Rx
   extern const double spread =0.1; // relative spread of cluster
   extern const std::size_t nn =20; // # of antennae, assumed same in Tx and Rx
   extern const std::size_t kk =4; // # of subcarrier
   extern const std::size_t s_blk_kk  =2; // # of subcarriers that share a precoder
   extern const std::size_t nn_blk_kk  =kk /s_blk_kk;
   extern const std::size_t uu =6; // # of user
   extern const std::size_t s_blk_uu  =2; // # of users that share a precoder
   extern const std::size_t nn_blk_uu  =kk /s_blk_uu;
   extern const std::size_t ff =(kk/s_blk_kk) *(uu/s_blk_uu); // # of different precoders
   extern const std::size_t nn_c =6; // # of cluster
   extern const std::size_t nn_s =6; // # of scatterer
   extern const std::size_t t_record =10; // interval between subsequent record
   extern const std::size_t t_change_temp =10; // interval between subsequent record
   extern const double delta =0.3; // simulation step
   extern const double temp_init =1; // initial temperature
}
