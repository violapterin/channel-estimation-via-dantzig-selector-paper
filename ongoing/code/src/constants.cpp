#include <map>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "constants.h"

namespace cst
{
   extern const double pi =3.14159265358979;
   extern const std::complex<double> ii =std::sqrt( static_cast<Comp>(-1) );
   extern const double a_t_phase =0.5; // a phase factor in uniform linear array of Tx
   extern const double a_r_phase =0.6; // a phase factor in uniform linear array of Rx
   extern const double spread_cluster =0.1; // relative spread of cluster
   extern const std::size_t nn =7; // # of antennae, assumed same in Tx and Rx
   extern const std::size_t nn_c =6; // # of cluster
   extern const std::size_t nn_s =6; // # of scatterer
   extern const double max_iter =800;
   extern const double step_size =1.3;
}
