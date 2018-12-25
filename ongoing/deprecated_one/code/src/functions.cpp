#include <iostream>
#include <cstdio>      /* printf, scanf, puts, NULL */
#include <cstdlib>     /* srand, rand */
#include <ctime>       /* time */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "constants.h"

// stackoverflow.com/questions/29877760/boostublas-how-to-get-determinant-of-int-matrix?rq=1
Comp determinant(const Mat& mm)
{
    Mat m =mm;
    ublas::permutation_matrix<std::size_t> pivots( m.size1() );
    auto be_singular = ublas::lu_factorize(m, pivots);
    if (be_singular) {return 0;}

    Comp det = 1;
    for (std::size_t i = 0; i < pivots.size(); ++i) 
    {
        if (pivots(i) != i) {det *= -1;}
        det *= m(i,i);
    //std::cout << "det= " << det << std::endl; //XXX
    }

    return det;
}

double frob_norm(const Mat& m)
{
   double ret =0;
   for(std::size_t i=0; i<=cst::nn-1; i++)
      for(std::size_t j=0; j<=cst::nn-1; j++)
      {
         double hold_abs =std::abs( m(i,j) );
         ret =ret +hold_abs *hold_abs; // hold
      }
   ret =std::sqrt(ret);
   return ret;
}

Mat get_array_response(double psi, std::size_t len)
{
   Mat ret(len,1);
   for(std::size_t n=0; n <=len -1; n++)
   {
      ret(n,0) =std::exp(cst::ii *static_cast<double>(n) *psi);
   }
   return ret;
}

void set_rand_channel( Mat& h )
{
   h =ublas::zero_matrix <Comp> (cst::nn, cst::nn);
   for(std::size_t n_c=0; n_c <cst::nn_c-1; n_c++)
   {
      double mu =static_cast<double> (rand()) /RAND_MAX;
      mu =2*mu -1;// [-1,1]
      for(std::size_t n_s=0; n_s <cst::nn_s-1; n_s++)
      {
         double delta =static_cast<double> (rand()) /RAND_MAX -0.5;
         delta =cst::spread *delta *delta *delta;// [0, cst::spread *0.5]
         double psi= mu+delta;// [-1-cst::spread *0.5, 1+cst::spread *0.5]
         if( psi > 1 ){ psi =psi-2; }
         if( psi < -1 ){ psi =psi+2; }
         psi =cst::pi*(psi+1);
         Mat a_t =get_array_response(cst::a_t_phase *psi,cst::nn);
         Mat tr_a_r =ublas::trans(get_array_response(cst::a_r_phase *psi,cst::nn));
         h= h +ublas::prod( a_t, tr_a_r );
      }
   }
}

void accept_random_step(Map_mat* p_map_f, Map_mat* p_map_d)
{
   for(std::size_t i=0; i<=cst::ff-1; i++)
      { (*p_map_f)[i] +=(*p_map_d)[i]; }
}

void set_random_step(Map_mat* p_map_d)
{
   for(std::size_t f=0; f<=cst::ff-1; f++)
   {
      Mat& m =(*p_map_d)[f];
      for(std::size_t i=0; i<=cst::nn-1; i++)
         for(std::size_t j=0; j<=cst::nn-1; j++)
            m(i,j) =static_cast<double> (rand()) /RAND_MAX;
      m =cst::step_size *m /frob_norm(m);
   }
}

double find_new_sum_rate(Map_mat* p_map_h, Map_mat* p_map_f, Map_mat* p_map_d)
{
   double ret =0;
   // sum rate
   for(std::size_t k=0; k<=cst::kk-1; k++)
   {
      double sq_norm[cst::uu];
      double sum_sq_norm =0;
      for(std::size_t u=0; u<=cst::uu-1; u++)
      {
         std::size_t f =k /cst::s_blk_kk +(u /cst::s_blk_uu) *cst::nn_blk_kk;
         Mat hold_eff_h =ublas::prod( (*p_map_h)[k], ((*p_map_f)[f] +(*p_map_d)[f]) );
         sq_norm[u] =std::abs( determinant(hold_eff_h) );// hold
         sq_norm[u] =sq_norm[u] * sq_norm[u];
         sum_sq_norm +=sq_norm[u];
         ret =ret +std::log( 1 +(sq_norm[u]) /(1 +sum_sq_norm -sq_norm[u]) );
      }
   }

   // punishment
   double hold_pp_punish =1;
   for(std::size_t f=0; f<=cst::ff-1; f++)
   {
      double hold =cst::pp -frob_norm( (*p_map_f)[f] +(*p_map_d)[f] );
      hold =1 /(hold *hold);
      hold_pp_punish /=(1 +std::exp(hold));
   }
   ret *= hold_pp_punish;
   return ret;
}

void write_idx( std::ostream* p_ios )
{
   for(std::size_t idx=0; idx <=cst::max_iter-1; idx++)
   {
      if( idx %cst::t_record ==0 ){ (*p_ios) << idx << ' '; }
   }
   (*p_ios) << std::endl;
}

void simulated_annealing( double cooling_param, std::ostream* p_ios )
{
   double temp =cst::temp_init;
   Map_mat map_f;
   for(std::size_t f=0; f<=cst::ff-1; f++)
      { map_f[f] =ublas::zero_matrix <Comp> (cst::nn, cst::nn); }

   Map_mat map_h;
   for(std::size_t i=0; i<=cst::kk-1; i++)
      { set_rand_channel(map_h[i]); }

   double old_sum_rate =0;
   double new_sum_rate =0;
   Map_mat map_d;
   for(std::size_t i=0; i<=cst::ff-1; i++)
      { map_d[i] =ublas::zero_matrix <Comp> (cst::nn,cst::nn); }

   for(std::size_t idx=0; idx <=cst::max_iter-1; idx++)
   {
      set_random_step(&map_d);
      new_sum_rate =find_new_sum_rate(&map_h, &map_f, &map_d);

      double bar =new_sum_rate -old_sum_rate;
      if( static_cast<double> (rand()) /RAND_MAX >std::exp(-bar/temp) )
      {
         accept_random_step(&map_f, &map_d);
         old_sum_rate =new_sum_rate;
      }

      if( idx %cst::t_record ==0 ){ (*p_ios) << old_sum_rate << ' '; }
      if( idx %cst::t_change_temp ==0 ){ temp *=cooling_param; }
   }
}
