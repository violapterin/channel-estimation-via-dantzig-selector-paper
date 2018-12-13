#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "constants.h"



void newton_method(std::vector<double> list_ang_dep, std::vector<double> list_ang_arr, std::vector<Comp> list_amp)
{
   
}

void walk_one_step(Parameters& param)
{
   
}

void write_idx( std::ostream* p_ios )
{
   for(std::size_t idx=0; idx <=cst::max_iter-1; idx++)
   {
      if( idx %cst::t_record ==0 ){ (*p_ios) << idx << ' '; }
   }
   (*p_ios) << std::endl;
}

double infty_norm(const Mat& m)
{
   double ret =0;
   for(std::size_t i=0; i<=cst::nn-1; i++)
      for(std::size_t j=0; j<=cst::nn-1; j++)
      {
         double hold_abs =std::abs( m(i,j) );
         if(hold_abs >ret){ ret =hold_abs; }
      }
   return ret;
}

double get_normalized_uniform()
{
   std::random_device dev;
   std::mt19937 gen(dev());
   std::uniform_real_distribution<> dist(0.0, 1.0);
   return dist(gen);
}

double get_normalized_gaussian()
{
   std::default_random_engine gen;
   std::normal_distribution<double> dist(0,1);
   return x+dist(gen);
}

double mod_2_pi(double x)
{
   double ret =0;
   bool sgn =false;

   if(x>0){sgn =true;}
   else{sgn =false;}
   if(!sgn){ x =-x; }

   ret =x -2 *cst::pi *static_cast<int>(x/2/cst::pi);
   if(!sgn){ ret =-ret; }
   return ret;
}

Mat get_array_response(double psi)
{
   Mat ret(len,1);
   for(std::size_t n=0; n <=cst::nn -1; n++)
   {
      ret(n,0) =std::exp(cst::ii *static_cast<double>(n) *psi);
   }
   return ret;
}
