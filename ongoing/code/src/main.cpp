#include <Python.h>
#include <cmath>
#include <cstdio>
#include <sstream> // stringstream
#include <fstream> // ofstream
#include <string> // string
#include <cstdlib>     /* srand, rand */
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "functions.h"
#include "constants.h"

int main( int argc, char** argv )
{
   srand(time(0));
   std::ofstream* p_file =new std::ofstream;
   p_file -> open (argv[1]);

   for(std::size_t i=1; i<=10; i++)
   {
      std::ostringstream stream_tot;
      stream_tot << i << ' ';
      (*p_file) << stream_tot.str();
   }

   (*p_file) << std::endl;

   for(std::size_t i=1; i<=10; i++)
   {
      std::ostringstream stream_tot;
      stream_tot << 3*i*i +1 << ' ';
      (*p_file) << stream_tot.str();
   }

   (*p_file) << std::endl;

   for(std::size_t i=1; i<=10; i++)
   {
      std::ostringstream stream_tot;
      stream_tot << 4*i*i +2 << ' ';
      (*p_file) << stream_tot.str();
   }

/*
   Map_mat map_f;
   for(std::size_t i=0; i<=num_ff-1; i++)
      { map_f[i] =ublas::zero_matrix <std::complex> (cst::nn,cst::nn); }

   Map_mat map_h;
   for(std::size_t i=0; i<=kk-1; i++)
      { set_rand_channel(map_h[i]); }

   // Simulated annealing
   double old_sum_rate =0;
   double new_sum_rate =0;
   Map_mat map_d;
   for(std::size_t i=0; i<=num_ff-1; i++)
      { map_d[i] =ublas::zero_matrix <std::complex> (cst::nn,cst::nn); }

   for(std::size_t idx=0; idx <=cst::max_iter-1; idx++)
   {
      set_random_step(&map_d);
      new_sum_rate =find_new_sum_rate(map_h, map_f, &map_d);

      double bar =new_sum_rate -old_sum_rate;
      if( rand() /std::RAND_MAX >std::exp(-bar/temp) )
      {
         accept_random_step(map_f, &map_d);
         old_sum_rate =new_sum_rate;
      }

      if( idx %cst::t_record ==0 ){ (*ios) << old_sum_rate << ' '; }
      if( idx %cst::t_change_temp ==0 ){ temp *=0.8; }
   }
*/

   delete p_file;
}
