#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <cassert>

#include <Eigen/Dense>

#include "const.h"
#include "functions.h"
#include "channel.h"
#include "ddss.h"
#include "oommpp.h"

void write_idx( std::ostream* p_ios )
{
   for(int ssnnrr=0; ssnnrr <=Const::max_idx_ssnnrr -1; ssnnrr++)
      { (*p_ios) << ssnnrr << ' '; }

   (*p_ios) << std::endl;
}

double get_uniform()
{
   std::random_device dev;
   std::mt19937 gen(dev());
   std::uniform_real_distribution<> distri(0.0, 1.0);
   return distri(gen);
}

double get_normal()
{
   std::default_random_engine gen;
   std::normal_distribution<double> distri(0,1);
   return distri(gen);
}

Vec get_array_response(double ang, int nn)
{
   assert(nn >=0);
   Vec ret(nn);
   for(int n=0; n <=nn -1; n++)
      { ret(n) =std::exp(Const::ii *static_cast<double>(n) *ang); }
   return ret;
}
