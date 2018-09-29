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

   write_idx(p_file);
   for(std::size_t i=0; i<=4; i++)
   {
      simulated_annealing( cst::arr_cooling_param[i], p_file);
      (*p_file) << std::endl;
   }

   delete p_file;

   /*
// Testing only
   for(std::size_t i=1; i<=10; i++) { (*p_file) << i << ' '; }
   (*p_file) << std::endl;
   for(std::size_t i=1; i<=10; i++) { (*p_file) << 2*i*i << ' '; }
   (*p_file) << std::endl;
   for(std::size_t i=1; i<=10; i++) { (*p_file) << 3*i*i+2*i << ' '; }
   (*p_file) << std::endl;

   delete p_file;
   */
}
