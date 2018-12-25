#include <cmath>
#include <cstdio>
#include <ctime>
#include <cstdlib>     /* srand, rand */
#include <sstream> // stringstream
#include <fstream> // ofstream
#include <string> // string
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "functions.h"
#include "constants.h"

int main( int argc, char** argv )
{
   srand(time(NULL));
   std::ofstream* p_file =new std::ofstream;
   p_file -> open (argv[1]);

   write_idx(p_file);
   for(std::size_t i=0; i<=3; i++)
   {
      simulated_annealing( cst::arr_cooling_param[i], p_file );
      (*p_file) << std::endl;
   }

   delete p_file;
}
