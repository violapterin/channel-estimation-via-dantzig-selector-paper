#include <cmath>
#include <cstdio>
#include <ctime>
#include <cstdlib>     /* srand, rand */
#include <sstream> // stringstream
#include <fstream> // ofstream
#include <string> // string

#include <Eigen/Dense>

#include "functions.h"
#include "const.h"
#include "channel.h"
#include "oommpp.h"
#include "ddss.h"

int main( int argc, char** argv )
{
   Eigen::VectorXd v(2);
   Eigen::MatrixXd m(2,2);

   srand(time(NULL));
   std::ofstream* p_file =new std::ofstream;
   p_file -> open (argv[1]);


   //write_idx(p_file);
   (*p_file) << v.lpNorm<Eigen::Infinity>() << std::endl;

   delete p_file;
}
