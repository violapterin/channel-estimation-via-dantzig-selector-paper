#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>

int main () {
   typedef boost::numeric::ublas::matrix<double> Mat;
   Mat m (3, 3);
   for (unsigned i = 0; i < m.size1 (); ++ i)
      for (unsigned j = 0; j < m.size2 (); ++ j)
         m (i, j) = 3 * i + j;
   std::cout << prod(m,m) << std::endl;
}
