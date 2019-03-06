#ifndef CLASS_H
#define CLASS_H

#include <vector>
#include "const.h"

//class Comp;
//class Mat;
//class Vec;

class Channel
{
public:
   Channel();
   Channel(int, int);
   ~Channel();
   void realize();
   Mat get_matrix();

private:
   int nn_path;
   std::vector<double> list_ang_dep;
   std::vector<double> list_ang_arr;
   std::vector<Comp> list_amp;

};


#endif
