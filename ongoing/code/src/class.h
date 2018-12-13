
class Parameters
{
public:
   Parameters();
   Parameters(int, int);
   ~Parameters();
   void generate_new_instance();
   Mat get_channel();
   
private:
   int nn_path;
   std::vector<double> list_ang_dep;
   std::vector<double> list_ang_arr;
   std::vector<Comp> list_amp;
   
}
