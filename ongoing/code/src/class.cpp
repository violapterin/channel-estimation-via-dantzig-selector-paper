
void Parameters::Parameters(int num_path) :
      num_path(num_path),
      list_ang_dep(cst::nn, 0),
      list_ang_arr(cst::nn, 0),
      list_amp(cst::nn, 0)
{
}

void Parameters::generate_new_instance()
{
   for(int i=0; i<=this->cst::nn-1; i++)
   {
      this ->list_ang_dep[i] =2 *cst::pi *get_normalized_uniform();
      this ->list_ang_arr[i] =2 *cst::pi *get_normalized_uniform();
      this ->list_amp.[i] +=get_normalized_gausian() +cst::ii *get_normalized_gausian();
   }
}

Mat Parameters::get_channel()
{
   Mat ret;
   ret =ublas::zero_matrix <Comp> (this->cst::nn, this->cst::nn);
   for(std::int i_path=0; i_path <num_path-1; i_path++)
   {
      Mat a_t =get_array_response(cst::a_t_phase *psi);
      Mat tr_a_r =ublas::trans(get_array_response(cst::a_r_phase *psi));
      ret =ret +ublas::prod( a_t, tr_a_r );
   }
   return ret;
}


