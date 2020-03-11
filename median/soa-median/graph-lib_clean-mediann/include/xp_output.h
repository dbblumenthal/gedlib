#ifndef __XP_OUTPUT_H__
#define __XP_OUTPUT_H__

#include <sstream>

std::stringstream _xp_out_;

void init_xp_output(){
  _xp_out_.precision(6);
  _xp_out_ << std::fixed;
  _xp_out_.str("");
}


#ifdef _OPENMP
  int* _distances_;

  void init_distances_omp(int p){
    _distances_ = new int[p];
  }

  void destroy_distances_omp(){
    delete[] _distances_;
  }
#endif

// count the number of insertion, deletions and substitutions 
int _n_ins_;
int _n_del_;
int _n_sub_;

#endif
