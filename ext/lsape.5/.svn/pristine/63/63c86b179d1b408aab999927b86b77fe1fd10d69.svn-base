/* -----------------------------------------------------------
  MEX file for computing an approximate solution to the LSAP with greedy algorithms

  This is a MEX-file for MATLAB.
  
  This file is part of LSAPE.
  LSAPE is free software: you can redistribute it and/or modify it
  under the terms of the CeCILL-C License. See README for more details.

     Copyright 2015
      authors: Sebastien Bougleux
  institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, France
   last modif: March 2018
 
  execute matlab file 'compile_mex' to compile this function and use it in matlab
 *  -----------------------------------------------------------
*/

#include <mex.h>
#include <lsap.h>

using namespace lsape;

template <typename DT>
void greedy_helper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], const mxClassID &mxC)
{
  // 1st input : edit cost matrix
  const DT *C = (DT*)mxGetPr(prhs[0]);
  // get dimensions
  const size_t *nrc = mxGetDimensions(prhs[0]);
  unsigned int nrows = (unsigned int)nrc[0];
  unsigned int ncols = (unsigned int)nrc[1];
  //if (nrows < ncols) mexErrMsgTxt("matrix dimension must be n>=m");
  const size_t dimr2[2] = { nrc[0], 1 }, dimc2[2] = { 1, nrc[1] };
  unsigned short init_type = 2;
  if (nrhs > 1)
  {
    init_type = (unsigned short)mxGetScalar(prhs[1]);
  }
  //------------------------------------------------------------------
  plhs[0] = mxCreateNumericArray(2, dimr2, mxUINT32_CLASS, mxREAL);
  unsigned int *rho = (unsigned int*)mxGetData(plhs[0]);
  unsigned int *varrho = NULL;
  if (nlhs == 3)
  {
    plhs[1] = mxCreateNumericArray(2, dimc2, mxUINT32_CLASS, mxREAL);
    varrho = (unsigned int*)mxGetData(plhs[1]);
  }
  DT cost = greedyLSAP<DT>(C,nrows,ncols,rho,varrho,init_type);
  plhs[nlhs-1] = mxCreateDoubleScalar((double)cost);
  for (unsigned int i = 0; i < nrows; i++) rho[i]++;
  if (nlhs == 3) for (unsigned int j = 0; j < ncols; j++) varrho[j]++;
}

//------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  if (nrhs < 1 || nrhs > 2 || nlhs < 2 || nlhs > 3) mexErrMsgTxt("USAGE: [rho,cost] = greedyLSAP(C,greedy_type)\nUSAGE: [rho,varrho,cost] = greedyLSAP(C,greedy_type)\n");
  mxClassID tclass = mxGetClassID(prhs[0]);
  switch(tclass)
  {
  case mxINT16_CLASS: greedy_helper<short>(nlhs,plhs,nrhs,prhs,tclass); break;
  case mxINT32_CLASS: greedy_helper<int>(nlhs,plhs,nrhs,prhs,tclass); break;
  case mxINT64_CLASS: greedy_helper<int64_T>(nlhs,plhs,nrhs,prhs,tclass); break;
  case mxSINGLE_CLASS: greedy_helper<float>(nlhs,plhs,nrhs,prhs,tclass); break;
  default: greedy_helper<double>(nlhs,plhs,nrhs,prhs,tclass);
  }
}
