/** -----------------------------------------------------------
  MEX file for computing a solution to the LSAP with the Hungarian algorithm
  
  see hungarianLSAP.m for more details
  compute this file with compile_mex.m or compile_octave.m
 *  -----------------------------------------------------------
*/

#include <mex.h>
#include <cstdio>
#include <string>
#include <lsap.h>

using namespace lsape;

template <typename DT>
void hungarian_helper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], const mxClassID &mxC)
{
  // 1st input : edit cost matrix
  const DT *C = (DT*)mxGetPr(prhs[0]);
  // get dimensions
  const size_t *nrc = mxGetDimensions(prhs[0]);
  const unsigned int nrows = (unsigned int)nrc[0], ncols = (unsigned int)nrc[1];
  size_t dimr2[2] = { nrc[0], 1 }, dimc2[2] = { 1, nrc[1] };
  bool forbassign = false;
  unsigned short init_type = 1;
  if (nrhs > 1) init_type = (unsigned short)mxGetScalar(prhs[1]);
  if (nrhs == 3) forbassign = (bool)mxGetScalar(prhs[2]);
  //------------------------------------------------------------------
  int i = 0;
  plhs[0] = mxCreateNumericArray(2, dimr2, mxUINT32_CLASS, mxREAL);
  unsigned int *rho = (unsigned int*)mxGetData(plhs[0]);
  unsigned int *varrho = NULL;
  if (nlhs == 2 || nlhs == 4)
  {
    plhs[1] = mxCreateNumericArray(2, dimc2, mxUINT32_CLASS, mxREAL);
    varrho = (unsigned int*)mxGetData(plhs[1]);
  }
  DT *u = NULL, *v = NULL;
  if (nlhs > 2)
  {
    plhs[nlhs-2] = mxCreateNumericArray(2, dimr2, mxC, mxREAL);
    u = (DT*)mxGetPr(plhs[nlhs-2]);
    plhs[nlhs-1] = mxCreateNumericArray(2, dimc2, mxC, mxREAL);
    v = (DT*)mxGetPr(plhs[nlhs-1]);
  }
  else { u = new DT[nrows]; v = new DT[ncols];}
  hungarianLSAP<DT>(C,nrows,ncols,rho,u,v,varrho,init_type,forbassign);
  for (unsigned int i = 0; i < nrows; i++) rho[i]++;
  if (nlhs == 2 || nlhs == 4) for (unsigned int j = 0; j < ncols; j++) varrho[j]++;
  if (nlhs < 3) { delete[] u; delete[] v; }
}

//------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  if (nrhs < 1 || nrhs > 3 || nlhs < 1 || nlhs > 4) mexErrMsgTxt("USAGE: [rho] = hungarianLSAP(C,init_type,forbassign)\nUSAGE: [rho,varrho] = hungarianLSAP(C,init_type,forbassign)\nUSAGE: [rho,u,v] = hungarianLSAP(C,init_type,forbassign)\nUSAGE: [rho,varrho,u,v] = hungarianLSAP(C,init_type,forbassign)\n");
  mxClassID tclass = mxGetClassID(prhs[0]);
  switch(tclass)
  {
    case mxINT16_CLASS: hungarian_helper<short>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxINT32_CLASS: hungarian_helper<int>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxINT64_CLASS: hungarian_helper<int64_T>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxSINGLE_CLASS: hungarian_helper<float>(nlhs,plhs,nrhs,prhs,tclass); break;
    default: hungarian_helper<double>(nlhs,plhs,nrhs,prhs,tclass);
  }
}
