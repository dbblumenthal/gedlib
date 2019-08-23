
#include <mex.h>
#include <cstdio>
#include <string>
#include <JacobiAuctionLSAP.h>
#include <utils.hh>

//------------------------------------------------------------------
template <typename DT>
void jacobiAuctionLSAP_helper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], const mxClassID &mxC)
{
  const DT *C = (DT*)mxGetPr(prhs[0]);
  int nrows = (int)mxGetDimensions(prhs[0])[0];
  int ncols = (int)mxGetDimensions(prhs[0])[1];
  int dimr[2] = { nrows, 1 }, dimc[2] = { 1, ncols };
  bool forb_assign = false;
  if (nrhs >= 2) forb_assign = (bool)mxGetScalar(prhs[1]);

  int k = -1;
  if (nrhs == 3) k = (int)mxGetScalar(prhs[2]);

  if (nrows != ncols) return;
  int n = nrows;

  int* rho12 = new int[n];
  DT* u = new DT[n];
  DT* v = new DT[n];

  JacobiAuctionLSAP<DT,int> auctions(C, n);
  auctions(C,n,u,v,rho12);

  //*
  plhs[0] = mxCreateNumericMatrix(n, 1, mxINT32_CLASS, mxREAL);

  int* tab = (int*)mxGetData(plhs[0]);
  for (int i=0; i<n; i++){
    tab[i] = rho12[i]+1;
  }

  if (nlhs == 3){
    plhs[1] = mxCreateNumericMatrix(n, 1, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(n, 1, mxDOUBLE_CLASS, mxREAL);
    double* mx_u = (double*)mxGetData(plhs[1]);
    double* mx_v = (double*)mxGetData(plhs[2]);

    for (int i=0; i<n; i++){
      mx_u[i] = u[i];
      mx_v[i] = v[i];
    }
  }

  delete[] rho12;
  delete[] u;
  delete[] v;
  //*/
}

//------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  if (nrhs < 1 || nrhs > 3)// || nlhs < 1 || nlhs > 3)
    mexErrMsgTxt("Auction algorithm for LSAP with square cost matrix C\n USE: [rho] = jacobiAuctionLSAP(C,forb)");

  // get the type of cost values
  mxClassID tclass = mxGetClassID(prhs[0]);

  switch(tclass)
  {
    case mxINT16_CLASS: jacobiAuctionLSAP_helper<short>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxINT32_CLASS: jacobiAuctionLSAP_helper<int>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxINT64_CLASS: jacobiAuctionLSAP_helper<int64_T>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxSINGLE_CLASS: jacobiAuctionLSAP_helper<float>(nlhs,plhs,nrhs,prhs,tclass); break;
    default: jacobiAuctionLSAP_helper<double>(nlhs,plhs,nrhs,prhs,tclass);
  }
}
