/*
 * @file RandomWalksGraphEditDistance.cpp
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Wed Feb  8 2017
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 * TODO : Optimiser avec eigen
 */

#include "RandomWalksGraphEditDistance.h"
int * RandomWalksGraphEditDistance::labeledKron(int *m1, int nb_rows_m1,int nb_cols_m1,
						int * m2, int nb_rows_m2, int nb_cols_m2,
						int sizeWx[2]){
  int i, j, ia, ib, k, l, ja, jb;
  //std::vector<std::pair<int,int> > mab;
  //------------------------------------------------------------------
  // 1st output : labeled Kronecker product
  // for (i = 0; i < nb_rows_m1; i++) // for each node of Wa
  //   {
  //     ia = nb_rows_m1*i + i;
  //     for (j = 0; j < nb_rows_m2; j++) // for each node of m2
  // 	{
  //         ib = nb_rows_m2*j+j;
  //         if (m1[ia] == m2[ib]) // same label
  // 	    {    mab.push_back(std::make_pair(i,j)); }
  // 	}
  //   }
  // sizeWx = mab.size();
  sizeWx[0] = nb_rows_m1 * nb_rows_m2;
  sizeWx[1] = nb_cols_m1 * nb_cols_m2;
  
  int *Wx = new int[sizeWx[0]*sizeWx[1]]; //mab.size()*mab.size()
  for (i = 0; i < nb_rows_m1; i++) // for each node of m1
    {
      ia = nb_rows_m1*i + i;
      for (j = 0; j < nb_cols_m1; j++)
	{
	  ja = nb_rows_m1*j + j;
	  for (k = 0; k < nb_rows_m2; k++) // for each node of m2
	    {
	      ib = nb_rows_m2*k+k;
	      for (l = 0; l < nb_cols_m2; l++)
		{
		  jb = nb_rows_m2*l+l;
		  if ((m1[ia] == m2[ib]) && (m1[ja] == m2[jb])) // same node label
		    {
		      if ((m1[nb_rows_m1*j+i] == m2[nb_rows_m2*l+k]) && (m1[nb_rows_m1*j+i] > 0) &&
			  (m2[nb_rows_m2*l+k] > 0)) // same edge label
			Wx[nb_rows_m1*nb_rows_m2*(j*nb_cols_m2+l)+(i*nb_rows_m2+k)] = m1[nb_rows_m1*j+i];
		      else Wx[nb_rows_m1*nb_rows_m2*(j*nb_cols_m2+l) +(i*nb_rows_m2+k)] = 0.0;
		    }
		  else Wx[nb_rows_m1*nb_rows_m2*(j*nb_cols_m2+l)+(i*nb_rows_m2+k)] = 0.0;
		}
	    }
	}
    }
  return Wx;
}

void RandomWalksGraphEditDistance::computeCostMatrix(Graph<int,int> * g1,
						     Graph<int,int> * g2){
  int * am_g1 = ((SymbolicGraph *)(g1))->getLabeledAdjacencyMatrix();
  int * am_g2 = ((SymbolicGraph *)(g2))->getLabeledAdjacencyMatrix();
  int sizeWx[2] = {-1,-1};
  int * Wx = labeledKron(am_g1,g1->Size(),g1->Size(),
			 am_g2,g2->Size(),g2->Size(),sizeWx);
  
  int n=g1->Size();
  int m=g2->Size();
  
  Map<MatrixXi> m1(am_g1,n,n);
  Map<MatrixXi> m2(am_g2,m,m);
  Map<MatrixXi> mX(Wx,sizeWx[0],sizeWx[1]);
  // histogram of similar walks
  RowVectorXi ILx = mX.diagonal();
  int nbLab = ILx.maxCoeff();
  MatrixXi mtX(sizeWx[0],sizeWx[1]);
  
  for(int i=0;i<sizeWx[0];i++)
    for(int j=0;j<sizeWx[1];j++)
      mtX(i,j) = (mX(i,j) > 0) && (i != j);

  MatrixXi mtX_pow = mtX;
  for (int k=1 ; k<_k;k++)
    mtX_pow = mtX_pow*mtX;
  
  MatrixXi Hx = histoLab(nbLab,ILx,mtX_pow);

  RowVectorXi IL1 = m1.diagonal();
  MatrixXi mt1(n,n);
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      mt1(i,j) = (m1(i,j) > 0) && (i != j);

  MatrixXi mt1_pow = mt1; 
  for (int k=1 ; k<_k;k++)
    mt1_pow *= mt1;
  MatrixXi H1 = histoLab(nbLab,IL1,mt1_pow);

  RowVectorXi IL2 = m2.diagonal();
  MatrixXi mt2(m,m);
  for(int i=0;i<m;i++)
    for(int j=0;j<m;j++)
      mt2(i,j) = (m2(i,j) > 0) && (i != j);

  MatrixXi mt2_pow = mt2; 
  for (int k=1 ; k<_k;k++)
    mt2_pow *= mt2;
  MatrixXi H2 = histoLab(nbLab,IL2,mt2_pow);
  Hx = Hx.array().sqrt().floor();
  MatrixXi Hi(Hx.rows(), Hx.cols());
  Hi = MatrixXi::Zero(Hx.rows(), Hx.cols());
  int k =0;
  for (int i=0;i<H1.cols();i++){
    for (int j=0;j<H2.cols();j++)
      Hi.col(k+j) = H1.col(i).array() - (H1.col(i).array().min(H2.col(j).array())).min(Hx.col(k+j).array());
    k = k + H2.cols();
  }

  MatrixXi Hj = Hx;
  k = 0;
  for (int i=0;i<H1.cols();i++){
    for (int j=0;j<H2.cols();j++)
      Hj.col(k+j) = H2.col(j).array() - (H1.col(i).array().min(H2.col(j).array())).min(Hx.col(k+j).array());
    k = k + H2.cols();
  }
  
  RowVectorXi S = (Hi.array().min(Hj.array())).colwise().sum();
  RowVectorXi Ri = Hi.colwise().sum() - S;
  RowVectorXi Rj = Hj.colwise().sum()-S;
  RowVectorXi ILxt(ILx.size());
  for(int i = 0;i<ILx.size();i++)
    ILxt(i) = (ILx(i) > 0);

  this->C = new double[(n+1)*(m+1)];
  assert(this->C != 0);
  
  Map<MatrixXd> matrixC(this->C,n+1,m+1);
			
  RowVectorXd C_sub = (((_k-ILxt.array())*cf->cns()+_k*cf->ces()).array()*S.array() +
  		       ((1+_k-ILxt.array())*cf->cns()+_k*cf->ces()).array()*Ri.array().min(Rj.array()) +
  		       ((1+_k-ILxt.array())*cf->cnd()+_k*cf->ced()).array() *(Ri-Rj).array().abs()).cast<double>();
  
  matrixC.block(0,0,n,m) = Map<MatrixXd>(C_sub.data(),m,n).transpose();
  
  RowVectorXd Cie = (((_k+1)*cf->cnd()+_k*cf->ced())* (H1.colwise().sum().array())).cast<double>();
  RowVectorXd Cej = ((_k+1)*cf->cnd()+_k*cf->ced())*H2.colwise().sum().array().cast<double>();

  matrixC.block(0,m,n,1) = Map<MatrixXd>(Cie.data(),n,1);  
  matrixC.block(n,0,1,m) = Map<MatrixXd>(Cej.data(),1,m);  

  matrixC(n,m) = 0.0;


  
  delete [] Wx;
  delete [] am_g1;
  delete [] am_g2;

  
}

MatrixXi RandomWalksGraphEditDistance::histoLab(int nbLab,  RowVectorXi IL, MatrixXi W){
  MatrixXi L(nbLab,W.rows());
  for (int i=0;i<nbLab;i++)
    for(int j=0;j<W.rows();j++)
      L(i,j) = (IL(j)-1== i);
#if DEBUG
  IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
  std::cout << "IL" << std::endl;
  std::cout << IL.format(OctaveFmt) << std::endl;
  std::cout << "L" << std::endl;
  std::cout << L.format(OctaveFmt) << std::endl;
  std::cout << "W" << std::endl;
  std::cout << W.format(OctaveFmt) << std::endl;
  std::cout << "L*W" << std::endl;
  std::cout << (L*W).format(OctaveFmt) << std::endl;
#endif
  return L * W;
    
}
