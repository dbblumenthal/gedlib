/**
 * @file MultipleIPFPGraphEditDistance.h
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version     Jun  9 2017
 *
 */

#ifndef __MULTISTARTMAPPINGREFINEMENT_H__
#define __MULTISTARTMAPPINGREFINEMENT_H__

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <sys/time.h>
#include <list>
#include "MappingRefinement.h"
#include "MappingGenerator.h"


/**
 * @brief A MappingRefinement method which uses a multistart approach
 *
 *   The multistart approach consists in refine several initial mappings
 *   and keep the best one. To this end, the initializations are here lists of
 *   mappings
 */



void printMatrix(double * M, int n, int m);
double * recenteredPsi(double * psi, int n, int m);
double * invertedPsi(double * psi, int n, int m);
double * invertedCenteredPsi(double * psi, int n, int m);
double * copyPsi(double * psi, int n, int m);
std::list<int*> getRandomMappingsPsi( int n, int m, int k, double * psi, double alpha, bool v2  );
std::list<int*> getRandomMappingsInvertedPsi( int n, int m, int k, double * psi, double alpha, bool v2,  std::list<int*> visitedMappings, double explorationFactor=0  );
std::list<int*> getRandomMappingsInvertedCenteredPsi( int n, int m, int k, double * psi, double alpha, bool v2, std::list<int*> visitedMappings);
std::list<int*> getKLSAPEMappings( int n,  int m, double* C, const int& k);
void printMapping(int* mapping, int n);
void printMappings(std::list<int*> mappings, int n);
std::list<int*> getRandomMappingsMixtPsi( int n, int m, int k, double * psiPre, double * psiPost, double alpha, bool v2  );
bool mappingIsEqual(int* map1, int * map2, int n);
bool mappingIsVisited(int* map, std::list<int*> visitedMappings, int n, int m);


bool mappingIsVisited(int* map, std::list<int*> visitedMappings, int n, int m){
	bool isVisited = 0;
	bool isEqual =1;
	int i;
	std::list<int*>::iterator it = visitedMappings.begin();
	while ((isVisited==0)&&(it!=visitedMappings.end())){
             int* localMap=*it;
	     i=0;
	     isEqual=1;
	     while((isEqual==1)&&(i<n)){
		  if ((map[i]!=localMap[i]) && !((map[i]>=m) && (localMap[i]>=m))) isEqual = 0;
		  i++;
	     }
	     if (isEqual == 1) isVisited =1;
	     it++;
	     }
	return isVisited;
}

bool mappingIsEqual(int* map1, int * map2, int n){
bool result =1;
int i=0;
while((result==1)&&(i<n)){
if (map1[i]!=map2[i]) result = 0;
i++;
}
return result;

}

std::list<int*> getRandomMappingsMixtPsi( int n, int m, int k, double * psiPre, double * psiPost, double alpha, bool v2  )
{
 // this method generates k random mapping where pairwise mappings that have a high value in psi are
 // more likely to be picked

  if (k < 0) k = 100;
  //int numEpsMax = (n>m) ? n-m : 0;
  //int numEpsAssigned;
  //bool allEpsAssigned;

  double * localPsi(0);

  double * incrLinePsi = new double[(n+1)*(m+1)];
  double * incrColumnPsi = new double[(n+1)*(m+1)];
  double X; // Variable aléatoire
  int maxPsi;
  std::list<int*> mappings;
  for (int i=0; i<k; i++){
    //allEpsAssigned = 0;
    //numEpsAssigned = 0;
    if (i<=k/2){
		localPsi = invertedPsi(psiPre,n,m); }
        else{
                localPsi = copyPsi(psiPost,n,m);}
    maxPsi=0;
    for (int j=0;j<sub2ind(n,m,n);j++)if (localPsi[j]>maxPsi) maxPsi = localPsi[j];
    int* _map = new int[n];

          // initialisation localPsi
	    	for (int j=0;j<n;j++){
		    //if ((allEpsAssigned == 0) && (numEpsAssigned==numEpsMax)){
		//	allEpsAssigned =1;
		//	for(int l=j;l<n;l++) localPsi[sub2ind(l,m,n+1)]=-maxPsi;
		        //std::cout<<"All Epsilon assigned"<< std::endl;
		  //  }
		    for(int h=0;h<=m;h++){
		         if (h==0) incrLinePsi[sub2ind(j,h,n+1)]=localPsi[sub2ind(j,h,n+1)]+maxPsi;
		         else incrLinePsi[sub2ind(j,h,n+1)]=incrLinePsi[sub2ind(j,h-1,n+1)]+localPsi[sub2ind(j,h,n+1)]+maxPsi;
		    }
                    incrLinePsi[sub2ind(j,m,n+1)]*=1+alpha;
		    X=((double)std::rand()/(double)RAND_MAX)*incrLinePsi[sub2ind(j,m,n+1)];
		    int h=0;
		    if (X!=0){

		    	while((X>incrLinePsi[sub2ind(j,h,n+1)])&&(h<m)) h++;
		    	_map[j]=h;
		    	if ((h<m)&&(j<n-1)) for(int l=j+1;l<n;l++) localPsi[sub2ind(l,h,n+1)]=-maxPsi;
		    	}
		      //  if (h==m) numEpsAssigned++;
	    	    else {
			_map[j]=m;
			}

		    //std::cout<<"Matrice IncrLine avec X pour Ligne"<<j << " = " << X <<  std::endl;
		    //std::cout<<"LocalPsi=" << std::endl;
		    //printMatrix(localPsi,n+1,m+1);
		    //std::cout<<"incrLinePsi=" << std::endl;
		    //printMatrix(incrLinePsi,n+1,m+1);
		    //std::cout<<"ligne: sommet "<<j<< " sur " << n << " mappé à " << h << std::endl;
    }
   /*
    else{
       for (int h=0;h<m;h++){
            for(int j=0;j<=n;j++){
                 if (j==0) incrColumnPsi[sub2ind(j,h,n)]=localPsi[sub2ind(j,h,n)];
                 else incrColumnPsi[sub2ind(j,h,n)]=incrColumnPsi[sub2ind(j-1,h,n)]+localPsi[sub2ind(j,h,n)];
            }
           X=((double)std::rand()/(double)RAND_MAX)*incrColumnPsi[sub2ind(n,h,n)];
           int j=0;
             if (X!=0){
             	 while((X>incrColumnPsi[sub2ind(j,h,n)])&&(j<n)) j++;
            	_map[j]=h;
            	if ((j<n)&&(h<m-1))for(int c=h+1;c<m;c++) localPsi[sub2ind(j,c,n)]=0;
                }
	     else{
		j=n;
   		 }
            //std::cout<<"colonne: sommet "<<h<< " sur " << m << " mappé à " << j << std::endl;
            std::cout<<"Matrice IncrColmun avec X pour Colonne "<< h <<" = "<< X <<  std::endl;
            printMatrix(incrColumnPsi,n,m);

         }


    }
*/
    //complete the n _map to n+m rhoperm

     int* rhoperm = new int[n+m];
    bool* epsAssign = new bool[n]; // is eps[i] assigned
    bool* epsAssignG2 = new bool[m];
    for (int j=0; j<n; j++) epsAssign[j] = false;
    for (int j=0; j<m; j++) epsAssignG2[j] = true;
    for (int jj=0; jj<n; jj++){
	 if (_map[jj] < m){
		rhoperm[jj] = _map[jj];
                epsAssignG2[_map[jj]] = false;}
	 else{
	         rhoperm[jj] = jj+n;
		 epsAssign[jj] = true;
	 }
     }



     int firstEpsNonAssign = 0;
     for (int j=0; j<m; j++){
          if (epsAssignG2[j] == true)
               rhoperm[j+n] = j;

          else{ // find the first epsilon not assigned
              while (firstEpsNonAssign < n && epsAssign[firstEpsNonAssign]) firstEpsNonAssign++;
              rhoperm[j+n] = firstEpsNonAssign + n;
              epsAssign[firstEpsNonAssign] = true;
          }
      }

   // printMapping(rhoperm,n);
    mappings.push_back(rhoperm);

    //std::cout<<"mapping "<< i << " ok, list-mappings- contains " << mappings.size() << "mappings."<< std::endl;
    delete[] _map;
    delete[] epsAssign;
    delete[] epsAssignG2;
    delete[] localPsi;
  }

//

   //delete [] epsAssign;

  delete[] incrLinePsi;
  delete[] incrColumnPsi;

  return mappings;
}

std::list<int*> getRandomMappingsInvertedCenteredPsi( int n, int m, int k, double * psi, double alpha, bool v2, std::list<int*> visitedMappings )
{
 // this method generates k random mapping where pairwise mappings that have a high value in psi are
 // more likely to be picked
  if (alpha == 1) alpha = 0.999 ;
  if (k < 0) k = 100;
   int maxIt = n*m;
   double deletionsProportion = 0;
   if (n>m) deletionsProportion = (n-m)/n ; else deletionsProportion = 0;
  double * localPsi(0);

  double * incrLinePsi = new double[(n+1)*(m+1)];
  double * incrColumnPsi = new double[(n+1)*(m+1)];
  double X; // Variable aléatoire
  int maxPsi;
  int iprec=0;
  int numVisitedComputed = 0;
  std::list<int*> mappings;
  for (int i=0; i<k; i++){
    if (iprec!=i) numVisitedComputed = 0;
    iprec=i;
    localPsi = invertedCenteredPsi(psi,n,m); // initialisation localPsi
    maxPsi=0;
    for (int j=0;j<sub2ind(n,m,n);j++)if (localPsi[j]>maxPsi) maxPsi = localPsi[j];
    int* _map = new int[n];
    	for (int j=0;j<n;j++){
            for(int h=0;h<=m;h++){
                 if (h==0) incrLinePsi[sub2ind(j,h,n+1)]=localPsi[sub2ind(j,h,n+1)]+maxPsi*(numVisitedComputed/maxIt)*(localPsi[sub2ind(j,h,n+1)]!=0);
                 else incrLinePsi[sub2ind(j,h,n+1)]=incrLinePsi[sub2ind(j,h-1,n+1)]+localPsi[sub2ind(j,h,n+1)]+maxPsi*(numVisitedComputed/maxIt)*(localPsi[sub2ind(j,h,n+1)]!=0);
            }

            incrLinePsi[sub2ind(j,m,n+1)]*=1/(1-alpha*(1-numVisitedComputed/maxIt));
            //for(int h=0;h<=m-1;h++){
            //    incrLinePsi[sub2ind(j,h,n+1)]=incrLinePsi[sub2ind(j,h,n+1)]*(1-alpha)*((maxIt-numVisitedComputed)/maxIt)*(numdeletions/k);
            //}
            //changes probabilities such that the probability to delete a vertex is exactly alpha.
	    X=((double)std::rand()/(double)RAND_MAX)*incrLinePsi[sub2ind(j,m,n+1)];
            int h=0;
            if (X!=0){

            	while((X>incrLinePsi[sub2ind(j,h,n+1)])&&(h<m)) h++;
            	_map[j]=h;
            	if ((h<m)&&(j<n-1)) for(int l=j+1;l<n;l++) localPsi[sub2ind(l,h,n+1)]=0;
            	}
            //    if (h==m) numEpsAssigned++;
    	    else {
		_map[j]=m;
		}
	   /*
	   std::cout<<"Matrice IncrLine avec X pour Ligne"<<j << " = " << X <<  std::endl;
           std::cout<<"LocalPsi=" << std::endl;
           printMatrix(localPsi,n+1,m+1);
           std::cout<<"incrLinePsi=" << std::endl;
	   printMatrix(incrLinePsi,n+1,m+1);
           std::cout<<"ligne: sommet "<<j<< " sur " << n << " mappé à " << h << std::endl;
        */}
    /*
    else{
       for (int h=0;h<m;h++){
            for(int j=0;j<=n;j++){
                 if (j==0) incrColumnPsi[sub2ind(j,h,n)]=localPsi[sub2ind(j,h,n)];
                 else incrColumnPsi[sub2ind(j,h,n)]=incrColumnPsi[sub2ind(j-1,h,n)]+localPsi[sub2ind(j,h,n)];
            }
           X=((double)std::rand()/(double)RAND_MAX)*incrColumnPsi[sub2ind(n,h,n)];
           int j=0;
             if (X!=0){
             	 while((X>incrColumnPsi[sub2ind(j,h,n)])&&(j<n)) j++;
            	_map[j]=h;
            	if ((j<n)&&(h<m-1))for(int c=h+1;c<m;c++) localPsi[sub2ind(j,c,n)]=0;
                }
	     else{
		j=n;
   		 }
            //std::cout<<"colonne: sommet "<<h<< " sur " << m << " mappé à " << j << std::endl;
            std::cout<<"Matrice IncrColmun avec X pour Colonne "<< h <<" = "<< X <<  std::endl;
            printMatrix(incrColumnPsi,n,m);

         }


    }
*/

 if((v2==1) && mappingIsVisited(_map,visitedMappings,n,m) && (numVisitedComputed <maxIt)){


    i--;
    numVisitedComputed++;
    //std::cout << "dejavu" << std::endl;

 }
 else{

    //complete the n _map to n+m rhoperm
    int* rhoperm = new int[n+m];
    bool* epsAssign = new bool[n]; // is eps[i] assigned
    bool* epsAssignG2 = new bool[m];
    for (int j=0; j<n; j++) epsAssign[j] = false;
    for (int j=0; j<m; j++) epsAssignG2[j] = true;
    for (int jj=0; jj<n; jj++){
	 if (_map[jj] < m){
		rhoperm[jj] = _map[jj];
                epsAssignG2[_map[jj]] = false;}
	 else{
	         rhoperm[jj] = jj+n;
		 epsAssign[jj] = true;
	 }
     }



     int firstEpsNonAssign = 0;
     for (int j=0; j<m; j++){
          if (epsAssignG2[j] == true)
               rhoperm[j+n] = j;

          else{ // find the first epsilon not assigned
              while (firstEpsNonAssign < n && epsAssign[firstEpsNonAssign]) firstEpsNonAssign++;
              rhoperm[j+n] = firstEpsNonAssign + n;
              epsAssign[firstEpsNonAssign] = true;
          }
      }

   // printMapping(rhoperm,n);
    mappings.push_back(rhoperm);
    visitedMappings.push_back(rhoperm);

    //std::cout<<"mapping "<< i << " ok, list-mappings- contains " << mappings.size() << "mappings."<< std::endl;

    delete[] epsAssign;
    delete[] epsAssignG2;
  }


    delete[] localPsi;
    delete[] _map;
  }

//

   //delete [] epsAssign;

  delete[] incrLinePsi;
  delete[] incrColumnPsi;

  return mappings;
}

std::list<int*> getRandomMappingsInvertedPsi( int n, int m, int k, double * psi, double alpha, bool v2, std::list<int*> visitedMappings, double exploFactor )
{
 // this method generates k random mapping where pairwise mappings that have a high value in psi are
 // more likely to be picked
  if (alpha == 1) alpha = 0.999 ;
  if (k < 0) k = 100;
   int maxIt = n*m;
   double deletionsProportion = 0;
   if (n>m) deletionsProportion = (n-m)/n ; else deletionsProportion = 0;

  double * localPsi(0);

  double * incrLinePsi = new double[(n+1)*(m+1)];
  double * incrColumnPsi = new double[(n+1)*(m+1)];
  double X; // Variable aléatoire
  int maxPsi;
  int iprec=0;
  int numVisitedComputed = 0;
  std::list<int*> mappings;
  for (int i=0; i<k; i++){
    if (iprec!=i) numVisitedComputed = 0;
    iprec=i;
    localPsi = copyPsi(psi,n,m); // initialisation localPsi
    maxPsi=0;
    for (int j=0;j<sub2ind(n,m,n);j++)if (localPsi[j]>maxPsi) maxPsi = localPsi[j];
    int* _map = new int[n];

    	for (int j=0;j<n;j++){

            for(int h=0;h<=m;h++){
                 if (h==0) incrLinePsi[sub2ind(j,h,n+1)]=maxPsi*exploFactor+localPsi[sub2ind(j,h,n+1)]*(1-exploFactor)+maxPsi*numVisitedComputed/maxIt*(localPsi[sub2ind(j,h,n+1)]!=0);
                 else incrLinePsi[sub2ind(j,h,n+1)]=incrLinePsi[sub2ind(j,h-1,n+1)]+maxPsi*exploFactor+localPsi[sub2ind(j,h,n+1)]*(1-exploFactor)+maxPsi*numVisitedComputed/maxIt*(localPsi[sub2ind(j,h,n+1)]!=0);
            }
             incrLinePsi[sub2ind(j,m,n+1)]*=1/(1-alpha*(1-numVisitedComputed/maxIt));
            //for(int h=0;h<=m-1;h++){
            //    incrLinePsi[sub2ind(j,h,n+1)]=incrLinePsi[sub2ind(j,h,n+1)]*(1-alpha)*((maxIt-numVisitedComputed)/maxIt)*(numdeletions/k);
            //}
            //changes probabilities such that the probability to delete a vertex is exactly alpha.
	    X=((double)std::rand()/(double)RAND_MAX)*incrLinePsi[sub2ind(j,m,n+1)];
            int h=0;
            if (X!=0){

            	while((X>incrLinePsi[sub2ind(j,h,n+1)])&&(h<m)) h++;
            	_map[j]=h;
            	if ((h<m)&&(j<n-1)) for(int l=j+1;l<n;l++) localPsi[sub2ind(l,h,n+1)]=0;
            	}
            //    if (h==m) numEpsAssigned++;
    	    else {
		_map[j]=m;
		}
	   /*
	   std::cout<<"Matrice IncrLine avec X pour Ligne"<<j << " = " << X <<  std::endl;
           std::cout<<"LocalPsi=" << std::endl;
           printMatrix(localPsi,n+1,m+1);
           std::cout<<"incrLinePsi=" << std::endl;
	   printMatrix(incrLinePsi,n+1,m+1);
           std::cout<<"ligne: sommet "<<j<< " sur " << n << " mappé à " << h << std::endl;
        */}
    /*
    else{
       for (int h=0;h<m;h++){
            for(int j=0;j<=n;j++){
                 if (j==0) incrColumnPsi[sub2ind(j,h,n)]=localPsi[sub2ind(j,h,n)];
                 else incrColumnPsi[sub2ind(j,h,n)]=incrColumnPsi[sub2ind(j-1,h,n)]+localPsi[sub2ind(j,h,n)];
            }
           X=((double)std::rand()/(double)RAND_MAX)*incrColumnPsi[sub2ind(n,h,n)];
           int j=0;
             if (X!=0){
             	 while((X>incrColumnPsi[sub2ind(j,h,n)])&&(j<n)) j++;
            	_map[j]=h;
            	if ((j<n)&&(h<m-1))for(int c=h+1;c<m;c++) localPsi[sub2ind(j,c,n)]=0;
                }
	     else{
		j=n;
   		 }
            //std::cout<<"colonne: sommet "<<h<< " sur " << m << " mappé à " << j << std::endl;
            std::cout<<"Matrice IncrColmun avec X pour Colonne "<< h <<" = "<< X <<  std::endl;
            printMatrix(incrColumnPsi,n,m);

         }


    }
*/

 if((v2==1) && mappingIsVisited(_map,visitedMappings,n,m) && (numVisitedComputed <maxIt)){


    i--;
    numVisitedComputed++;
    //std::cout << "dejavu" << std::endl;

 }
 else{

    //complete the n _map to n+m rhoperm
    int* rhoperm = new int[n+m];
    bool* epsAssign = new bool[n]; // is eps[i] assigned
    bool* epsAssignG2 = new bool[m];
    for (int j=0; j<n; j++) epsAssign[j] = false;
    for (int j=0; j<m; j++) epsAssignG2[j] = true;
    for (int jj=0; jj<n; jj++){
	 if (_map[jj] < m){
		rhoperm[jj] = _map[jj];
                epsAssignG2[_map[jj]] = false;}
	 else{
	         rhoperm[jj] = jj+n;
		 epsAssign[jj] = true;
	 }
     }



     int firstEpsNonAssign = 0;
     for (int j=0; j<m; j++){
          if (epsAssignG2[j] == true)
               rhoperm[j+n] = j;

          else{ // find the first epsilon not assigned
              while (firstEpsNonAssign < n && epsAssign[firstEpsNonAssign]) firstEpsNonAssign++;
              rhoperm[j+n] = firstEpsNonAssign + n;
              epsAssign[firstEpsNonAssign] = true;
          }
      }

   // printMapping(rhoperm,n);
    mappings.push_back(rhoperm);
    visitedMappings.push_back(rhoperm);

    //std::cout<<"mapping "<< i << " ok, list-mappings- contains " << mappings.size() << "mappings."<< std::endl;

    delete[] epsAssign;
    delete[] epsAssignG2;
  }


    delete[] localPsi;
    delete[] _map;
  }

//

   //delete [] epsAssign;

  delete[] incrLinePsi;
  delete[] incrColumnPsi;

  return mappings;
}

std::list<int*> getRandomMappingsPsi( int n, int m, int k, double * psi, double alpha, bool v2 )
{
 // this method generates k random mapping where pairwise mappings that have a high value in psi are
 // less likely to be picked

  if (k < 0) k = 100;
  //int numEpsMax = (n>m) ? n-m : 0;
  //int numEpsAssigned;
  //bool allEpsAssigned;

  double * localPsi(0);

  double * incrLinePsi = new double[(n+1)*(m+1)];
  double * incrColumnPsi = new double[(n+1)*(m+1)];
  double X; // Variable aléatoire
  int maxPsi;
  std::list<int*> mappings;
  for (int i=0; i<k; i++){
    //allEpsAssigned = 0;
    //numEpsAssigned = 0;
    localPsi = invertedPsi(psi,n,m); // initialisation localPsi
    maxPsi=0;
    for (int j=0;j<sub2ind(n,m,n);j++)if (localPsi[j]>maxPsi) maxPsi = localPsi[j];
    int* _map = new int[n];

    	for (int j=0;j<n;j++){
	//    if ((allEpsAssigned == 0) && (numEpsAssigned==numEpsMax)){
	//	allEpsAssigned =1;
	//	for(int l=j;l<n;l++) localPsi[sub2ind(l,m,n+1)]=-maxPsi;
                //std::cout<<"All Epsilon assigned"<< std::endl;
         //   }
            for(int h=0;h<=m;h++){
                 if (h==0) incrLinePsi[sub2ind(j,h,n+1)]=localPsi[sub2ind(j,h,n+1)]+maxPsi;
                 else incrLinePsi[sub2ind(j,h,n+1)]=incrLinePsi[sub2ind(j,h-1,n+1)]+localPsi[sub2ind(j,h,n+1)]+maxPsi;
            }
            incrLinePsi[sub2ind(j,m,n+1)]*=1+alpha;
	    X=((double)std::rand()/(double)RAND_MAX)*incrLinePsi[sub2ind(j,m,n+1)];
            int h=0;
            if (X!=0){

            	while((X>incrLinePsi[sub2ind(j,h,n+1)])&&(h<m)) h++;
            	_map[j]=h;
            	if ((h<m)&&(j<n-1)) for(int l=j+1;l<n;l++) localPsi[sub2ind(l,h,n+1)]=-maxPsi;
            	}
           //     if (h==m) numEpsAssigned++;
    	    else {
		_map[j]=m;
		}
/*
	     std::cout<<"n= " << n << " , m= " << m << std::endl << "Matrix psi" << std::endl;
             printMatrix(psi, n+1, m+1);
	    std::cout<<"Matrice IncrLine avec X pour Ligne"<<j << " = " << X <<  std::endl;
            std::cout<<"LocalPsi=" << std::endl;
            printMatrix(localPsi,n+1,m+1);
           std::cout<<"incrLinePsi=" << std::endl;
	    printMatrix(incrLinePsi,n+1,m+1);
            std::cout<<"ligne: sommet "<<j<< " sur " << n << " mappé à " << h << std::endl;
*/
        }
    /*
    else{
       for (int h=0;h<m;h++){
            for(int j=0;j<=n;j++){
                 if (j==0) incrColumnPsi[sub2ind(j,h,n)]=localPsi[sub2ind(j,h,n)];
                 else incrColumnPsi[sub2ind(j,h,n)]=incrColumnPsi[sub2ind(j-1,h,n)]+localPsi[sub2ind(j,h,n)];
            }
           X=((double)std::rand()/(double)RAND_MAX)*incrColumnPsi[sub2ind(n,h,n)];
           int j=0;
             if (X!=0){
             	 while((X>incrColumnPsi[sub2ind(j,h,n)])&&(j<n)) j++;
            	_map[j]=h;
            	if ((j<n)&&(h<m-1))for(int c=h+1;c<m;c++) localPsi[sub2ind(j,c,n)]=0;
                }
	     else{
		j=n;
   		 }
            //std::cout<<"colonne: sommet "<<h<< " sur " << m << " mappé à " << j << std::endl;
            std::cout<<"Matrice IncrColmun avec X pour Colonne "<< h <<" = "<< X <<  std::endl;
            printMatrix(incrColumnPsi,n,m);

         }


    }
*/
    //complete the n _map to n+m rhoperm

    int* rhoperm = new int[n+m];
    bool* epsAssign = new bool[n]; // is eps[i] assigned
    bool* epsAssignG2 = new bool[m];
    for (int j=0; j<n; j++) epsAssign[j] = false;
    for (int j=0; j<m; j++) epsAssignG2[j] = true;
    for (int jj=0; jj<n; jj++){
	 if (_map[jj] < m){
		rhoperm[jj] = _map[jj];
                epsAssignG2[_map[jj]] = false;}
	 else{
	         rhoperm[jj] = jj+n;
		 epsAssign[jj] = true;
	 }
     }



     int firstEpsNonAssign = 0;
     for (int j=0; j<m; j++){
          if (epsAssignG2[j] == true)
               rhoperm[j+n] = j;

          else{ // find the first epsilon not assigned
              while (firstEpsNonAssign < n && epsAssign[firstEpsNonAssign]) firstEpsNonAssign++;
              rhoperm[j+n] = firstEpsNonAssign + n;
              epsAssign[firstEpsNonAssign] = true;
          }
      }

   // printMapping(rhoperm,n);
    mappings.push_back(rhoperm);

    //std::cout<<"mapping "<< i << " ok, list-mappings- contains " << mappings.size() << "mappings."<< std::endl;
    delete[] _map;
    delete[] epsAssign;
    delete[] epsAssignG2;
    delete[] localPsi;
  }

//

   //delete [] epsAssign;

  delete[] incrLinePsi;
  delete[] incrColumnPsi;

  return mappings;
}



void printMapping(int* mapping, int n){
for (int i=0; i<n;i++) std::cout<<i<< "-"<< mapping[i] << " , " ;
std::cout<< std::endl;
}

void printMappings(std::list<int*> mappings, int n){
for(std::list<int*>::iterator iter = mappings.begin(); iter != mappings.end(); iter++){
   printMapping(*iter,n);
}
std::cout<< std::endl;
}

void printMatrix(double * M, int n, int m){
for (int i=0; i<n;i++){
   	for (int j=0; j<m;j++) std::cout << M[sub2ind(i,j,n)] << " " ;
	std::cout<< std::endl;
}
}




double * recenteredPsi(double * psi, int n, int m){
double * rPsi = new double [(n+1)*(m+1)];
double psimax = -1;
for ( int i=0;i<(n+1)*(m+1);i++)
    if (psi[i]>psimax) psimax = psi[i];
for (int i=0;i<(n+1)*(m+1);i++)
    rPsi[i] = psimax/2 - std::abs(psimax/2 - psi[i]);
return rPsi;

}

double * invertedCenteredPsi(double * psi, int n, int m){
double * rPsi = new double [(n+1)*(m+1)];
double psimax = -1;
for ( int i=0;i<(n+1)*(m+1);i++)
    if (psi[i]>psimax) psimax = psi[i];
for (int i=0;i<(n+1)*(m+1);i++)
    rPsi[i] = std::abs(psimax/2 - psi[i]);
return rPsi;

}

double * invertedPsi(double * psi, int n, int m){
double * rPsi = new double [(n+1)*(m+1)];
double psimax = -1;
for ( int i=0;i<(n+1)*(m+1);i++)
    if (psi[i]>psimax) psimax = psi[i];
for (int i=0;i<(n+1)*(m+1);i++)
    rPsi[i] = psimax - psi[i];
return rPsi;
}



double * copyPsi(double * psi, int n, int m){
double * rPsi = new double [(n+1)*(m+1)];
for (int i=0;i<(n+1)*(m+1);i++)
    rPsi[i] = psi[i];
return rPsi;
}

std::list<int*> getKLSAPEMappings( int n,
                     int m,
                     double* C,     const int& k)
{



// Compute cost Matrix that LSAPE -> LSAP

double* Clsap = new  double[(n+m)*(n+m)];
  //memset(Clsap, -1.0, sizeof(double)*(n+m)*(n+m)); // inf costs to all non-possible mappings
  for (int j=0; j<m+n; j++)
    for (int i=0; i<m+n; i++)
      if (i>=n && j>=m) Clsap[sub2ind(i,j,n+m)] = 0;
      else Clsap[sub2ind(i,j,n+m)] = -1.0;

  //XXX changer column first
  for(int i=0;i<n;i++)
    for(int j=0;j<m;j++)
      Clsap[sub2ind(i,j,n+m)] = C[sub2ind(i,j,n+1)];

   for(int i=0;i<n;i++)
     Clsap[sub2ind(i,m+i,n+m)] = C[sub2ind(i,m,n+1)];
   for(int j=0;j<m;j++)
     Clsap[sub2ind(n+j,j,n+m)] = C[sub2ind(n,j,n+1)];


// Compute an optimal assignement
  double *u = new double[n+1];
  double *v = new double[m+1];
  int* G1_to_G2 = new int[n+1];
  int* G2_to_G1 = new int[m+1];



  //hungarianLSAP<double,int>(this->_Clsap, n+m, m+n, G1_to_G2, u, v, true);
  hungarianLSAPE(C, n+1, m+1, G1_to_G2, G2_to_G1, u, v, false);



  // Compute LSAP solution from LSAPE
  int* rhoperm = new int[n+m];
  bool* epsAssign = new bool[n]; // is eps[i] assigned
  for (int j=0; j<n; j++) epsAssign[j] = false;
  for (int i=0; i<n; i++){
    if (G1_to_G2[i] < m)
      rhoperm[i] = G1_to_G2[i];
    else{
      rhoperm[i] = i+m;
      epsAssign[i] = true;
    }
  }
  int firstEpsNonAssign = 0;
  for (int j=0; j<m; j++){
    if (G2_to_G1[j] == n)
      rhoperm[j+n] = j;

    else{ // find the first epsilon not assigned
      while (firstEpsNonAssign < n && epsAssign[firstEpsNonAssign]) firstEpsNonAssign++;
      rhoperm[j+n] = firstEpsNonAssign + m;
      epsAssign[firstEpsNonAssign] = true;
    }
  }

  // Compute LSAP u and v
  double *lu = new double[n+m];
  double *lv = new double[n+m];
  for (int i=0; i<n; i++) lu[i] = u[i];
  for (int i=n; i<n+m; i++) lu[i] = 0;
  for (int j=0; j<m; j++) lv[j] = v[j];
  for (int j=m; j<n+m; j++) lv[j] = 0;



  // Compute the k optimal mappings
  cDigraph<int> edg = equalityDigraph<double,int> (Clsap, n+m, n+m, rhoperm, lu, lv);
  AllPerfectMatchingsEC<int> apm(edg, n, m);
  apm.enumPerfectMatchings(edg,k);


  delete [] epsAssign;
  delete [] u;
  delete [] v;
  delete [] lu;
  delete [] lv;
  delete [] G2_to_G1;
  delete [] G1_to_G2;
  delete [] Clsap;

  return apm.getPerfectMatchings();


}

template<class NodeAttribute, class EdgeAttribute>
class MultistartMappingRefinement
{

protected:

  MappingGenerator<NodeAttribute, EdgeAttribute> * initGen; //!< Generator of initializations
  // MultiGed<NodeAttribute, EdgeAttribute> * secondGen; // Generator of mappins that uses _psi to explore new areas of the polytope.
  int _nbIterations; // number of times the matrix _exploredMappingsSum will be updated and used to generate new initial solutions using secondGen.
  int k; //!< Number of initial mapping to generate from \ref initGen
  int kr;//number of initial mappings to return when kr have converged among the k launched
  double adap; //adaptative factor 0=not adaptative, 1=strongly adaptative
  std::list<int*> refinedMappings; //!< The last set of refined mappings
  int _numBestMap;
  double* _psiPre;
  double* _psiPost;
  int _psiMethod;
  double _alpha;
  std::default_random_engine randGen;

public:


virtual std::list<int*> getInitMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
	                                 Graph<NodeAttribute,EdgeAttribute> * g2,
	                                 int k );

  /**
   * @brief Outputs the best mapping refined by \ref algorithm from initializations given by \ref initGen
   *
   *  During the process, only the best mapping is kept in memory
   *  allowing the procedure to require less memory than \ref getBetterMappings.
   *
   * @param  algorithm   the refinement method
   * @param  g1          First graph
   * @param  g2          Second graph
   * @param  G1_to_G2    forward output mapping
   * @param  G2_to_G1    reverse output mapping, useful for the graph edit distance
   * @see getBestMappingFromSet
   */
  virtual void getBestMapping( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                               Graph<NodeAttribute,EdgeAttribute> * g1,
                               Graph<NodeAttribute,EdgeAttribute> * g2,
                               int * G1_to_G2, int * G2_to_G1 );

  /**
   * @brief Returns the list of refined mappings from initializations given by \reg initGen
   *
   * @param  algorithm   the refinement method
   * @param  g1          First graph
   * @param  g2          Second graph
   * @see getBetterMappingsFromSet
   */
  virtual const std::list<int*>&
  getBetterMappings( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                     Graph<NodeAttribute,EdgeAttribute> * g1,
                     Graph<NodeAttribute,EdgeAttribute> * g2 );

  /**
   * @brief Returns the best mappings refined by \ref algorithm from the given \ref mappings
   */
  virtual void getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                                      Graph<NodeAttribute,EdgeAttribute> * g1,
                                      Graph<NodeAttribute,EdgeAttribute> * g2,
                                      int * G1_to_G2, int * G2_to_G1,
                                      std::list<int*>& mappings, double * _psiPre, double * _psiPost, 	double * cost);

  virtual void getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                                      Graph<NodeAttribute,EdgeAttribute> * g1,
                                      Graph<NodeAttribute,EdgeAttribute> * g2,
                                      int * G1_to_G2, int * G2_to_G1,
                                      std::list<int*>& mappings, double * _psiPre, double * _psiPost, 	double * cost, int kr);

  virtual void getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                                      Graph<NodeAttribute,EdgeAttribute> * g1,
                                      Graph<NodeAttribute,EdgeAttribute> * g2,
                                      int * G1_to_G2, int * G2_to_G1,
                                      std::list<int*>& mappings, double * _psiPre, double * _psiPost, 	double * cost, int kr,double adap);


  /**
   * @brief Returns the list of refined mappings from the given \ref mappings
   *
   * @param  algorithm   the refinement method
   * @param  g1          First graph
   * @param  g2          Second graph
   * @param  mapping     a list of arrays representing initial mappings
   * @note   mappings are allocated on the heap and memory management is left to the user
   * @see getBetterMappings
   */
  virtual const std::list<int*>&
  getBetterMappingsFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                            Graph<NodeAttribute,EdgeAttribute> * g1,
                            Graph<NodeAttribute,EdgeAttribute> * g2,
                            std::list<int*>& mappings );


// constructeur avec la meme signature que la version base de la classe
MultistartMappingRefinement( MappingGenerator<NodeAttribute, EdgeAttribute> * gen,
                               int nSol
                             ):
    initGen(gen),
    k(nSol),
     kr(nSol),
     adap(0),
    _numBestMap(-1),
    _psiPre(0),
    _psiPost(0)
   {
        randGen.seed(123);
        //std::cout << "Enter number of iterations: ";
       //std::cin >> _nbIterations;
       //std::cout << "psi generation: 0=flat 1=recentered ";
       //std::cin >> _psiMethod;
       //valeur de alpha, proportion de noeuds supprimés
       //std::cin >> _alpha;
	_nbIterations = 0;
	_psiMethod = 12;
	_alpha = 0.0;


  }


  MultistartMappingRefinement( MappingGenerator<NodeAttribute, EdgeAttribute> * gen,
                               int nSol, int nReturned
                             ):
  initGen(gen),
    k(nSol),
    kr(nReturned),
    adap(0),
    _numBestMap(-1),
    _psiPre(0),
    _psiPost(0)


  {
       randGen.seed(123);
        //std::cout << "Enter number of iterations: ";
       //std::cin >> _nbIterations;
       //std::cout << "psi generation: 0=flat 1=recentered ";
       //std::cin >> _psiMethod;
       //valeur de alpha, proportion de noeuds supprimés
       //std::cin >> _alpha;
	_nbIterations = 0;
	_psiMethod = 12;
	_alpha = 0.0;

  }
   MultistartMappingRefinement( MappingGenerator<NodeAttribute, EdgeAttribute> * gen,
                               int nSol, int nReturned, double adapFactor
                             ):
  initGen(gen),
    k(nSol),
    kr(nReturned),
    adap(adapFactor),
    _numBestMap(-1),
    _psiPre(0),
    _psiPost(0)


  {
       randGen.seed(123);
        //std::cout << "Enter number of iterations: ";
       //std::cin >> _nbIterations;
       //std::cout << "psi generation: 0=flat 1=recentered ";
       //std::cin >> _psiMethod;
       //valeur de alpha, proportion de noeuds supprimés
       //std::cin >> _alpha;
	_nbIterations = 0;
	_psiMethod = 12;
	_alpha = 0.0;

  }

  virtual ~MultistartMappingRefinement(){


  }


};

//---
template<class NodeAttribute, class EdgeAttribute>
std::list<int*> MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getInitMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
	     Graph<NodeAttribute,EdgeAttribute> * g2,
	     int k )
{

double alpha = this->_alpha;

   if (k < 0) k = 100;
     int n = g1->Size();
     int m = g2->Size();
     std::list<int*> mappings;
  if (alpha==-1){




	  for (int i=0; i<k; i++){
	    int* _map = new int[n+m];
	    for (int a=0; a<n+m; a++) _map[a] = a;

	    std::shuffle(&_map[0], &_map[n+m], this->randGen);
	    mappings.push_back(_map);
	  }
  }
  else{

	  //int mx = std::max(n,m);
          int mx = n;
          int alphaNM=alpha*(n+m-mx);

	  for (int i=0; i<k; i++){
	    int* _map = new int[mx];
	    for (int a=0; a<mx; a++)  _map[a] = a;

	    std::shuffle(&_map[0], &_map[mx], this->randGen);



	    int* mapping = new int[n+m];
            for(int j=0;j<alphaNM;j++) mapping[j]=mx+j; // les alphaNM premiers sommets sont supprimés
	    for (int j=0;j<mx;j++) mapping[j+alphaNM]=_map[j];// les mx sommets suivants sont ordonnés comme dans _map
	    for (int j=mx+alphaNM;j<n+m;j++) mapping[j]=j; // les sommets restants sont supprimés
	    std::shuffle(&mapping[0], &mapping[mx], this->randGen); // on reshuffle les mx premiers
	    mappings.push_back(mapping);
	    delete [] _map;

            }

  }

  return mappings;
}





template<class NodeAttribute, class EdgeAttribute>
void MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getBestMapping( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                Graph<NodeAttribute,EdgeAttribute> * g1,
                Graph<NodeAttribute,EdgeAttribute> * g2,
                int * G1_to_G2, int * G2_to_G1 )
{

int k2 = k/2+1;

double * bestCost (0);
bestCost = new double;
*bestCost = -1;
double bestCostInitLoop=*bestCost;
double exploFactor=0;

 std::list<int*> mappings2;
 std::list<int*> visitedMappings;



  //Compute Mapping init
  struct timeval  tv1, tv2;
 gettimeofday(&tv1, NULL);
  std::list<int*> mappings = initGen->getMappings(g1, g2, k);
 // std::cout << "mappings created" << std::endl;
  int kk = 0;
  /*
  for (std::list<int*>::iterator it = mappings.begin(); it!=mappings.end(); it++){
      std::cout <<" mapping " << kk << "[";
      for (int i = 0; i<g1->Size();i++){
        std::cout << (*it)[i] << " ";
      }
      std::cout << "]" << std::endl;
     delete[] *it;
     //
     kk++;
}
*/
 // std::list<int*> mappings = getInitMappings(g1, g2, k);
  visitedMappings=mappings;
  //std::cout << "mappings copied" << std::endl;
  //for (std::list<int*>::iterator it = visitedMappings.begin(); it!=visitedMappings.end(); it++){
  //      *it = new int[g1->Size()];
  //      for (int i=0;i<g1->Size();i++){
  //          *it[i]= mappings[]
  //       }

  //}

  gettimeofday(&tv2, NULL);

  // initialisation de _psi



  _psiPre=new double[(g1->Size()+1)*( g2->Size()+1)];
  std::memset(_psiPre,0,sizeof(double)*(g1->Size()+1)*( g2->Size()+1));
  _psiPost=new double[(g1->Size()+1)*( g2->Size()+1)];
  std::memset(_psiPost,0,sizeof(double)*(g1->Size()+1)*( g2->Size()+1));

  this->getBestMappingFromSet(algorithm, g1, g2, G1_to_G2, G2_to_G1, mappings, _psiPre,_psiPost, bestCost, kr);

        //  std::cout << std::endl;
         // std::cout << "0  " << _nbIterations<< " " << g1->Size() << " " <<  g2->Size() << std::endl;
         // printMatrix(_psiPost,g1->Size()+1,g2->Size()+1);  // decommenter ces 3 lignes pour v10


  // AJOUT SAUVEGARDE BEST INIT MAPPING
  // int num = 0;
  // for ( std::list<int*>::iterator it = mappings.begin(), ite = mappings.end(); num != _numBestMap && it != ite; ++it, num++)
    // COMMENT SAUVEGARDER *it QUI DONNE UN int* (qui correspond à l'assignement initial qui a conduit à la meilleure solution)
    ;

  //nico Compute, improve and select explorative mappings







   for (int i=0; i<_nbIterations ;i++){
          bestCostInitLoop=*bestCost;
	  //mappings.clear();
	  switch(_psiMethod){
		  case 1:{
			  double * rPsi = new double[(g1->Size()+1)*( g2->Size()+1)];
			  rPsi = recenteredPsi(_psiPre,g1->Size(), g2->Size());

			  mappings = getKLSAPEMappings(g1->Size(), g2->Size(),rPsi, k);
			  delete[] rPsi;}
			  break;
		  case 2:{
			  double * rPsi = new double[(g1->Size()+1)*( g2->Size()+1)];
			  rPsi = invertedPsi(_psiPost,g1->Size(), g2->Size());

			  mappings = getKLSAPEMappings(g1->Size(), g2->Size(),rPsi, k);
			  delete[] rPsi;}
			  break;
		  case 3:{
			  double * iPsi = new double[(g1->Size()+1)*( g2->Size()+1)];
			  iPsi = invertedPsi(_psiPost,g1->Size(), g2->Size());
			  mappings2= getKLSAPEMappings(g1->Size(), g2->Size(),iPsi, k2);
			  mappings = getKLSAPEMappings(g1->Size(), g2->Size(),_psiPre, k2);
			  mappings.splice(mappings.end(),mappings2,mappings2.begin(),mappings2.end());
			  mappings2.clear();}
			  break;
		  case 4:{
    			  mappings = getRandomMappingsPsi(g1->Size(), g2->Size(), k,_psiPre,0,0);}
                          break;
                  case 5:{
    			  mappings = getRandomMappingsPsi(g1->Size(), g2->Size(), k,_psiPre,_alpha,0);}
                          break;
   		  case 6:{
    			  mappings = getRandomMappingsInvertedPsi(g1->Size(), g2->Size(), k,_psiPost,0,0,visitedMappings);}
                          break;
		  case 7:{
    			  mappings = getRandomMappingsInvertedPsi(g1->Size(), g2->Size(), k,_psiPost,_alpha,0,visitedMappings);}
                          break;
                  case 8:{
    			  mappings = getRandomMappingsMixtPsi(g1->Size(), g2->Size(), k,_psiPre,_psiPost,0,0);}
                          break;
                  case 9:{
    			  mappings = getRandomMappingsMixtPsi(g1->Size(), g2->Size(), k,_psiPre,_psiPost,_alpha,0);}
                          break;
                  case 10:{
			  mappings = getRandomMappingsInvertedPsi(g1->Size(), g2->Size(), k,_psiPost,2*_alpha,1,visitedMappings);}
                          break;
                  case 11:{
    			   mappings = getRandomMappingsInvertedPsi(g1->Size(), g2->Size(), k,_psiPost,3*_alpha,1,visitedMappings);}
                          break;
                  case 12:{
    			  mappings = getRandomMappingsInvertedPsi(g1->Size(), g2->Size(), k,_psiPost,0,1,visitedMappings,exploFactor);}
                          break;
                  case 13:{
    			  mappings = getRandomMappingsInvertedPsi(g1->Size(), g2->Size(), k,_psiPost,_alpha,1,visitedMappings);}
                          break;
                  case 14:{
    			  mappings = getRandomMappingsInvertedCenteredPsi(g1->Size(), g2->Size(), k,_psiPost,0,0,visitedMappings);}
                          break;
                  case 15:{
    			  mappings = getRandomMappingsInvertedCenteredPsi(g1->Size(), g2->Size(), k,_psiPost,_alpha,0,visitedMappings);}
                          break;

		  default:{
		  	  mappings = getKLSAPEMappings(g1->Size(), g2->Size(),_psiPre, k);}
                          break;

	   	  }


		  //printMappings(mappings, g1->Size());
		  //printMapping(G1_to_G2, g1->Size());
  		  //for (int j=0;j<g1->Size();j++)
	          //bestMapping[j] = G1_to_G2[j];
                  //printMapping(bestMapping, g1->Size());
		  //mappings.push_back(bestMapping);
   		this->getBestMappingFromSet(algorithm, g1, g2, G1_to_G2, G2_to_G1, mappings, _psiPre, _psiPost, bestCost, kr);
                if (bestCostInitLoop != 0)
                exploFactor=(((bestCostInitLoop-*bestCost)/bestCostInitLoop)>0.05) ? adap*(1-(bestCostInitLoop-*bestCost)/bestCostInitLoop) : 0;
                else
                exploFactor=0;
  //               std::cout << std::endl << i +1 << " " << _nbIterations<< " " << g1->Size() << " " <<  g2->Size() << std::endl;
//          printMatrix(_psiPost,g1->Size()+1,g2->Size()+1); //décommenter cette ligne et celle d'avant pour v10



 //if (*bestCost < bestCostInitLoop)
	         //std::cout<< "improvement of Best Cost by " << * bestCost - bestCostInitLoop << " in loop number " << i << ". New Best Cost = " << *bestCost << std::endl;


}




  // Memoy cleaning
// for (std::list<int*>::iterator it = mappings.begin(); it!=mappings.end(); it++){
//    delete[] *it;
//}
kk = 0;

for (std::list<int*>::iterator it = visitedMappings.begin(); it!=visitedMappings.end(); it++){
     delete[] *it;
}

//for (std::list<int*>::iterator it = visitedMappings.begin(); it!=visitedMappings.end(); it++){
//      std::cout << "deleting mapping " << kk << ":" << std::endl << "[";
//      for (int ii = 0; ii<g1->Size();ii++){
//        std::cout << (*it)[ii] << " ";
//      }
//      std::cout << "]" << std::endl;
//     delete[] *it;
     //
//     kk++;
//}
//visitedMappings.clear();
delete bestCost;
}



template<class NodeAttribute, class EdgeAttribute>
void MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                Graph<NodeAttribute,EdgeAttribute> * g1,
                Graph<NodeAttribute,EdgeAttribute> * g2,
                int * G1_to_G2, int * G2_to_G1,
                std::list<int*>& mappings, double * _psiPre,double * _psiPost, double* pcost)
{
this->kr=this->k;
this->getBestMappingFromSet(algorithm, g1, g2, G1_to_G2, G2_to_G1, mappings, _psiPre, _psiPost, pcost, kr);
}

template<class NodeAttribute, class EdgeAttribute>
void MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                Graph<NodeAttribute,EdgeAttribute> * g1,
                Graph<NodeAttribute,EdgeAttribute> * g2,
                int * G1_to_G2, int * G2_to_G1,
                std::list<int*>& mappings, double * _psiPre,double * _psiPost, double* pcost, int kr)
{
this->adap=0;
this->getBestMappingFromSet(algorithm, g1, g2, G1_to_G2, G2_to_G1, mappings, _psiPre, _psiPost, pcost, kr,adap);
}


template<class NodeAttribute, class EdgeAttribute>
void MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                Graph<NodeAttribute,EdgeAttribute> * g1,
                Graph<NodeAttribute,EdgeAttribute> * g2,
                int * G1_to_G2, int * G2_to_G1,
                std::list<int*>& mappings, double * _psiPre,double * _psiPost, double* pcost, int kr, double adap)
{

//debug
//std::cout << "init getBestMappingFromSet" << std::endl;

  struct timeval  tv1, tv2;
  int n = g1->Size();
  int m = g2->Size();
int numConv=0;

  typename std::list<int*>::const_iterator it;
  double cost = -1;
  double ncost;


  // Multithread
  #ifdef _OPENMP
    gettimeofday(&tv1, NULL);
    int** arrayMappings = new int*[mappings.size()];
    int* arrayCosts = new int[mappings.size()];
    int* arrayLocal_G1_to_G2 = new int[n * mappings.size()];
    int* arrayLocal_G2_to_G1 = NULL;
    if (G2_to_G1 != NULL)
      arrayLocal_G2_to_G1 = new int[m * mappings.size()];

    int i=0; for (it=mappings.begin(); it!=mappings.end(); it++){
      arrayMappings[i] = *it;
      i++;
    }

    //omp_set_dynamic(0);
    omp_set_num_threads(6);
    #pragma omp parallel for schedule(dynamic) reduction(+:numConv) //private(tid, i, j, ncost, ipfpGed )
    for (unsigned int tid=0; tid<mappings.size(); tid++){
      int* lsapMapping = arrayMappings[tid];
      int* local_G1_to_G2 = &(arrayLocal_G1_to_G2[tid*n]);

      int* local_G2_to_G1 = NULL;
      if (G2_to_G1 != NULL)
        local_G2_to_G1 = &(arrayLocal_G2_to_G1[tid*m]);

  // Sequential
  #else

  //debug
std::cout << "init creation localmappings" << std::endl;
    int* local_G1_to_G2 = new int[n];
    int* local_G2_to_G1 = NULL;
    if (G2_to_G1 != NULL)
      local_G2_to_G1 = new int[m];

    double t_acc = 0; // accumulated time
    int numMap = 0;
    for (it=mappings.begin(); it!=mappings.end(); it++, numMap++){
      gettimeofday(&tv1, NULL);
      int* lsapMapping = *it;

  #endif

    // Copy the mapping into the local array
    for (int i=0; i<n; i++)
      local_G1_to_G2[i] = lsapMapping[i];

    if (local_G2_to_G1 != NULL){
      for (int j=0; j<m; j++) local_G2_to_G1[j] = -1;
      for (int i=0; i<n; i++)
        if (local_G1_to_G2[i] >= 0)
          local_G2_to_G1[local_G1_to_G2[i]] = i;
    }
    // incrementation of psi before the descent
//debug
std::cout << " init incrementation of psi before the descent" << std::endl;

for (int ii=0; ii<n; ii++)
     _psiPre[sub2ind(ii,local_G1_to_G2[ii],n+1)]++;
if (G2_to_G1 != NULL)
    for (int jj=0; jj<m; jj++)
        if (local_G2_to_G1[jj]==n)
            _psiPre[sub2ind(n,jj,n+1)]++;
//debug
std::cout << " end incrementation of psi before the descent" << std::endl;

    MappingRefinement<NodeAttribute, EdgeAttribute> * local_method;

    #ifdef _OPENMP
      local_method = algorithm->clone();
    #else
      local_method = algorithm;
    #endif

    local_method->getBetterMapping(kr,&numConv,g1, g2, local_G1_to_G2, local_G2_to_G1, true);
    ncost = local_method->mappingCost(g1, g2, local_G1_to_G2, local_G2_to_G1);


    // Multithread
    #ifdef _OPENMP
      // save the approx cost
      arrayCosts[tid] = ncost;
     // _distances_[tid] = ncost;



    // Sequential
    #else
 //debug
std::cout << " init incrementation of psi after the descent" << std::endl;


      //  incrementation of _psi after the descent

for (int ii=0; ii<n; ii++)
     _psiPost[sub2ind(ii,local_G1_to_G2[ii],n+1)]++;
if (G2_to_G1 != NULL)
    for (int jj=0; jj<m; jj++)
        if (local_G2_to_G1[jj]==n)
            _psiPost[sub2ind(n,jj,n+1)]++;


    //debug
std::cout << " end incrementation of psi after the descent" << std::endl;



      // if ncost is better : save the mapping and the cost
      if (cost > ncost || cost == -1){
        cost = ncost;
        std::cout << "improving solution found" << std::endl;
	//_numBestMap = numMap;
        for (int i=0; i<n; i++) G1_to_G2[i] = local_G1_to_G2[i];
        if (G2_to_G1 != NULL)
          for (int j=0; j<m; j++) G2_to_G1[j] = local_G2_to_G1[j];
      }
      gettimeofday(&tv2, NULL);
      t_acc += ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));

    #endif

  } //end for

  // Multithread : Reduction
  #ifdef _OPENMP
    gettimeofday(&tv2, NULL);

    gettimeofday(&tv1, NULL);

    int i_optim=0;
    for (unsigned int i=0; i<mappings.size(); i++){
      if (cost > arrayCosts[i] || cost == -1){
         cost = arrayCosts[i];
         i_optim = i;
      }
    }
    for (int i=0; i<n; i++) G1_to_G2[i] = arrayLocal_G1_to_G2[i_optim*n + i];

    if (G2_to_G1 != NULL)
      for (int j=0; j<m; j++) G2_to_G1[j] = arrayLocal_G2_to_G1[i_optim*m + j];

    // To match the output format size in XPs
    //for (int i=mappings.size(); i<k; i++) _distances_[i] = 9999;

    gettimeofday(&tv2, NULL);
    //_xp_out_ <<  ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << ", ";
    //_xp_out_ << ((float)t) / CLOCKS_PER_SEC << ", ";

    delete[] arrayLocal_G1_to_G2;
    delete[] arrayCosts;
    delete[] arrayMappings;
    if (arrayLocal_G2_to_G1 != NULL)
      delete[] arrayLocal_G2_to_G1;

  // Sequential : deletes
  #else

    delete [] local_G1_to_G2;
    if (G2_to_G1 != NULL)
      delete [] local_G2_to_G1;

  #endif


}


template<class NodeAttribute, class EdgeAttribute>
const std::list<int*>& MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getBetterMappings( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                   Graph<NodeAttribute,EdgeAttribute> * g1,
                   Graph<NodeAttribute,EdgeAttribute> * g2 )
{
  //Compute Mapping init
  struct timeval  tv1, tv2;
  gettimeofday(&tv1, NULL);

  std::list<int*> mappings = initGen->getMappings(g1, g2, k);
  gettimeofday(&tv2, NULL);

  return this->getBetterMappingsFromSet(algorithm, g1, g2, mappings);
}



template<class NodeAttribute, class EdgeAttribute>
const std::list<int*>& MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getBetterMappingsFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                          Graph<NodeAttribute,EdgeAttribute> * g1,
                          Graph<NodeAttribute,EdgeAttribute> * g2,
                          std::list<int*>& mappings )
{
  int n = g1->Size();

  typename std::list<int*>::const_iterator it;

  int* arrayLocal_G1_to_G2 = new int[n * mappings.size()];

  // Multithread
  #ifdef _OPENMP
    int** arrayMappings = new int*[mappings.size()];

    int i=0; for (it=mappings.begin(); it!=mappings.end(); it++){
      arrayMappings[i] = *it;
      i++;
    }

    //omp_set_dynamic(0);
    omp_set_num_threads(6);
    // Set a dynamic scheduler : algorithms can have different execution time given different initial mappings
    #pragma omp parallel for schedule(dynamic) //private(tid, i, j, ncost, ipfpGed )
    for (unsigned int tid=0; tid<mappings.size(); tid++){
      int* lsapMapping = arrayMappings[tid];

  // Sequential
  #else
    unsigned int tid=0;
    for (it=mappings.begin(); it!=mappings.end(); it++){
      int* lsapMapping = *it;
  #endif


    int* local_G1_to_G2 = &(arrayLocal_G1_to_G2[tid*n]);

    // Copy the mapping into the local array
    for (int i=0; i<n; i++)
      local_G1_to_G2[i] = lsapMapping[i];

    MappingRefinement<NodeAttribute, EdgeAttribute> * local_method;

    #ifdef _OPENMP
      local_method = algorithm->clone();
    #else
      local_method = algorithm;
      tid++;
    #endif

    local_method->getBetterMapping(g1, g2, local_G1_to_G2, NULL, true);

  } //end for

  // Reduction - indexation of the list
    refinedMappings.clear();
    int _i_=0;
    for (it=mappings.begin(); it!=mappings.end(); it++){
      refinedMappings.push_back(arrayLocal_G1_to_G2 + _i_*n);
      _i_++;
    }

   #ifdef _OPENMP
    delete[] arrayMappings;
   #endif

   return refinedMappings;
}




#endif // __MULTIPLEIPFPGRAPHEDITDISTANCE_H__elete[] rPsi;
