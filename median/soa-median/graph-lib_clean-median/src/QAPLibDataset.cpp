#include "QAPLibDataset.h"
#include "utils.h"

QAPLibDataset::QAPLibDataset(const char* filename)
{
  const char * ext = strrchr(filename,'.');
  if (strcmp(ext,".ds") == 0){
    loadDS(filename);
  }
}


void QAPLibDataset::loadDS(const char* filename)
{

  std::ifstream f_tmp (filename);
  char * unconst_filename = new char[strlen(filename)+1];
  unconst_filename = strcpy(unconst_filename, filename);
  char * path = dirname(unconst_filename);
  if (f_tmp.is_open()){
    std::string s;
    while (getline(f_tmp, s)){
      if (s[0] != '#'){
	      std::string path_ctfile(path);
	      path_ctfile += std::string("/");
	      std::istringstream liness(s);
	      std::string ctfile;
	      liness >> ctfile;
	      int y;
	      liness >> y;
	      std::string full_ctfile = path_ctfile;
	      full_ctfile += ctfile;
	      loadQAP(full_ctfile.c_str());
      }
    }
  }
  f_tmp.close();

  delete[] unconst_filename;
}


void QAPLibDataset::loadQAP(const char* filename)
{
  int input;
  std::ifstream file (filename);

  int n;
  file >> n;

  if (file.fail()){
    std::cerr << "[E] Failed reading file " << filename << std::endl;
    return;
  }

  int* matA = new int[n*n];
  int* matB = new int[n*n];

  // Extract Matrix A
  int i=0;  int j=0;
  while (i < n){
    j=0;
    while (!file.eof() && j<n){
      file >> input;
      matA[sub2ind(i,j,n)] = input;
      j++;
    }
    i++;
  }
  if (i < n){
    std::cerr << "[E] Failed reading file " << filename << std::endl;
    return;
  }


  // Extract Matrix B
  i=0; j=0;
  while (i<n){
    j=0;
    while (!file.eof() && j<n){
      file >> input;
      matB[sub2ind(i,j,n)] = input;
      j++;
    }
    i++;
  }
  if (i < n){
    std::cerr << "[E] Failed reading file " << filename << std::endl;
    return;
  }


  QAPLibGraph * gA = new QAPLibGraph(matA, n);
  QAPLibGraph * gB = new QAPLibGraph(matB, n);
  this->add(gA,0);
  this->add(gB,0);

  delete[] matA;
  delete[] matB;
}
