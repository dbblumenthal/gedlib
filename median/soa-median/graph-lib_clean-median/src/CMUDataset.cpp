#include "CMUDataset.h"


CMUDataset::CMUDataset(const char* filename) 
{
  const char * ext = strrchr(filename,'.'); 
  if (strcmp(ext,".ds") == 0){
    loadDS(filename);
  }
}


void CMUDataset::loadDS(const char* filename)
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
	      CMUGraph * g = new CMUGraph(full_ctfile.c_str());
	      this->add(g,y);
      }
    }
  }
  f_tmp.close();

  delete[] unconst_filename;
}
