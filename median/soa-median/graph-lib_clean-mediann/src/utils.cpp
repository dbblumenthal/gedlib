/*
 * @file utils.cpp
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>>
 * @version     0.0.1 - Thu Feb  2 2017
 *
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *
 * Description of the program objectives.
 * All necessary references.
 *
 */
#include <cstring>
#include <iostream>

#include "utils.h"

std::vector<char*> split (const char* chaine, const char* sep){
  std::vector<char*> v;

  char* s = strtok ((char*)chaine, (char*)sep);

  while (s != NULL)
    {
      v.push_back (s);
      s = strtok (NULL, (char*)sep);
    }

  return v;
}

//int sub2ind(int i, int j, int n){return i + j*n;};

void print_array(double* arr, int n, int m){
  for (int i=0; i<n; i++){
    for (int j=0; j<m; j++)
      std::cout << arr[i + j*n] << "  ";
    std::cout << std::endl;
  }
}

//double abs(double x) { return (x >= 0) ? x : -x; }
