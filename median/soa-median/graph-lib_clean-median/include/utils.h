/**
 * @file utils.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Thu Feb  2 2017
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>


#define sub2ind(i, j, n)    (i + (j) * (n))

std::vector<char*> split (const char* chaine, const char* sep);

//int sub2ind(int i, int j, int n);

template<typename T>
T mean(T * tab, int size){
  T sum = 0;
  for (int i =0;i<size;i++)
    sum += tab[i];
  return sum/size;
}

#endif // __UTILS_H__
