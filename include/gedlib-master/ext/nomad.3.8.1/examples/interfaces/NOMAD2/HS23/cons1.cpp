// cons1.C -- PB23 general constraint #1 treated as a black box.

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const double evaluate(const double coord[]);

int main()
{
  unsigned int dimension = 2;
  double coord[2];
  ifstream input("input.txt");
  for (unsigned int i = 0; i < dimension; i++)
    input >> coord[i];
  input.close();

  double value = evaluate(coord);
  cout << value;

  return 0;
}

const double evaluate(const double coord[])
{
  return (-coord[0] - coord[1] + 1);
}
