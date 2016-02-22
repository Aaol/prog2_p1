#include <iostream>
#include <iomanip>

#include <string>
#include <cmath>
#include <cassert>
#include <cstdlib>

using namespace std;

const int MAX_SIZE = 15;
const int STORAGE_SIZE = MAX_SIZE * MAX_SIZE;

/*******************************************************/

class Matrix {

private:
  int my_size_i, my_size_j;	// Les tailles en i et j
  double contents[STORAGE_SIZE];     // Les éléments de la matrice

/**************************/

// Une méthode privée pour la projection d'une matrice 2D sur un tableau 1D
  int embedding(int i, int j); 
 
/**************************/

public:
  Matrix(int n, int p);	// Constructeur
  int num_lines(void);
  int num_columns(void);
  void set(int i, int j, double x);
  double get(int i, int j);

  void print();
};

/*******************************************************/

int Matrix::embedding(int i, int j) {
  return (my_size_j * i) + j;
}

Matrix::Matrix(int n, int p) {
  assert(1 <= n &&  n < MAX_SIZE);	// Pour détecter les cas aberrants
  assert(1 <= p &&  p < MAX_SIZE);
  my_size_i = n;
  my_size_j = p;

  for (int i = 0; i < my_size_i; i++)
    for (int j = 0; j < my_size_j; j++) {
      int a = embedding(i, j);
      assert(0 <= a && a < STORAGE_SIZE);
      contents[a] = 0.0;
    }
}

int Matrix::num_lines(void) {
  return my_size_i;
}

int Matrix::num_columns(void) {
  return my_size_j;
}

void Matrix::set(int i, int j, double x) {
  assert(0 <= i && i < num_lines());
  assert(0 <= j && j < num_columns());
  int n = embedding(i, j);
  assert(0 <= n && n < STORAGE_SIZE);
  contents[n] = x;
}

double Matrix::get(int i, int j) {
  // TODO
}

void Matrix::print() {
  // TODO
}

/*******************************************************/



int main() {

  int n = 3;

  cout << "Setting n = " << n << endl;
  
  Matrix M(n, n+2);

  for (int i = 0; i < M.num_lines(); i++)
    for (int j = 0; j < M.num_columns(); j++)
      M.set(i, j, i+j);
  
  M.print();

  return 0;
}
