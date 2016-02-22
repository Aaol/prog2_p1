#include <iostream>
#include <iomanip>

#include <string>
#include <cmath>
#include <cassert>
#include <cstdlib>

#include <vector>

#include "matrix.hh"

using namespace std;

Matrix::Matrix(unsigned n, unsigned p) {
  size_i = n;
  size_j = p;
  assert(1 <= size_i);
  assert(1 <= size_j);

  contents = vector<vector<scalar_t>/**/>(size_i);
  for (unsigned i = 0; i < size_i; i++)
    contents.at(i) = vector<scalar_t>(size_j, 0.0);	// Initialization to 0.0
}

unsigned Matrix::get_size_i() const {
  return size_i;
}

unsigned Matrix::get_size_j() const {
  return size_j;
}

void Matrix::set(unsigned i, unsigned j, scalar_t x) {
  assert(0 <= i && i < size_i);
  assert(0 <= j && j < size_j);
  contents.at(i).at(j) = x;
}

Matrix::scalar_t Matrix::get(unsigned i, unsigned j) const {
  assert(0 <= i && i < size_i);
  assert(0 <= j && j < size_j);
  return contents.at(i).at(j);
}

void Matrix::print() const {
  for (unsigned i = 0; i < size_i; i++) {
    for (unsigned j = 0; j < size_j; j++) {
      cout << setprecision(2) << setw(8) << contents.at(i).at(j);
    }
    cout << endl;
  }
  cout << "___________________________" << endl;
}

/*****************************************************/

Matrix operator+(const Matrix& M1, const Matrix& M2) {
  assert(M1.get_size_i() == M2.get_size_i());
  assert(M1.get_size_j() == M2.get_size_j());
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  Matrix M(size_i, size_j);

  for (unsigned i = 0; i < size_i; i++) {
    for (unsigned j = 0; j < size_j; j++) {
      Matrix::scalar_t x =  M1.get(i, j) + M2.get(i, j);
      M.set(i, j, x);
    }
  }

  return M;
}

Matrix operator-(const Matrix& M1, const Matrix& M2) {
  assert(M1.get_size_i() == M2.get_size_i());
  assert(M1.get_size_j() == M2.get_size_j());
  // ...
  return Id(1);
}

Matrix operator*(Matrix::scalar_t a, const Matrix& M1) {
  // ...
  return Id(1);
}

Matrix operator*(const Matrix& M1, const Matrix& M2) {
  assert(M1.get_size_j() == M2.get_size_i());
  // ...
  return Id(1);
}

Matrix Id(unsigned n) {
  unsigned size = n;
  Matrix M(size, size);
  for (unsigned i = 0; i < size; i++) {
    M.set(i, i, 1.0);
  }
  return M;
}

Matrix::scalar_t norm(const Matrix& M1) {
  // ...
  return 0;
}

static Matrix submatrix(const Matrix& M1, unsigned a, unsigned b) {	//Note it is static!
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  assert (0 <= a && a < size_i);
  assert (0 <= b && b < size_j);
  assert(size_i >= 2);
  assert(size_j >= 2);
  // ...
  return Id(1);
}

static int toggle(unsigned i) {		//Note it is static!
  // ...
  return 0;
}

Matrix::scalar_t determinant(const Matrix& M1) {
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  assert(size_i == size_j);
  // ...
  return 0;
}

Matrix transpose(const Matrix& M1) {
  // ...
  return Id(1);
}

Matrix inverse(const Matrix& M1) {
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  assert(size_i == size_j);

  Matrix::scalar_t det = determinant(M1);
  if (det == 0) {
    cerr << "Cannot invert a matrix with null determinant" << endl;
    M1.print();
    exit(1);
  }

  unsigned size = size_i;
  Matrix cofactors(size, size);
  Matrix M2 = transpose(M1);	// Do not forget it!

  for (unsigned i = 0; i < size; i++) {
    for (unsigned j = 0; j < size; j++) {
      Matrix M = submatrix(M2, i, j);

      Matrix::scalar_t x = determinant(M);
      Matrix::scalar_t y =  toggle(i+j) * (x / det);
      cofactors.set(i, j, y);
    }
  }

  return cofactors;
}
