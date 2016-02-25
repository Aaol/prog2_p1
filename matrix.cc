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
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  Matrix M(size_i,size_j);
  
  for ( unsigned i = 0; i < size_i; i++)
    {
      for(unsigned j = 0; j< size_j; j++)
	{
	  Matrix::scalar_t x = M1.get(i, j) - M2.get(i, j);
	  M.set(i, j, x);
	}
    }
  return M;
}

Matrix operator*(Matrix::scalar_t a, const Matrix& M1) {
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  Matrix M(size_i, size_j);
  for ( unsigned i = 0; i < size_i; i++)
    {
      for(unsigned j = 0; j< size_j; j++)
	{
	  Matrix::scalar_t x = a * M1.get(i ,j);
	  M.set(i, j, x);
	}
    }
  return M;
}
/* M1 size (p,q), M2 size(q,r)*/
Matrix operator*(const Matrix& M1, const Matrix& M2) {
  assert(M1.get_size_j() == M2.get_size_i());
  unsigned p = M1.get_size_i();
  unsigned q = M2.get_size_i();
  unsigned r = M2.get_size_j();
  Matrix M(p,r);

  for(unsigned i = 0; i < p; i++)
    {
      for( unsigned j = 0; j < r; j++)
	{
	  Matrix::scalar_t x = 0;
	  for (unsigned k = 0; k < q; k++)
	    {
	      x+= M1.get(i,k) * M2.get(k,j);
	    }
	  M.set( i, j, x);
	}
    }
  return M;
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
  Matrix::scalar_t result = 0.0;
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  for( unsigned i = 0 ; i < size_i; i++)
    {
      for(unsigned j = 0; j < size_j ; j++)
	{
	  if(abs(M1.get(i,j))> result )
	    result = M1.get(i,j);
	}
    }
  return result;
}

static Matrix submatrix(const Matrix& M1, unsigned a, unsigned b) {	//Note it is static!
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  assert (0 <= a && a < size_i);
  assert (0 <= b && b < size_j);
  assert(size_i >= 2);
  assert(size_j >= 2);
  Matrix M(size_i-1, size_j-1);
  for (unsigned i = 0; i < size_i; i++)
    {
      for(unsigned j = 0 ; j < size_j; j++)
	{
	  if(i < a)
	    {
	      if(j < b)
		{
		  M.set(i, j, M1.get(i,j));
		}
	      else if(j >b)
		{
		  M.set(i, j-1, M1.get(i,j) );
		}
	    }
	  else if(i > a)
	    {
	      if(j < b)
		{
		  M.set(i-1, j, M1.get(i,j) );
		}
	      else if (j > b)
		{
		  M.set( i-1, j-1, M1.get(i,j ) );
		}
	    }
	}
    }
  return M;
}

static int toggle(unsigned i) {		//Note it is static!
  if( i % 2 == 0)
    return 1;
  else 
    return -1;
}

Matrix::scalar_t determinant(const Matrix& M1) {
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  assert(size_i == size_j);
  if(size_i == 1)
    return M1.get(0,0);
  else
    {
      Matrix::scalar_t result = 0.0;
      for (unsigned i = 0; i < size_i; i++)
	{
	  result += toggle(i)*M1.get(0,i)*determinant(submatrix ( M1, 0, i));
	}
      return result;
    } 
}

Matrix transpose(const Matrix& M1) {
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();

  Matrix M(size_j,size_i);
  for(unsigned i = 0; i < size_i; i++)
    {
      for(unsigned j = 0; j <size_j; j++)
	{
	   M.set( j, i, M1.get(j,i));
	}
    }
  return M;
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
  Matrix M2 = transpose(M1);

  for (unsigned i = 0; i < size; i++) {
    for (unsigned j = 0; j < size; j++) {
      Matrix M = submatrix(M2, i, j);

      Matrix::scalar_t x = determinant(M);
      Matrix::scalar_t y =  toggle(i+j) * (x / det);
      cofactors.set(i, j, y);
    }
  }

  return transpose(cofactors);
}
Matrix::scalar_t verification (const Matrix& M1)
{
  assert(M1.get_size_i() == M1.get_size_j());
  unsigned i = M1.get_size_i();
  Matrix MM = inverse(M1);

  Matrix N1 = M1 * MM - Id(i);
  Matrix N2 = MM * M1 - Id(i);
  return norm(N1) - norm(N2);
}
/* Return an hilbert matrix of size n */
Matrix hilbert(unsigned n)
{
  Matrix M(n,n);
  for(unsigned i = 0; i < n ; i++)
    {
      for(unsigned j = 0; j < n; j++)
	{
	  //this is i + j +1 because i and j start at 0 and not at 1
	  M.set(i,j, 1.0/(i + j +1.0));
	}
    }
  return M;
}

