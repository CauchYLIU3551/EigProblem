/*
 * @file   mintrace.h
 * @author Chengyu Liu
 * @date   Fri Jan 22 2020
 *
 * @ brief trace minimization method to get the p smallest eigenvalues.
 *
 */

#ifndef _mintrace_h_
#define _mintrace_h_

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <list>
#include <base/exceptions.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>

#include <AFEPack/Miscellaneous.h>

class TraceSolver{
 public:
  typedef SparseMatrix<double> Matrix;
  TraceSolver();
  TraceSolver(const Matrix& a, const Matrix& m);

  void solve(std::vector<double>& x,
	     double tol=0.0,
	     u_int step=20); 
 private:
  Matrix A;
  Matrix M;
};

#endif

