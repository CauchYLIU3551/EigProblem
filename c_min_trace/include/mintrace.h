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

#include <CG/CGSolver.h>
#include <AFEPack/Miscellaneous.h>
#include <AFEPack/AMGSolver.h>

class TraceSolver:public CGSolver
{
 public:
  typedef SparseMatrix<double> Matrix;
  TraceSolver();
  TraceSolver(const Matrix& a, const Matrix& m);

  void rand_V(int p, std::vector<std::vector<double>>& V);
  void Householder(std::vector<std::vector<double>> &a);
  void QR(std::vector<std::vector<double>> &a, std::vector<std::vector<double>>& Q);
  void get_VtAV(std::vector<std::vector<double>> a);
  void QRSolver(std::vector<std::vector<double>> a);
  double get_residual();

  virtual void get_Ap(std::vector<double> x);
  
  void solve(int p,
	     double tol=0.0,
	     u_int step=20);
 private:
  const Matrix* A;
  const Matrix* M;
  std::vector<std::vector<double>> X; // solution matrix, which contains the eigenvectors of the eigenvalues.
  std::vector<double> theta;// storing the diagonal entries of the QR factorization and maybe used
  // as shift in future advancing function.
  std::vector<double> lambda; // storing the eigenvalues;

};

#endif

