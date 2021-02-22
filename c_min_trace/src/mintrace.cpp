#include<iostream>
#include<assert.h>
#include<lac/sparsity_pattern.h>
//#include<AFEPack/AMGSolver.h>
#include <trace/mintrace.h>

TraceSolver::TraceSolver(){};
TraceSolver::TraceSolver(const Matrix& a, const Matrix& m){
  A=&a;
  M=&m;
  assert(A->m()==A->n());
  assert(M->m()==M->n());
  assert(A->m()==M->m());

  // initializing the vector of theta;
  theta.resize(A->m());
}

void TraceSolver::rand_V(int p, std::vector<std::vector<double>>& V){};
void TraceSolver::get_VtAV(std::vector<std::vector<double>> V){};
void TraceSolver::QRSolver(std::vector<std::vector<double>> V){};
void TraceSolver::get_Ap(std::vector<double> x){};

double TraceSolver::get_residual()
{
  double res=0;
  return res;
}

/*void TraceSolver::solve(std::vector<double>& x,double tol, u_int step){}*/

// solve function computes the smallest p eigenvalues of the Matrix A and M;
// and the eigenvectors will be stored in the matrix X;
void TraceSolver::solve(int p, double tol, u_int max_iter)
{
  std::vector<std::vector<double>> V(A->n());
  // using rand_V function to get the initial matrix V which is orthogonal corresponding with M;
  rand_V(p, V);
  int iter = 0;
  while(iter<max_iter)
    {
      get_VtAV(V);
      QRSolver(V);// return the diagonal entries of VtAV in theta and the V_k* U_k stored in X;
      //Householder(V);//this function will return the tridiagonal function in V and the Householder matrices stored in X;
      //QR(V);
      double residual=get_residual();
      if (residual<tol)
	{
	  break;
	}
      for (int i=0;i<p;i++)
	{
	  std::vector<double> delta(A->m(),0), rhs(A->m(),0);// using a n x 1 dimension initial vector to compute the Conjugate gradient process.
	  //computing rhs vector!!
	  //
	  //
	  //////////////////////////////
	  CGSolver::solve(delta,rhs,1.0e-3,20);
	}
      iter++;
    }
  
}



