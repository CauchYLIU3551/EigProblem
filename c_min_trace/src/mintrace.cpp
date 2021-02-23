#include<iostream>
#include<assert.h>
#include<lac/sparsity_pattern.h>
//#include<AFEPack/AMGSolver.h>
#include <trace/mintrace.h>
/////////////////////////////////////////////
// These functions can be gathered into one .cpp file after testing;

// this function will turn the matrix Q into the nxn identity matrix;
void identitymatrix(std::vector<std::vector<double>> &Q, int n)
{
  Q.clear();
  std::vector<double> temp(n,0);
  for (int i=0;i<n;i++)
    {
      Q.push_back(temp);
      Q[i][i]=1;
    }
}

//This function computes the matrix A times matrix B, and its result will be stored in A.
// It can be improved if reducing the use of the temp matrix to store the original values which
// is used in computing.
void multiply(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> B)
{
    int n = A.size();
    int m = B.size();
    int col = B[0].size();
    std::vector<double> temp(m);
    std::vector<std::vector<double>> A0(n, temp);

    if (col < n)
      {
        for (int i = 0; i < n; i++)
	  {
            for (int j = 0; j < col; j++)
	      {
                double tempValue = 0;
                for (int t = 0; t < m; t++)
		  {
                    tempValue += A[i][t] * B[t][j];
		  }
                A0[i][j] = tempValue;
	      }
            for (int j = col; j < m; j++)
	      {
                A0[i].pop_back();
	      }
	  }

        A = A0;
      }
    else
      {
        for (int i = 0; i <n; i++)
	  {
            for (int j = 0; j < m; j++)
	      {
                double tempValue = 0;
                for (int t = 0; t < m; t++)
		  {
                    tempValue += A[i][t] * B[t][j];
		  }
                A0[i][j] = tempValue;
	      }
            for (int j = m; j < col; j++)
	      {
                double tempValue = 0;
                for (int t = 0; t < m; t++)
		  {
                    tempValue += A[i][t] * B[t][j];
		  }
                A0[i].push_back(tempValue);
	      }
	  }
        A = A0;
      }
}



/////////////////////////////////////////////
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
//void TraceSolver::Householder(std::vector<std::vector<double>> &a){};
void TraceSolver::Householder(std::vector<std::vector<double>> &a)
{
  //initialize the matrix X
  std::vector<double> temp(a.size(),0);
  X.resize(a.size(),temp);// just for test this function, will be deleted after debug;
  for(int j=0;j<X.size();j++)
    {
      X[j][j]=1;
    }
  int n=a.size();
  std::vector<std::vector<double>> Pk;

  for(int k=0;k<n-2;k++)
    {
      identitymatrix(Pk, n);
      double q=0;
      double alpha=0;
      // Computing the value of q in the Householder function;
      for(int i=k+1;i<n;i++)
	{
	  q+=a[i][k]*(a[i][k]);
	}

      // Computing the value of alpha
      if(a[k+1][k]==0)
	{
	  alpha=-sqrt(q);
	}
      else
	{
	  alpha=-sqrt(q)*a[k+1][k]/abs(a[k+1][k]);
	}
      
      double RSQ=alpha*alpha-alpha*a[k+1][k]; // RSQ = 2r^2;

      std::vector<double> v(n-k,0), u(n-k,0), z(n-k,0);
      v[0]=0;
      v[1]=a[k+1][k]-alpha;

      for(int j=2;j<v.size();j++)
	{
	  v[j]=a[k+j][k];
	}
      
      for(int j=0;j<n-k;j++)
	{
	  for(int i=1;i<n-k;i++)
	    {
	      u[j]+=a[k+j][k+i]*v[i];
	    }
	  u[j]/=RSQ;
	}
      
      double PROD=0;
      for(int i=1;i<n-k;i++)
	{
	  PROD+=v[i]*u[i];
	}
      
      for(int j=0;j<n-k;j++)
	{
	  z[j]=u[j]-PROD/(2*RSQ)*v[j];
	}

      for(int l=k+1;l<n-1;l++)
	{
	  for(int j=l+1;j<n;j++)
	    {
	      a[j][l]=a[j][l]-v[l-k]*z[j-k]-v[j-k]*z[l-k];
	      a[l][j]=a[j][l];
	    }
	  a[l][l]=a[l][l]-2*v[l-k]*z[l-k];
	}
      
      a[n-1][n-1]=a[n-1][n-1]-2*v[n-k-1]*z[n-k-1];
	
      for(int j=k+2;j<n;j++)
	{
	  a[k][j]=0;
	  a[j][k]=0;
	}

      a[k+1][k]=a[k+1][k]-v[1]*z[0];
      a[k][k+1]=a[k+1][k];


      for (int j=k+1;j<n;j++)
	{
	  std::vector<double> temp(j-k,0);
	  for(int i=k+1;i<=j;i++)
	    {
	      for(int t=k+1;t<n;t++)
		{
		  temp[i-k-1]+=v[i-k]*v[t-k]*Pk[t][j]/(2*RSQ);
		}
	    }
	  for (int i=k+1;i<=j;i++)
	    {
	      Pk[i][j]-=2*temp[i-k-1];
	      Pk[j][i]=Pk[i][j];
	    }
	}
      multiply(X,Pk);
    }

}

//void TraceSolver::QR(std::vector<std::vector<double>> &a){};
///////////////
// This function computes the QR factorization of the tridiagonal matrix A.
// Input: tridiagonal matrix A, identity matrix Q
// Output: uptriangle matrix A_, the orthonormal matrix Q; A=Q*A_;
//
// tips: Maybe I can improve the function by replace the matrix A by two vectors an bn-1;
// an contains the diagonal entries of A; bn-1 contains the sub-diagonal entries of A;

void TraceSolver::QR(std::vector<std::vector<double>> &a,   std::vector<std::vector<double>>& Q)
{
  int n=a.size();

  for (int i=0;i<a.size()-1;i++)
    {
      double theta=0;
      double tempa1=a[i][i],tempa2=a[i+1][i+1],tempb=a[i+1][i],tempc=a[i][i+1];
      theta=atan(a[i+1][i]/a[i][i]);
      double c=cos(theta), s=sin(theta);

      a[i][i]=c*tempa1+s*tempb;
      a[i][i+1]=c*tempc+s*tempa2;
      a[i+1][i+1]=-s*tempc+c*tempa2;
      a[i+1][i]=-s*tempa1+c*tempb;
      if (i<a.size()-2)
	{
	  a[i][i+2]=s*a[i+1][i+2];
	  a[i+1][i+2]=c*a[i+1][i+2];
	}


      // update the orthogonal matrix Q by Q_k=Q_k-1*Q;
      
      for(int k=0;k<n;k++)
	{
	  double tmp1=Q[k][i],tmp2=Q[k][i+1];
	  Q[k][i]=Q[k][i]*c+Q[k][i+1]*s;
	  Q[k][i+1]=tmp1*(-s)+tmp2*c;
	}
      
    }
  //multiply(X,Q); this update the X=V*Qk but in QR factorization, it does not contain this command;
  //multiply(a,Q); this update a=RQ, because after this QR function, a=R, after this command, a=RQ;
}

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
  rand_V(p, V);// at this time it will send X = V; cause X=VU;
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



