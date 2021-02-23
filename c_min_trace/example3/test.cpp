#include <trace/mintrace.h>
#include <CG/CGSolver.h>

int main()
{
  //CGSolver aaa;
  SparseMatrix<double> A;
  TraceSolver B;
  std::vector<double> temp(4,0);
  std::vector<std::vector<double>>a(4,temp),Q(4,temp);
  Q[0][0]=1;
  Q[1][1]=1;
  Q[2][2]=1;
  Q[3][3]=1;
  a[0][0]=4;
  a[0][1]=1;
  a[0][2]=-2;
  a[0][3]=2;
  a[1][0]=1;
  a[1][1]=2;
  a[1][3]=1;
  a[2][0]=-2;
  a[2][2]=3;
  a[2][3]=-2;
  a[3][0]=2;
  a[3][1]=1;
  a[3][2]=-2;
  a[3][3]=-1;
  B.Householder(a);
  for(int i=0;i<a.size();i++)
    {
      for(int j=0;j<a.size();j++)
	{
	  std::cout<<a[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  //identitymatrix(Q,a.size());
  B.QR(a, Q);
  for(int i=0;i<a.size();i++)
    {
      for(int j=0;j<a.size();j++)
	{
	  std::cout<<a[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"\n";
    for(int i=0;i<a.size();i++)
    {
      for(int j=0;j<a.size();j++)
	{
	  std::cout<<Q[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"Hello world!\n";
  return 0;
}
