#include <trace/mintrace.h>
#include <CG/CGSolver.h>

int main()
{
  SparseMatrix<double> A, M;
  // std::ofstream sparsematrix1 ("original_matrix.1");
  //A.print(sparsematrix1);
  std::vector<unsigned int> row_length(10,10);
  unsigned int t=3;
  SparsityPattern sparsity_pattern(10,10,{2,3,3,3,3,3,3,3,3,2});
  SparsityPattern sp2(10,10,1);
  //sparsity_pattern.add(1,2);
  sparsity_pattern.add(0,1);
  sparsity_pattern.add(9,8);
  for (int i=1;i<9;i++)
    {
      sparsity_pattern.add(i,i+1);
      sparsity_pattern.add(i,i-1);
    }
  
  A.reinit(sparsity_pattern);
  M.reinit(sp2);

  // In this way, I can construct a SparseMatrix to be the test data for the algorithm.
  for (int k=0;k<A.m();k++)
    {
      SparseMatrix<double>::iterator i=A.begin(k);
      i->value()=2;
      while(i!=A.end(k))
	{
	  i->value()+=1;
	  ++i;
	}
    }

  SparseMatrix<double>::iterator i=M.begin();
  double num=0;
  while(i!=M.end())
    {
      num=num+1;
      i->value()=num;
      ++i;
    }
  
  //CGSolver aaa;
  //SparseMatrix<double> A;
  TraceSolver B;
  TraceSolver C(A,M);
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

  for(int i=0;i<a.size();i++)
    {
      for(int j=0;j<a.size();j++)
	{
	  std::cout<<a[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"\n";
  
  B.QRSolver(a);
  for(int i=0;i<a.size();i++)
    {
      for(int j=0;j<a.size();j++)
	{
	  std::cout<<a[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"\n";

  std::vector<double> tempc(3,0);
  std::vector<std::vector<double>> tempC(10,tempc);
  tempC[0][0]=1;
  tempC[1][1]=sqrt(2)/2;
  tempC[2][2]=sqrt(3)/3;
  C.X=tempC;
  for(int i=0;i<C.X.size();i++)
    {
      for(int j=0;j<C.X[0].size();j++)
	{
	  std::cout<<C.X[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"\n";
  C.get_MX();
  
  for(int i=0;i<C.MX.size();i++)
    {
      for(int j=0;j<C.MX[0].size();j++)
	{
	  std::cout<<C.MX[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"\n";
  /*
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
    */
  std::cout<<"Hello world!\n";
  std::ofstream sparsematrix1 ("original_matrix.1");
  A.print(sparsematrix1);
  std::ofstream sparsematrix2 ("original_matrix.2");
  M.print(sparsematrix2);
  return 0;
}
