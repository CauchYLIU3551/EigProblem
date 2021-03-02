#include<cmath>
#include<lac/sparsity_pattern.h>
#include<lac/sparse_matrix.h>
// Multiplication of SparseMatrix A and vector x0;
std::vector<double> multiply(const SparseMatrix<double>* A, std::vector<double> x0)
{
  std::vector<double> x(x0.size(),0);
  //  std::cout<<"Flag 1!!\n";
  //std::cout<<A.n()<<"\n";
  
  if (A->n()!=x0.size())
    {
      std::cout<<"Function multiply() ERROR: The sizes of matrix and vector are not same! Please check it!\n";
    }
  else
    {
      // std::cout<<"TEST POINT 2 \n";
      for(int k=0;k<A->m();k++)
	{
	  //std::cout<<"this is the "<<k<<"th iterations \n";
	  SparseMatrix<double>::const_iterator i=A->begin(k);
	  
	  while(i!=A->end(k))
	    {
	      x[k]+=i->value()*x0[i->column()];
	      ++i;
	    }
	}
      return x;
    }
}
