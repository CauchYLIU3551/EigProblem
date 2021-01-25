#include <iostream>
#include <fstream>
#include <iterator>
#include <mintrace.h>
// Following header files are used in the step-3 in the examples of deal.II;
// Here just using these to get a large scale matrix as the test data.
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
// Finally, this is for output to a file and to the console:
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

using namespace dealii;

class Step3
{
public:
  Step3 ();

  void run ();
  void make_grid ();
  void setup_system ();
  void assemble_system ();

  //void multiply(std::vector<double>& x);
  std::vector<double> multiply(std::vector<double> x0);

  
  //void solve ();
  //
  void output_results () const;

  Triangulation<2>     triangulation;
  FE_Q<2>              fe;
  DoFHandler<2>        dof_handler;
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;
  
  //std::vector<double> x;

};

Step3::Step3 ()
  :
  fe (1),
  dof_handler (triangulation)
{}

void Step3::make_grid ()
{
  // First create the grid and refine all cells five times. Since the initial
  // grid (which is the square [-1,1]x[-1,1]) consists of only one cell, the
  // final grid has 32 times 32 cells, for a total of 1024.
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (5);
  // Unsure that 1024 is the correct number?  Let's see: n_active_cells
  // returns the number of active cells:
  std::cout << "Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl;
  
  std::cout << "Total number of cells: "
            << triangulation.n_cells()
            << std::endl;
  // Note the distinction between n_active_cells() and n_cells().
}

void Step3::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  std::cout << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;
  
  CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
  sparsity_pattern.copy_from(c_sparsity);

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}

void Step3::assemble_system ()
{

  QGauss<2>  quadrature_formula(2);

  FEValues<2> fe_values (fe, quadrature_formula,
                         update_values | update_gradients | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  DoFHandler<2>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      cell_matrix = 0;
      cell_rhs = 0;

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                 fe_values.shape_grad (j, q_point) *
                                 fe_values.JxW (q_point));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                          1 *
                          fe_values.JxW (q_point));

      cell->get_dof_indices (local_dof_indices);
      ///////////////////////
      // Following commands compute the interior of the matrix and store them into the matrix.
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             cell_matrix(i,j));

      // And again, we do the same thing for the right hand side vector.
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }

  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<2>(),
                                            boundary_values);

  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
}

//void Step3::multiply(std::vector<double>& x)
std::vector<double> Step3::multiply(std::vector<double> x0)
{
  std::vector<double> x(x0.size(),0);
  
  if (system_matrix.n()!=x0.size())
    {
      std::cout<<"Function multiply() ERROR: The sizes of matrix and vector are not same! Please check it!\n";
    }
  else
    {
      // std::cout<<"TEST POINT 2 \n";
      for(int k=0;k<system_matrix.m();k++)
	{
	  //std::cout<<"this is the "<<k<<"th iterations \n";
	  SparseMatrix<double>::const_iterator i=system_matrix.begin(k);
	  
	  while(i!=system_matrix.end(k))
	    {
	      x[k]+=i->value()*x0[i->column()];
	      ++i;
	    }
	}
      return x;
    }
}

//void Householder(std::vector<std::vector<double>> *A, std::vector<std::vector<double>> *P)
//void Householder(std::vector<std::vector<double>> *A, std::vector<std::vector<double>> *P)
void Householder(std::vector<std::vector<double>> A)
{
  int n=A.size();
  std::cout<<"This is the size of the matrix A"<<n<<" \n";
  for(int k=0;k<n-2;k++)
    {
      //      std::cout<<"Check Point 1!!! \n";
      double q=0;
      double alpha=0;
      // Computing the value of q in the Householder function;
      for(int i=k+1;i<n;i++)
	{
	  q+=A[i][k]*(A[i][k]);
	}

      // Computing the value of alpha
      if(A[k+1][k]==0)
	{
	  alpha=-sqrt(q);
	}
      else
	{
	  alpha=-sqrt(q)*A[k+1][k]/abs(A[k+1][k]);
	  std::cout<<"the abs A[k+1][k] is "<<A[k+1][k]<<std::endl;
	}
      
      double RSQ=alpha*alpha-alpha*A[k+1][k]; // RSQ = 2r^2;

      //std::cout<<"Check Point 2, This is RSQ:::"<<RSQ<<"\n";
	
      std::vector<double> v(n-k), u(n-k,0), z(n-k);
      v[0]=0;
      v[1]=A[k+1][k]-alpha;

      //std::cout<<"Check Point 3::: "<<v[1]<<" \n";
      // w=(1/sqrt(2*RSQ)*v=1/2r*v);
      for(int j=2;j<v.size();j++)
	{
	  v[j]=A[j][k];
	}
      
      // std::cout<<"CheckPPPPPoint555 \n";
      
      for(int j=0;j<n-k;j++)
	{
	  for(int i=1;i<n-k;i++)
	    {
	      u[j]+=A[j][i]*v[i];
	    }
	  u[j]/=RSQ;
	}

      //std::cout<<"Check POOOIONIOHIHO\n";
      
      double PROD=0;
      for(int i=1;i<n-k;i++)
	{
	  PROD+=v[i]*u[i];
	}

      //std::cout<<"CheckPoint 4!!!!!!!!!\n";
      
      for(int j=0;j<n-k;j++)
	{
	  z[j]=u[j]-PROD/(2*RSQ)*v[j];
	}

      for(int l=k+1;l<n-1;l++)
	{
	  for(int j=l+1;j<n;j++)
	    {
	      A[j][l]=A[j][l]-v[j-k]*z[j-k]-v[j-k]*z[l-k];
	      A[l][j]=A[j][l];
	    }
	  A[l][l]=A[l][l]-2*v[l-k]*z[l-k];
	}
      for(int j=k+2;j<n;j++)
	{
	  A[k][j]=0;
	  A[j][k]=0;
	}

      A[k+1][k]=A[k+1][k]-v[1]*z[0];
      A[k][k+1]=A[k+1][k];

      std::cout<<"This is the result of "<<k<<"th iteration \n";
      std::cout<<"The alpha is "<<alpha<<"\n";
      std::cout<<"The r is :: "<<sqrt(RSQ/2)<<"\n";
      for(int i=0;i<A.size();i++)
	{
	  for(int j=0; j<A.size();j++)
	    {
	      std::cout<<A[i][j]<<" ";
	    }
	  std::cout<<"\n";
	}
      std::cout<<"The vector of partial w is \n";
      for (int i=0;i<v.size();i++)
	{
	  std::cout<<v[i]/sqrt(2*RSQ)<<" ";
	}
      std::cout<<"\n";

      std::cout<<"The vector of u is \n";
      for (int i=0;i<u.size();i++)
	{
	  std::cout<<u[i]<<" ";
	}
      std::cout<<"\n";
      
    }



}


void Step3::run ()
{
  make_grid ();
  setup_system();
  assemble_system ();
  std::cout<<system_matrix.m()<<"\n";
  std::cout<<system_matrix.n()<<"\n";

  //std::cout<<system_matrix.begin()->value()<<"\n";
  //std::cout<<system_matrix.begin()->column()<<"\n";

  std::vector<double> x(system_matrix.n(),1);
  //
  // Due to the system_matrix is a member of the class step3, but x is not a member of the class,
  // So we can not directly define the multiply function out of the class and using the
  // system_matrix directly.
  // we need to get a variable A equals to the system_matrix but A is not a member of the class;
  
  std::cout << "CheckPoint 1 \n";

  // multiply(x);
  x=multiply(x);
  //std::cout<<"This is the value of the vector x"<<x[0]<<std::endl;
  std::ofstream out ("sparse_matrix");
  system_matrix.print(out);

  std::ofstream output_file ("solution");
  for (const auto &e : x) output_file << e << "\n";
  
}



int main()
{
  Step3 laplace_problem;
  /////laplace_problem.run ();

  
  //std::stringstream result;
  std::vector<double> b(4,0);
  std::vector<std::vector<double>> A(4, b);
  A[0][0]=4;
  A[0][1]=1;
  A[0][2]=-2;
  A[0][3]=2;
  A[1][0]=1;
  A[1][1]=2;
  A[1][2]=0;
  A[1][3]=1;
  A[2][0]=-2;
  A[2][1]=0;
  A[2][2]=3;
  A[2][3]=-2;
  A[3][0]=2;
  A[3][1]=1;
  A[3][2]=-2;
  A[3][3]=-1;
  Householder(A);


  std::cout<<"hello world!"<<std::endl;
  //SparseMatrix<double> A;
  return 0;
}
