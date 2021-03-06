////////////////////////////////////////////////////////////////////////////////////////////
// main1.cpp :
//

#include <iostream>
#include <fstream>

#include <AFEPack/AMGSolver.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>
#include <AFEPack/EasyMesh.h>

#include <CG/CGSolver.h>

#define PI (4.0*atan(1.0))

double u(const double *);
double f(const double *);

clock_t start, end, start2, end2;

int main(int argc, char * argv[])
{
  EasyMesh mesh;
  mesh.readData(argv[1]);

  TemplateGeometry<2>	triangle_template_geometry;
  triangle_template_geometry.readData("triangle.tmp_geo");
  CoordTransform<2,2>	triangle_coord_transform;
  triangle_coord_transform.readData("triangle.crd_trs");
  TemplateDOF<2>	triangle_template_dof(triangle_template_geometry);
  triangle_template_dof.readData("triangle.1.tmp_dof");
  BasisFunctionAdmin<double,2,2> triangle_basis_function(triangle_template_dof);
  triangle_basis_function.readData("triangle.1.bas_fun");

  std::vector<TemplateElement<double,2,2> > template_element(1);
  template_element[0].reinit(triangle_template_geometry,
			     triangle_template_dof,
			     triangle_coord_transform,
			     triangle_basis_function);

  FEMSpace<double,2> fem_space(mesh, template_element);
	
  int n_element = mesh.n_geometry(2);
  fem_space.element().resize(n_element);
  for (int i = 0;i < n_element;i ++)
    fem_space.element(i).reinit(fem_space,i,0);

  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  StiffMatrix<2,double> stiff_matrix(fem_space);
  stiff_matrix.algebricAccuracy() = 4;
  stiff_matrix.build();

  FEMFunction<double,2> solution(fem_space), solution2(fem_space);
  Vector<double> right_hand_side;
  Operator::L2Discretize(&f, fem_space, right_hand_side, 4);

  BoundaryFunction<double,2> boundary(BoundaryConditionInfo::DIRICHLET, 1, &u);
  BoundaryConditionAdmin<double,2> boundary_admin(fem_space);
  boundary_admin.add(boundary);
  boundary_admin.apply(stiff_matrix, solution, right_hand_side);
  boundary_admin.apply(stiff_matrix, solution2, right_hand_side);

  AMGSolver solver(stiff_matrix);
  start=clock();
  solver.solve(solution, right_hand_side, 1.0e-08, 200);
  end=clock();

  CGSolver sol2(stiff_matrix);
  start2=clock();
  sol2.solve(solution2, right_hand_side, 1.0e-08, 200);
  end2=clock();
  
  solution.writeOpenDXData("u.dx");
  double error = Functional::L2Error(solution, FunctionFunction<double>(&u), 3);
  std::cerr << "\nL2 error = " << error << std::endl;

  double error2 = Functional::L2Error(solution2, FunctionFunction<double>(&u), 3);
  std::cerr << "\nL2 error2 from CG Solver = " << error2 << std::endl;

  double endtime=(double)(end-start)/CLOCKS_PER_SEC;
  double endtime2=(double)(end2-start2)/CLOCKS_PER_SEC;
  std::cout<<"The computing time of AMGSolver is :::"<<1000*endtime<<" ms \n";
  std::cout<<"The computing time of CGSolver is :::"<<1000*endtime2<<" ms \n"; 
  
  return 0;
};

double u(const double * p)
{
  return sin(PI*p[0]) * sin(2*PI*p[1]);
};

double f(const double * p)
{
  return 5*PI*PI*u(p);
};

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
