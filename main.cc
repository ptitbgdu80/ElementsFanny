#include "matricesEF.h"

int main()
{
  Eigen::SparseMatrix<double> M = createMavecCL(2,2);

  Eigen::VectorXd F = createFpourMavecCL(2, 2);

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;

  solver.compute(M);

  Eigen::VectorXd U = solver.solve(F);

  std::cout << U << std::endl;
}
