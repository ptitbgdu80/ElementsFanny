#include "matricesEF.h"

int main()
{
  int choix = 1;
  int choixCL = 1;
  double epsilon = 0.1;
  int Nk = 40;

  Eigen::SparseMatrix<double> M = createMavecCL(choix, Nk, epsilon);

  Eigen::VectorXd F = createFpourMavecCL(choix, Nk);

  // Eigen::SparseMatrix<double> M = createMsansCL(choix, Nk);
  //
  // Eigen::VectorXd F = createFbasique(choix, Nk);

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;

  solver.compute(M);

  Eigen::VectorXd U = solver.solve(F);

  saveSol("resultat", choix, Nk, U);
}
