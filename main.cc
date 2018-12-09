#include "matricesEF.h"

int main()
{
  int choix = 1;
  int Nk = 20;

  Eigen::SparseMatrix<double> M = createMavecCL(choix, Nk);

  Eigen::VectorXd F = createFpourMavecCL(choix, Nk);

  // Eigen::SparseMatrix<double> M = createMsansCL(choix, Nk);
  //
  // Eigen::VectorXd F = createFbasique(choix, Nk);

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;

  solver.compute(M);

  Eigen::VectorXd U = solver.solve(F);

  createVTK("resultat", choix, Nk, U);
}
