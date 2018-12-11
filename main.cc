#include "matricesEF.h"

int main()
{
  int choix = 1;
  int choixCL = 2;
  double epsilon = 0.0001;
  int Nk = 20;

  Eigen::SparseMatrix<double> M = createM(choix, Nk, epsilon);

  Eigen::VectorXd F = createF(choix, choixCL, Nk);

  // Eigen::SparseMatrix<double> M = createMsansCL(choix, Nk);
  //
  // Eigen::VectorXd F = createFbasique(choix, Nk);

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;

  solver.compute(M);

  Eigen::VectorXd U = solver.solve(F);

  saveSol("resultat", choix, Nk, U);
}
