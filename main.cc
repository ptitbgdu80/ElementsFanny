#include "matricesEF.h"

int main()
{
  int choix = 2;
  int Nk = 20;

  Eigen::SparseMatrix<double> M = createMsansCL(choix, Nk);

  Eigen::VectorXd F = createFpourMavecCL(choix, Nk);

  std::vector<double> Fk=createFK(getQ2PolVect());

  insertSource(Fk,Nk,F);

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;

  solver.compute(M);

  Eigen::VectorXd U = solver.solve(F);

  createVTK("resultat", choix, Nk, U);
}
