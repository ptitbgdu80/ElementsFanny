#include "matricesEF.h"

int main()
{
  int choix = 2; // choix 1= (P0,Q1)  et choix 2= (Q1,Q2=)
  int choixCL = 2; //  choix CL=1 correspond à p=0 et u=(1,1) et choixCL=2 correspond p=x-y+c et u = (x,-y)
  double epsilon = 0.0001; // pour la matrice Espilon identité
  int Nk = 20; //nombre de mailles dans chaque direction

  Eigen::SparseMatrix<double> M = createM(choix, Nk, epsilon); // matrice M

  Eigen::VectorXd F = createF(choix, choixCL, Nk); //second membre

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;

  solver.compute(M);

  Eigen::VectorXd U = solver.solve(F); //vecteur solution U

  saveSol("resultat", choix, Nk, U);
}
