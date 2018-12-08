#include "polynomes.h"
#include <vector>
#include <iostream>
#include <ostream>
#include <math.h>
#include "Sparse"
#include "Dense"

std::vector<Polynome2D> getP0PolVect();

std::vector<Polynome2D> getQ1PolVect();

std::vector<Polynome2D> getQ2PolVect();

std::vector<std::vector<double> > createAK(std::vector<Polynome2D> polVect);

std::vector<std::vector<double> > createB1K(std::vector<Polynome2D> polVect1, std::vector<Polynome2D> polVect2);

std::vector<std::vector<double> > createB2K(std::vector<Polynome2D> polVect1, std::vector<Polynome2D> polVect2);

std::vector<double> createFK(std::vector<Polynome2D> polVect);

void insertSource(int Nk, Eigen::VectorXd &F);

std::vector<double> CLvitesse (double x, double y);

double CLpression (double x, double y);

Eigen::VectorXd createFpourMavecCL(int choix, int Nk);

void insertAsansCL(std::vector<std::vector<double> > Ak, int Nk, Eigen::SparseMatrix<double> &M);

void insertB1B2sansCL(std::vector<std::vector<double> > B1k, std::vector<std::vector<double> > B2k, int Nk,Eigen::SparseMatrix<double> &M);

Eigen::SparseMatrix<double> createMsansCL(int choix, int Nk);

void insertAavecCL(std::vector<std::vector<double> > Ak, int Nk, Eigen::SparseMatrix<double> &M);

void insertB1B2avecCL(std::vector<std::vector<double> > B1k, std::vector<std::vector<double> > B2k, int Nk,Eigen::SparseMatrix<double> &M);

Eigen::SparseMatrix<double> createMavecCL(int choix, int Nk);

int localToGlobalQ1(int elementK, int numeroSommet, int Nk);

int localToGlobalQ2(int elementK, int numeroSommet, int Nk);
