#include "polynomes.h"
#include <vector>
#include <iostream>
#include <ostream>
#include <math.h>

std::vector<Polynome2D> getP0PolVect();

std::vector<Polynome2D> getQ1PolVect();

std::vector<Polynome2D> getQ2PolVect();

std::vector<std::vector<double> > createAK(std::vector<Polynome2D> polVect);

std::vector<std::vector<double> > createB1K(std::vector<Polynome2D> polVect1, std::vector<Polynome2D> polVect2);

std::vector<std::vector<double> > createB2K(std::vector<Polynome2D> polVect1, std::vector<Polynome2D> polVect2);

std::vector<std::vector<double> > createA(std::vector<std::vector<double> > Ak, int Nk);

std::vector<std::vector<double> > createB1(std::vector<std::vector<double> > B1k, int Nk);

int localToGlobalQ1(int elementK, int numeroSommet, int Nk);

int localToGlobalQ2(int elementK, int numeroSommet, int Nk);
