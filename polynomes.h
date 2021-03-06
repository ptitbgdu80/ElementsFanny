#include <vector>
#include <iostream>
#include <ostream>
#include <math.h>

struct Polynome2D
{
  int degreX, degreY;
  std::vector<std::vector<double> > coeffs;
};

int max(int a, int b);

Polynome2D createZeroPol();

Polynome2D createPol(std::vector<std::vector<double> > coeffs);

Polynome2D dx(Polynome2D pol);

Polynome2D dy(Polynome2D pol);

Polynome2D primxy(Polynome2D pol);

std::vector<Polynome2D> gradient(Polynome2D pol);

Polynome2D dotProduct(std::vector<Polynome2D> v1, std::vector<Polynome2D> v2);

double integraleSurUnCarreUnitaire(Polynome2D pol);

void checkDegre(Polynome2D &pol);

Polynome2D operator+(Polynome2D pol1, Polynome2D pol2);

Polynome2D operator+(double lambda, Polynome2D pol);

Polynome2D operator+(Polynome2D pol, double lambda);

Polynome2D operator+(Polynome2D pol);

Polynome2D operator-(Polynome2D pol1, Polynome2D pol2);

Polynome2D operator-(double lambda, Polynome2D pol);

Polynome2D operator-(Polynome2D pol, double lambda);

Polynome2D operator-(Polynome2D pol);

Polynome2D operator*(Polynome2D pol1, Polynome2D pol2);

Polynome2D operator*(double lambda, Polynome2D pol);

Polynome2D operator*(Polynome2D pol, double lambda);

Polynome2D pow(Polynome2D pol, int exposant);

std::ostream &operator<<(std::ostream &out, Polynome2D pol);

std::ostream &operator<<(std::ostream &out, std::vector<Polynome2D> vect);
