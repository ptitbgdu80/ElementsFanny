#include "matricesEF.h"

std::vector<Polynome2D> getP0PolVect()
{
  Polynome2D psi = createPol({{1}});
  return {psi};
}

std::vector<Polynome2D> getQ1PolVect()
{
  std::vector<std::vector<double> > v1 {{1,-1},{-1,1}}, v2 {{0,0},{1,-1}}, v3 {{0,1},{0,-1}}, v4 {{0,0},{0,1}};

  Polynome2D phi1 = createPol(v1);
  Polynome2D phi2 = createPol(v2);
  Polynome2D phi3 = createPol(v3);
  Polynome2D phi4 = createPol(v4);

  std::vector<Polynome2D> result {phi1, phi2, phi3, phi4};

  return result;
}

std::vector<Polynome2D> getQ2PolVect()
{
  Polynome2D x = createPol({{0},{1}});
  Polynome2D y = createPol({{0,1}});
  Polynome2D phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9;

  phi1 = 4*(x-1)*(x-0.5)*(y-0.5)*(y-1);
  phi2 = -8*x*(x-1)*(y-0.5)*(y-1);
  phi3 = 4*x*(x-0.5)*(y-0.5)*(y-1);
  phi4 = -8*(x-1)*(x-0.5)*y*(y-1);
  phi5 = 16*(x-1)*x*y*(y-1);
  phi6 = -8*x*(x-0.5)*y*(y-1);
  phi7 = 4*(x-1)*(x-0.5)*(y-0.5)*y;
  phi8 = -8*(x-1)*x*(y-0.5)*y;
  phi9 = 4*x*(x-0.5)*y*(y-0.5);

  std::vector<Polynome2D> result {phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9};

  return result;
}

std::vector<std::vector<double> > createAK(std::vector<Polynome2D> polVect)
{
  std::vector<std::vector<double> > result;
  result.resize(polVect.size());
  for (int i = 0; i < polVect.size(); i++)
  {
    result[i].resize(polVect.size());
    for (int j = 0; j < polVect.size(); j++)
    {
      result[i][j] = integraleSurUnCarreUnitaire(dotProduct(gradient(polVect[i]),gradient(polVect[j])));
    }
  }
  return result;
}

std::vector<std::vector<double> > createB1K(std::vector<Polynome2D> polVect1, std::vector<Polynome2D> polVect2)
{
  std::vector<std::vector<double> > result;
  result.resize(polVect2.size());
  for (int i = 0; i < polVect2.size(); i++)
  {
    result[i].resize(polVect1.size());
    for (int j = 0; j < polVect1.size(); j++)
    {
      result[i][j] = - integraleSurUnCarreUnitaire(dx(polVect2[i])*polVect1[j]);
    }
  }
  return result;
}

std::vector<std::vector<double> > createB2K(std::vector<Polynome2D> polVect1, std::vector<Polynome2D> polVect2)
{
  std::vector<std::vector<double> > result;
  result.resize(polVect2.size());
  for (int i = 0; i < polVect2.size(); i++)
  {
    result[i].resize(polVect1.size());
    for (int j = 0; j < polVect1.size(); j++)
    {
      result[i][j] = - integraleSurUnCarreUnitaire(dy(polVect2[i])*polVect1[j]);
    }
  }
  return result;
}

std::vector<std::vector<double> > createA(std::vector<std::vector<double> > Ak, int Nk) //Nk nombre d'éléments (ou de mailles) par ligne
{
  if (Nk < 1)
  {
    std::cout << "Le nombre d'éléments par ligne doit être positif pour créer A" << std::endl;
    exit(1);
  }
  int dim = Ak.size();

  if (Ak[0].size() != dim)
  {
    std::cout << "La matrice Ak n'est pas carrée" << std::endl;
    exit(1);
  }

  int Nx;

  if (dim == 4) //cas Q1
  {
    Nx = Nk +1;
  }
  else if (dim == 9) //cas Q2
  {
    Nx = 2*Nk + 1;
  }

  std::vector<std::vector<double> > MatA;

  MatA.resize(Nx*Nx);
  for (int i = 0; i < Nx*Nx; i++)
  {
    MatA[i].resize(Nx*Nx);
    for (int j = 0; j < Nx*Nx; j++)
    {
      MatA[i][j] = 0;
    }
  }

  if (dim == 4) //cas Q1
  {
    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim; i++)
      {
        for (int j = 0; j < dim; j++)
        {
          MatA[localToGlobalQ1(elementK,i+1,Nk)][localToGlobalQ1(elementK,j+1,Nk)] += Ak[i][j];
        }
      }
    }
  }

  else if (dim == 9) //cas Q2
  {
    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim; i++)
      {
        for (int j = 0; j < dim; j++)
        {
          MatA[localToGlobalQ2(elementK,i+1,Nk)][localToGlobalQ2(elementK,j+1,Nk)] += Ak[i][j];
        }
      }
    }
  }
  return MatA;
}

std::vector<std::vector<double> > createB1(std::vector<std::vector<double> > B1k, int Nk)
{
  if (Nk < 1)
  {
    std::cout << "Le nombre d'éléments par ligne doit être positif pour créer B1" << std::endl;
    exit(1);
  }

  int dim1 = B1k.size();
  int dim2 = B1k[0].size();

  int Nx1, Nx2;

  if (dim1 == 4 and dim2 == 1) //cas (P0,Q1)
  {
    Nx1 = Nk +1;
    Nx2 = Nk;
  }
  else if (dim1 == 9 and dim2 == 4) //cas (Q1,Q2)
  {
    Nx1 = 2*Nk + 1;
    Nx2 = Nk + 1;
  }

  std::vector<std::vector<double> > MatB1;

  MatB1.resize(Nx1*Nx1);
  for (int i = 0; i < Nx1*Nx1; i++)
  {
    MatB1[i].resize(Nx2*Nx2);
    for (int j = 0; j < Nx2*Nx2; j++)
    {
      MatB1[i][j] = 0;
    }
  }

  if (dim1 == 4 and dim2 ==1) //cas (P0,Q1)
  {
    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim1; i++)
      {
        for (int j = 0; j < dim2; j++)
        {
          MatB1[localToGlobalQ1(elementK,i+1,Nk)][j] += B1k[i][j];
        }
      }
    }
  }

  else if (dim1 == 9 and dim2 == 4) //cas Q2
  {
    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim1; i++)
      {
        for (int j = 0; j < dim2; j++)
        {
          MatB1[localToGlobalQ2(elementK,i+1,Nk)][localToGlobalQ1(elementK,j+1,Nk)] += B1k[i][j];
        }
      }
    }
  }
  return MatB1;
}

int localToGlobalQ1(int elementK, int numeroSommet, int Nk) //Donne l'indice dans le maillage du sommet numeroSommet appartenant à l'élément elementK, Nk nombre d'éléments par ligne
{
  if (elementK < 0 or elementK > Nk*Nk-1)
  {
    std::cout << "L'élément " << elementK << " n'appartient pas au domaine" << std::endl;
    std::cout << "Il doit être compris entre 0 et " << Nk-1 << std::endl;
    exit(1);
  }

  if (numeroSommet < 1 or numeroSommet > 4)
  {
    std::cout << "Le numéro du sommet doit être compris entre 1 et 4, numeroSommet = " << numeroSommet << std::endl;
    exit(1);
  }

  int Nx = Nk+1;
  int L = elementK/Nk; //numéro de ligne de l'élément K
  int C = elementK%Nk; //numéro de la colonne de l'élément K

  return L*Nx + C + (numeroSommet - 1)/2*Nx + (numeroSommet-1)%2;
}

int localToGlobalQ2(int elementK, int numeroSommet, int Nk) //Donne l'indice dans le maillage du sommet numeroSommet appartenant à l'élément elementK, Nk nombre d'éléments par ligne
{
  if (elementK < 0 or elementK > Nk*Nk-1)
  {
    std::cout << "L'élément " << elementK << " n'appartient pas au domaine" << std::endl;
    std::cout << "Il doit être compris entre 0 et " << Nk*Nk-1 << std::endl;
    exit(1);
  }

  if (numeroSommet < 1 or numeroSommet > 9)
  {
    std::cout << "Le numéro du sommet doit être compris entre 1 et 9, numeroSommet = " << numeroSommet << std::endl;
    exit(1);
  }

  int Nx = 2*Nk+1;
  int L = elementK/Nk; //numéro de ligne de l'élément K
  int C = elementK%Nk; //numéro de la colonne de l'élément K

  return 2*L*Nx + 2*C + (numeroSommet - 1)/3*Nx + (numeroSommet-1)%3;
}
