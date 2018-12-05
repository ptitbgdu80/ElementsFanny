#include "polynomes.h"

int max(int a, int b)
{
  if (a > b)
  {
    return a;
  }
  else
  {
    return b;
  }
}

Polynome2D createZeroPol()
{
  Polynome2D result;
  result.degreX = 0;
  result.degreY = 0;
  result.coeffs.resize(1);
  result.coeffs[0].resize(1);
  result.coeffs[0][0] = 0;
  return result;
}

Polynome2D createPol(std::vector<std::vector<double> > coeffs)
{
  Polynome2D result;
  result.degreX = coeffs.size() - 1;
  result.degreY = coeffs[0].size() - 1;
  for (int i = 1; i <= result.degreX; i++)
  {
    if (coeffs[i].size() - 1 != result.degreY)
    {
      std::cout << "Les lignes de coeffs doivent être de même taille" << std::endl;
      exit(1);
    }
  }
  result.coeffs = coeffs;
  return result;
}

Polynome2D dx(Polynome2D pol)
{
  Polynome2D result;
  result.degreX = pol.degreX - 1;
  result.degreY = pol.degreY;
  if (pol.degreX == 0)
  {
    result.degreX = 0;
    result.degreY = 0;
  }

  result.coeffs.resize(result.degreX+1);
  for (int i = 0; i <= result.degreX; i++)
  {
    result.coeffs[i].resize(result.degreY+1);
  }

  for (int i = 0; i < pol.degreX; i++)
  {
    for (int j = 0; j <= pol.degreY; j++)
    {
      result.coeffs[i][j] = pol.coeffs[i+1][j]*(i+1);
    }
  }

  if (pol.degreX == 0)
  {
    result.coeffs[0][0] = 0;
  }

  return result;
}

Polynome2D dy(Polynome2D pol)
{
  Polynome2D result;
  result.degreY = pol.degreY - 1;
  result.degreX = pol.degreX;
  if (pol.degreY == 0)
  {
    result.degreX = 0;
    result.degreY = 0;
  }


  result.coeffs.resize(result.degreX + 1);
  for (int i = 0; i <= result.degreX; i++)
  {
    result.coeffs[i].resize(result.degreY + 1);
  }

  for (int i = 0; i <= pol.degreX; i++)
  {
    for (int j = 0; j < pol.degreY; j++)
    {
      result.coeffs[i][j] = pol.coeffs[i][j+1]*(j+1);
    }
  }

  if (pol.degreY == 0)
  {
    result.coeffs[0][0] = 0;
  }

  return result;
}

Polynome2D primxy(Polynome2D pol)
{
  Polynome2D result = createZeroPol();
  if (pol.degreX != 0 or pol.degreY != 0)
  {
    result.degreX = pol.degreX+1;
    result.degreY = pol.degreY+1;
    result.coeffs.resize(result.degreX + 1);
    for (int i = 0; i <= result.degreX; i++)
    {
      result.coeffs[i].resize(result.degreY + 1);
      for (int j = 0; j <= result.degreY; j++)
      {
        if (i == 0 or j == 0)
        {
          result.coeffs[i][j] = 0;
        }
        else
        {
          result.coeffs[i][j] = pol.coeffs[i-1][j-1]/(i*j);
        }
      }
    }
  }
  return result;
}

std::vector<Polynome2D> gradient(Polynome2D pol)
{
  std::vector<Polynome2D> result;
  result.resize(2);

  result[0] = dx(pol);
  result[1] = dy(pol);

  return result;
}

Polynome2D dotProduct(std::vector<Polynome2D> v1, std::vector<Polynome2D> v2)
{
  if (v1.size() != v2.size())
  {
    std::cout << "Le produit scalaire doit se faire avec deux vecteurs de même taille" << std::endl;
    exit(1);
  }
  Polynome2D result = createZeroPol();
  for (int i = 0; i < v1.size(); i++)
  {
    result = result + v1[i]*v2[i];
  }
  return result;
}

double integraleSurUnCarreUnitaire(Polynome2D pol)
{
  Polynome2D prim = primxy(pol);
  double result = 0;
  for (int i = 0; i <= prim.degreX; i++)
  {
    for (int j = 0; j <= prim.degreY; j++)
    {
      result += prim.coeffs[i][j];
    }
  }
  return result;
}

void checkDegre(Polynome2D &pol)
{
  bool baisseDeDegreX = true;
  while (baisseDeDegreX and pol.degreX != 0)
  {
    for (int j = 0; j <= pol.degreY; j++)
    {
      baisseDeDegreX = baisseDeDegreX and (pol.coeffs[pol.degreX][j] == 0);
    }
    if (baisseDeDegreX)
    {
      pol.degreX -= 1;
      pol.coeffs.resize(pol.degreX);
    }
  }

  bool baisseDeDegreY = true;
  while (baisseDeDegreY and pol.degreY != 0)
  {
    for (int i = 0; i <= pol.degreX; i++)
    {
      baisseDeDegreY = baisseDeDegreY and (pol.coeffs[i][pol.degreY] == 0);
    }
    if (baisseDeDegreY)
    {
      pol.degreY -= 1;
      for (int i = 0; i <= pol.degreX; i++)
      {
        pol.coeffs[i].resize(pol.degreY);
      }
    }
  }
}

Polynome2D operator+(Polynome2D pol1, Polynome2D pol2)
{
  Polynome2D result;
  result.degreX = max(pol1.degreX,pol2.degreX);
  result.degreY = max(pol1.degreY,pol2.degreY);

  result.coeffs.resize(result.degreX+1);
  for (int i = 0; i <= result.degreX; i++)
  {
    result.coeffs[i].resize(result.degreY+1);
    for (int j = 0; j <= result.degreY; j++)
    {
      result.coeffs[i][j] = 0;
    }
  }

  for (int i = 0; i <= pol1.degreX; i++)
  {
    for (int j = 0; j <= pol1.degreY; j++)
    {
      result.coeffs[i][j] += pol1.coeffs[i][j];
    }
  }

  for (int i = 0; i <= pol2.degreX; i++)
  {
    for (int j = 0; j <= pol2.degreY; j++)
    {
      result.coeffs[i][j] += pol2.coeffs[i][j];
    }
  }

  checkDegre(result);

  return result;
}

Polynome2D operator+(double lambda, Polynome2D pol)
{
  Polynome2D result = pol;
  result.coeffs[0][0] += lambda;
  return result;
}

Polynome2D operator+(Polynome2D pol, double lambda)
{
  Polynome2D result = pol;
  result.coeffs[0][0] += lambda;
  return result;
}

Polynome2D operator+(Polynome2D pol)
{
  return pol;
}

Polynome2D operator-(Polynome2D pol1, Polynome2D pol2)
{
  Polynome2D result;
  result.degreX = max(pol1.degreX,pol2.degreX);
  result.degreY = max(pol1.degreY,pol2.degreY);

  result.coeffs.resize(result.degreX+1);
  for (int i = 0; i <= result.degreX; i++)
  {
    result.coeffs[i].resize(result.degreY+1);
    for (int j = 0; j <= result.degreY; j++)
    {
      result.coeffs[i][j] = 0;
    }
  }

  for (int i = 0; i <= pol1.degreX; i++)
  {
    for (int j = 0; j <= pol1.degreY; j++)
    {
      result.coeffs[i][j] += pol1.coeffs[i][j];
    }
  }

  for (int i = 0; i <= pol2.degreX; i++)
  {
    for (int j = 0; j <= pol2.degreY; j++)
    {
      result.coeffs[i][j] -= pol2.coeffs[i][j];
    }
  }

  checkDegre(result);

  return result;
}

Polynome2D operator-(double lambda, Polynome2D pol)
{
  Polynome2D result = pol;
  for (int i = 0; i <= pol.degreX; i++)
  {
    for (int j = 0; j <= pol.degreY; j++)
    {
      result.coeffs[i][j] *= -1;
    }
  }
  result.coeffs[0][0] += lambda;
  return result;
}

Polynome2D operator-(Polynome2D pol, double lambda)
{
  Polynome2D result = pol;
  result.coeffs[0][0] -= lambda;
  return result;
}

Polynome2D operator-(Polynome2D pol)
{
  Polynome2D result = pol;
  for (int i = 0; i <= pol.degreX; i++)
  {
    for (int j = 0; j <= pol.degreY; j++)
    {
      result.coeffs[i][j] *= -1;
    }
  }
  return result;
}

Polynome2D operator*(Polynome2D pol1, Polynome2D pol2)
{
  Polynome2D result = createZeroPol();

  if ((pol1.degreX != 0 or pol1.degreY != 0) and (pol2.degreX != 0 or pol2.degreY != 0))
  {
    result.degreX = pol1.degreX + pol2.degreX;
    result.degreY = pol1.degreY + pol2.degreY;

    result.coeffs.resize(result.degreX+1);
    for (int i = 0; i <= result.degreX; i++)
    {
      result.coeffs[i].resize(result.degreY+1);
      for (int j = 0; j <= result.degreY; j++)
      {
        result.coeffs[i][j] = 0;
      }
    }

    for (int i1 = 0; i1 <= pol1.degreX; i1++)
    {
      for (int j1 = 0; j1 <= pol1.degreY; j1++)
      {
        for (int i2 = 0; i2 <= pol2.degreX; i2++)
        {
          for (int j2 = 0; j2 <= pol2.degreY; j2++)
          {
            result.coeffs[i1+i2][j1+j2] += pol1.coeffs[i1][j1]*pol2.coeffs[i2][j2];
          }
        }
      }
    }
  }

  return result;
}

Polynome2D operator*(double lambda, Polynome2D pol)
{
  Polynome2D result;
  if (lambda == 0)
  {
    result = createZeroPol();
  }
  else
  {
    result = pol;
    for (int i = 0; i <= pol.degreX; i++)
    {
      for (int j = 0; j <= pol.degreY; j++)
      {
        result.coeffs[i][j] *= lambda;
      }
    }
  }
  return result;
}

Polynome2D operator*(Polynome2D pol, double lambda)
{
  Polynome2D result;
  if (lambda == 0)
  {
    result = createZeroPol();
  }
  else
  {
    result = pol;
    for (int i = 0; i <= pol.degreX; i++)
    {
      for (int j = 0; j <= pol.degreY; j++)
      {
        result.coeffs[i][j] *= lambda;
      }
    }
  }
  return result;
}

std::ostream &operator<<(std::ostream &out, Polynome2D pol)
{
  bool premier = true;
  for (int i = 0; i <= pol.degreX; i++)
  {
    for (int j = 0; j <= pol.degreY; j++)
    {
      double coeff = floor(pol.coeffs[i][j]*100)/100.;
      if (coeff != 0 || (pol.degreX == 0 and pol.degreY ==0))
      {
        if (!premier and coeff > 0)
        {
          out << " + ";
        }
        else if(coeff < 0)
        {
          out << " - ";
          coeff = -coeff;
        }
        premier = false;
        if (coeff != 1 || (i == 0 and j == 0))
        {
          out << coeff;
        }

        if (i != 0)
        {
          out << "x";
          if (i != 1)
          {
            out << "^" << i;
          }
        }
        if (j != 0)
        {
          out << "y";
          if (j != 1)
          {
            out << "^" << j;
          }
        }
      }
    }
  }
  return out;
}

std::ostream &operator<<(std::ostream &out, std::vector<Polynome2D> vect)
{
  out << "(";
  for (int i = 0; i < vect.size(); i++)
  {
    out << vect[i];
    if (i == vect.size() - 1)
    {
      out << ")";
    }
    else
    {
      out << ";" << std::endl;
    }
  }
  return out;
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
      result[i][j] = - integraleSurUnCarreUnitaire(gradient(polVect2[i])[0]*polVect1[j]);
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
      result[i][j] = - integraleSurUnCarreUnitaire(gradient(polVect2[i])[1]*polVect1[j]);
    }
  }
  return result;
}

std::vector<std::vector<double> > createA(std::vector<std::vector<double> > Ak, int Nx)
{
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

  int dim = Ak.size();
  if (dim == 4)
  {
    for (int elementK = 0; elementK < (Nx-1)*(Nx-1); elementK++)
    {
      for(int i = 0; i < dim; i++)
      {
        for (int j = 0; j < dim; j++)
        {
          MatA[localToGlobalQ1(elementK,i+1,Nx)][localToGlobalQ1(elementK,j+1,Nx)] += Ak[i][j];
        }
      }
    }
  }

  if (dim == 9)
  {
    for (int elementK = 0; elementK < (Nx-1)*(Nx-1)/4; elementK++)
    {
      for(int i = 0; i < dim; i++)
      {
        for (int j = 0; j < dim; j++)
        {
          MatA[localToGlobalQ2(elementK,i+1,Nx)][localToGlobalQ2(elementK,j+1,Nx)] += Ak[i][j];
        }
      }
    }
  }
  return MatA;
}

int localToGlobalQ1(int elementK, int numeroSommet, int Nx) //Donne l'indice dans le maillage du sommet numeroSommet appartenant à l'élément elementK
{
  if (elementK < 0 or elementK > Nx*Nx-1)
  {
    std::cout << "L'élément " << elementK << " n'appartient pas au domaine" << std::endl;
    std::cout << "Il doit être compris entre 0 et " << Nx*Nx-1 << std::endl;
    exit(1);
  }

  if (numeroSommet < 1 or numeroSommet > 4)
  {
    std::cout << "Le numéro du sommet doit être compris entre 1 et 4, numeroSommet = " << numeroSommet << std::endl;
    exit(1);
  }

  int L = elementK/(Nx-1); //numéro de ligne de l'élément K
  int C = elementK%(Nx-1); //numéro de la colonne de l'élément K

  return L*Nx + C + (numeroSommet - 1)/2*Nx + (numeroSommet-1)%2;
}

int localToGlobalQ2(int elementK, int numeroSommet, int Nx) //Donne l'indice dans le maillage du sommet numeroSommet appartenant à l'élément elementK
{
  if (elementK < 0 or elementK > Nx*Nx-1)
  {
    std::cout << "L'élément " << elementK << " n'appartient pas au domaine" << std::endl;
    std::cout << "Il doit être compris entre 0 et " << Nx*Nx-1 << std::endl;
    exit(1);
  }

  if (numeroSommet < 1 or numeroSommet > 9)
  {
    std::cout << "Le numéro du sommet doit être compris entre 1 et 9, numeroSommet = " << numeroSommet << std::endl;
    exit(1);
  }

  int L = elementK/((Nx-1)/2); //numéro de ligne de l'élément K
  int C = elementK%((Nx-1)/2); //numéro de la colonne de l'élément K

  return 2*L*Nx + 2*C + (numeroSommet - 1)/3*Nx + (numeroSommet-1)%3;
}
