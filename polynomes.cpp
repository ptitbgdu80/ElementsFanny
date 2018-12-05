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

Polynome2D pow(Polynome2D pol, int exposant)
{
  Polynome2D result;
  if (exposant < 0)
  {
    std::cout << "L'exposant d'un polynôme doit être supérieur à 0" << std::endl;
    exit(1);
  }
  else
  {
    result.degreX = exposant*pol.degreX;
    result.degreY = exposant*pol.degreY;
    result.coeffs.resize(result.degreX + 1);
    for (int i = 0; i <= result.degreX; i++)
    {
      result.coeffs[i].resize(result.degreY + 1);
      for (int j = 0; j <= result.degreY; j++)
      {
        result.coeffs[i][j] = 0;
      }
    }
    for (int i = 0; i <= pol.degreX; i++)
    {
      for (int j = 0; j <= pol.degreY; j++)
      {
        result.coeffs[exposant*i][exposant*j] = pow(pol.coeffs[i][j],exposant);
      }
    }

    checkDegre(result);
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
