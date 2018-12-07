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

std::vector<double> createFK(std::vector<Polynome2D> polVect)
{
  std::vector<double> result;
  result.resize(polVect.size());
  for (int i = 0; i < polVect.size(); i++)
  {
        result[i]= integraleSurUnCarreUnitaire(polVect[i]);

  }
  return result;
}

std::vector<double> createF(std::vector<double> Fk, int Nk) //Nk nombre d'éléments (ou de mailles) par ligne
{
  std::vector<double> F;
  if (Nk < 1)
  {
    std::cout << "Le nombre d'éléments par ligne doit être positif pour créer F" << std::endl;
    exit(1);
  }
  int dim = Fk.size();

  // if (Ak[0].size() != dim)
  // {
  //   std::cout << "La matrice Ak n'est pas carrée" << std::endl;
  //   exit(1);
  // }

  int Nx;

  if (dim == 4) //cas Q1
  {
    Nx = Nk +1;
  }
  else if (dim == 9) //cas Q2
  {
    Nx = 2*Nk + 1;
  }
  else
  {
    std::cout << "La matrice Fk n'a pas une dimension correspondant à Q1 ou Q2" << std::endl;
    exit(1);
  }

  if (dim == 4) //cas Q1
  {
    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim; i++)
      {
        F[localToGlobalQ1(elementK,i+1,Nk)]+=Fk[i];
      }
    }
  }

  else if (dim == 9) //cas (Q1,Q2)
  {
    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim; i++)
      {
        F[localToGlobalQ2(elementK,i+1,Nk)]+=Fk[i];
      }
    }
  }
  return F;
}

void insertA(std::vector<std::vector<double> > Ak, int Nk, Eigen::SparseMatrix<double> &M) //Nk nombre d'éléments (ou de mailles) par ligne
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
  else
  {
    std::cout << "La matrice Ak n'a pas une dimension correspondant à Q1 ou Q2" << std::endl;
    exit(1);
  }

  if (dim == 4) //cas Q1
  {
    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim; i++)
      {
        for (int j = 0; j < dim; j++)
        {
          M.coeffRef(localToGlobalQ1(elementK,i+1,Nk),localToGlobalQ1(elementK,j+1,Nk))+=Ak[i][j];
          M.coeffRef(localToGlobalQ1(elementK,i+1,Nk)+Nx*Nx,localToGlobalQ1(elementK,j+1,Nk)+Nx*Nx)+=Ak[i][j];
        }
      }
    }
  }

  else if (dim == 9) //cas (Q1,Q2)
  {
    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim; i++)
      {
        for (int j = 0; j < dim; j++)
        {
          M.coeffRef(localToGlobalQ2(elementK,i+1,Nk),localToGlobalQ2(elementK,j+1,Nk))+=Ak[i][j];
          M.coeffRef(localToGlobalQ2(elementK,i+1,Nk)+Nx*Nx,localToGlobalQ2(elementK,j+1,Nk)+Nx*Nx)+=Ak[i][j];
        }
      }
    }
  }
}

void insertB1B2(std::vector<std::vector<double> > B1k, std::vector<std::vector<double> > B2k, int Nk,Eigen::SparseMatrix<double> &M)
{
  if (Nk < 1)
  {
    std::cout << "Le nombre d'éléments par ligne doit être positif pour créer B1" << std::endl;
    exit(1);
  }

  int dim1 = B1k.size();
  int dim2 = B1k[0].size();

  if (B2k.size() != dim1 or B2k[0].size() != dim2)
  {
    std::cout << "Les matrices B1k et B2k n'ont pas les mêmes dimensions" << std::endl;
    exit(1);
  }

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
  else
  {
    std::cout << "La matrice Bk n'a pas une dimension correspondant à (P0,Q1) ou (Q1,Q2)" << std::endl;
    exit(1);
  }


  if (dim1 == 4 and dim2 ==1) //cas (P0,Q1)
  {
    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim1; i++)
      {
        for (int j = 0; j < dim2; j++)
        {
          M.coeffRef(localToGlobalQ1(elementK,i+1,Nk),elementK+2*Nx1*Nx1) += B1k[i][j];
          M.coeffRef(elementK+2*Nx1*Nx1,localToGlobalQ1(elementK,i+1,Nk)) += B1k[i][j];
          M.coeffRef(localToGlobalQ1(elementK,i+1,Nk)+Nx1*Nx1,elementK+2*Nx1*Nx1) += B2k[i][j];
          M.coeffRef(elementK+2*Nx1*Nx1,localToGlobalQ1(elementK,i+1,Nk)+Nx1*Nx1) += B2k[i][j];
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
          M.coeffRef(localToGlobalQ2(elementK,i+1,Nk),localToGlobalQ1(elementK,j+1,Nk)+2*Nx1*Nx1) += B1k[i][j];
          M.coeffRef(localToGlobalQ1(elementK,j+1,Nk)+2*Nx1*Nx1,localToGlobalQ2(elementK,i+1,Nk)) += B1k[i][j];
          M.coeffRef(localToGlobalQ2(elementK,i+1,Nk)+Nx1*Nx1,localToGlobalQ1(elementK,j+1,Nk)+2*Nx1*Nx1) += B2k[i][j];
          M.coeffRef(localToGlobalQ1(elementK,j+1,Nk)+2*Nx1*Nx1,localToGlobalQ2(elementK,i+1,Nk)+Nx1*Nx1) += B2k[i][j];
        }
      }
    }
  }
}

Eigen::SparseMatrix<double> createM(int choix, int Nk)
{
  int Nx1,Nx2;
  std::vector<std::vector<double> > AK, B1K, B2K;
  Eigen::SparseMatrix<double> M;

  switch (choix) {
    case 1:
    AK = createAK(getQ1PolVect());
    B1K = createB1K(getP0PolVect(),getQ1PolVect());
    B2K = createB2K(getP0PolVect(),getQ1PolVect());
    Nx1=Nk+1;
    Nx2=Nk;
    break;
    case 2:
    AK = createAK(getQ2PolVect());
    B1K = createB1K(getQ1PolVect(),getQ2PolVect());
    B2K = createB2K(getQ1PolVect(),getQ2PolVect());
    Nx1=2*Nk+1;
    Nx2=Nk+1;
    break;
    default:
    std::cout<<"Le choix doit être 1 ou 2"<<std::endl;
    exit(1);
  }

  M.resize(2*(Nx1*Nx1)+Nx2*Nx2,2*(Nx1*Nx1)+Nx2*Nx2);

  insertA (AK, Nk, M);
  insertB1B2(B1K,B2K,Nk,M);

  return M;
}

void insertAavecCL(std::vector<std::vector<double> > Ak, int Nk, Eigen::SparseMatrix<double> &M) //Nk nombre d'éléments (ou de mailles) par ligne
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

  if (dim == 4) //cas (P0,Q1)
  {
    Nx = Nk +1;
  }
  else if (dim == 9) //cas (Q1,Q2)
  {
    Nx = 2*Nk + 1;
  }
  else
  {
    std::cout << "La matrice Ak n'a pas une dimension correspondant à Q1 ou Q2" << std::endl;
    exit(1);
  }

  for (int elementK = 0; elementK < Nk*Nk; elementK++)
  {
    for(int i = 0; i < dim; i++)
    {
      int iA;
      if (dim == 4) //cas (P0,Q1)
      {
        iA = localToGlobalQ1(elementK,i+1,Nk);
      }
      else if (dim == 9) //cas (Q1,Q2)
      {
        iA = localToGlobalQ2(elementK,i+1,Nk);
      }

      if (iA < Nx or iA%Nx == 0 or iA%Nx == Nx-1 or iA > Nx*(Nx-1) - 1) //On est sur un bord, il faut imposer la valeur de la vitesse
      {
        M.coeffRef(iA,iA) = 1;
        M.coeffRef(iA+Nx*Nx,iA+Nx*Nx) = 1;
      }
      else
      {
        for (int j = 0; j < dim; j++)
        {
          int jA;
          if (dim == 4) //cas (P0,Q1)
          {
            jA = localToGlobalQ1(elementK,j+1,Nk);
          }
          else if (dim == 9) //cas (Q1,Q2)
          {
            jA = localToGlobalQ2(elementK,j+1,Nk);
          }
          M.coeffRef(iA,jA)+=Ak[i][j];
          M.coeffRef(iA+Nx*Nx,jA+Nx*Nx)+=Ak[i][j];
        }
      }
    }
  }
}

void insertB1B2avecCL(std::vector<std::vector<double> > B1k, std::vector<std::vector<double> > B2k, int Nk,Eigen::SparseMatrix<double> &M)
{
  if (Nk < 1)
  {
    std::cout << "Le nombre d'éléments par ligne doit être positif pour créer B1" << std::endl;
    exit(1);
  }

  int dim1 = B1k.size();
  int dim2 = B1k[0].size();

  if (B2k.size() != dim1 or B2k[0].size() != dim2)
  {
    std::cout << "Les matrices B1k et B2k n'ont pas les mêmes dimensions" << std::endl;
    exit(1);
  }

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
  else
  {
    std::cout << "La matrice Bk n'a pas une dimension correspondant à (P0,Q1) ou (Q1,Q2)" << std::endl;
    exit(1);
  }

  for (int elementK = 0; elementK < Nk*Nk; elementK++)
  {
    for(int i = 0; i < dim1; i++)
    {
      for (int j = 0; j < dim2; j++)
      {
        int iB,jB;
        if (dim1 == 4 and dim2 == 1) //cas (P0,Q1)
        {
          iB = localToGlobalQ1(elementK,i+1,Nk);
          jB = elementK;
        }
        else if (dim1 == 9 and dim2 == 4) //cas (Q1,Q2)
        {
          iB = localToGlobalQ2(elementK,i+1,Nk);
          jB = localToGlobalQ1(elementK,j+1,Nk);
        }

        if (iB >= Nx1 and iB%Nx1 != 0 and iB%Nx1 != Nx1-1 and iB < Nx1*(Nx1-1))
        {
          M.coeffRef(iB,jB+2*Nx1*Nx1) += B1k[i][j];
          M.coeffRef(iB+Nx1*Nx1,jB+2*Nx1*Nx1) += B2k[i][j];
        }

        if (jB >= Nx2 and jB%Nx2 != 0 and jB%Nx2 != Nx2-1 and jB < Nx2*(Nx2-1)) //remplissage des B^T
        {
          M.coeffRef(jB+2*Nx1*Nx1,iB) += B1k[i][j];
          M.coeffRef(jB+2*Nx1*Nx1,iB+Nx1*Nx1) += B2k[i][j];
        }
        else
        {
          M.coeffRef(jB+2*Nx1*Nx1,jB+2*Nx1*Nx1) = 1;
        }
      }
    }
  }
}

Eigen::SparseMatrix<double> createMavecCL(int choix, int Nk)
{
  int Nx1,Nx2;
  std::vector<std::vector<double> > AK, B1K, B2K;
  Eigen::SparseMatrix<double> M;

  switch (choix) {
    case 1:
    AK = createAK(getQ1PolVect());
    B1K = createB1K(getP0PolVect(),getQ1PolVect());
    B2K = createB2K(getP0PolVect(),getQ1PolVect());
    Nx1=Nk+1;
    Nx2=Nk;
    break;
    case 2:
    AK = createAK(getQ2PolVect());
    B1K = createB1K(getQ1PolVect(),getQ2PolVect());
    B2K = createB2K(getQ1PolVect(),getQ2PolVect());
    Nx1=2*Nk+1;
    Nx2=Nk+1;
    break;
    default:
    std::cout<<"Le choix doit être 1 ou 2"<<std::endl;
    exit(1);
  }

  M.resize(2*(Nx1*Nx1)+Nx2*Nx2,2*(Nx1*Nx1)+Nx2*Nx2);

  insertAavecCL(AK, Nk, M);
  insertB1B2avecCL(B1K,B2K,Nk,M);

  return M;
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
