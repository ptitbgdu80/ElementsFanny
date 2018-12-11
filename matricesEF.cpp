#include "matricesEF.h"

std::vector<Polynome2D> getP0PolVect() //creation des fonctions de forme de P0
{
  Polynome2D psi = createPol({{1}});
  return {psi};
}

std::vector<Polynome2D> getQ1PolVect() //creation des fonctions de forme de Q1
{
  std::vector<std::vector<double> > v1 {{1,-1},{-1,1}}, v2 {{0,0},{1,-1}}, v3 {{0,1},{0,-1}}, v4 {{0,0},{0,1}};

  Polynome2D phi1 = createPol(v1);
  Polynome2D phi2 = createPol(v2);
  Polynome2D phi3 = createPol(v3);
  Polynome2D phi4 = createPol(v4);

  std::vector<Polynome2D> result {phi1, phi2, phi3, phi4};

  return result;
}

std::vector<Polynome2D> getQ2PolVect() //creation des fonctions de forme de Q1
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

std::vector<std::vector<double> > createAk(std::vector<Polynome2D> polVect) // creation de la matrice Ak pour l'element K (en pratique toutes les memes)
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

std::vector<std::vector<double> > createB1k(std::vector<Polynome2D> polVect1, std::vector<Polynome2D> polVect2) // creation de la matrice B1k pour l'element K (en pratique toutes les memes)
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

std::vector<std::vector<double> > createB2k(std::vector<Polynome2D> polVect1, std::vector<Polynome2D> polVect2)// creation de la matrice B2k pour l'element K (en pratique toutes les memes)
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

std::vector<double> createF1K(int choixCL, std::vector<Polynome2D> polVect)// creation du vecteur Fk pour u1  pour l'element K (en pratique toutes les memes)
{
  std::vector<double> result;
  result.resize(polVect.size());

  if(choixCL==2)
  {
    for (int i = 0; i < polVect.size(); i++)
    {
      result[i]= integraleSurUnCarreUnitaire(polVect[i]);
    }
  }
  return result;
}

std::vector<double> createF2K(int choixCL, std::vector<Polynome2D> polVect)// creation du vecteur Fk pour u2  pour l'element K (en pratique toutes les memes)
{
  std::vector<double> result;
  result.resize(polVect.size());

  if(choixCL==2)
  {
    for (int i = 0; i < polVect.size(); i++)
    {
      result[i]= -integraleSurUnCarreUnitaire(polVect[i]);
    }
  }
  return result;
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

void insertSource(std::vector<double> F1k, std::vector<double> F2k, int Nk, Eigen::VectorXd &F) // remplissage du vecteur F sans CL pour les vitesses u1 et u2
{
  if (Nk < 1)
  {
    std::cout << "Le nombre d'éléments par ligne doit être positif pour créer F" << std::endl;
    exit(1);
  }

  int dim = F1k.size();

  if (dim == 4) //cas (P0,Q1)
  {
    int Nx1 = Nk +1;
    int Nx2 = Nk;
    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim; i++)
      {
        F[localToGlobalQ1(elementK,i+1,Nk)]+=F1k[i];
        F[localToGlobalQ1(elementK,i+1,Nk)+Nx1*Nx1]+=F2k[i];
      }
    }
  }

  else if (dim == 9) //cas (Q1,Q2)
  {
    int Nx1 = 2*Nk+1;
    int Nx2 = Nk +1;
    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim; i++)
      {
        F[localToGlobalQ2(elementK,i+1,Nk)]+=F1k[i];
        F[localToGlobalQ2(elementK,i+1,Nk)+Nx1*Nx1]+=F2k[i];
      }
    }
  }
}

std::vector<double> CLvitesse (int choixCL, double x, double y) //définition des CL pour les vitesses u1 et u2
{
  std::vector<double> u;
  switch (choixCL)
  {
    case 1:
    u = {1,1};
    break;

    case 2:
    u = {x,-y};
    break;

    default:
    std::cout << "Le choix de CL ne correspond pas à un cas implémenté" << std::endl;
    exit(1);
  }
  return u;
}

Eigen::VectorXd createF(int choix, int choixCL, int Nk) // creation du second membre F
{
  double dx1, dx2;
  int Nx1;//nombre de points par ligne pour les vitesses
  int Nx2;//nombre de points par ligne pour les pressions

  std::vector<double> F1k, F2k;
  Eigen::VectorXd F;

  switch (choix)
  {
    case 1: //cas (P0,Q1)
    Nx1 = Nk +1;
    Nx2 = Nk;
    dx1 = 1./Nk;
    dx2 = 1./Nk;

    if (choixCL != 1)
    {
      F1k = createF1K(choixCL, getQ1PolVect());
      F2k = createF2K(choixCL, getQ1PolVect());
    }
    break;

    case 2: //cas (Q1,Q2)
    Nx1 = 2*Nk + 1;
    Nx2 = Nk + 1;
    dx1 = 1./(2*Nk);
    dx2 = 1./Nk;

    if (choixCL != 1)
    {
      F1k = createF1K(choixCL, getQ2PolVect());
      F2k = createF2K(choixCL, getQ2PolVect());
    }
    break;

    default:
    std::cout << "Le choix ne correspond ni au cas (P0,Q1), ni au cas (Q1,Q2)" << std::endl;
    exit(1);
  }

  F.resize(2*Nx1*Nx1+Nx2*Nx2);
  for (int i = 0; i < 2*Nx1*Nx1+Nx2*Nx2; i++)
  {
    F[i] = 0.; //initialisation du vecteur F
  }

  if (choixCL != 1)
  {
    insertSource(F1k, F2k, Nk, F);
  }

  for (int i=0; i<Nx1; i++) //vitesses de bord
  {
    //premiere ligne
    std::vector<double> Ubord = CLvitesse(choixCL,i*dx1,0);
    F[i] = Ubord[0]; //u1
    F[i+Nx1*Nx1] = Ubord[1]; //u2

    //derniere ligne
    Ubord = CLvitesse(choixCL,i*dx1,1);
    F[i+Nx1*(Nx1-1)] = Ubord[0]; //u1
    F[i+Nx1*(Nx1-1)+Nx1*Nx1] = Ubord[1]; //u2

    //première colonne
    Ubord = CLvitesse(choixCL,0,i*dx1);
    F[i*Nx1] = Ubord[0]; //u1
    F[i*Nx1+Nx1*Nx1] = Ubord[1]; //u2

    //dernière colonne
    Ubord = CLvitesse(choixCL,1,i*dx1);
    F[(i+1)*Nx1-1] = Ubord[0]; //u1
    F[(i+1)*Nx1-1+Nx1*Nx1] = Ubord[1]; //u2
  }

  return F;
}

void insertA(std::vector<std::vector<double> > Ak, int Nk, Eigen::SparseMatrix<double> &M) //remplissage de la matrice M avec les éléments de A
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

    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim; i++)
      {
        int iA;
        iA = localToGlobalQ1(elementK,i+1,Nk);

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
            jA = localToGlobalQ1(elementK,j+1,Nk);

            M.coeffRef(iA,jA)+=Ak[i][j];//chaque élément K apporte sa contribution à Aij
            M.coeffRef(iA+Nx*Nx,jA+Nx*Nx)+=Ak[i][j];//remplissage de la deuxieme matrice A
          }
        }
      }
    }
  }

  else if (dim == 9) //cas (Q1,Q2)
  {
    Nx = 2*Nk + 1;

    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim; i++)
      {
        int iA;

        iA = localToGlobalQ2(elementK,i+1,Nk);

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

            jA = localToGlobalQ2(elementK,j+1,Nk);

            M.coeffRef(iA,jA)+=Ak[i][j];//chaque élément K apporte sa contribution à Aij
            M.coeffRef(iA+Nx*Nx,jA+Nx*Nx)+=Ak[i][j];//remplissage de la deuxieme matrice A
          }
        }
      }
    }
  }

  else
  {
    std::cout << "La matrice Ak n'a pas une dimension correspondant à Q1 ou Q2" << std::endl;
    exit(1);
  }
}

void insertB1B2(std::vector<std::vector<double> > B1k, std::vector<std::vector<double> > B2k, int Nk,Eigen::SparseMatrix<double> &M) //remplissage de la matrice M avec les éléments de B1 et B2
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

    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim1; i++)
      {
        for (int j = 0; j < dim2; j++)
        {
          int iB,jB;

          iB = localToGlobalQ1(elementK,i+1,Nk);
          jB = elementK;

          if (iB >= Nx1 and iB%Nx1 != 0 and iB%Nx1 != Nx1-1 and iB < Nx1*(Nx1-1)) //conditions au bords pour les vitesses
          {
            M.coeffRef(iB,jB+2*Nx1*Nx1) += B1k[i][j];
            M.coeffRef(iB+Nx1*Nx1,jB+2*Nx1*Nx1) += B2k[i][j];
          }

          M.coeffRef(jB+2*Nx1*Nx1,iB) += B1k[i][j];
          M.coeffRef(jB+2*Nx1*Nx1,iB+Nx1*Nx1) += B2k[i][j];
        }
      }
    }
  }

  else if (dim1 == 9 and dim2 == 4) //cas (Q1,Q2)
  {
    Nx1 = 2*Nk + 1;
    Nx2 = Nk + 1;

    for (int elementK = 0; elementK < Nk*Nk; elementK++)
    {
      for(int i = 0; i < dim1; i++)
      {
        for (int j = 0; j < dim2; j++)
        {
          int iB,jB;

          iB = localToGlobalQ2(elementK,i+1,Nk);
          jB = localToGlobalQ1(elementK,j+1,Nk);

          if (iB >= Nx1 and iB%Nx1 != 0 and iB%Nx1 != Nx1-1 and iB < Nx1*(Nx1-1)) //conditions aux bords pour les vitesses
          {
            M.coeffRef(iB,jB+2*Nx1*Nx1) += B1k[i][j];
            M.coeffRef(iB+Nx1*Nx1,jB+2*Nx1*Nx1) += B2k[i][j];
          }

          M.coeffRef(jB+2*Nx1*Nx1,iB) += B1k[i][j];
          M.coeffRef(jB+2*Nx1*Nx1,iB+Nx1*Nx1) += B2k[i][j];
        }
      }
    }
  }

  else
  {
    std::cout << "La matrice Bk n'a pas une dimension correspondant à (P0,Q1) ou (Q1,Q2)" << std::endl;
    exit(1);
  }

}

Eigen::SparseMatrix<double> createM(int choix, int Nk, double epsilon) //assemblage de la matrice M
{
  int Nx1,Nx2;
  std::vector<std::vector<double> > Ak, B1k, B2k;
  Eigen::SparseMatrix<double> M;

  switch (choix) {
    case 1:
    Ak = createAk(getQ1PolVect());
    B1k = createB1k(getP0PolVect(),getQ1PolVect());
    B2k = createB2k(getP0PolVect(),getQ1PolVect());
    Nx1=Nk+1;
    Nx2=Nk;
    break;
    case 2:
    Ak = createAk(getQ2PolVect());
    B1k = createB1k(getQ1PolVect(),getQ2PolVect());
    B2k = createB2k(getQ1PolVect(),getQ2PolVect());
    Nx1=2*Nk+1;
    Nx2=Nk+1;
    break;
    default:
    std::cout<<"Le choix doit être 1 ou 2"<<std::endl;
    exit(1);
  }

  M.resize(2*(Nx1*Nx1)+Nx2*Nx2,2*(Nx1*Nx1)+Nx2*Nx2);

  insertA(Ak, Nk, M);
  insertB1B2(B1k,B2k,Nk,M);
  insertEpsId(choix, epsilon, Nk, M);

  return M;
}

void insertEpsId(int choix, double epsilon, int Nk,Eigen::SparseMatrix<double> &M)//création de la matrice Epsilon Identité
{
  if (Nk < 1)
  {
    std::cout << "Le nombre d'éléments par ligne doit être positif pour créer Epsilon Id" << std::endl;
    exit(1);
  }

  int Nx1, Nx2;

  if (choix == 1) //cas (P0,Q1)
  {
    Nx1 = Nk +1;
    Nx2 = Nk;
  }
  else if (choix == 2) //cas (Q1,Q2)
  {
    Nx1 = 2*Nk + 1;
    Nx2 = Nk + 1;
  }
  else
  {
    std::cout << "La matrice Epsilon Id n'a pas une dimension correspondant à (P0,Q1) ou (Q1,Q2)" << std::endl;
    exit(1);
  }

  for (int i=0; i<Nx2*Nx2; i++)
  {
    M.coeffRef(2*Nx1*Nx1+i,2*Nx1*Nx1+i)=epsilon;
  }
}

void saveSol(std::string fichier, int choix, int Nk, Eigen::VectorXd U) //creation d'un fichier .txt pour sauver les solutions
{
  int Nx1, Nx2;
  double dx1, dx2, decalage;
  switch (choix)
  {
    case 1:
    Nx1 = Nk+1;
    Nx2 = Nk;
    dx1 = 1./Nk;
    dx2 = 1./Nk;
    decalage = dx2/2.;
    break;

    case 2:
    Nx1 = 2*Nk+1;
    Nx2 = Nk+1;
    dx1 = 1./(2*Nk);
    dx2 = 1./Nk;
    decalage = 0;
    break;

    default:
    std::cout << "Le choix ne correspond ni au cas (P0,Q1), ni au cas (Q1,Q2)" << std::endl;
    exit(1);
  }

  std::ofstream mon_fluxP;
  mon_fluxP.open(fichier + "_p.txt", std::ios::out);
  mon_fluxP << "# champ de pressions sur un maillage carré" << std::endl;

  for (int i = 0; i < Nx2; i++)
  {
    for (int j = 0; j < Nx2; j++)
    {
      mon_fluxP << j*dx2+decalage << " " << i*dx2 + decalage << " " << U[2*Nx1*Nx1 + j + i*Nx2] << std::endl;
    }
  }

  mon_fluxP.close();

  std::ofstream mon_fluxU;
  mon_fluxU.open(fichier + "_u.txt", std::ios::out);
  mon_fluxU << "# champ de vitesses sur un maillage carré" << std::endl;

  std::vector<double> norme;
  norme.resize(Nx1*Nx1);

  double normeMax = 0;

  for (int i = 0; i < Nx1; i++)
  {
    for (int j = 0; j < Nx1; j++)
    {
      double u1 = U[j + i*Nx1];
      double u2 = U[Nx1*Nx1 + j + i*Nx1];
      norme[j + i*Nx1] = sqrt(u1*u1 + u2*u2);
      if (norme[j+i*Nx1] > normeMax)
      {
        normeMax = norme[j + i*Nx1];
      }
    }
  }

  for (int i = 0; i < Nx1; i++)
  {
    for (int j = 0; j < Nx1; j++)
    {
      double u1 = U[j + i*Nx1];
      double u2 = U[Nx1*Nx1 + j + i*Nx1];
      double coeff = dx1/normeMax;
      mon_fluxU << j*dx1 << " " << i*dx1 << " " << u1*coeff << " " << u2*coeff << " " << norme[j + i*Nx1] << std::endl;
    }
  }

  mon_fluxU.close();
}
