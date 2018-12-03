#include "polynomes.h"

int main()
{
  // std::vector<std::vector<double> > AK = createAK(getQ2PolVect());
  //
  // std::cout << "AK = 1/90*" << std::endl;
  // for (int i = 0; i < AK.size(); i++)
  // {
  //   for (int j = 0; j < AK[0].size(); j++)
  //   {
  //     std::cout << floor(900*AK[i][j]+0.1)/10 << " ";
  //   }
  //   std::cout << std::endl;
  // }

  // std::vector<std::vector<double> > B1K = createB1K(getQ1PolVect(),getQ2PolVect());
  //
  // std::cout << "B1K = 1/360*" << std::endl;
  // for (int i = 0; i < B1K.size(); i++)
  // {
  //   for (int j = 0; j < B1K[0].size(); j++)
  //   {
  //     std::cout << floor(36000*B1K[i][j] + 0.1)/100 << "   ";
  //   }
  //   std::cout << std::endl;
  // }

  std::vector<std::vector<double> > B2K = createB2K(getQ1PolVect(),getQ2PolVect());

  std::cout << "B2K = 1/360*" << std::endl;
  for (int i = 0; i < B2K.size(); i++)
  {
    for (int j = 0; j < B2K[0].size(); j++)
    {
      std::cout << floor(36000*B2K[i][j] + 0.1)/100 << "   ";
    }
    std::cout << std::endl;
  }
}
