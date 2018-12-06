#include "matricesEF.h"

int main()
{
  Eigen::SparseMatrix<double> M = createM(2,2);

  std::cout <<M<<std::endl;

}
