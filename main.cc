#include "matricesEF.h"

int main()
{
  Eigen::SparseMatrix<double> M = createMavecCL(2,1);

  std::cout << M << std::endl;

}
