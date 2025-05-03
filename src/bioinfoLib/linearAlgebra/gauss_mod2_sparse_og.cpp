#include <linbox/linbox-config.h>
#include <linbox/solutions/solve.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/solve/solve-wiedemann.h>
#include <linbox/field/gf2.h>
#include <linbox/blackbox/zo-gf2.h>
#include <givaro/modular.h>
#include <iostream>

using namespace LinBox;
int main(int argc, char **argv) {
  typedef Givaro::Modular<double> Field;
  typedef SparseMatrix<Field> Mat;
  Field F2(2);
  Mat A(F2, 4, 4);

  // A = [
  //   1 0 0 0
  //   1 0 0 1
  //   0 1 1 0
  //   0 1 1 1
  // ]
  int nelem = 8;
  int row[nelem] = {0, 1, 1, 2, 2, 3, 3, 3};
  int col[nelem] = {0, 0, 3, 1, 2, 1, 2, 3};
  for (int i = 0; i < nelem; ++i) {
    A.setEntry(
      row[i],
      col[i],
      F2.one
    );
  }

  // solvable: b1 = [1, 1, 0, 0]
  // not solvable: b2 = [0, 0, 1, 0]
  DenseVector<Field> b1(F2, 4), b2(F2, 4), x1(F2, 4), x2(F2, 4), b1s(F2, 4), b2s(F2, 4);
  F2.assign(b1[0], F2.one);
  F2.assign(b1[1], F2.one);
  F2.assign(b2[2], F2.one);

  Method::SparseElimination W;
  // BW.blockingFactor = 4;
  try {
    solve(
      x1,
      A,
      b1,
      RingCategories::ModularTag(),
      W);
  } catch (const LinboxError&) {    
    std::cout << "b1 not solved" << "\n";
  }

  try {
    solve(
      x2,
      A,
      b2,
      RingCategories::ModularTag(),
      W);
  } catch (const LinboxError&) {
    std::cout << "b2 not solved" << "\n";
  }

  A.apply(b1s, x1);
  A.apply(b2s, x2);
  VectorDomain<Field> VD(F2);

  if (VD.areEqual(b1s, b1)) {
    std::cout << "b1 solved" << "\n";
  }
  if (VD.areEqual(b2s, b2)) {
    std::cout << "b2 solved" << "\n";
  }

  for (std::size_t i = 0; i < x1.size(); ++i) {
    std::cout << "b1[" << i << "] = " << b1[i] << "\n";
  }
  for (std::size_t i = 0; i < x1.size(); ++i) {
    std::cout << "x1[" << i << "] = " << x1[i] << "\n";
  }
  for (std::size_t i = 0; i < x2.size(); ++i) {
    std::cout << "b2[" << i << "] = " << b2[i] << "\n";
  }
  for (std::size_t i = 0; i < x2.size(); ++i) {
    std::cout << "x2[" << i << "] = " << x2[i] << "\n";
  }
}
