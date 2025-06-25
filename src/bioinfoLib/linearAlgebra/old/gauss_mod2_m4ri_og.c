#include "m4ri/m4ri.h"

int main() {
  mzd_t *mat = mzd_init(2, 2);
  mzd_write_bit(mat, 0, 0, 0);
  mzd_write_bit(mat, 1, 0, 1);
  mzd_write_bit(mat, 0, 1, 1);
  mzd_write_bit(mat, 1, 1, 1);
  mzd_print_row(mat, 0);
  mzd_print_row(mat, 1);

  rci_t rmat = mzd_echelonize_m4ri(mat, 1, 0);
  printf("Value: %d\n", rmat);

  mzd_print_row(mat, 0);
  mzd_print_row(mat, 1);
  mzd_free(mat);
}
