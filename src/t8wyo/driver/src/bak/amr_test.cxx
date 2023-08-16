/**
 * \file    amr_test.cxx
 * \ingroup amr_group
 * \author  akirby, mbrazell
 *
 * \brief Initialization functions for the AMR code module.
 */

/* header files */
#include "amr_test.h"

int test(int a){
  printf("Value of A: %d\n",a);

  wyo::memory<int> my_int_array(100);

  printf("Size of array: %d\n",my_int_array.length());

  my_int_array.resize(200);
  return 0;
}
