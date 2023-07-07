/**
 * File:   t8wyo_interface.h
 * Author: akirby
 *
 * Created on July 6, 2023, 2:35 PM
 */

#ifndef T8WYO_INTERFACE_H
#define T8WYO_INTERFACE_H

/* header files */
#include "t8wyo_globals.h"
#include "t8wyo_initialize.h"

/* system header files */
#include <mpi.h>

T8WYO_EXTERN_C_BEGIN();
void t8wyo_interface_init_(MPI_Comm comm);
void t8wyo_interface_fortran_init_(int *fcomm);
T8WYO_EXTERN_C_END();

#endif /* T8WYO_INTERFACE_H */