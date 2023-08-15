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
#include "t8wyo_adapt.h"

/* system header files */
#include <mpi.h>

T8WYO_EXTERN_C_BEGIN();

/* functions */
void t8wyo_interface_init_(MPI_Comm comm);

void t8wyo_interface_fortran_init_(int *fcomm);

void t8wyo_build_cmesh_mcell_(int *level_cmesh,
                              int *ntetra,int *npyr,int *nprizm,int *nhex,
                              int *ntetra_ng,int *npyr_ng,int *nprizm_ng,int *nhex_ng,
                              int *nnode,int *ngpt,
                              int *nbface3,int *nbface4,
                              int *nface3,int *nface4,
                              int *ndc4,int *ndc5,int *ndc6,int *ndc8,
                              int *nbf3,int *nbf4,
                              int *ifpat3,int *ifpat4,
                              Real *xgeom,
                              int *ndf3,int *ndf4);

void t8wyo_build_forest_(int *level_forest,
                         int *ntetra,int *npyr,int *nprizm,int *nhex,
                         int *ntetra_ng,int *npyr_ng,int *nprizm_ng,int *nhex_ng);

void t8wyo_build_lists_(int *ncell_real,int *ncell,int *nface,
                        int **face2cellptr,int **ifacetypeptr,
                        int **cellinfoptr,
                        Real **cellvolptr,Real **facenormptr);

void t8wyo_exchange_ghost_data_(void *data,
                                size_t *bytes_per_element,
                                int *barrier_flag,
                                int *mess_flag);

void t8wyo_allocate_solution_(int *nvar,int *ncell,Real **wvalues_new);

void t8wyo_adapt_(tag_callback_t *tag_function,
                  int *nvar,int *ncell,int *ncell_real,
                  Real *wvalues,Real **wvalues_new);

void t8wyo_write_vtk_(int *vtk_counter,Real *wvalues_in);

T8WYO_EXTERN_C_END();

#endif /* T8WYO_INTERFACE_H */