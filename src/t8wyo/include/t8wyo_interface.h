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
#include <iostream>

/* ======================= */
/* T8WYO LIBRARY FUNCTIONS */
/* ======================= */
T8WYO_EXTERN_C_BEGIN();

/** Initialize t8wyo library from C/C++ MPI communicator.
 *
 * @param [in] comm     C/C++ MPI communicator.
 */
void t8wyo_interface_init_(MPI_Comm comm);

/** Initialize t8wyo library from Fortran MPI communicator.
 *
 * @param [in] fcomm    Fortran MPI communicator
 */
void t8wyo_interface_fortran_init_(int *fcomm);

/** Construct the coarse mesh (cmesh) from mcell data structures.
 *
 * @param [in] level_cmesh  level to refine the initial coarse mesh
 * @param [in] ntetra       number of tets in mcell mesh
 * @param [in] npyr         number of pyramids in mcell mesh
 * @param [in] nprizm       number of prisms in mcell mesh
 * @param [in] nhex         number of hexs in mcell mesh
 * @param [in] ntetra_ng    number of ghost tets in mcell mesh
 * @param [in] npyr_ng      number of ghost pyramids in mcell mesh
 * @param [in] nprizm_ng    number of ghost prisms in mcell mesh
 * @param [in] nhex_ng      number of ghost hexs in mcell mesh
 * @param [in] nnode        number of nodes in cmell mesh
 * @param [in] ngpt         number of ghost nodes in cmell mesh
 * @param [in] nbface3      
 * @param [in] nbface4      
 * @param [in] nface3       
 * @param [in] nface4       
 * @param [in] ndc4         
 * @param [in] ndc5         
 * @param [in] ndc6         
 * @param [in] ndc8         
 * @param [in] nbf3         
 * @param [in] nbf4         
 * @param [in] ifpat3       
 * @param [in] ifpat4       
 * @param [in] xgeom        
 * @param [in] ndf3         
 * @param [in] ndf4         
 */
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

/** Construct the forest.
 *
 * @param [out] ntetra      number of tets in forest
 * @param [out] npyr        number of pyramids in forest
 * @param [out] nprizm      number of prisms in forest
 * @param [out] nhex        number of hexs in forest
 * @param [out] ntetra_ng   number of ghost tets in forest (set to 0)
 * @param [out] npyr_ng     number of ghost pyramids in forest (set to 0)
 * @param [out] nprizm_ng   number of ghost prisms in forest (set to 0)
 * @param [out] nhex_ng     number of ghost hexs in forest (set to 0)
 */
void t8wyo_build_forest_(int *ntetra,
                         int *npyr,
                         int *nprizm,
                         int *nhex,
                         int *ntetra_ng,
                         int *npyr_ng,
                         int *nprizm_ng,
                         int *nhex_ng);

/** Construct the forest connectivity information.
 *
 * @param [out] ncell_real      number of real cells on rank
 * @param [out] ncell           number of total cells on rank (real+ghost)
 * @param [out] nface           number of mesh faces
 * @param [out] nnodes          number of mesh nodes (unique)
 * @param [out] ntet            number of real tetrahedra
 * @param [out] npyr            number of real pyramids
 * @param [out] nprism          number of real prisms
 * @param [out] nhex            number of real hexahedra
 * @param [out] face2cellptr    array pointer to face2cell data structure
 * @param [out] ifacetypeptr    array pointer to ifacetype data structure
 * @param [out] cellinfoptr     array pointer to cellinfo data structure
 * @param [out] ndc4ptr         array pointer to ndc4 data structure: tet connectivity, size=ntet
 * @param [out] ndc5ptr         array pointer to ndc5 data structure: pyramid connectivity, size=npyr
 * @param [out] ndc6ptr         array pointer to ndc6 data structure: prism connectivity, size=nprism
 * @param [out] ndc8ptr         array pointer to ndc8 data structure: hex connectivity, size=nhex
 * @param [out] xgeomptr        array pointer to xgeom data structure: node vertices, size=nnodes
 * @param [out] cellvolptr      array pointer to cellvol data structure: cell volumes, size=ncell_real
 * @param [out] facenormptr     array pointer to facenorm data structure: face normal vectors, size=nface
 */
void t8wyo_build_lists_(int *ncell_real,int *ncell,
                        int *nface,int *nnodes,
                        int *ntet,int *npyr,int *nprism,int *nhex,
                        int **face2cellptr,
                        int **ifacetypeptr,
                        int **cellinfoptr,
                        int **ndc4ptr,
                        int **ndc5ptr,
                        int **ndc6ptr,
                        int **ndc8ptr,
                        Real **xgeomptr,
                        Real **cellvolptr,
                        Real **facenormptr);

/** Ghost exchange abstraction function.
 *
 * @param [inout] data              pointer to solution array for mpi exchange (must be length=ncell+nghost)
 * @param [in] bytes_per_element    number of bytes per entry in the solution array
 * @param [in] barrier_flag         flag to call mpi_barrier
 * @param [in] mess_flag            flag to display mpi communication time
 */
void t8wyo_exchange_ghost_data_(void *data,
                                size_t *bytes_per_element,
                                int *barrier_flag,
                                int *mess_flag);

/** Allocate solution vector as it can be recursively adapted
 *
 * @param [in] nvar         number of variables per cell
 * @param [in] ncell        number of total cells
 * @param [out] wvalues_new solution vector array (allocated)
 */
void t8wyo_allocate_solution_(int *nvar,
                              int *ncell,
                              Real **wvalues_new);

/** Adaption forest based on tag_function callback function.
 *
 * @param [in] tag_function tagging callback function supplied by user (see tag_callback_t for function signature)
 * @param [in] nvar         number of variables per cell
 * @param [out] ncell       number of total cells (real+ghost) after adaption
 * @param [out] ncell_real  number of real cells after adaption
 * @param [in] wvalues      solution vector to be adapted
 * @param [out] wvalues_new new solution vector
 */
void t8wyo_adapt_(tag_callback_t *tag_function,
                  int *nvar,
                  int *ncell,
                  int *ncell_real,
                  Real *wvalues,
                  Real **wvalues_new);

/** Output VTK visualization files (ascii format only).
 *
 * @param [in] vtk_counter  index of output files
 * @param [in] wvalues_in   solution vector to be visualized
 */
void t8wyo_write_vtk_(int *vtk_counter,
                      Real *wvalues_in);

T8WYO_EXTERN_C_END();

#endif /* T8WYO_INTERFACE_H */
