/**
 * \file   amr_initialize.h
 * \author akirby, mbrazell
 */

#ifndef AMR_INITIALIZE_H
#define AMR_INITIALIZE_H

/* header files */
#include "t8wyo_solver.hxx"
#include "t8wyo_cmesh_from_mcell.h"

/* 3PL header files */
#include <sc.h>
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_general.h>

#ifdef __cplusplus
extern "C" {
#endif

/** MPI and p4est initialization function wrapper
 *
 * @param [in]    argc      number of command line arguments
 * @param [in]    argv      command line arguments
 * @param [inout] ctx       context data
 */
int initialize_libs(int argc,char **argv,ctx_t *ctx);

/** MPI and p4est initialization function wrapper given a mpi communicator
 *
 * @param [in]    comm              mpi communicator
 * @param [inout] ctx               context data
 * @param [in]    log_threshold     logging threshold flag
 */
void initialize_libs_from_comm(MPI_Comm *comm,ctx_t *ctx,int log_threshold);

/** MPI finalization function wrapper
 *
 * @return mpi finalize success
 */
int mpi_finalize();

/** t8code coarse mesh construction.
 *
 * @param [in] mode                 mesh mode: MCELL,GMESH,CART
 * @param [in] mshfile              gmsh file name
 * @param [in] level                level to refine coarse mesh
 * @param [in] dim                  dimension of mesh
 * @param [in] use_occ_geometry     set geometry to use open cascade cad
 * @param [in] comm                 mpi communicator
 * @return     cmesh
 */
t8_cmesh_t
t8wyo_create_cmesh(char mode,const char *mshfile,
                   int level,int dim,int use_occ_geometry,
                   MPI_Comm comm);

/** t8code uniform forest construction.
 *
 * @param [in] cmesh    coarse mesh to build uniform forest
 * @param [in] level    level to uniformly refine forest
 * @param [in] comm     mpi communicator
 * @return     forest
 */
t8_forest_t
t8wyo_build_forest(t8_cmesh_t cmesh,int level,sc_MPI_Comm comm);


void t8wyo_build_lists(t8_forest_t forest,external_t *ext);

#ifdef __cplusplus
}
#endif
#endif /* AMR_INITIALIZE_H */