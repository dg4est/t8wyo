/**
 * \file   amr_initialize.h
 * \author akirby, mbrazell
 */

#ifndef AMR_INITIALIZE_H
#define AMR_INITIALIZE_H

/* header files */
#include "t8wyo_solver.hxx"

/* 3PL header files */
#include <sc.h>
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

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
 * @param comm
 * @param cube_type
 * @param mshfile
 * @param level
 * @param dim
 * @param use_occ_geometry
 * @return
 */
t8_cmesh_t t8wyo_create_cmesh(MPI_Comm comm,
                              const char *mshfile, int level, int dim,
                              int use_occ_geometry);

#ifdef __cplusplus
}
#endif
#endif /* AMR_INITIALIZE_H */