/** Top-level library code-interface.
 * \file    t8wyo_interface.cxx
 * \author  akirby
 *
 * \brief   Interface to t8wyo.
 */

/* header files */
#include "t8wyo_interface.h"

void t8wyo_interface_init_(MPI_Comm comm){
    t8wyo_initialize_libs_from_comm(&comm,&t8wyo.ctx);
}

void t8wyo_interface_fortran_init_(int *fcomm){
    MPI_Comm comm = MPI_Comm_f2c(*fcomm);
    t8wyo_initialize_libs_from_comm(&comm,&t8wyo.ctx);
}

void t8wyo_build_cmesh_forest(){
    Real cmesh_time,forest_time;

    /* build coarse mesh */
    TIMER(cmesh_time,
        cmesh = t8wyo_create_cmesh(MESH_MODE_MCELL,mshfile,
                                   level_cmesh,dim,use_occ_geometry,
                                   t8wyo.ctx.comm);
    );

    /* build forest */
    TIMER(forest_time,
        forest = t8wyo_build_forest(cmesh,level_forest,t8wyo.ctx.comm);
    );
}