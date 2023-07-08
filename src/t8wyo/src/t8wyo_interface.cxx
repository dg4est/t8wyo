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
                              int *ndf3,int *ndf4){
    mcell_t mcell;
    mcell.ntetra = *ntetra;
    mcell.npyr = *npyr;
    mcell.nprizm = *nprizm;
    mcell.nhex = *nhex;
    mcell.ntetra_ng = *ntetra_ng;
    mcell.npyr_ng = *npyr_ng;
    mcell.nprizm_ng = *nprizm_ng;
    mcell.nhex_ng = *nhex_ng;
    mcell.nnode = *nnode;
    mcell.ngpt = *ngpt;
    mcell.nbface3 = *nbface3;
    mcell.nbface4 = *nbface4;
    mcell.nface3 = *nface3;
    mcell.nface4 = *nface4;
    mcell.ndc4 = ndc4;
    mcell.ndc5 = ndc5;
    mcell.ndc6 = ndc6;
    mcell.ndc8 = ndc8;
    mcell.nbf3 = nbf3;
    mcell.nbf4 = nbf4;
    mcell.ifpat3 = ifpat3;
    mcell.ifpat4 = ifpat4;
    mcell.xgeom = xgeom;
    mcell.ndf3 = ndf3;
    mcell.ndf4 = ndf4;

    /* build coarse mesh */
    int use_occ_geometry = 0;   /* open cascade geometry flag */
    int dim = DIM;              /* grid dimension */

    /* build coarse mesh */
    Real cmesh_time;
    TIMER(cmesh_time,
        cmesh = t8wyo_create_cmesh(MESH_MODE_MCELL,NULL,&mcell,
                                   *level_cmesh,dim,use_occ_geometry,
                                   t8wyo.ctx.comm);
    );

    if(t8wyo.ctx.rank==0) printf("[t8wyo] COARSE MESH CONSTRUCTION: %f (sec)\n",cmesh_time);
}

void t8wyo_build_forest_(int *level_forest){
    /* build forest */
    Real forest_time;
    TIMER(forest_time,
        forest = t8wyo_build_forest(cmesh,*level_forest,t8wyo.ctx.comm);
    );

    /* send back new grid info */
//    ntetra;
//    npyr;
//    nprizm;
//    nhex;
//    ntetra_ng;
//    npyr_ng;
//    nprizm_ng;
//    nhex_ng;
//    nnode;
//    ngpt;
    if(t8wyo.ctx.rank==0) printf("[t8wyo] FOREST CONSTRUCTION: %f (sec)\n",forest_time);
}