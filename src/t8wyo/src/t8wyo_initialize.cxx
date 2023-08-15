/**
 * \file    t8wyo_initialize.c
 * \author  akirby
 *
 * \brief Initialization functions for the AMR code module.
 */

/* header files */
#include "t8wyo_initialize.h"

int t8wyo_initialize_libs(int argc,char **argv,ctx_t *ctx){
    /* ============== */
    /* initialize MPI */
    /* ============== */
    int mpi_return = MPI_Init(&argc, &argv);
    ctx->comm = MPI_COMM_WORLD;

    /* ======================= */
    /* initialize sc and p4est */
    /* ======================= */
    sc_init(ctx->comm,1,1,NULL,SC_LP_DEFAULT);
    p4est_init(NULL,SC_LP_ERROR);
    t8_init(SC_LP_PRODUCTION);

    /* ======================== */
    /* assign MPI rank and size */
    /* ======================== */
    MPI_Comm_rank(ctx->comm,&ctx->rank);
    MPI_Comm_size(ctx->comm,&ctx->nranks);

    DEBUG_MESG(fprintf(ctx->log_io,"[t8wyo] My mpi rank: %d\n",ctx->rank))
    return mpi_return;
}

void t8wyo_initialize_libs_from_comm(MPI_Comm *comm,ctx_t *ctx){
    /* ============================================ */
    /* initialize MPI information from communicator */
    /* ============================================ */
    ctx->comm = *comm;

    /* ======================= */
    /* initialize sc and p4est */
    /* ======================= */
    sc_init(ctx->comm,1,1,NULL,SC_LP_ALWAYS);
    p4est_init(NULL,SC_LP_ERROR);
    t8_init(SC_LP_PRODUCTION);

    /* ======================== */
    /* assign MPI rank and size */
    /* ======================== */
    MPI_Comm_rank(ctx->comm,&ctx->rank);
    MPI_Comm_size(ctx->comm,&ctx->nranks);

    DEBUG_MESG(DGOUT(ctx->log_io,"[t8wyo] My mpi rank: %d\n",ctx->rank))
}

t8_cmesh_t
t8wyo_create_cmesh(char mode,const char *mshfile,mcell_t *mcell,
                   int level,int dim,int use_occ_geometry,
                   MPI_Comm comm){

    /* load from gmsh or mcell file and partition */
    if (mode == MESH_MODE_GMSH || mode == MESH_MODE_MCELL) {
        t8_cmesh_t cmesh = (mode == MESH_MODE_GMSH) ? t8_cmesh_from_msh_file(mshfile,0,comm,dim,0,use_occ_geometry):
                           (mode == MESH_MODE_MCELL)? t8_cmesh_from_mcell(mcell,1,comm,dim,use_occ_geometry):
                                                      NULL;

        /* partitioning of the occ geometry is not yet available */
        if (use_occ_geometry) {
            t8_productionf("cmesh was not partitioned. Partitioning is not yet "
                           "available with the curved geometry\n");
            return cmesh;
        }

        /* partition this cmesh according to the initial refinement level */
        t8_cmesh_t cmesh_partition;
        t8_cmesh_init(&cmesh_partition);
        t8_cmesh_set_partition_uniform(cmesh_partition,level,t8_scheme_new_default_cxx());
        t8_cmesh_set_derive(cmesh_partition,cmesh);
        t8_cmesh_commit(cmesh_partition,comm);
        return cmesh_partition;
    } else
    if(mode == MESH_MODE_CART) {
        /* t8code example meshes */
        if(dim == 2) return t8_cmesh_new_disjoint_bricks(10,10, 0,0,0,0,comm); // 2D
        if(dim == 3) return t8_cmesh_new_disjoint_bricks(10,10,10,0,0,0,comm); // 3D
        return NULL;
    } else {
        printf("UNKNOWN MESH TYPE!\n");
        exit(1);
    }
}

t8_forest_t
t8wyo_build_forest(t8_cmesh_t cmesh,int level,sc_MPI_Comm comm){
    t8_scheme_cxx_t *scheme = t8_scheme_new_default_cxx();
    int ghost_flag = 1;

    /* return uniform forest with ghost faces */
    return t8_forest_new_uniform(cmesh,scheme,level,ghost_flag,comm);
}
