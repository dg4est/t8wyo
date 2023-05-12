/**
 * \file    amr_initialize.c
 * \ingroup amr_group
 * \author  akirby, mbrazell
 *
 * \brief Initialization functions for the AMR code module.
 */

/* header files */
#include "amr_initialize.h"

int initialize_libs(int argc,char **argv,ctx_t *ctx){
    /* ============== */
    /* initialize MPI */
    /* ============== */
    int mpi_return = MPI_Init(&argc, &argv);
    ctx->comm = MPI_COMM_WORLD;

    /* ======================= */
    /* initialize sc and p4est */
    /* ======================= */
    sc_init(ctx->comm,1,1,NULL,SC_LP_DEFAULT);
    p4est_init(NULL,SC_LP_DEFAULT);
    t8_init(SC_LP_DEFAULT);

    /* ======================== */
    /* assign MPI rank and size */
    /* ======================== */
    MPI_Comm_rank(ctx->comm,&ctx->rank);
    MPI_Comm_size(ctx->comm,&ctx->nranks);

    DEBUG_MESG(fprintf(ctx->log_io,"[ amr ] My mpi rank: %d\n",ctx->rank))
    return mpi_return;
}

void initialize_libs_from_comm(MPI_Comm *comm,ctx_t *ctx,int log_threshold){
    (void) log_threshold;
    /* ============================================ */
    /* initialize MPI information from communicator */
    /* ============================================ */
    ctx->comm = *comm;

    /* ======================= */
    /* initialize sc and p4est */
    /* ======================= */
    sc_init(ctx->comm,1,1,NULL,SC_LP_ALWAYS);
    t8_init(SC_LP_DEFAULT);

    /* ======================== */
    /* assign MPI rank and size */
    /* ======================== */
    MPI_Comm_rank(ctx->comm,&ctx->rank);
    MPI_Comm_size(ctx->comm,&ctx->nranks);

    DEBUG_MESG(DGOUT(ctx->log_io,"[ amr ] My mpi rank: %d\n",ctx->rank))
}

int mpi_finalize(){
    return MPI_Finalize();
}

t8_cmesh_t t8wyo_create_cmesh(MPI_Comm comm,
                              const char *mshfile,int level,int dim,
                              int use_occ_geometry){
    if (mshfile != NULL) {
        /* Load from .msh file and partition */
        t8_cmesh_t cmesh, cmesh_partition;
        T8_ASSERT(mshfile != NULL);

        cmesh = t8_cmesh_from_msh_file(mshfile,0,comm,dim,0,use_occ_geometry);

        /* The partitioning of the occ geometry is not yet available */
        if (use_occ_geometry) {
          t8_productionf ("cmesh was not partitioned. Partitioning is not yet "
                          "available with the curved geometry\n");
          return cmesh;
        }
        /* partition this cmesh according to the initial refinement level */
        t8_cmesh_init(&cmesh_partition);
        t8_cmesh_set_partition_uniform(cmesh_partition,level,t8_scheme_new_default_cxx());
        t8_cmesh_set_derive(cmesh_partition,cmesh);
        t8_cmesh_commit(cmesh_partition, comm);
        return cmesh_partition;
    } else {
        /* t8code example meshes */
        if(dim == 2) return t8_cmesh_new_disjoint_bricks(10,10, 0,0,0,0,comm); // 2D
        if(dim == 3) return t8_cmesh_new_disjoint_bricks(10,10,10,0,0,0,comm); // 3D
  }
}