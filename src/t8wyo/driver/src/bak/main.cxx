/** Top-level main function.
 * \file    main.cxx
 * \author  akirby
 *
 * \brief   Main function for t8wyo.
 */

/* header files */
#include "main.h"

#ifdef _2D_
#  define DIM 2
#else
#  define DIM 3
#endif

/** Main
 * @param [in] argc     number of command line arguments
 * @param [in] argv     command line arguments
 * @return              program status
 */
int main(int argc, char **argv){
    t8_forest_t forest;
    t8_cmesh_t cmesh;
    t8wyo_t t8wyo;

    /* timers */
    Real cmesh_time;
    Real forest_time;
    Real build_time;

    /* initialize */
    initialize_libs(argc,argv,&t8wyo.ctx);
    cell3d_initialize();

    /* build coarse mesh */
    char *mshfile = NULL;       /* gmsh (ascii) file to read */
    int use_occ_geometry = 0;   /* open cascade geometry flag */
    int level_cmesh = 0;        /* cmesh refinement level */
    int level_forest = 0;       /* forest refinement level */
    int dim = DIM;              /* grid dimension */
    char outfile[BUFF_SIZE];

    (DIM==2) ? strcpy(outfile,"cmesh_demo_2D"):
               strcpy(outfile,"cmesh_demo_3D");

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

    /* build cell/face data structures */
    TIMER(build_time,
        t8wyo_build_lists(forest,t8wyo.external.ptr());
    );

    printf("[%4d] Timers: cmesh=%f, forest=%f, lists=%f\n",
            t8wyo.ctx.rank,cmesh_time,forest_time,build_time);

    /* visualize cmesh */
    t8_cmesh_vtk_write_file(cmesh,outfile,1.0);

    /* shut down */
    t8_forest_unref(&forest);
    t8_cmesh_destroy(&cmesh);
    sc_finalize();
    return mpi_finalize();
}