/** Top-level main function.
 * \file    main.cxx
 * \author  akirby
 *
 * \brief   Main function for t8wyo.
 */

/* header files */
#include "main.h"

/** Main
 * @param [in] argc     number of command line arguments
 * @param [in] argv     command line arguments
 * @return              program status
 */

#ifdef _2D_
# define DIM 2
#else
# define DIM 3
#endif

int main(int argc, char **argv){
    t8_cmesh_t cmesh;
    t8wyo_t t8wyo;

    /* initialize */
    initialize_libs(argc,argv,&t8wyo.ctx);

    /* build coarse mesh */
    char *mshfile = NULL;       /* gmsh (ascii) file to read */
    int use_occ_geometry = 0;   /* open cascade geometry flag */
    int level = 3;              /* refinement level */
    int dim = DIM;              /* grid dimension */
    char outfile[BUFF_SIZE];

    (DIM==2) ? strcpy(outfile,"cmesh_demo_2D"):
               strcpy(outfile,"cmesh_demo_3D");

    cmesh = t8wyo_create_cmesh(t8wyo.ctx.comm,mshfile,level,dim,use_occ_geometry);
    t8_cmesh_vtk_write_file(cmesh,outfile,1.0);

    /* shut down */
    t8_cmesh_destroy(&cmesh);
    sc_finalize();
    return mpi_finalize();
}