/** Top-level library code-interface.
 * \file    t8wyo_interface.cxx
 * \author  akirby
 *
 * \brief   Interface to t8wyo.
 */

/* header files */
#include "t8wyo_interface.h"
#include "memory_utilities.h"

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
        t8wyo_cmesh = t8wyo_create_cmesh(MESH_MODE_MCELL,NULL,&mcell,
                                        *level_cmesh,dim,use_occ_geometry,
                                         t8wyo.ctx.comm);
    );

    if(t8wyo.ctx.rank==0) printf("[t8wyo] COARSE MESH CONSTRUCTION: %f (sec)\n",cmesh_time);
}

void t8wyo_build_forest_(int *ntetra,int *npyr,int *nprizm,int *nhex,
                         int *ntetra_ng,int *npyr_ng,int *nprizm_ng,int *nhex_ng){
    /* build forest */
    Real forest_time;
    TIMER(forest_time,
        t8wyo_forest = t8wyo_build_forest(t8wyo_cmesh,0,t8wyo.ctx.comm);
    );

    /* count number of elements of each type*/
    int elem_counts[T8_ECLASS_COUNT] = {0};
    t8_locidx_t num_local_trees = t8_forest_get_num_local_trees(t8wyo_forest);
    for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
        t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements(t8wyo_forest,itree);
        t8_eclass_t tree_class = t8_forest_get_tree_class(t8wyo_forest,itree);
        t8_eclass_scheme_c *eclass_scheme = t8_forest_get_eclass_scheme(t8wyo_forest,tree_class);

        for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
            const t8_element_t *element = t8_forest_get_element_in_tree(t8wyo_forest,itree,ielement);
            const t8_element_shape_t element_shape = eclass_scheme->t8_element_shape(element);

            elem_counts[element_shape]++;
        }
    }

    /* set local element counts */
    *ntetra = elem_counts[T8_ECLASS_TET];
    *npyr   = elem_counts[T8_ECLASS_PYRAMID];
    *nprizm = elem_counts[T8_ECLASS_PRISM];
    *nhex   = elem_counts[T8_ECLASS_HEX];

    /* set local element ghost type counts */
    // TODO: get ghost element type counts
    *ntetra_ng = 0;
    *npyr_ng   = 0;
    *nprizm_ng = 0;
    *nhex_ng   = 0;

    //printf("T8WYO ELEM COUNTS: %d %d %d %d\n",*ntetra,*npyr,*nprizm,*nhex);
    if(t8wyo.ctx.rank==0) printf("[t8wyo] FOREST CONSTRUCTION: %f (sec)\n",forest_time);
}

void t8wyo_build_lists_(int *ncell_real,int *ncell,
                        int *nface,int *nnodes,int *nmortar,
                        int *ntet,int *npyr,int *nprism,int *nhex,
                        int **face2cellptr,
                        int **ifacetypeptr,
                        int **mortarptr,
                        int **cellinfoptr,
                        int **ndc4ptr,
                        int **ndc5ptr,
                        int **ndc6ptr,
                        int **ndc8ptr,
                        Real **xgeomptr,
                        Real **cellvolptr,
                        Real **facenormptr){
    /* construct face2cell,facetype,elem_info,elem_vol data structures */
    if(t8wyo.ctx.rank==0) std::cout << "[t8wyo] BUILDING CONNECTIVITY LISTS...\n"; std::cout.flush();
    Real lists_time;
    TIMER(lists_time,
        t8wyo_build_lists_ext(t8wyo_cmesh,t8wyo_forest,
                              face2cell,ifacetype,mortar_info,
                              elem_info,ndc4,ndc5,ndc6,ndc8,
                              forest_xgeom,elem_vol,face_norm);
    );

    t8_locidx_t num_elements = t8_forest_get_local_num_elements(t8wyo_forest);
    t8_locidx_t num_ghosts = t8_forest_get_num_ghosts(t8wyo_forest);

    /* set counts */
    *ncell_real = (int) num_elements;
    *ncell = (int) (num_elements+num_ghosts);
    *nface = face2cell.length()/2;     // storing two values per face
    *nnodes = forest_xgeom.length()/3; // storing three coordinates per node
    *nmortar = mortar_info.length()/10;// storing ten integers per mortar

    *ntet   = ndc4.length()/4;
    *npyr   = ndc5.length()/5;
    *nprism = ndc6.length()/6;
    *nhex   = ndc8.length()/8;

    /* fill Fortran data */
    *face2cellptr = face2cell.ptr();
    *ifacetypeptr = ifacetype.ptr();
    *mortarptr = mortar_info.ptr();
    *cellinfoptr = elem_info.ptr();
    *ndc4ptr = ndc4.ptr();
    *ndc5ptr = ndc5.ptr();
    *ndc6ptr = ndc6.ptr();
    *ndc8ptr = ndc8.ptr();
    *xgeomptr = forest_xgeom.ptr();
    *cellvolptr = elem_vol.ptr();
    *facenormptr = face_norm.ptr();

    if(t8wyo.ctx.rank==0) {std::cout << "done: " << lists_time << " (sec)" << std::endl;}
}

void t8wyo_allocate_solution_(int *nvar,int *ncell,Real **wvalues_new){
    /* allocate solution */
    t8wyo_wvalues.malloc((*nvar)*(*ncell));

    /* fill Fortran data */
    *wvalues_new = t8wyo_wvalues.ptr();
}

void t8wyo_exchange_ghost_data_(void *data,size_t *bytes_per_element,int *barrier_flag,int *mess_flag){
    t8_locidx_t num_elements = t8_forest_get_local_num_elements(t8wyo_forest);
    t8_locidx_t num_ghosts = t8_forest_get_num_ghosts(t8wyo_forest);
    Real exc_time;

    TIMER(exc_time,
        /* t8_forest_ghost_exchange_data expects an sc_array:
         *    length = num_local_elements + num_ghosts
         *    Wrap our data array to an sc_array.
         */
        sc_array *sc_array_wrapper = sc_array_new_data(data,
                                                      *bytes_per_element,
                                                       num_elements + num_ghosts);

        /* exchange data: entries with indices > num_local_elements get overwritten */
        t8_forest_ghost_exchange_data(t8wyo_forest,sc_array_wrapper);

        /* destroy the wrapper array: this will not free the data memory since we used sc_array_new_data. */
        sc_array_destroy(sc_array_wrapper);
    );
    if(*barrier_flag == 1) MPI_Barrier(t8wyo.ctx.comm);
    if(*mess_flag == 1 && t8wyo.ctx.rank==0){
        printf("[t8wyo] Exchange Time: %f (sec) w/ %zu bytes per elem\n",
               exc_time,*bytes_per_element);
    }
}

void t8wyo_adapt_(tag_callback_t *tag_function,
                  int *nvar,int *ncell,int *ncell_real,
                  Real *wvalues,Real **wvalues_new){
    /* adapt the grid*/
    t8wyo_forest = t8wyo_adapt_ext(t8wyo_forest,
                                   tag_function,nvar,ncell,ncell_real,
                                   wvalues,t8wyo_wvalues_new);

    /* swap pointers and clean the old element data */
    t8wyo_wvalues.swap(t8wyo_wvalues_new);
    t8wyo_wvalues_new.free();

    /* set external solution pointer */
    *wvalues_new = t8wyo_wvalues.ptr();
    if(t8wyo.ctx.rank==0) memory_usage(0,1);
}

void t8wyo_write_vtk_(int *vtk_counter,Real *wvalues_in){
    t8_locidx_t num_local_elements,ielem;
    t8_vtk_data_field_t vtk_data[5];
    double *wvalues[5];
    char fileprefix[BUFSIZ];

    if(t8wyo.ctx.rank==0) printf("[t8wyo] Outputting VTK...\n");

    Real vtk_time;
    TIMER(vtk_time,
        /* Allocate num_local_elements doubles to store flow */
        num_local_elements = t8_forest_get_local_num_elements(t8wyo_forest);

        /* density */
        wvalues[0] = T8_ALLOC_ZERO (Real,num_local_elements);
        /* u */
        wvalues[1] = T8_ALLOC_ZERO (Real,num_local_elements);
        /* v */
        wvalues[2] = T8_ALLOC_ZERO (Real,num_local_elements);
        /* w */
        wvalues[3] = T8_ALLOC_ZERO (Real,num_local_elements);
        /* E */
        wvalues[4] = T8_ALLOC_ZERO (Real,num_local_elements);

        /* fill flow variables */
        for (ielem = 0; ielem < num_local_elements; ielem++) {
            const Real *W = &wvalues_in[5*ielem];
            const Real iDensity = 1.0/W[0];

            wvalues[0][ielem] = W[0];
            wvalues[1][ielem] = W[1]*iDensity;
            wvalues[2][ielem] = W[2]*iDensity;
            wvalues[3][ielem] = W[3]*iDensity;
            wvalues[4][ielem] = W[4];
        }

        /* write meta data for vtk */
        snprintf(vtk_data[0].description, BUFSIZ, "Density");
        vtk_data[0].type = T8_VTK_SCALAR;
        vtk_data[0].data = wvalues[0];

        snprintf(vtk_data[1].description, BUFSIZ, "Flow-U");
        vtk_data[1].type = T8_VTK_SCALAR;
        vtk_data[1].data = wvalues[1];

        snprintf(vtk_data[2].description, BUFSIZ, "Flow-V");
        vtk_data[2].type = T8_VTK_SCALAR;
        vtk_data[2].data = wvalues[2];

        snprintf(vtk_data[3].description, BUFSIZ, "Flow-W");
        vtk_data[3].type = T8_VTK_SCALAR;
        vtk_data[3].data = wvalues[3];

        snprintf(vtk_data[4].description, BUFSIZ, "Energy");
        vtk_data[4].type = T8_VTK_SCALAR;
        vtk_data[4].data = wvalues[4];

        /* write filename */
        snprintf(fileprefix, BUFSIZ, "./solution/t8wyo_flow_%03i",*vtk_counter);

        /* write vtk files */
        t8_forest_write_vtk_ext(t8wyo_forest,fileprefix,
                                1, //write_treeid
                                1, //write_mpirank
                                1, //write_level
                                1, //write_element_id
                                0, //write_ghosts
                                0, //write_curved
                                0, //do_not_use_API
                                5, //num_data
                                vtk_data);
        /* clean-up */
        T8_FREE (wvalues[0]);
        T8_FREE (wvalues[1]);
        T8_FREE (wvalues[2]);
        T8_FREE (wvalues[3]);
        T8_FREE (wvalues[4]);
    );
    /* display vtk message */
    if(t8wyo.ctx.rank==0) printf("[t8wyo] VTK Complete: %f (sec)\n",vtk_time);
}