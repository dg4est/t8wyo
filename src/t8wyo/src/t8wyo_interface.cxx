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
        t8wyo_cmesh = t8wyo_create_cmesh(MESH_MODE_MCELL,NULL,&mcell,
                                        *level_cmesh,dim,use_occ_geometry,
                                         t8wyo.ctx.comm);
    );

    if(t8wyo.ctx.rank==0) printf("[t8wyo] COARSE MESH CONSTRUCTION: %f (sec)\n",cmesh_time);
}

void t8wyo_build_forest_(int *level_forest,
                         int *ntetra,int *npyr,int *nprizm,int *nhex,
                         int *ntetra_ng,int *npyr_ng,int *nprizm_ng,int *nhex_ng,
                         int *nnode,int *ngpt){
//                         int *nbface3,int *nbface4,
//                         int *nface3,int *nface4,
//                         int *ndc4,int *ndc5,int *ndc6,int *ndc8,
//                         int *nbf3,int *nbf4,
//                         int *ifpat3,int *ifpat4,
//                         Real *xgeom,
//                         int *ndf3,int *ndf4){
    /* build forest */
    Real forest_time;
    TIMER(forest_time,
        t8wyo_forest = t8wyo_build_forest(t8wyo_cmesh,*level_forest,t8wyo.ctx.comm);
    );

    /* send back new grid info */
    t8_locidx_t num_trees,num_elems_in_tree,itree;
    t8_eclass_scheme_c *ts;
    int nelem_type[T8_ECLASS_COUNT] = {0};

    num_trees = t8_forest_get_num_local_trees(t8wyo_forest);
    for (itree = 0; itree < num_trees; itree++) {
        ts = t8_forest_get_eclass_scheme(t8wyo_forest,t8_forest_get_tree_class(t8wyo_forest,itree));
        num_elems_in_tree = t8_forest_get_tree_num_elements(t8wyo_forest,itree);

        // add to element type count
        nelem_type[ts->eclass] += num_elems_in_tree;
    }

    /* set local element counts */
    *ntetra = nelem_type[T8_ECLASS_TET];
    *npyr = nelem_type[T8_ECLASS_PYRAMID];
    *nprizm = nelem_type[T8_ECLASS_PRISM];
    *nhex = nelem_type[T8_ECLASS_HEX];

    /* set local element ghost type counts */
    // TODO: get ghost element type counts
    *ntetra_ng = 0;
    *npyr_ng = 0;
    *nprizm_ng = 0;
    *nhex_ng = 0;

    //printf("T8WYO ELEM COUNTS: %d %d %d %d\n",*ntetra,*npyr,*nprizm,*nhex);
    if(t8wyo.ctx.rank==0) printf("[t8wyo] FOREST CONSTRUCTION: %f (sec)\n",forest_time);
}

void t8wyo_build_lists_(int *ncell_real,int *ncell,int *nface,
                        int **face2cellptr,int **ifacetypeptr,
                        int **cellinfoptr,
                        Real **cellvolptr,Real **facenormptr){

    /* construct face2cell,facetype,elem_info,elem_vol data structures */
    Real lists_time;
    TIMER(lists_time,
        t8wyo_build_lists_ext(t8wyo_cmesh,t8wyo_forest,face2cell,ifacetype,
                              elem_info,elem_vol,face_norm);
    );

    t8_locidx_t num_elements = t8_forest_get_local_num_elements(t8wyo_forest);
    t8_locidx_t num_ghosts = t8_forest_get_num_ghosts(t8wyo_forest);

    /* set counts */
    *ncell_real = (int) num_elements;
    *ncell = (int) (num_elements+num_ghosts);
    *nface = face2cell.length()/2;

    /* fill Fortran data */
    *face2cellptr = face2cell.ptr();
    *ifacetypeptr = ifacetype.ptr();
    *cellinfoptr = elem_info.ptr();
    *cellvolptr = elem_vol.ptr();
    *facenormptr = face_norm.ptr();

    if(t8wyo.ctx.rank==0) printf("[t8wyo] BUILD LISTS CONSTRUCTION: %f (sec)\n",lists_time);
}

void t8wyo_exchange_ghost_data_(void *data,size_t *bytes_per_element,int *barrier_flag){
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
    if(t8wyo.ctx.rank==0) printf("[t8wyo] Exchange Time: %f (sec) w/ %zu bytes per elem\n",
                                 exc_time,*bytes_per_element);
}