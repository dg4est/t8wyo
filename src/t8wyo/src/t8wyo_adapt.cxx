/** t8wyo adapt functions.
 * \file    t8wyo_adapt.cxx
 * \author  akirby
 *
 * \brief   t8wyo adaption functions.
 */

/* header files */
#include "t8wyo_adapt.h"
#include "memory_utilities.h"

/* Adapt the forest and interpolate the phi values to the new grid,
 * compute the new u values on the grid */
t8_forest_t t8wyo_adapt_ext(t8_forest_t forest,
                            tag_callback_t *tag_function,
                            int *nvar_in,int *ncell,int *ncell_real,
                            Real *wvalues,wyo::memory<Real> &wvalues_new){
    t8_forest_t forest_adapt;
    t8_locidx_t num_elems,num_elems_p_ghosts;
  //double adapt_time, balance_time, ghost_time, replace_time;
    int nvar;

    adapt_info_t adapt_info;
    adapt_info.tag_func = tag_function;
    adapt_info.wvalues = wvalues;
    adapt_info.nvar = nvar = *nvar_in;

    /* adapt the forest, but keep the old one */
    t8_forest_ref(forest);
    t8_forest_init(&forest_adapt);
    t8_forest_set_adapt(forest_adapt, forest, t8wyo_tag_callback, 0);
    t8_forest_set_ghost(forest_adapt, 1, T8_GHOST_FACES);
    t8_forest_set_balance(forest_adapt, NULL, 1);
    t8_forest_set_user_data(forest_adapt, &adapt_info);
    t8_forest_commit(forest_adapt);

    /* allocate new memory for the element_data of the adapted forest */
    num_elems = t8_forest_get_local_num_elements(forest_adapt);
    num_elems_p_ghosts = num_elems + t8_forest_get_num_ghosts(forest_adapt);

    /* set external cell counts and allocate new solution array */
    *ncell_real = num_elems;
    *ncell = num_elems_p_ghosts;
    wvalues_new.malloc(nvar*num_elems_p_ghosts);
    adapt_info.wvalues_new = wvalues_new.ptr();

    /* We now call iterate_replace in which we interpolate the new element data.
     * It is necessary that the old and new forest only differ by at most one level.
     * We guarantee this by calling adapt non-recursively and calling balance without
     * repartitioning.
     */
    //replace_time = -sc_MPI_Wtime();
    t8_forest_iterate_replace(forest_adapt,forest,t8wyo_adapt_replace);
    //replace_time += sc_MPI_Wtime();

    /* Free memory for the forest */
    t8_forest_unref(&forest);
    return forest_adapt;
}

int t8wyo_tag_callback(t8_forest_t forest,
                       t8_forest_t forest_from,
                       t8_locidx_t ltree_id,
                       t8_locidx_t lelement_id,
                       t8_eclass_scheme_c *ts,
                       const int is_family,
                       const int num_elements,
                       t8_element_t *elements[]){
    t8_locidx_t offset = t8_forest_get_tree_element_offset(forest_from,ltree_id);
    t8_locidx_t ielement = lelement_id + offset;
    int level = ts->t8_element_level(elements[0]);
    int tag;

    // retrieve external callback tagging function
    adapt_info_t *adapt_info = (adapt_info_t *) t8_forest_get_user_data(forest);
    int nvar =  adapt_info->nvar;

    // set solution pointer for cell
    Real *Wvalues_cell = &adapt_info->wvalues[nvar*ielement];

    // set tag value based on callback function
    // tag: (1) refine this element
    //      (0) do not change this element
    //     (-1) coarsen this element (only is_family=1)
    ielement += FBASE;
    (*(adapt_info->tag_func))(&ielement,&level,&nvar,Wvalues_cell,&tag);

    // if coarsen, check if is_family, since returning < 0 is illegal
    if(tag == -1) return (is_family) ? -1:0;

    // return tag otherwise
    return tag;
}

/* Replace callback to decide how to interpolate a refined or coarsened element.
 * If an element is refined, each child gets the phi value of its parent.
 * If elements are coarsened, the parent gets the average phi value of the children.
 * >> outgoing are the old elements and incoming the new ones
 */
void t8wyo_adapt_replace(t8_forest_t forest_old,
                         t8_forest_t forest_new,
                         t8_locidx_t which_tree,
                         t8_eclass_scheme_c *ts,
                         int refine,
                         int num_outgoing,t8_locidx_t first_outgoing,
                         int num_incoming,t8_locidx_t first_incoming){
    adapt_info_t *adapt_info = (adapt_info_t *) t8_forest_get_user_data(forest_new);
    const int nvar = adapt_info->nvar;

    t8_locidx_t first_incoming_data, first_outgoing_data;
    Real *elem_data_in, *elem_data_out;
    int i,n;

    /* Get pointers to the element data */
    first_incoming_data = first_incoming + t8_forest_get_tree_element_offset(forest_new,which_tree);
    first_outgoing_data = first_outgoing + t8_forest_get_tree_element_offset(forest_old,which_tree);

    /* fill new wvalues array */
    if (refine == 0) {
        /* The element is not changed, copy solution */
        T8_ASSERT (num_incoming == num_outgoing && num_incoming == 1);
        elem_data_out = &adapt_info->wvalues[nvar*first_outgoing_data];
        elem_data_in = &adapt_info->wvalues_new[nvar*first_incoming_data];

        memcpy(elem_data_in, elem_data_out, nvar*sizeof(Real));
    } else
    if (refine == 1) {
        /* The old element is refined:
         *   we copy the solution values
         */
        elem_data_out = &adapt_info->wvalues[nvar*first_outgoing_data];
        for (i = 0; i < num_incoming; i++) {
            elem_data_in = &adapt_info->wvalues_new[nvar*(first_incoming_data+i)];
            memcpy(elem_data_in, elem_data_out, nvar*sizeof(Real));
        }
    } else {
        /* The old elements form a family which is coarsened:
         *   we compute the average solution and set it as the new value
         */
        T8_ASSERT(refine = -1);

        /* Compute average of solution */
        Real wtotal[nvar]; memset(wtotal, 0, nvar*sizeof(Real));
        Real voltotal = {0.0};
        for (i = 0; i < num_outgoing; i++) {
            const t8_element_t *element = t8_forest_get_element_in_tree(forest_old,which_tree,first_outgoing+i);
            Real vol = t8_forest_element_volume(forest_old,first_outgoing+i,element);
            voltotal += vol;

            elem_data_out = &adapt_info->wvalues[nvar*(first_outgoing_data+i)];
            for (n = 0; n < nvar; n++) {
                wtotal[n] += elem_data_out[n]*vol;
            }
        }
        for(n = 0; n < nvar; n++) wtotal[n] /= (voltotal);

        elem_data_in = &adapt_info->wvalues_new[nvar*first_incoming_data];
        memcpy(elem_data_in, wtotal, nvar*sizeof(Real));
    }
}