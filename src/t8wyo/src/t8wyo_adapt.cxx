/** t8wyo adapt functions.
 * \file    t8wyo_adapt.cxx
 * \author  akirby
 *
 * \brief   t8wyo adaption functions.
 */

/* header files */
#include "t8wyo_adapt.h"

/* Adapt the forest and interpolate the phi values to the new grid,
 * compute the new u values on the grid */
void t8wyo_adapt_ext(t8_forest_t forest,t8_forest_t forest_adapt,
                     Real *wvalues,Real **wvalues_new){
    t8_locidx_t num_elems,num_elems_p_ghosts;
    double adapt_time, balance_time, ghost_time, replace_time;

    /* Adapt the forest, but keep the old one */
    t8_forest_ref(forest);
    t8_forest_init(&forest_adapt);

    /* Enable profiling to measure the runtime */
    //t8_forest_set_profiling(forest_adapt, 1);

    /* Set the user data pointer of the new forest */
//    t8_forest_set_user_data(forest_adapt,tag_callback);

    /* Set the adapt function */
    t8_forest_set_adapt(forest_adapt,forest,t8wyo_tag_callback, 0);

    /* balance forest */
    t8_forest_set_balance(forest_adapt, NULL, 1);

    /* We also want ghost elements in the new forest */
    t8_forest_set_ghost(forest_adapt, 1, T8_GHOST_FACES);

    /* Commit the forest, adaptation and balance happens here */
    t8_forest_commit(forest_adapt);

    /* Allocate new memory for the element_data of the adapted forest */
    num_elems = t8_forest_get_local_num_elements(forest_adapt);
    num_elems_p_ghosts = num_elems + t8_forest_get_num_ghosts(forest_adapt);

//    problem->element_data_adapt = sc_array_new_count(sizeof (t8_advect_element_data_t), num_elems);
//    problem->phi_values_adapt = sc_array_new_count((problem->dummy_op ? 2 : 1) * sizeof (double),
//                                                    num_elems_p_ghosts);

    /* We now call iterate_replace in which we interpolate the new element data.
     * It is necessary that the old and new forest only differ by at most one level.
     * We guarantee this by calling adapt non-recursively and calling balance without
     * repartitioning.
     */
    replace_time = -sc_MPI_Wtime();
    t8_forest_iterate_replace(forest_adapt,
                              forest,
                              t8wyo_adapt_replace);
    replace_time += sc_MPI_Wtime();

    /* clean the old element data */
//    t8_advect_problem_elements_destroy(problem);
//    sc_array_destroy(problem->element_data);
//    sc_array_destroy(problem->phi_values);

    /* Free memory for the forest */
    t8_forest_unref(&forest);

    /* Set the forest to the adapted one */
    forest = forest_adapt;
    forest_adapt = NULL;
}

int t8wyo_tag_callback(t8_forest_t forest, t8_forest_t forest_from,
                         t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                         t8_eclass_scheme_c *ts, const int is_family,
                         const int num_elements, t8_element_t *elements[]){
    t8_locidx_t offset;
    int tag;

    offset = t8_forest_get_tree_element_offset(forest_from,ltree_id);
    t8_locidx_t ielement = lelement_id + offset;

    //user_data = (t8_advect_problem_t *) t8_forest_get_user_data(forest);
    // TODO: set tag value based on callback function
    // amr_adapt_function(&ielement,&tag);
    // return tag;
    return 1;
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
                         int num_outgoing,
                         t8_locidx_t first_outgoing,
                         int num_incoming, t8_locidx_t first_incoming){
    t8_locidx_t first_incoming_data, first_outgoing_data;
//    t8_advect_problem_t *problem;
//    t8_advect_element_data_t *elem_data_in, *elem_data_out;
//    t8_element_t *element;
//    int i, iface;

    /* Get pointers to the element data */
    first_incoming_data = first_incoming + t8_forest_get_tree_element_offset(forest_new,
                                                                             which_tree);
    first_outgoing_data = first_outgoing + t8_forest_get_tree_element_offset(forest_old,
                                                                             which_tree);

//    elem_data_out = (t8_advect_element_data_t *)
//        t8_sc_array_index_locidx(problem->element_data,first_outgoing_data);
//
//    elem_data_in = (t8_advect_element_data_t *)
//        t8_sc_array_index_locidx(problem->element_data_adapt,
//                                 first_incoming_data);

    /* Get the old phi value (used in the cases with num_outgoing = 1) */
//    phi_old = t8_advect_element_get_phi(problem, first_outgoing_data);
    if (refine == 0) {
        /* The element is not changed, copy solution */
        T8_ASSERT (num_incoming == num_outgoing && num_incoming == 1);

//        memcpy(elem_data_in, elem_data_out, sizeof (t8_advect_element_data_t));
//        t8_advect_element_set_phi_adapt(problem, first_incoming_data, phi_old);
//
//        /* Get a pointer to the new element */
//        element = t8_forest_get_element_in_tree(problem->forest_adapt, which_tree,
//                                                first_incoming);
//        /* Set the neighbor entries to uninitialized */
//        T8_ASSERT (elem_data_in->num_faces == ts->t8_element_num_faces (element));
//        for (iface = 0; iface < elem_data_in->num_faces; iface++) {
//            elem_data_in->num_neighbors[iface] = 0;
//            elem_data_in->flux_valid[iface] = -1;
//            elem_data_in->dual_faces[iface] = NULL;
//            elem_data_in->fluxes[iface] = NULL;
//            elem_data_in->neighs[iface] = NULL;
//        }
    } else
    if (refine == 1) {
        /* The old element is refined:
         *   we copy the solution values and compute the new midpoints
         */
//        for (i = 0; i < num_incoming; i++) {
//            /* Get a pointer to the new element */
//            element =
//            t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
//                                           first_incoming + i);
//            /* Compute midpoint and vol of the new element */
//            t8_advect_compute_element_data (problem, elem_data_in + i, element,
//                                            which_tree, ts);
//            t8_advect_element_set_phi_adapt(problem, first_incoming_data + i,
//                                            phi_old);
//
//            /* Set the neighbor entries to uninitialized */
//            elem_data_in[i].num_faces = elem_data_out->num_faces;
//
//            for (iface = 0; iface < elem_data_in[i].num_faces; iface++) {
//                elem_data_in[i].num_neighbors[iface] = 0;
//                elem_data_in[i].flux_valid[iface] = -1;
//                elem_data_in[i].dual_faces[iface] = NULL;
//                elem_data_in[i].fluxes[iface] = NULL;
//                elem_data_in[i].neighs[iface] = NULL;
//            }
//            /* Update the level */
//            elem_data_in[i].level = elem_data_out->level + 1;
//        }
    } else {
        /* The old elements form a family which is coarsened:
         *   we compute the average phi value and set it as the new value
         */
        T8_ASSERT(refine = -1);

//        element = t8_forest_get_element_in_tree(problem->forest_adapt,which_tree,
//                                                first_incoming);
//
//        /* Compute midpoint and vol of the new element */
//        t8_advect_compute_element_data(problem, elem_data_in, element,
//                                       which_tree, ts);
//
//        /* Compute average of phi */
//        double phi = 0;
//        for (i = 0; i < num_outgoing; i++) {
//            phi += t8_advect_element_get_phi(problem, first_outgoing_data + i);
//        }
//        phi /= num_outgoing;
//
//        t8_advect_element_set_phi_adapt(problem, first_incoming_data, phi);
//
//        /* Set the neighbor entries to uninitialized */
//        elem_data_in->num_faces = elem_data_out[0].num_faces;
//        T8_ASSERT (elem_data_in->num_faces == ts->t8_element_num_faces (element));
//        for (iface = 0; iface < elem_data_in->num_faces; iface++) {
//            elem_data_in->num_neighbors[iface] = 0;
//            elem_data_in->flux_valid[iface] = -1;
//            elem_data_in->dual_faces[iface] = NULL;
//            elem_data_in->fluxes[iface] = NULL;
//            elem_data_in->neighs[iface] = NULL;
//        }
//
//        /* update the level */
//        elem_data_in->level = elem_data_out->level - 1;
    }
}