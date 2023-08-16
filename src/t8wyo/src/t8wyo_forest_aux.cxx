/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* header files */
#include "t8wyo_forest_aux.h"
//#include "t8_element.h"

static t8_locidx_t
t8_forest_bin_search_lower (t8_element_array_t *elements,
                            t8_linearidx_t element_id, int maxlevel)
{
  t8_element_t       *query;
  t8_linearidx_t      query_id;
  t8_locidx_t         low, high, guess;
  t8_eclass_scheme_c *ts;

  ts = t8_element_array_get_scheme (elements);
  /* At first, we check whether any element has smaller id than the
   * given one. */
  query = t8_element_array_index_int (elements, 0);
  query_id = ts->t8_element_get_linear_id (query, maxlevel);
  if (query_id > element_id) {
    /* No element has id smaller than the given one */
    return -1;
  }

  /* We now perform the binary search */
  low = 0;
  high = t8_element_array_get_count (elements) - 1;
  while (low < high) {
    guess = (low + high + 1) / 2;
    query = t8_element_array_index_int (elements, guess);
    query_id = ts->t8_element_get_linear_id (query, maxlevel);
    if (query_id == element_id) {
      /* we are done */
      return guess;
    }
    else if (query_id > element_id) {
      /* look further left */
      high = guess - 1;
    }
    else {
      /* look further right, but keep guess in the search range */
      low = guess;
    }
  }
  T8_ASSERT (low == high);
  return low;
}

void t8wyo_forest_leaf_face_neighbors(t8_forest_t forest, t8_locidx_t ltreeid,
                                      const t8_element_t *leaf,
                                      int face,
                                      int *num_neighbors,
                                      t8_locidx_t *pelement_indices,
                                      int forest_is_balanced){
  t8_eclass_t         neigh_class, eclass;
  t8_gloidx_t         gneigh_treeid;
  t8_locidx_t         lneigh_treeid = -1;
  t8_locidx_t         lghost_treeid = -1, *element_indices, element_index;
  t8_eclass_scheme_c *ts, *neigh_scheme;
  t8_element_array_t *element_array;
  t8_element_t       *ancestor, *neighbor_leafs[4];
  t8_linearidx_t      neigh_id;
  int                 num_children_at_face, at_maxlevel;
  int                 ineigh, different_owners, have_ghosts;
  int                 owners[4],dual_faces[4];

  /* TODO: implement is_leaf check to apply to leaf */
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (!forest_is_balanced || t8_forest_is_balanced (forest));

    if (forest_is_balanced) {
        /* In a balanced forest, the leaf neighbor of a leaf is either the neighbor element itself,
         * its parent or its children at the face.
         */
        eclass = t8_forest_get_tree_class (forest, ltreeid);
        ts = t8_forest_get_eclass_scheme (forest, eclass);

        /* At first we compute these children of the face neighbor elements of leaf. For this, we need the
         * neighbor tree's eclass, scheme, and tree id */
        neigh_class = t8_forest_element_neighbor_eclass (forest, ltreeid, leaf, face);
        neigh_scheme = t8_forest_get_eclass_scheme (forest, neigh_class);

        /* If we are at the maximum refinement level, we compute the neighbor instead */
        at_maxlevel = ts->t8_element_level (leaf) == t8_forest_get_maxlevel (forest);

        if (at_maxlevel) {
            num_children_at_face = 1;
            neigh_scheme->t8_element_new (num_children_at_face, neighbor_leafs);

            /* Compute neighbor element and global treeid of the neighbor */
            gneigh_treeid =
                t8_forest_element_face_neighbor (forest, ltreeid, leaf,
                                                 neighbor_leafs[0], neigh_scheme,
                                                 face, dual_faces);
        } else {
            /* Allocate neighbor element */
            num_children_at_face = ts->t8_element_num_face_children (leaf, face);
            neigh_scheme->t8_element_new (num_children_at_face, neighbor_leafs);

            /* Compute neighbor elements and global treeid of the neighbor */
            gneigh_treeid =
                t8_forest_element_half_face_neighbors (forest, ltreeid, leaf,
                                                       neighbor_leafs,
                                                       neigh_scheme, face,
                                                       num_children_at_face,
                                                       dual_faces);
        }

        if (gneigh_treeid < 0) {
            /* There exists no face neighbor across this face, we return with this info */
            neigh_scheme->t8_element_destroy (num_children_at_face, neighbor_leafs);
            *num_neighbors = 0;
            return;
        }
        T8_ASSERT (gneigh_treeid >= 0
                   && gneigh_treeid < forest->global_num_trees);

        /* We have computed the half face neighbor elements, we now compute their owners,
         * if they differ, we know that the half face neighbors are the neighbor leafs.
         * If the owners do not differ, we have to check if the neighbor leaf is their
         * parent or grandparent. */
        different_owners = 0;
        have_ghosts = 0;
        for (ineigh = 0; ineigh < num_children_at_face; ineigh++) {
          /* At first, we check whether the current rank owns the neighbor, since
           * this is a constant time check and it is the most common case */
            if (t8_forest_element_check_owner (forest, neighbor_leafs[ineigh],
                                               gneigh_treeid, neigh_class,
                                               forest->mpirank, at_maxlevel)) {
                owners[ineigh] = forest->mpirank;

                /* The neighbor tree is also a local tree. we store its local treeid */
                lneigh_treeid = t8_forest_get_local_id (forest, gneigh_treeid);
            } else {
                owners[ineigh] =
                  t8_forest_element_find_owner (forest, gneigh_treeid,
                                                neighbor_leafs[ineigh], neigh_class);
                /* Store that at least one neighbor is a ghost */
                have_ghosts = 1;
            }
            if (ineigh > 0) {
                /* Check if all owners are the same for all neighbors or not */
                different_owners = different_owners || (owners[ineigh] != owners[ineigh - 1]);
            }
        }

        if (have_ghosts) {
            /* At least one neighbor is a ghost, we compute the ghost treeid of the neighbor tree. */
            lghost_treeid = t8_forest_ghost_get_ghost_treeid (forest, gneigh_treeid);
            T8_ASSERT (lghost_treeid >= 0);
        }

        /* TODO: Maybe we do not need to compute the owners. It suffices to know
         *       whether the neighbor is owned by mpirank or not.
         */
        if (!different_owners) {
            /* The face neighbors belong to the same process, we thus need to determine
             * if they are leafs or their parent or grandparent.
             */
            neigh_id = neigh_scheme->t8_element_get_linear_id(neighbor_leafs[0],forest->maxlevel);

            if (owners[0] != forest->mpirank) {
                /* The elements are ghost elements of the same owner */
                element_array = t8_forest_ghost_get_tree_elements(forest, lghost_treeid);

                /* Find the index in element_array of the leaf ancestor of the first neighbor.
                 * This is either the neighbor itself or its parent, or its grandparent */
                element_index = t8_forest_bin_search_lower(element_array, neigh_id, forest->maxlevel);

                /* Get the element */
                ancestor = t8_forest_ghost_get_element(forest, lghost_treeid, element_index);

                /* Add the number of ghost elements on previous ghost trees and the number
                 * of local elements. */
                element_index += t8_forest_ghost_get_tree_element_offset (forest, lghost_treeid);
                element_index += t8_forest_get_local_num_elements (forest);

                T8_ASSERT (forest->local_num_elements <= element_index
                           && element_index <
                           forest->local_num_elements +
                           t8_forest_get_num_ghosts (forest));
            } else {
                /* the elements are local elements */
                element_array = t8_forest_get_tree_element_array (forest, lneigh_treeid);

                /* Find the index in element_array of the leaf ancestor of the first neighbor.
                 * This is either the neighbor itself or its parent, or its grandparent */
                element_index = t8_forest_bin_search_lower(element_array, neigh_id,
                                                           forest->maxlevel);
                /* Get the element */
                ancestor = t8_forest_get_tree_element(t8_forest_get_tree(forest, lneigh_treeid),
                                                      element_index);

                /* Add the element offset of this tree to the index */
                element_index += t8_forest_get_tree_element_offset(forest, lneigh_treeid);
            }

            if (neigh_scheme->t8_element_compare(ancestor,neighbor_leafs[0]) < 0) {
                /* set return values */
                *num_neighbors = 1;
                pelement_indices[0] = element_index;

                /* clean up memory */
                neigh_scheme->t8_element_destroy (num_children_at_face, neighbor_leafs);
                return;
            }
        }

        /* The leafs are the face neighbors that we are looking for. */
        /* The face neighbors either belong to different processes and thus must be leafs
         * in the forest, or the ancestor leaf of the first half neighbor is the half
         * neighbor itself and thus all half neighbors must be leafs.
         * Since the forest is balanced, we found all neighbor leafs.
         * It remains to compute their local ids */
        *num_neighbors = num_children_at_face;
        element_indices = pelement_indices;
        for (ineigh = 0; ineigh < num_children_at_face; ineigh++) {
            /* Compute the linear id at maxlevel of the neighbor leaf */
            neigh_id = neigh_scheme->t8_element_get_linear_id (neighbor_leafs[ineigh], forest->maxlevel);

            /* Get a pointer to the element array in which the neighbor lies and search
             * for the element's index in this array.
             * This is either the local leaf array of the local tree or the corresponding leaf array
             * in the ghost structure
             */
            if (owners[ineigh] == forest->mpirank) {
                /* The neighbor is a local leaf */
                element_array = t8_forest_get_tree_element_array (forest, lneigh_treeid);

                /* Find the index of the neighbor in the array */
                element_indices[ineigh] = t8_forest_bin_search_lower (element_array, neigh_id, forest->maxlevel);
                T8_ASSERT(element_indices[ineigh] >= 0);

                /* We have to add the tree's element offset to the index found to get the actual local element id */
                element_indices[ineigh] += t8_forest_get_tree_element_offset (forest, lneigh_treeid);
            } else {
                /* The neighbor is a ghost */
                element_array = t8_forest_ghost_get_tree_elements(forest, lghost_treeid);

                /* Find the index of the neighbor in the array */
                element_indices[ineigh] = t8_forest_bin_search_lower(element_array, neigh_id, forest->maxlevel);

                /* Add the element offset of previous ghosts to this index */
                element_indices[ineigh] += t8_forest_ghost_get_tree_element_offset(forest, lghost_treeid);

                /* Add the number of all local elements to this index */
                element_indices[ineigh] += t8_forest_get_local_num_elements(forest);
            }
        }
        /* clean up memory */
        neigh_scheme->t8_element_destroy (num_children_at_face, neighbor_leafs);
    } else {
        /* TODO: implement unbalanced version */
        SC_ABORT_NOT_REACHED ();
    }
}