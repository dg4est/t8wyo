/**
 * File:   t8wyo_forest_aux.h
 * Author: akirby
 *
 * Created on August 10, 2023, 1:13 PM
 */

#ifndef T8WYO_FOREST_AUX_H
#define T8WYO_FOREST_AUX_H

/* 3PL header files */
#include <t8_cmesh.h>
#include <t8_element.h>
#include <t8_data/t8_containers.h>

#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_balance.h>
#include <t8_element_cxx.hxx>

#ifdef __cplusplus
extern "C" {
#endif

void t8wyo_forest_leaf_face_neighbors(t8_forest_t forest, t8_locidx_t ltreeid,
                                      const t8_element_t *leaf,
                                      int face,
                                      int *num_neighbors,
                                      t8_locidx_t *pelement_indices,
                                      int forest_is_balanced);

#ifdef __cplusplus
}
#endif
#endif /* T8WYO_FOREST_AUX_H */