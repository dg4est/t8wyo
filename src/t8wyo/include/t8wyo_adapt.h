/**
 * File:   t8wyo_adapt.h
 * Author: akirby
 *
 * Created on August 2, 2023, 9:20 PM
 */

#ifndef T8WYO_ADAPT_H
#define T8WYO_ADAPT_H

/* header files */
#include "precision_types.h"

/* 3PL header files */
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_iterate.h>

#ifdef __cplusplus
extern "C" {
#endif

void t8wyo_adapt_ext(t8_forest_t forest,t8_forest_t forest_adapt,
                     Real *wvalues,Real **wvalues_new);

int t8wyo_tag_callback(t8_forest_t forest, t8_forest_t forest_from,
                       t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                       t8_eclass_scheme_c *ts, const int is_family,
                       const int num_elements, t8_element_t *elements[]);

void t8wyo_adapt_replace(t8_forest_t forest_old,
                         t8_forest_t forest_new,
                         t8_locidx_t which_tree,
                         t8_eclass_scheme_c *ts,
                         int refine,
                         int num_outgoing,
                         t8_locidx_t first_outgoing,
                         int num_incoming, t8_locidx_t first_incoming);

#ifdef __cplusplus
}
#endif
#endif /* T8WYO_ADAPT_H */