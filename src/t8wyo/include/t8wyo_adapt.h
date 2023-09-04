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
#include "memory.hxx"

/* 3PL header files */
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

#ifdef __cplusplus
extern "C" {
#endif

/** AMR Tagging Callback Function Template.
 *
 * @param [in] elem_id      element id
 * @param [in] level        level number
 * @param [in] nvar         number of variables per cell
 * @param [in] Wvalues_cell solution values stored in cell
 * @param [out] tag         adaption tag: (-1) coarsen, (0) unchanged, (1) refine
 */
typedef void (*tag_callback_t)(int *elem_id,
                               int *level,
                               int *nvar,
                               Real *Wvalues_cell,
                               int *tag);

/** Adaption Information Data Structure */
typedef struct {
    tag_callback_t *tag_func; /**< tagging callback function */
    Real *wvalues;            /**< solution values */
    Real *wvalues_new;        /**< new solution values */
    int nvar;                 /**< number of variables per cell */
}
adapt_info_t;

t8_forest_t t8wyo_adapt_ext(t8_forest_t forest,
                            tag_callback_t *tag_function,
                            int *nvar_in,int *ncell,int *ncell_real,
                            Real *wvalues,wyo::memory<Real> &wvalues_new);

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