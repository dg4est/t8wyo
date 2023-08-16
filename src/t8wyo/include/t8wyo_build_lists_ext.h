/**
 * File:   t8wyo_build_lists_ext.h
 * Author: akirby
 *
 * Created on August 1, 2023, 4:57 PM
 */

#ifndef T8WYO_BUILD_LISTS_EXT_H
#define T8WYO_BUILD_LISTS_EXT_H

/* header files */
#include "t8wyo_solver.hxx"
#include "t8wyo_forest_aux.h"

/* 3PL header files */
#include <sc.h>
#include <t8.h>
#include <t8_element_c_interface.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ======= */
/* Defines */
/* ======= */
/* cmesh geometry offset info */
#define T8WYO_CMESH_OFFSET_KEY  T8_CMESH_NEXT_POSSIBLE_KEY
#define T8WYO_FID_KEY 0
#define T8WYO_FTYPE_KEY 1

/* =============== */
/* Data Structures */
/* =============== */
/* struct stores all information associated to a elements face */
typedef struct {
    t8_locidx_t tree_id;    /* local tree id this face belongs to */
    t8_locidx_t lelem_id;   /* local element id in tree */
    t8_locidx_t elem_id;    /* local element id this face belongs to */
    int face_number;     /* face number within the element */
    int num_neighbors;      /* number of elements on face */
    t8_locidx_t neighids[4];  /* neighbor element indices */
}
t8wyo_face_full_t;

typedef struct {
    t8_locidx_t e1; /* element 1 on face */
    t8_locidx_t e2; /* element 2 on face */
    t8_locidx_t face_index; /* local face index */
    t8_locidx_t nvert;   /* number of face vertices */
    Real normal[3];
    Real area;
}
t8wyo_face_t;

/* ========= */
/* functions */
/* ========= */

/** t8wyo connectivity construction.
 *
 * @param [in] forest   t8code forest data structure
 */
void t8wyo_build_lists_ext(t8_cmesh_t cmesh,
                           t8_forest_t forest,
                           wyo::memory<int> &face2cell,
                           wyo::memory<int> &facetype,
                           wyo::memory<int> &elem_info,
                           wyo::memory<Real> &elem_vol,
                           wyo::memory<Real> &face_norm);

#ifdef __cplusplus
}
#endif
#endif /* T8WYO_BUILD_LISTS_EXT_H */