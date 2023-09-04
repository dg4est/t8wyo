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

/* system header files */
#include <unordered_set>

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
    t8_element_t *elem;     /* local element */
    t8_locidx_t ltree_id;    /* local tree id this face belongs to */
    t8_locidx_t lelem_id;   /* local element id in tree */
    t8_locidx_t elem_id;    /* local element id this face belongs to */
    int face_number;        /* face number within the element */
}
t8wyo_face_full_t;

typedef struct {
    t8_locidx_t e1; /* element 1 on face */
    t8_locidx_t e2; /* element 2 on face */
    t8_locidx_t face_index; /* local face index */
    t8_locidx_t nvert;   /* number of face vertices */
    Real normal[3];
    Real area;
    Real sign;
}
t8wyo_face_t;

#define TOL 1.0E-12
struct Node {
    int id;
    double x;
    double y;
    double z;

    Node() { }
    Node(int id,double *geo){
        this->id = id;
        this->x = geo[0];
        this->y = geo[1];
        this->z = geo[2];
    }

    bool operator==(const Node& otherNode) const {
        return (otherNode.id == id) || ((abs(this->x - otherNode.x) <= TOL) &&
                                        (abs(this->y - otherNode.y) <= TOL) &&
                                        (abs(this->z - otherNode.z) <= TOL));
    }

    struct HashFunction {
        size_t operator()(const Node& node) const {
            size_t xHash = std::hash<int>()(int(node.x));
            size_t yHash = std::hash<int>()(int(node.y)) << 1;
            size_t zHash = std::hash<int>()(int(node.z)) << 2;
            return xHash ^ yHash ^ zHash;
        }
    };
};

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
                           wyo::memory<int> &ndc4,
                           wyo::memory<int> &ndc5,
                           wyo::memory<int> &ndc6,
                           wyo::memory<int> &ndc8,
                           wyo::memory<Real> &xgeom,
                           wyo::memory<Real> &elem_vol,
                           wyo::memory<Real> &face_norm);

#ifdef __cplusplus
}
#endif
#endif /* T8WYO_BUILD_LISTS_EXT_H */