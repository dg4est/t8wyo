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
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>

/* system header files */
#include <assert.h>
#include <unordered_set>
#include <set>

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

#define MAX_FACES 6
#define MAX_SUBFACES 4
#define MAX_SUBELEMENTS 10

/* =============== */
/* Data Structures */
/* =============== */
/* struct stores all information associated to a elements face */
typedef struct {
    const t8_element_t *elem; /* local element */
    t8_locidx_t ltree_id;   /* local tree id this face belongs to */
    t8_locidx_t lelem_id;   /* local element id in tree */
    t8_locidx_t elem_id;    /* local element id this face belongs to */
    int face_number;        /* face number within the element */
    int nvert;
    int orientation;
}
t8wyo_face_full_t;

typedef struct {
    t8_locidx_t e1; /* element 1 on face */
    t8_locidx_t e2; /* element 2 on face */
    int iface1;
    int iface2;
    int orientation;
    t8_locidx_t face_index; /* local face index */
    t8_locidx_t nvert;   /* number of face vertices */
    Real normal[3];
    Real area;
}
t8wyo_face_t;

struct Mortar {
    t8_locidx_t elemidx_minus; /**< Local index of the element corresponding to minus (big element) */
    int iface_minus; /**< Local face index inside the big element */
    mutable int cnt; // counter for hanging faces filled
    mutable t8_locidx_t elemidx_plus[MAX_SUBFACES]; /**< Local index of the elements corresponding to plus (hanging elements) */
    mutable int iface_plus[MAX_SUBFACES]; /**< Local face indices inside the hanging elements */

    Mortar(int elemidx,int iface) : cnt{0} {
        elemidx_minus = elemidx;
        iface_minus = iface;

        // initialize hanging elements
        elemidx_plus[0]=elemidx_plus[1]=elemidx_plus[2]=elemidx_plus[3] = -1;
        iface_plus[0]=iface_plus[1]=iface_plus[2]=iface_plus[3] = -1;
    }

    bool operator==(const Mortar& otherMortar) const {
        return (otherMortar.elemidx_minus == this->elemidx_minus) &&
               (otherMortar.iface_minus == this->iface_minus);
    }

    struct HashFunction {
        size_t operator()(const Mortar& mortar) const {
            return std::hash<t8_locidx_t>()(mortar.elemidx_minus) ^ std::hash<int>()(mortar.iface_minus);
        }
    };

    struct KeyEqual {
        bool operator()(const Mortar& lhs, const Mortar& rhs) const {
            return lhs == rhs;
        }
    };
};

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
        return (otherNode.id == id) || ((std::abs(this->x - otherNode.x) <= TOL) &&
                                        (std::abs(this->y - otherNode.y) <= TOL) &&
                                        (std::abs(this->z - otherNode.z) <= TOL));
    }

    struct HashFunction {
        size_t operator()(const Node& node) const {
            size_t xHash = std::hash<int>()(int(1000.0*node.x));
            size_t yHash = std::hash<int>()(int(1000.0*node.y)) << 1;
            size_t zHash = std::hash<int>()(int(1000.0*node.z)) << 2;
            return xHash ^ yHash ^ zHash;
        }
    };

    struct CompareFunction {
        bool operator()(const Node &n1, const Node &n2) const {
            if(n1.x < n2.x - TOL) return true;
            if(n1.x > n2.x + TOL) return false;
            if(n1.y < n2.y - TOL) return true;
            if(n1.y > n2.y + TOL) return false;
            if(n1.z < n2.z - TOL) return true;
            return false;
        }
    };
};

static void correct_neg_volume(t8_eclass eclass,double *tree_vertices){
   /* The volume described is negative-- change vertices
    *   tets:     switch 0 and 3
    *   prisms:   switch 0 and 3, 1 and 4, 2 and 5
    *   pyramids: switch 0 and 4
    *   hexs:     switch 0 and 4, 1 and 5, 2 and 6, 3 and 7
    */
    double temp;
    int i;
    int iswitch;
    int num_switches = 0;
    int switch_indices[4] = { 0 };
    T8_ASSERT(t8_eclass_to_dimension[eclass] == 3);

    switch (eclass) {
        case T8_ECLASS_TET:
            /* We switch vertex 0 and vertex 3 */
            num_switches = 1;
            switch_indices[0] = 3;
            break;
        case T8_ECLASS_PRISM:
            num_switches = 3;
            switch_indices[0] = 3;
            switch_indices[1] = 4;
            switch_indices[2] = 5;
            break;
        case T8_ECLASS_HEX:
            num_switches = 4;
            switch_indices[0] = 4;
            switch_indices[1] = 5;
            switch_indices[2] = 6;
            switch_indices[3] = 7;
            break;
        case T8_ECLASS_PYRAMID:
            num_switches = 1;
            switch_indices[0] = 4;
            break;
        default:
            SC_ABORT_NOT_REACHED();
    }

    /* switch vertex 0 + iswitch and vertex switch_indices[iswitch] */
    for (iswitch = 0; iswitch < num_switches; ++iswitch) {
        for (i = 0; i < 3; i++) {
            temp = tree_vertices[3 * iswitch + i];
            tree_vertices[3 * iswitch + i] = tree_vertices[3 * switch_indices[iswitch] + i];
            tree_vertices[3 * switch_indices[iswitch] + i] = temp;
        }
    }
}

static void correct_node_ordering(t8_eclass eclass,int *nodes){
   /* The volume described is negative-- change vertices
    *   tets:     switch 0 and 3
    *   prisms:   switch 0 and 3, 1 and 4, 2 and 5
    *   pyramids: switch 0 and 4
    *   hexs:     switch 0 and 4, 1 and 5, 2 and 6, 3 and 7
    */
    int temp;
    int iswitch;
    int num_switches = 0;
    int switch_indices[4] = { 0 };
    T8_ASSERT(t8_eclass_to_dimension[eclass] == 3);

    switch (eclass) {
        case T8_ECLASS_TET:
            /* We switch vertex 0 and vertex 3 */
            num_switches = 1;
            switch_indices[0] = 3;
            break;
        case T8_ECLASS_PRISM:
            num_switches = 3;
            switch_indices[0] = 3;
            switch_indices[1] = 4;
            switch_indices[2] = 5;
            break;
        case T8_ECLASS_HEX:
            num_switches = 4;
            switch_indices[0] = 4;
            switch_indices[1] = 5;
            switch_indices[2] = 6;
            switch_indices[3] = 7;
            break;
        case T8_ECLASS_PYRAMID:
            num_switches = 1;
            switch_indices[0] = 4;
            break;
        default:
            SC_ABORT_NOT_REACHED();
    }

    /* switch vertex 0 + iswitch and vertex switch_indices[iswitch] */
    for (iswitch = 0; iswitch < num_switches; ++iswitch) {
        temp = nodes[iswitch];
        nodes[iswitch] = nodes[switch_indices[iswitch]];
        nodes[switch_indices[iswitch]] = temp;
    }
}

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
                           wyo::memory<int> &mortar_info,
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
