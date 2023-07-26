/**
 * File:   t8wyo_cmesh_from_mcell.h
 * Author: akirby
 *
 * Created on June 24, 2023, 11:37 AM
 */

#ifndef T8WYO_CMESH_FROM_MCELL_FILE_H
#define T8WYO_CMESH_FROM_MCELL_FILE_H

/* header files */
#include "precision_types.h"
#include "memory.hxx"
#include "t8wyo_cell3d.h"

/* 3PL header files */
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_eclass.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_cmesh/t8_cmesh_types.h> // INTERNAL HACK: copy to 3PL/t8code
#include <t8_cmesh/t8_cmesh_stash.h> // INTERNAL HACK: copy to 3PL/t8code

//#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>
//#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.h>

#ifdef __cplusplus
extern "C" {
#endif

/* offset for cmesh geometry info */
#define T8WYO_CMESH_OFFSET_KEY  T8_CMESH_NEXT_POSSIBLE_KEY
#define T8WYO_FID_KEY 0
#define T8WYO_FTYPE_KEY 1

/* =============== */
/* data structures */
/* =============== */
typedef struct {
    t8_locidx_t index;
    double coordinates[3];
}
t8wyo_mcell_node_t;

typedef struct {
    t8_locidx_t index;
    double coordinates[3];
    double parameters[2];
    int parametric;
    int entity_dim;
    t8_locidx_t entity_tag;
}
t8wyo_mcell_node_parametric_t;

/* ========= */
/* functions */
/* ========= */
t8_cmesh_t
t8_cmesh_from_mcell(mcell_t *mcell,
                    int do_bcast,MPI_Comm comm,
                    int dim,int use_occ_geometry);

#ifdef __cplusplus
}
#endif
#endif /* T8WYO_CMESH_FROM_MCELL_FILE_H */
