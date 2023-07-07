/**
 * \file    t8wyo_solver.hxx
 * \author  akirby
 */

#ifndef T8WYO_SOLVER_H
#define T8WYO_SOLVER_H

/* header files */
#include "precision_types.h"
#include "memory.hxx"

/* system header files */
#ifndef DOXYGEN_IGNORE
#  include <mpi.h>
#  include <stdio.h>
#  include <stdlib.h>
#  include <stddef.h>
#  include <sys/stat.h>
#  include <memory>
#endif



#ifdef __cplusplus
extern "C" {
#endif

#define T8WYO_VERSION 101   /**< @brief t8wyo_solver.hxx version number for compatibility */

#define BUFF_SIZE 1024      /**< @brief string buffer size */
#define MAX_LEVELS 29       /**< @brief number of maximum levels allowed */

/* offsets for \elem_info: NOTES end of file */
#define ETYPE_IND 0     /**< @brief Element type     offset into \elem_info */
#define ELEVL_IND 1     /**< @brief Element level    offset into \elem_info */
#define INFO_ELEM_SIZE 2

/* offsets for \face_info: NOTES end of file */
#define FTYPE_IND  0    /**< @brief            face type                   offset into \face_info */
#define FFELEM_IND 1    /**< @brief       this quad index in \elem_info    offset into \face_info */
#define FFSIDE_IND 2    /**< @brief       this quad side                   offset into \face_info */
#define F1ELEM_IND 3    /**< @brief neighbor 1 quad index in \elem_info    offset into \face_info */
#define F1SIDE_IND 4    /**< @brief neighbor 1 quad side                   offset into \face_info */
#define F2ELEM_IND 5    /**< @brief neighbor 2 quad index in \elem_info    offset into \face_info */
#define F2SIDE_IND 6    /**< @brief neighbor 2 quad side                   offset into \face_info */
#define F3ELEM_IND 7    /**< @brief neighbor 3 quad index in \elem_info    offset into \face_info */
#define F3SIDE_IND 8    /**< @brief neighbor 3 quad side                   offset into \face_info */
#define F4ELEM_IND 9    /**< @brief neighbor 4 quad index in \elem_info    offset into \face_info */
#define F4SIDE_IND 10   /**< @brief neighbor 4 quad side                   offset into \face_info */
#define FBCID_IND  3    /**< @brief physical ID (only for boundary sides!) offset into \face_info */
#define INFO_FACE_SIZE 11

#define BC_TYPE   1
#define FULL_TYPE 2
#define HANG_TYPE 3

/* offset for geometry \face_nodes_info */
#define EFENT_IND  0    /**< @brief Element face entity tag */
#define EFTAG_IND  1    /**< @brief Element face tag */
#define EFNDS_IND  2    /**< @brief Element face node list */

/* patch id or tag not assigned */
#define NO_PATCH   -1
#define NO_TAG     -1

/* face indicator bit masks */
#define xlo_mask 0b00000001 /**< xlo face intersection flag */
#define xhi_mask 0b00000010 /**< xhi face intersection flag */
#define ylo_mask 0b00000100 /**< ylo face intersection flag */
#define yhi_mask 0b00001000 /**< yhi face intersection flag */
#define zlo_mask 0b00010000 /**< zlo face intersection flag */
#define zhi_mask 0b00100000 /**< zhi face intersection flag */

#define MESH_MODE_MCELL   0
#define MESH_MODE_GMSH    1
#define MESH_MODE_CART    2
#define MESH_MODE_COUNT   3
#define MESH_MODE_INVALID 4

#define TIMER(time,...)      \
  do {                       \
      double t_start,t_end;  \
      t_start = MPI_Wtime(); \
      __VA_ARGS__            \
      t_end = MPI_Wtime();   \
      time = t_end-t_start;  \
  } while(0)

/**
 * @brief Contains all grid data.
 */
typedef struct {
    int nlevels;             /**< Number of AMR levels */
    int min_level;           /**< Minumum level of grid refinement */
    int max_level;           /**< Maximum level of grid refinement */
    int max_level_pmax;      /**< Maximum level allowing pdegree_max */
    int max_level_global;    /**< Maximum level override if hbox level > max_level */
    int nelem[3];            /**< Level-0 number of elements in mesh */
    int periodic[3];         /**< Periodic boundary condition indicator */
    int construct_grid;      /**< Initialize grid by building up levels */
    int total_quads;         /**< Total Mesh Quadrants */
    unsigned int total_dofs; /**< Total Degrees of Freedom */
    Real xlo[3];             /**< Lower coordinates of the grid domain */
    Real xhi[3];             /**< Upper coordinates of the grid domain */
    Real min_dx[3];          /**< Size of the finest element */
    Real max_dx[3];          /**< Size of the coarsest element */
    Real domain_volume;      /**< Grid domain volume */
}
grid_t;

/**
 * @brief Contains all MPI communication data.
 */
typedef struct {
    int nsend;                 /**< MPI send count */
    int nrecv;                 /**< MPI receive count */
    int send_info_size;        /**< Number of entries in the mpi send */
    int recv_info_size;        /**< Number of entries in the mpi recv */

    int *data_counts;          /**< Number of data counts on each rank */
    Real *sendbuf;             /**< MPI data communication buffer */
    MPI_Request *request;      /**< MPI request buffer */
    wyo::memory<int> send_info; /**< MPI send info: [0] = send quad id */
    wyo::memory<int> recv_info; /**< MPI receive info */
}
mpi_t;

/**
 * @brief Contains all external solver data.
 */
class external_t {
  public:
    Real dof_nopad_total;       /**< Number of total real degrees of freedom */
    int dof_nopad;              /**< Number of actual real degrees of freedom */
    int dof;                    /**< Number of real degrees of freedom */
    int nquad;                  /**< Number of real quads in array */
    int nface_full;             /**< Number of full faces in array */
    int nface_hang;             /**< Number of hang faces in array */
    int nface_bdry;             /**< Number of boundary condition faces in array */
    int nface_full_ghost;       /**< Number of ghost full faces in array */
    int nface_hang_ghost;       /**< Number of ghost hang faces in array */
    int nface;                  /**< Number of faces in array */
    int nedge;                  /**< Number of edges in array */
    int nnode;                  /**< Number of nodes in array */
    int nghost;                 /**< Number of ghost quads in array */
    int quad_info_size;         /**< Number of entries per quad in quad_info */
    int face_info_size;         /**< Number of entries per face in face_info */
    int edge_info_size;         /**< Number of entries per edge in egde_info */
    int node_info_size;         /**< Number of entries per node in node_info */
    int elem2face_info_size;    /**< Number of entries per quad in q2face_info */
    int elem2node_info_size;    /**< Number of entries per quad in q2node_info */

    wyo::memory<int> elem_info;  /**< Element information: notes at bottom of file */
    wyo::memory<int> face_info;  /**< Face information: notes at bottom of file */
    wyo::memory<int> node_info;  /**< Node information: notes at bottom of file */
    wyo::memory<int> elem2face_info;       /**< Element to face information */
    wyo::memory<int> elem2node_info;       /**< Element to node information */
    wyo::memory<int> elem_data_sizes;      /**< Element data size used */
    wyo::memory<int> full_face_info;       /**< Full face indices into face_info */
    wyo::memory<int> hang_face_info;       /**< Hang face indices into face_info */
    wyo::memory<int> bdry_face_info;       /**< Boundary condition face indices into face_info */
    wyo::memory<int> full_face_ghost_info; /**< Ghost full face indices into face_info */
    wyo::memory<int> hang_face_ghost_info; /**< Ghost hang face indices into face_info */

    wyo::memory<Real> soln;          /**< External solution buffer */
    wyo::memory<Real> soln_interp;   /**< External interpolated solution buffer */
    wyo::memory<Real> dsoln;         /**< External solution correction buffer */
    wyo::memory<Real> resid;         /**< External total residual buffer */
    wyo::memory<Real> src;           /**< External source buffer */
    wyo::memory<Real> R0;            /**< External Acoustic storage: time (n) */
    wyo::memory<Real> R1;            /**< External Acoustic storage: time (n+1) */
    wyo::memory<Real> Q0;            /**< External BDF storage: time (n) */
    wyo::memory<Real> Q1;            /**< External BDF storage: time (n-1) */
    wyo::memory<Real> Q2;            /**< External BDF storage: time (n-2) */

    wyo::memory<Real> QBC;           /**< External solution BC buffer */
    wyo::memory<Real> elem_geom;     /**< Grid element coordinates */
    wyo::memory<Real> elem_volume;   /**< Grid element volume */

    wyo::memory<Real> elem_vol2surf; /**< Grid element volume to surface area ratio */
    wyo::memory<Real> vgeo_co;       /**< Volume collocation geometry Jacobian factors */
    wyo::memory<Real> vgeo;          /**< Volume  geometry Jacobian factors */
    wyo::memory<Real> sgeo;          /**< Surface geometry Jacobian factors */
    wyo::memory<Real> svgeo;         /**< Surface volume geometry Jacobian factors */
    wyo::memory<Real> vjac;          /**< Volume  Jacobian determinants */
    wyo::memory<Real> sjac;          /**< Surface Jacobian determinants */
    wyo::memory<Real> ivjac;         /**< Inverse volume  Jacobian determinants */
    wyo::memory<int> vgeo_co_index;  /**< Volume collocation geometry Jacobian factors index pointer */
    wyo::memory<int> vgeo_index;     /**< Volume  geometry Jacobian factors index pointer */
    wyo::memory<int> sgeo_index;     /**< Surface geometry Jacobian factors index pointer */
    wyo::memory<int> svgeo_index;    /**< Surface volume geometry Jacobian factors index pointer */
    wyo::memory<int> vjac_index;     /**< Volume  Jacobian determinants index pointer */
    wyo::memory<int> sjac_index;     /**< Surface Jacobian determinants index pointer */

    mpi_t d_mpi;                     /**< MPI communication data */

    /* constructors */
    external_t()=default;
   ~external_t()=default;
};

/**
 * @brief Contains boundary condition data.
 */
typedef struct {
    int  npatch;                          /**< Number of the physical patches defined */
    Uint nentity_face;                    /**< Number of boundary geometry entities (edges/faces) */
    Uint nentity_vol;                     /**< Number of volume geometry entities (faces/volumes) */
    Uint nelem_face;                      /**< Number of element faces in the grid_file */
    char wake3d_wbc[BUFF_SIZE];           /**< Char containing the patches name of wall boundaries */
    char wake3d_obc[BUFF_SIZE];           /**< Char containing the patches name of outern boundaries */
    char **patch_name;                    /**< List of physical patches' names (size npatch)*/
    wyo::memory<int> patch_dim;            /**< List of physical patches's dimension 2: 1: 1D, 2D, 3: 3D (size npatch) */
    wyo::memory<int> entity_face_patchtag; /**< List of boundary geometry entities' patch tag (edges/faces) */
    wyo::memory<int> entity_vol_patchtag;  /**< List of volume geometry entities' patch tag (faces/volumes) */
    wyo::memory<Uint> entity_face_tag;     /**< List of boundary geometry entities' tag (edges/faces) */
    wyo::memory<Uint> entity_vol_tag;      /**< List of volume geometry entities' tag (faces/volumes) */
    wyo::memory<Uint> face_nodes_info;     /**< List of element face (edges/faces) info */
}
geometry_t;

/**
 * @brief Contains gridfile line number of elements (used to retrieve HiO info faster)
 */
typedef struct {
    char gridfile_name[BUFF_SIZE]; /**< Unstructured grid file name */
    int binary;                    /**< Binary file (1: yes, 0: no) */
    int swap;                      /**< Binary swap needed (1:yes, 0: no) */
    int file_format;               /**< File format (41 or 2.0) */
    size_t entities_bytes;         /**< number of bytes for gmsh $ENTITIES section (only for binary) */
    size_t nodes_bytes;            /**< number of bytes for gmsh $NODES section (only for binary) */
    int nvertices;                 /**< number of vertices per high-order volume element (only for binary) */
    wyo::memory<Uint> elem_tag;     /**< Element volume tag */
    wyo::memory<Uint> p4est_part;   /**< p4est partition mapping */
}
gridfile_t;

/**
 * @brief Contains all context data of T8WYO.
 *        This contains all the other data structures.
 */
class ctx_t {
  public:
    int t8wyo_version;         /**< Header file version flag */
    int rank;                   /**< MPI rank */
    int nranks;                 /**< Number of ranks in this communicator */
    MPI_Comm comm;              /**< MPI communicator */
    FILE *log_io;               /**< Output file stream */
    int log_info;               /**< Logging threshold flag:
                                 *      *Values:
                                 *      * 0     log everything
                                 *      * 1     prefix file and line number
                                 *      * 2     information on the internal state
                                 *      * 3     information on conditions, decisions
                                 *      * 4     main information a function is doing
                                 *      * 5     important consistency/performance information
                                 *      * 6     few lines for a major API function
                                 *      * 7     logs a few lines max per program
                                 *      * 8     logs errors only
                                 *      * 9     disable logging
                                 */
    grid_t d_grid;              /**< Grid data */
    gridfile_t d_gridfile;      /**< Gridfile data */

    /* constructors */
    ctx_t()=default;
   ~ctx_t()=default;
};

/**
 * @brief Contains all data of T8WYO.
 */
class t8wyo_t {
  public:
    ctx_t ctx;                        /**< Context data */
    wyo::memory<external_t> external; /**< External solver data */

    /* constructors */
    t8wyo_t()=default;
   ~t8wyo_t()=default;
};

#ifdef __cplusplus
}
#endif
#endif /* T8WYO_SOLVER_H */

/* ========================================================================== *
 * NOTES: Data Structures                                                     *
 * -------------------------------------------------------------------------- *
 * int *elem_info: local information on element (4 fields per element)        *
 *      elem_info[0] = element type                                           *
 *      elem_info[1] = solution index                                         *
 *      elem_info[2] = geometry index                                         *
 *      elem_info[3] = grid level                                             *
 * -------------------------------------------------------------------------- *
 * int *face_info: Face information on element (7/11 fields per face)         *
 *      face_info[0] = face_type                                              *
 *                    *face_type == 1: boundary                               *
 *                    *face_type == 2: full-full                              *
 *                    *face_type == 3: hanging: check side (face_info[2])     *
 *      face_info[1] = this elem index in elem_info                           *
 *      face_info[2] = side of element                                        *
 *                    *side == 0: xlo face                                    *
 *                      | Left  | Right |                                     *
 *                      | ----: | :---- |                                     *
 *                      | hang  | full  |                                     *
 *                                                                            *
 *                    *side == 1: xhi face                                    *
 *                      | Left  | Right |                                     *
 *                      | ----: | :---- |                                     *
 *                      | full  | hang  |                                     *
 *                                                                            *
 *                    *side == 2: ylo face                                    *
 *                      | Left  | Right |                                     *
 *                      | ----: | :---- |                                     *
 *                      | hang  | full  |                                     *
 *                                                                            *
 *                    *side == 3: yhi face                                    *
 *                      | Left  | Right |                                     *
 *                      | ----: | :---- |                                     *
 *                      | full  | hang  |                                     *
 *                                                                            *
 *                    *side == 4: zlo face                                    *
 *                      | Left  | Right |                                     *
 *                      | ----: | :---- |                                     *
 *                      | hang  | full  |                                     *
 *                                                                            *
 *                    *side == 5: zhi face                                    *
 *                      | Left  | Right |                                     *
 *                      | ----: | :---- |                                     *
 *                      | full  | hang  |                                     *
 *                                                                            *
 *      face_info[3] = neighbor 1 quad index in quad_info (if present)        *
 *      face_info[4] = neighbor 1 side (if present)                           *
 *      face_info[5] = neighbor 2 quad index in quad_info (if present)        *
 *      face_info[6] = neighbor 2 side (if present)                           *
 *      face_info[7] = neighbor 3 quad index in quad_info (if present)        *
 *      face_info[8] = neighbor 3 side (if present)                           *
 *      face_info[9] = neighbor 4 quad index in quad_info (if present)        *
 *      face_info[10]= neighbor 4 side (if present)                           *
 * -------------------------------------------------------------------------- *
 * int *edge_info: Edge information on element (? fields per face)            *
 *      edge_info[0] = {NOT SET}                                              *
 * -------------------------------------------------------------------------- *
 * int *node_info: Node information on element (9/17 fields per face)         *
 *      node_info[0] = number of quadrants touching node                      *
 *      node_info[1] = quad 0 index in quad_info                              *
 *      node_info[2] = quad 0 corner side                                     *
 *      node_info[3] = quad 1 index in quad_info (if present)                 *
 *      node_info[4] = quad 1 corner side (if present)                        *
 *      node_info[5] = quad 2 index in quad_info (if present)                 *
 *      node_info[6] = quad 2 corner side (if present)                        *
 *      node_info[7] = quad 3 index in quad_info (if present)                 *
 *      node_info[8] = quad 3 corner side (if present)                        *
 *      ----- 3D Only -----                                                   *
 *      node_info[9]  = quad 4 index in quad_info (if present)                *
 *      node_info[10] = quad 4 corner side (if present)                       *
 *      node_info[11] = quad 5 index in quad_info (if present)                *
 *      node_info[12] = quad 5 corner side (if present)                       *
 *      node_info[13] = quad 6 index in quad_info (if present)                *
 *      node_info[14] = quad 6 corner side (if present)                       *
 *      node_info[15] = quad 7 index in quad_info (if present)                *
 *      node_info[16] = quad 7 corner side (if present)                       *
 * ========================================================================== */