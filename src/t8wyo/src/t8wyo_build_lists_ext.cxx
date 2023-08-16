/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* header files */
#include "t8wyo_build_lists_ext.h"

/* Hash a face: value is the sum of its element indices */
unsigned
t8wyo_face_hash(const void *face, const void *data){
    t8wyo_face_t *Face = (t8wyo_face_t *) face;
    return (unsigned) (abs(Face->e1) + abs(Face->e2));
}

/* two faces are considered equal if they have the same elements on each side */
static int
t8wyo_face_equal(const void *facea,
                 const void *faceb,
                 const void *data){
    t8wyo_face_t *Face_a = (t8wyo_face_t *) facea;
    t8wyo_face_t *Face_b = (t8wyo_face_t *) faceb;

    /* if both have different number of vertices they can't be equal */
    if(Face_a->e2 >= 0 && Face_b->e2 >= 0){
        return ((Face_a->e1 == Face_b->e1) && (Face_a->e2 == Face_b->e2)) ||
               ((Face_a->e1 == Face_b->e2) && (Face_a->e2 == Face_b->e1));
    } else {
        /* boundary face: uniquely assembled during full grid traversal */
        return 0;
    }
}

static int
t8wyo_face2cell_fill(void **face,const void *data){
    t8wyo_face_t *Face = *(t8wyo_face_t **) face;
    int *face2cell = (int *) data;

    face2cell[2*Face->face_index+0] = Face->e1;
    face2cell[2*Face->face_index+1] = Face->e2;
    return 1;
}

static int
t8wyo_facetype_fill(void **face,const void *data){
    t8wyo_face_t *Face = *(t8wyo_face_t **) face;
    int *face_info = (int *) data;

    face_info[Face->face_index] = (int) Face->nvert;
    return 1;
}

static int
t8wyo_facenormal_fill(void **face,const void *data){
    t8wyo_face_t *Face = *(t8wyo_face_t **) face;
    Real *fnorm = (Real *) data;

    fnorm[3*Face->face_index+0] = -Face->normal[0]*Face->area;
    fnorm[3*Face->face_index+1] = -Face->normal[1]*Face->area;
    fnorm[3*Face->face_index+2] = -Face->normal[2]*Face->area;
    return 1;
}

void t8wyo_build_lists_ext(t8_cmesh_t cmesh,
                           t8_forest_t forest,
                           wyo::memory<int> &face2cell,
                           wyo::memory<int> &facetype,
                           wyo::memory<int> &elem_info,
                           wyo::memory<Real> &elem_vol,
                           wyo::memory<Real> &face_norm){
    T8_ASSERT(t8_forest_is_committed(forest));
    static int elem_class[] = { 0, // T8_ECLASS_ZERO,T8_ECLASS_VERTEX
                                1, // T8_ECLASS_LINE,
                                2, // T8_ECLASS_QUAD,
                                3, // T8_ECLASS_TRIANGLE,
                                8, // T8_ECLASS_HEX,
                                4, // T8_ECLASS_TET,
                                6, // T8_ECLASS_PRISM,
                                5, // T8_ECLASS_PYRAMID,
                               -1, // T8_ECLASS_COUNT
                               -1};// T8_ECLASS_INVALID

    /* clean up memory to avoid excessive usage */
    face2cell.free();
    facetype.free();
    elem_info.free();
    elem_vol.free();
    face_norm.free();

    /* setup memory pools */
    sc_mempool_t *face_mempool = sc_mempool_new(sizeof(t8wyo_face_full_t));
    sc_list_t *faces = sc_list_new(NULL);

    sc_mempool_t *face_unique_mempool = sc_mempool_new(sizeof(t8wyo_face_t));
    sc_hash_t *faces_unique = sc_hash_new(t8wyo_face_hash,
                                          t8wyo_face_equal,
                                          NULL,NULL);

    /* number of trees that have elements of this process */
    t8_locidx_t num_local_trees = t8_forest_get_num_local_trees(forest);

    t8_locidx_t ghost_count = 0;
    t8_locidx_t elem_count = 0;
    t8_locidx_t face_count = 0;

    /* count number of local elements; allocate elem_info () */
    elem_count = t8_forest_get_local_num_elements(forest);
    ghost_count = t8_forest_get_num_ghosts(forest);
    elem_info.malloc(INFO_ELEM_SIZE*elem_count,-1); // auto-deletes if pre-allocated
    elem_vol.malloc(elem_count+ghost_count,-1.0);

    /* count max faces in forest; stash all faces into list */
    for (t8_locidx_t itree = 0,elem_index = 0; itree < num_local_trees; ++itree) {
        t8_eclass_t tree_class = t8_forest_get_tree_class(forest,itree);
        t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements(forest,itree);
        t8_eclass_scheme_c *eclass_scheme = t8_forest_get_eclass_scheme(forest,tree_class);

        /* loop over all local elements in the tree. */
        for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement, ++elem_index) {
            t8_element_t *element = t8_forest_get_element_in_tree(forest,itree,ielement);
            int num_faces = eclass_scheme->t8_element_num_faces(element);

            /* set local element information */
            int *einfo = &elem_info[INFO_ELEM_SIZE*elem_index];
            einfo[ETYPE_IND] = elem_class[tree_class];
            einfo[ELEVL_IND] = eclass_scheme->t8_element_level(element);

            /* initialize element's volume */
            elem_vol[elem_index] = t8_forest_element_volume(forest,itree,element);

            /* loop over all faces of an element: count total neighboring elements */
            for (int iface = 0; iface < num_faces; iface++) {

                /* allocate face and set element id */
                t8wyo_face_full_t *Face = (t8wyo_face_full_t *) sc_mempool_alloc(face_mempool);
                Face->tree_id = itree;
                Face->lelem_id = ielement;
                Face->elem_id = elem_index;
                Face->face_number = iface;

                /* collect neighbors info at current face */
                t8wyo_forest_leaf_face_neighbors(forest,itree,element,iface,
                                                &Face->num_neighbors, /*< Number of neighbors for each face */
                                                 Face->neighids,      /*< Indices of the neighbor elements */
                                                 1);

                /* update small face count */
                face_count += Face->num_neighbors;

                /* add face to faces list */
                (void) sc_list_append(faces,Face);
            }
        }
    }

    /* make unique list of faces */
    int face_unique_count = 0;
    int bface_count = 0;
    while(faces->elem_count > 0){
        t8wyo_face_full_t *Face_full = (t8wyo_face_full_t *) sc_list_pop(faces);

        if(Face_full->num_neighbors > 0){
            /* interior face:
             *  loop of elements on large face and
             *  construct neighbor info for each local small face
             */
            for (int ielem = 0; ielem < Face_full->num_neighbors; ielem++) {
                t8wyo_face_t *Face = (t8wyo_face_t *) sc_mempool_alloc(face_unique_mempool);

                /* set owner and neighbor element info */
                Face->e1 = Face_full->elem_id + FBASE;
                Face->e2 = Face_full->neighids[ielem] + FBASE;

                /* try to insert the face into the hash */
                if(sc_hash_insert_unique(faces_unique,Face,NULL)){
                    /* compute geometry info */
                    t8_element_t *element = t8_forest_get_element_in_tree(forest,
                                                                          Face_full->tree_id,
                                                                          Face_full->lelem_id);
                    /* face normal */
                    t8_forest_element_face_normal(forest, Face_full->tree_id,
                                                  element,Face_full->face_number,
                                                  Face->normal);
                    /* face area */
                    Face->area = t8_forest_element_face_area(forest, Face_full->tree_id,
                                                             element,Face_full->face_number);
                    Face->face_index = face_unique_count;
                    face_unique_count++;
                } else {
                    /* face already existed: clean up memory */
                    sc_mempool_free(face_unique_mempool,Face);
                }
            }
        } else {
            /* boundary face:
             *  compute base tree local information and
             *  retrieve saved boundary face information
             *   -- contains original BC-FACE ID info from input mesh
             */
            t8_eclass_t tree_class = t8_forest_get_tree_class(forest,Face_full->tree_id);
            t8_locidx_t cmesh_ltreeid = t8_forest_ltreeid_to_cmesh_ltreeid(forest,Face_full->tree_id);
            t8_eclass_scheme_c *eclass_scheme = t8_forest_get_eclass_scheme(forest,tree_class);
            t8_element_t *element = t8_forest_get_element_in_tree(forest,Face_full->tree_id,Face_full->lelem_id);

            /* compute original tree face number */
            int tree_face = t8_element_tree_face(eclass_scheme,element,Face_full->face_number);

            t8_locidx_t *face_info = (t8_locidx_t *)
                t8_cmesh_get_attribute(cmesh,t8_get_package_id(),
                                       T8WYO_CMESH_OFFSET_KEY+tree_face, // use original tree face number!
                                       cmesh_ltreeid);

            /* set owner element info */
            t8wyo_face_t *Face = (t8wyo_face_t *) sc_mempool_alloc(face_unique_mempool);
            Face->e1 = Face_full->elem_id + FBASE;
            Face->e2 = face_info[T8WYO_FID_KEY]; // face id in cmesh (has FBASE)

            /* try to insert the face into the hash */
            if(sc_hash_insert_unique(faces_unique,Face,NULL)) {
                /* compute geometry info */
                element = t8_forest_get_element_in_tree(forest,
                                                        Face_full->tree_id,
                                                        Face_full->lelem_id);
                /* face normal */
                t8_forest_element_face_normal(forest, Face_full->tree_id,
                                              element,Face_full->face_number,
                                              Face->normal);
                /* face area */
                Face->area = t8_forest_element_face_area(forest, Face_full->tree_id,
                                                         element,Face_full->face_number);

                /* set face index and increment counter */
                Face->nvert = face_info[T8WYO_FTYPE_KEY]; // lookup value from cmesh
                Face->face_index = face_unique_count;
                face_unique_count++;
                bface_count++;
            } else {
                /* face already existed: cleanup memory */
                sc_mempool_free(face_unique_mempool,Face);
            }
        }
    }

    //printf("NUMBER OF ELEMENTS: %d\n",elem_count);
    //printf("NUMBER OF MAX FACES %d, %d unique, %d boundary\n",
    //        face_count,face_unique_count,bface_count);

    /* free non-unique face mempool/list */
    sc_mempool_destroy(face_mempool);
    sc_list_destroy(faces);

    /* allocate and fill face2cell */
    face2cell.malloc(2*face_unique_count,-1);
    faces_unique->user_data = face2cell.ptr(); // set user_data for loop
    sc_hash_foreach(faces_unique,t8wyo_face2cell_fill);

    /* allocate and fill facetype */
    facetype.malloc(face_unique_count,-1);
    faces_unique->user_data = facetype.ptr();
    sc_hash_foreach(faces_unique,t8wyo_facetype_fill);

    /* allocate and fill face normals (non-normalized) */
    face_norm.malloc(3*face_unique_count,-1.0);
    faces_unique->user_data = face_norm.ptr(); // set user_data for loop
    sc_hash_foreach(faces_unique,t8wyo_facenormal_fill);

    /* free unique face mempool/hash */
    sc_mempool_destroy(face_unique_mempool);
    sc_hash_destroy(faces_unique);
}