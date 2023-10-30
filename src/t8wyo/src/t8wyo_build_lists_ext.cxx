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
    int *iface2cell = &face2cell[2*Face->face_index];

    iface2cell[0] = Face->e1;
    iface2cell[1] = Face->e2;
    return 1;
}

static int
t8wyo_facetype_fill(void **face,const void *data){
    t8wyo_face_t *Face = *(t8wyo_face_t **) face;
    int *facetype = (int *) data;
    int *ifacetype = &facetype[4*Face->face_index];

    ifacetype[0] = (int) Face->iface1;
    ifacetype[1] = (int) Face->iface2;
    ifacetype[2] = (int) Face->orientation;
    ifacetype[3] = (int) Face->nvert;
    return 1;
}

static int
t8wyo_facenormal_fill(void **face,const void *data){
    t8wyo_face_t *Face = *(t8wyo_face_t **) face;
    Real *fnorm = (Real *) data;

    Real fac = Face->area;
    fnorm[3*Face->face_index+0] = fac*Face->normal[0];
    fnorm[3*Face->face_index+1] = fac*Face->normal[1];
    fnorm[3*Face->face_index+2] = fac*Face->normal[2];
    return 1;
}

static int faceswap_tet[] = {3,1,2,0};       // swap 0 & 3
static int faceswap_prism[] = {0,1,2,4,3,5}; // swap 3 & 4
static int faceswap_pyramid[] = {0,1,2,3,4}; // TODO: NOT CLEAR WHAT RIGHT ANS

static int swap_face(const int eclass,const int face_id){
   switch (eclass){
       case 4: //T8_ECLASS_TET
            return faceswap_tet[face_id];
            break;
       case 5://T8_ECLASS_PYRAMID:
            return faceswap_pyramid[face_id]; // TODO: fix pyramid swap
            break;
       case 6: //T8_ECLASS_PRISM:
            return faceswap_prism[face_id];
            break;
        default:
            return face_id;
            break;
    }
}

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
                           wyo::memory<Real> &face_norm){
    T8_ASSERT(t8_forest_is_committed(forest));

    const int F = t8_eclass_max_num_faces[cmesh->dimension];
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
    mortar_info.free();
    elem_info.free();
    elem_vol.free();
    face_norm.free();

    ndc4.free();
    ndc5.free();
    ndc6.free();
    ndc8.free();
    xgeom.free();

    /* setup memory pools */
    sc_mempool_t *face_mempool = sc_mempool_new(sizeof(t8wyo_face_full_t));
    sc_list_t *faces = sc_list_new(NULL);

    sc_mempool_t *face_unique_mempool = sc_mempool_new(sizeof(t8wyo_face_t));
    sc_hash_t *faces_unique = sc_hash_new(t8wyo_face_hash,
                                          t8wyo_face_equal,
                                          NULL,NULL);

    /* count number of local elements; allocate elem_info () */
    t8_locidx_t num_local_trees = t8_forest_get_num_local_trees(forest);
    t8_locidx_t elem_count = t8_forest_get_local_num_elements(forest);
    t8_locidx_t ghost_count = t8_forest_get_num_ghosts(forest);

    elem_info.malloc(INFO_ELEM_SIZE*elem_count,-1); // auto-deletes if pre-allocated
    elem_vol.malloc(elem_count+ghost_count,-1.0);

    /* count number of elements of each type*/
    int elem_counts[T8_ECLASS_COUNT] = {0};
    for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
        t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements(forest,itree);
        t8_eclass_t tree_class = t8_forest_get_tree_class(forest,itree);
        t8_eclass_scheme_c *eclass_scheme = t8_forest_get_eclass_scheme(forest,tree_class);

        for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
            t8_element_t *element = t8_forest_get_element_in_tree(forest,itree,ielement);
            const t8_element_shape_t element_shape = eclass_scheme->t8_element_shape(element);

            elem_counts[element_shape]++;
        }
    }
    int ntet   = elem_counts[T8_ECLASS_TET];
    int npyr   = elem_counts[T8_ECLASS_PYRAMID];
    int nprism = elem_counts[T8_ECLASS_PRISM];
    int nhex   = elem_counts[T8_ECLASS_HEX];

    /* allocate element connectivities */
    ndc4.malloc(4*ntet);
    ndc5.malloc(5*npyr);
    ndc6.malloc(6*nprism);
    ndc8.malloc(8*nhex);

    /* reset element counts */
    ntet = npyr = nprism = nhex = 0;

     // make unordered hash sets
    std::unordered_set<Mortar, Mortar::HashFunction, Mortar::KeyEqual> mortars;
    std::unordered_set<Node, Node::HashFunction> nodes;
//    wyo::memory<char> reorder(elem_count,0);
    double vertex_coords[3];
    double coords[3];
    int icorner;
    int iface;
    int tmp = -1000000;

    t8_locidx_t *faces1;
    int8_t *ttf;

    /* count max faces in forest; stash all faces into list */
    int node_count = 0;
    for (t8_locidx_t itree = 0,elem_index = 0; itree < num_local_trees; ++itree) {
        t8_eclass_t tree_class = t8_forest_get_tree_class(forest,itree);
        t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements(forest,itree);
        t8_eclass_scheme_c *eclass_scheme = t8_forest_get_eclass_scheme(forest,tree_class);
        t8_gloidx_t gtreeid = t8_forest_global_tree_id(forest,itree);
        t8_cmesh_t cmesh = t8_forest_get_cmesh(forest);

        /* loop over all local elements in the tree. */
        for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement, ++elem_index) {
            t8_element_t *element = t8_forest_get_element_in_tree(forest,itree,ielement);
            const t8_element_shape_t element_shape = eclass_scheme->t8_element_shape(element);
            const int num_faces = eclass_scheme->t8_element_num_faces(element);
            const int num_corners = t8_eclass_num_vertices[element_shape];

            int *elem_conn = (element_shape == T8_ECLASS_TET)     ? ndc4.ptr():
                             (element_shape == T8_ECLASS_PYRAMID) ? ndc5.ptr():
                             (element_shape == T8_ECLASS_PRISM)   ? ndc6.ptr():
                             (element_shape == T8_ECLASS_HEX)     ? ndc8.ptr():
                                                                    nullptr;

            int &eidx = (element_shape == T8_ECLASS_TET)     ? ntet:
                        (element_shape == T8_ECLASS_PYRAMID) ? npyr:
                        (element_shape == T8_ECLASS_PRISM)   ? nprism:
                        (element_shape == T8_ECLASS_HEX)     ? nhex:
                                                               tmp;
            /* set local element information */
            int *einfo = &elem_info[INFO_ELEM_SIZE*elem_index];
            einfo[ETYPE_IND] = elem_class[element_shape];
            einfo[ELEVL_IND] = eclass_scheme->t8_element_level(element);

            /* initialize element's volume */
            elem_vol[elem_index] = t8_forest_element_volume(forest,itree,element);

            /* fill nodes */
            /* TODO: check for a better way to get correct ordering insert vertices */
            int tree_node_ids[T8_ECLASS_MAX_CORNERS];
            double tree_vertices[T8_ECLASS_MAX_CORNERS * 3];
            for (icorner = 0; icorner < num_corners; icorner++) {
                const size_t num_coords = 1;
                eclass_scheme->t8_element_vertex_reference_coords(element,icorner,vertex_coords);
                t8_geometry_evaluate(cmesh,gtreeid,vertex_coords,num_coords,coords);
                memcpy(&tree_vertices[3*icorner],coords,3*sizeof(double));

                /* try to insert new node */
                auto ret = nodes.insert(Node(node_count,coords));

                /* set connectivity for node */
                /* calculate node id:
                 * if ret.second == true: (new node inserted)
                 *      set to node_count, post-increment node_count
                 * else: (node already existed)
                 *      set to the found node id (i.e., Node node = *ret.first;)
                 */
                const int node_id = (ret.second) ? node_count++:(*ret.first).id;
                tree_node_ids[icorner] = node_id+FBASE;
            }

            // check ordering
//            if(t8_cmesh_tree_vertices_negative_volume(element_shape,tree_vertices,num_corners)){
//                correct_node_ordering(element_shape,tree_node_ids);
//                reorder[elem_index] = 1; // set node connectivity reorder flag
//            }

            // fill element connectivity
            memcpy(&elem_conn[num_corners*eidx],tree_node_ids,num_corners*sizeof(int));

            // update element counter
            eidx++;

            // get the tree2face info for orientation
            t8_locidx_t lctreeid = t8_forest_ltreeid_to_cmesh_ltreeid(forest,itree);
            (void) t8_cmesh_trees_get_tree_ext(cmesh->trees,lctreeid,&faces1,&ttf);

            /* loop over all faces of an element: count total neighboring elements */
            for (iface = 0; iface < num_faces; iface++) {
                t8_element_shape_t face_shape = eclass_scheme->t8_element_face_shape(element,iface);

                /* allocate face and set element id */
                t8wyo_face_full_t *Face = (t8wyo_face_full_t *) sc_mempool_alloc(face_mempool);
                Face->elem = element;
                Face->ltree_id = itree;
                Face->lelem_id = ielement;
                Face->face_number = iface;
                Face->elem_id = elem_index;
                Face->nvert = t8_eclass_num_vertices[face_shape];
                Face->orientation = ttf[iface] / F; // see t8_cmesh_trees.c:L1008

                /* add face to faces list */
                (void) sc_list_append(faces,Face);
            }
        }
    }

    /* allocate and fill xgeom list */
    xgeom.malloc(3*node_count);

    // fill xgeom: must use vertex->id because nodes is an unordered set
    for (auto vertex = nodes.begin(); vertex != nodes.end(); ++vertex) {
        xgeom[3*vertex->id+0] = vertex->x;
        xgeom[3*vertex->id+1] = vertex->y;
        xgeom[3*vertex->id+2] = vertex->z;
    }
    /* free hash set memory */
    nodes.clear();

    /* ============================================== */
    /* make unique list of faces/hanging face mortars */
    /* ============================================== */
    int face_unique_count = 0;
    int bface_count = 0;
    int nmortar = 0;
    while(faces->elem_count > 0){
        t8wyo_face_full_t *Face_full = (t8wyo_face_full_t *) sc_list_pop(faces);
        const t8_locidx_t elem_idx = Face_full->elem_id;
        const int iface = Face_full->face_number;
        const int own_level = elem_info[INFO_ELEM_SIZE*elem_idx+ELEVL_IND];

        // local variables
        t8_eclass_scheme_c *neigh_scheme;
        t8_element_t **face_children;
        t8_element_t **neighs;
        t8_locidx_t *neighids;
        int *dual_faces;
        int num_neighbors;
        int num_face_children;

        /* gather face neighbor element information */
        t8_forest_leaf_face_neighbors(forest,
                                      Face_full->ltree_id,
                                      Face_full->elem,
                                     &neighs,iface,
                                     &dual_faces,
                                     &num_neighbors,
                                     &neighids,
                                     &neigh_scheme,1);

        /* ------------------------- */
        /* build hanging mortar list */
        /* ------------------------- */
        if (num_neighbors == 1 && (neigh_scheme->t8_element_level(neighs[0]) < own_level)) {
            /* hanging face: the neighbor element is the coarser one (lower level number) */

            // note: skip local big neighbor as it will filled in below (num_neighbors>1)
            if (neighids[0] >= elem_count) {
                /* neighbor element is a ghost: new mortar */
                auto ret = mortars.insert(Mortar(neighids[0],dual_faces[0]));
                if(ret.second) nmortar++;

                /* get hanging subface id */
                t8_eclass_scheme_c *ts;
                t8_eclass_scheme_c *boundary_scheme;
                t8_eclass_t boundary_class;
                t8_element_t *face_element;

                t8_tree_t tree;
                int tree_face;
                int isubface;

                /* get pointer to the tree to read its element class */
                tree = t8_forest_get_tree(forest,Face_full->ltree_id);
                ts = t8_forest_get_eclass_scheme(forest,tree->eclass);

                /* Get the scheme associated to the element class of the boundary element. */
                /* Compute the face of elem_tree at which the face connection is. */
                tree_face = ts->t8_element_tree_face(Face_full->elem,iface);

                /* get the eclass scheme for the boundary */
                boundary_class = (t8_eclass_t) t8_eclass_face_types[tree->eclass][tree_face];
                boundary_scheme = t8_forest_get_eclass_scheme(forest,boundary_class);
                boundary_scheme->t8_element_new(1,&face_element);
                ts->t8_element_boundary_face(Face_full->elem,iface,face_element,boundary_scheme);
                isubface = boundary_scheme->t8_element_child_id(face_element);

                /* retrieve mortar (may have existed) */
                auto &mortar = (*ret.first);
                assert(mortar.cnt < MAX_SUBFACES);
                assert(mortar.elemidx_plus[isubface] == -1);

                /* insert hang element info */
                mortar.elemidx_plus[isubface] = elem_idx;
                mortar.iface_plus[isubface] = iface;
                mortar.cnt++;

                /* free memory */
                boundary_scheme->t8_element_destroy(1,&face_element);
            }
        } else
        if (num_neighbors > 1) {
            /* hanging face: this element is the bigger one: try to insert new mortar */
            auto ret = mortars.insert(Mortar(elem_idx,iface));
            if(ret.second) nmortar++;

            /* insert hang element info */
            auto &mortar = (*ret.first); // retrieve mortar (may have existed)
            assert(mortar.cnt == 0);

            for(int ihang = 0; ihang < num_neighbors; ihang++){
                mortar.elemidx_plus[ihang] = neighids[ihang];
                mortar.iface_plus[ihang] = dual_faces[ihang];
            }
            mortar.cnt = MAX_SUBFACES;
        }

        /* ---------------------- */
        /* build unique face list */
        /* ---------------------- */
        if (num_neighbors > 0) {
            /* interior face:
             *  loop of elements on large face and
             *  construct neighbor info for each local small face
             */

            t8_eclass eclass = t8_forest_get_tree_class(forest,Face_full->ltree_id);
            t8_eclass_scheme_c *ts = t8_forest_get_eclass_scheme(forest,eclass);

            if (num_neighbors>1) {
                num_face_children = ts->t8_element_num_face_children(Face_full->elem,iface);

                face_children = T8_ALLOC(t8_element_t *,num_face_children);
                ts->t8_element_new(num_face_children,face_children);
                ts->t8_element_children_at_face(Face_full->elem,iface,face_children,num_face_children,NULL);
            }

            for (int ielem = 0; ielem < num_neighbors; ielem++) {
                t8wyo_face_t *Face = (t8wyo_face_t *) sc_mempool_alloc(face_unique_mempool);

                /* set owner and neighbor element info */
                Face->e1 = Face_full->elem_id + FBASE;
                Face->e2 = neighids[ielem] + FBASE;
                Face->nvert = Face_full->nvert;
                Face->iface1 = iface + FBASE;
                Face->iface2 = dual_faces[ielem] + FBASE;
                Face->orientation = Face_full->orientation;

                /* try to insert the face into the hash */
                if (sc_hash_insert_unique(faces_unique,Face,NULL)) {
                    /* compute geometry info */
                    char full = (num_neighbors == 1);
                    t8_element_t *element = (full) ? Face_full->elem:face_children[ielem];
                    int face_number = (full) ? Face_full->face_number:ts->t8_element_face_child_face(Face_full->elem,iface,ielem);

                    /* face information */
                    t8_forest_element_face_normal(forest,Face_full->ltree_id,element,face_number,Face->normal);
                    Face->area = t8_forest_element_face_area(forest,Face_full->ltree_id,element,face_number);
                    Face->face_index = face_unique_count++;
                } else {
                    /* face already existed: clean up memory */
                    sc_mempool_free(face_unique_mempool,Face);
                }
            }

            /* clean up memory */
            if(num_neighbors>1){
                ts->t8_element_destroy(num_face_children,face_children);
                T8_FREE(face_children);
            }

            neigh_scheme->t8_element_destroy(num_neighbors,neighs);
            T8_FREE(neighs);
            T8_FREE(dual_faces);
            T8_FREE(neighids);
        } else {
            /* boundary face:
             *  compute base tree local information and
             *  retrieve saved boundary face information
             *   -- contains original BC-FACE ID info from input mesh
             */
            t8_eclass_t tree_class = t8_forest_get_tree_class(forest,Face_full->ltree_id);
            t8_locidx_t cmesh_ltreeid = t8_forest_ltreeid_to_cmesh_ltreeid(forest,Face_full->ltree_id);
            t8_eclass_scheme_c *eclass_scheme = t8_forest_get_eclass_scheme(forest,tree_class);
            t8_element_t *element = t8_forest_get_element_in_tree(forest,Face_full->ltree_id,Face_full->lelem_id);

            /* compute original tree face number */
            int tree_face = t8_element_tree_face(eclass_scheme,element,Face_full->face_number);

            t8_locidx_t *face_info = (t8_locidx_t *)
                t8_cmesh_get_attribute(cmesh,t8_get_package_id(),
                                       T8WYO_CMESH_OFFSET_KEY+tree_face, // use original tree face number!
                                       cmesh_ltreeid);

            /* set owner element info */
            t8wyo_face_t *Face = (t8wyo_face_t *) sc_mempool_alloc(face_unique_mempool);
            Face->e1 = Face_full->elem_id + FBASE;
            Face->e2 = face_info[T8WYO_FID_KEY]; // face id in cmesh (FBASE included)

            /* try to insert the face into the hash */
            if (sc_hash_insert_unique(faces_unique,Face,NULL)) {
                /* compute geometry info */
                element = t8_forest_get_element_in_tree(forest,
                                                        Face_full->ltree_id,
                                                        Face_full->lelem_id);
                /* face normal */
                t8_forest_element_face_normal(forest, Face_full->ltree_id,
                                              element,Face_full->face_number,
                                              Face->normal);
                /* face area */
                Face->area = t8_forest_element_face_area(forest, Face_full->ltree_id,
                                                         element,Face_full->face_number);

                /* set face index and increment counter */
                Face->orientation = Face_full->orientation;
                Face->iface1 = Face_full->face_number + FBASE;
                Face->iface2 = -1; // TODO: lookup value from cmesh: patch number
                Face->nvert = face_info[T8WYO_FTYPE_KEY];   // lookup value from cmesh: face type
                Face->face_index = face_unique_count;
                face_unique_count++;
                bface_count++;
            } else {
                /* face already existed: cleanup memory */
                sc_mempool_free(face_unique_mempool,Face);
            }
        }
    }

    /* fill hanging mortar information */
    mortar_info.malloc(10*mortars.size());

    int mcount=0;
    for (auto mortar = mortars.begin(); mortar != mortars.end(); ++mortar,mcount++) {
        int * const minfo = &mortar_info[10*mcount];

        // big element mortar data (1st two entries)
        minfo[0] = mortar->elemidx_minus+FBASE;
        minfo[1] = mortar->iface_minus;
//        if(mortar->elemidx_minus < elem_count && reorder[mortar->elemidx_minus]){
//            const int eclass_retype = elem_info[INFO_ELEM_SIZE*mortar->elemidx_minus+ETYPE_IND];
//            minfo[1] = swap_face(eclass_retype,minfo[1]);
//        }

        // hang element mortar data (entries 2-9)
        int * const mhang = &minfo[2];
        for(int i=0; i<MAX_SUBFACES; i++){
            mhang[2*i] = mortar->elemidx_plus[i]+FBASE;
            mhang[2*i+1] = mortar->iface_plus[i];

            // swap faces
//            if(0 <= mortar->elemidx_plus[i] &&
//                    mortar->elemidx_plus[i] < elem_count &&
//                    reorder[mortar->elemidx_plus[i]]){
//                const int eclass_retype = elem_info[INFO_ELEM_SIZE*(mortar->elemidx_plus[i])+ETYPE_IND];
//                mhang[2*i+1] = swap_face(eclass_retype,mhang[2*i+1]);
//            }
        }
    }
    mortars.clear();

    /* free non-unique face mempool/list */
    sc_mempool_destroy(face_mempool);
    sc_list_destroy(faces);

    /* allocate and fill face2cell */
    face2cell.malloc(2*face_unique_count,-1); // two ints per face
    faces_unique->user_data = face2cell.ptr(); // set user_data for loop
    sc_hash_foreach(faces_unique,t8wyo_face2cell_fill);

    /* allocate and fill facetype */
    facetype.malloc(4*face_unique_count,-1);
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