/**
 * \file    t8wyo_cmesh_from_mcell.cxx
 * \author  akirby
 *
 * \brief   Construct cmesh from mcell data structures.
 */

/* header files */
#include "t8wyo_cmesh_from_mcell.h"

/* translate the mcell file vertex number to the t8code vertex number.
 * See also (1) https://github.com/DLR-AMR/t8code/wiki/Build-Cmesh
 *          (2) https://server.scientific-sims.com/cfdlab/scientific-sims/solver-file-formats.html#mcell-format-description
 */
const int t8wyo_mcell_tree_vertex_to_t8_vertex_num[T8_ECLASS_COUNT][8]
  = {
  {0},                          /* VERTEX */
  {0, 1},                       /* LINE */
  {0, 1, 3, 2},                 /* QUAD */
  {0, 1, 2},                    /* TRIANGLE */
  {0, 1, 3, 2, 4, 5, 7, 6},     /* HEX */
  {1, 2, 3, 0},                 /* TET */
  {0, 1, 2, 3, 4, 5, 6},        /* PRISM */
  {0, 1, 3, 2, 4}               /* PYRAMID */
};

/* translate the t8code vertex number to the .msh file vertex number.
 * See also (1) https://github.com/DLR-AMR/t8code/wiki/Build-Cmesh
 *          (2) https://server.scientific-sims.com/cfdlab/scientific-sims/solver-file-formats.html#mcell-format-description
 */
const int t8wyo_vertex_to_mcell_vertex_num[T8_ECLASS_COUNT][8]
  = {
  {0},                          /* VERTEX */
  {0, 1},                       /* LINE */
  {0, 1, 3, 2},                 /* QUAD */
  {0, 1, 2},                    /* TRIANGLE */
  {0, 1, 3, 2, 4, 5, 7, 6},     /* HEX */
  {1, 2, 3, 0},                 /* TET */
  {0, 1, 2, 3, 4, 5, 6},        /* PRISM */
  {0, 1, 3, 2, 4}               /* PYRAMID */
};

/* TODO: if partitioned then only add the needed face-connections to join faces
 *       maybe also only trees and ghosts to classes.
 *       Specifying all face-connections makes commit algorithm slow! */

/* Return the hash value of a node.
 * \param [in]  node        The node whose hash value should be computed.
 * \param [in]  num_nodes   A pointer to a locidx_t storing the total number of nodes.
 * \return   The hash value for a node. This is its index modulo the number of nodes.
 */
static unsigned
t8wyo_mcell_node_hash(const void *node,const void *num_nodes){
    t8wyo_mcell_node_t *Node;
    t8_locidx_t Num_nodes;

    T8_ASSERT(node != NULL);
    T8_ASSERT(num_nodes != NULL);

    /* the data parameter stores the total number of nodes */
    Num_nodes = *(t8_locidx_t *) num_nodes;

    /* the node parameter stores a node structure */
    Node = (t8wyo_mcell_node_t *) node;

    /* the hash value of the node is its index modulo the number of nodes */
    return Node->index % Num_nodes;
}

/* Returns true if two given nodes are the same.
 * False otherwise.
 * Two nodes are considered equal if their indices are the same.
 * u_data is not needed.
 */
static int
t8wyo_mcell_node_compare (const void *node_a,
                               const void *node_b,
                               const void *u_data){
    t8wyo_mcell_node_t *Node_a = (t8wyo_mcell_node_t *) node_a;
    t8wyo_mcell_node_t *Node_b = (t8wyo_mcell_node_t *) node_b;
    return Node_a->index == Node_b->index;
}

/* create node hash table from mcell data */
static sc_hash_t *
t8wyo_mcell_nodes(const mcell_t *mcell,sc_mempool_t **node_mempool){
    t8wyo_mcell_node_t *Node;
    t8_locidx_t n;
    int retval;

    /* create mempool for nodes */
    *node_mempool = sc_mempool_new(sizeof(t8wyo_mcell_node_t));

    /* create hash table */
    sc_hash_t *node_table = sc_hash_new(t8wyo_mcell_node_hash,
                                        t8wyo_mcell_node_compare,
                                        (void *) &mcell->nnode,
                                        NULL);

    /* set each node and add it to the hash table */
    for (n = 0; n < mcell->nnode; n++) {
        Node = (t8wyo_mcell_node_t *) sc_mempool_alloc(*node_mempool);

        /* fill node entries */
        Node->coordinates[0] = mcell->xgeom[3*n+0];
        Node->coordinates[1] = mcell->xgeom[3*n+1];
        Node->coordinates[2] = mcell->xgeom[3*n+2];
        Node->index = n;

        /* insert the node in the hash table */
        retval = sc_hash_insert_unique(node_table,Node,NULL);
    }

    printf("Successfully constructed Nodes.\n");
    return node_table;
}

static void
correct_neg_volume(t8_eclass eclass,double *tree_vertices,int tree_id){
   /* The volume described is negative-- change vertices
    *   tets: switch 0 and 3
    *   prisms: switch 0 and 3, 1 and 4, 2 and 5
    *   pyramids: switch 0 and 4
    *   hexs: switch 0 and 4, 1 and 5, 2 and 6, 3 and 7
    */
    double temp;
    int i;
    int iswitch;
    int num_switches = 0;
    int switch_indices[4] = { 0 };
    T8_ASSERT(t8_eclass_to_dimension[eclass] == 3);

    printf("Correcting negative volume of tree %i\n",tree_id);
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
    T8_ASSERT(!t8_cmesh_tree_vertices_negative_volume(eclass, tree_vertices, num_nodes));
}

static void
fill_element_class(t8_cmesh_t cmesh,sc_hash_t *vertices,sc_array_t **vertex_indices,
                   t8_eclass_t eclass,t8_locidx_t nelem,
                   int *ndc,t8_gloidx_t &tree_count){

    t8wyo_mcell_node_parametric_t Node;
    t8wyo_mcell_node_parametric_t **found_node;

    double tree_vertices[T8_ECLASS_MAX_CORNERS * 3];
    long node_indices[T8_ECLASS_MAX_CORNERS];
    long *stored_indices;
    t8_locidx_t tree_loop;
    int t8_vertex_num;
    int i;

    /* add element to cmesh and its nodes */
    for (tree_loop = 0; tree_loop < nelem; tree_loop++) {
        t8_cmesh_set_tree_class(cmesh,tree_count,eclass);

        /* set elements nodes */
        int num_nodes = t8_eclass_num_vertices[eclass];
        for (i = 0; i < num_nodes; i++) {
            Node.index = node_indices[i] = ndc[num_nodes*tree_loop + i] - FBASE; // subtract F90-BASE
            sc_hash_lookup(vertices,(void *) &Node,(void ***) &found_node);

            /* add node coordinates to the tree vertices (NOTE T8CODE ORDERING) */
            t8_vertex_num = t8wyo_mcell_tree_vertex_to_t8_vertex_num[eclass][i];
            tree_vertices[3*t8_vertex_num + 0] = (*found_node)->coordinates[0];
            tree_vertices[3*t8_vertex_num + 1] = (*found_node)->coordinates[1];
            tree_vertices[3*t8_vertex_num + 2] = (*found_node)->coordinates[2];
        }

        /* detect and correct negative volumes */
        if (t8_cmesh_tree_vertices_negative_volume(eclass,tree_vertices,num_nodes)) {
            correct_neg_volume(eclass,tree_vertices,tree_count);
        }

        /* set the vertices of this tree */
        t8_cmesh_set_tree_vertices(cmesh,tree_count,tree_vertices,num_nodes);

        /* store the vertex indices of tree */
        if (vertex_indices != NULL) {
            /* allocate memory for the indices */
            stored_indices = T8_ALLOC(long,t8_eclass_num_vertices[eclass]);
            for (i = 0; i < t8_eclass_num_vertices[eclass]; i++) {
                /* Get the i-th node index in t8code order and store it. */
                stored_indices[i] = node_indices[t8wyo_vertex_to_mcell_vertex_num[eclass][i]];
            }
            /* set the index array as a new entry in the array */
            *(long **) sc_array_push(*vertex_indices) = stored_indices;
        }

        /* advance the tree counter */
        tree_count++;
    }
}

/* If vertex_indices is not NULL, it is allocated and will store
 * the indices of each tree vertices; stored as array of long ints.
 *
 * If occ geometry is used, the geometry is passed as a pointer here.
 * We cannot access this geometry over the cmesh interface since the
 * cmesh is not committed yet.
 */
static int
t8_cmesh_mcell_elements(const mcell_t *mcell,t8_cmesh_t cmesh,
                        sc_hash_t *vertices,
                        sc_array_t **vertex_indices,
                        const t8_geometry_c *linear_geometry_base,
                        const int use_occ_geometry,
                        const t8_geometry_c *occ_geometry_base){

    /* allocate list of vertex indices for each element */
    if(vertex_indices != NULL) *vertex_indices = sc_array_new(sizeof(long *));

    /* reset index of next tree to insert */
    t8_gloidx_t tree_count = 0;

    /* ============================ */
    /* set elements for each eclass */
    /* ============================ */
    fill_element_class(cmesh,vertices,vertex_indices,T8_ECLASS_TET,    mcell->ntetra,mcell->ndc4,tree_count);
    fill_element_class(cmesh,vertices,vertex_indices,T8_ECLASS_PYRAMID,mcell->npyr,  mcell->ndc5,tree_count);
    fill_element_class(cmesh,vertices,vertex_indices,T8_ECLASS_PRISM,  mcell->nprizm,mcell->ndc6,tree_count);
    fill_element_class(cmesh,vertices,vertex_indices,T8_ECLASS_HEX,    mcell->nhex,  mcell->ndc8,tree_count);

    printf("Successfully constructed Elements.\n");
    return 0;
}

/* struct stores all information associated to a tree's face */
typedef struct {
    t8_locidx_t ltree_id; /* local id of the tree this face belongs */
    int8_t face_number;   /* number of that face within the tree */
    int num_vertices;     /* number of vertices of this face */
    long *vertices;       /* indices of these vertices */
}
t8_mcell_face_t;

/* Hash a face. The hash value is the sum of its vertex indices */
static unsigned
t8_mcell_face_hash(const void *face, const void *data){
  t8_mcell_face_t *Face;
  int iv;
  long sum = 0;

  Face = (t8_mcell_face_t *) face;
  for (iv = 0; iv < Face->num_vertices; iv++) {
    sum += Face->vertices[iv];
  }
  T8_ASSERT(sum >= 0);
  return (unsigned) sum;
}

/* two faces are considered equal if they have the same vertices up to renumeration. */
static int
t8_mcell_face_equal(const void *facea,
                    const void *faceb,
                    const void *data){
    int  iv, jv, ret;
    long vertex;
    t8_mcell_face_t *Face_a;
    t8_mcell_face_t *Face_b;

    Face_a = (t8_mcell_face_t *) facea;
    Face_b = (t8_mcell_face_t *) faceb;

    /* if both have different number of vertices they can't be equal */
    ret = Face_a->num_vertices == Face_b->num_vertices;
    if(!ret) return 0;

    /* check for each vertex of Face_a whether it is a vertex of Face_b */
    for (iv = 0; iv < Face_a->num_vertices; iv++) {
        vertex = Face_a->vertices[iv];
        ret = 0;
        for (jv = 0; jv < Face_b->num_vertices; jv++) {
            ret |= vertex == Face_b->vertices[jv];
        }
        /* if vertex was a vertex of Face_b then ret is true here */
        if(!ret) return 0;
    }
    return 1;
}

/* We use this function in a loop over all elements
 * in the hash table, to free the memory of the vertices array */
static int
t8_mcell_face_free(void **face,const void *data){
    t8_mcell_face_t *Face = *(t8_mcell_face_t **) face;
    T8_FREE(Face->vertices);
    return 1;
}

/* Given a face and a cmesh set the face as a domain boundary.
 * We use this function in a loop over all hashed faces.
 */
static int
t8_mcell_face_set_boundary(void **face, const void *data){
    t8_mcell_face_t *Face;
    t8_gloidx_t gtree_id;
    t8_cmesh_t cmesh;

    cmesh = (t8_cmesh_t) data;
    Face = *(t8_mcell_face_t **) face;

    /* get the global tree id */
    gtree_id = Face->ltree_id;

    /* Set the Face as a domain boundary */
    t8_cmesh_set_join(cmesh,gtree_id,gtree_id,
                      Face->face_number,
                      Face->face_number,
                      0);
    return 1;
}

/* Given two faces and the classes of their volume trees,
 * compute the orientation of the faces to each other */
static int
t8_mcell_face_orientation(t8_mcell_face_t *Face_a,
                          t8_mcell_face_t *Face_b,
                          t8_eclass_t tree_class_a,
                          t8_eclass_t tree_class_b){
    t8_mcell_face_t *smaller_Face,*bigger_Face;
    t8_eclass_t bigger_class;
    int orientation = -1;
    int compare,iv;

    compare = t8_eclass_compare(tree_class_a,tree_class_b);
    if (compare > 0) {
        /* tree_class_a is bigger than tree_class_b */
        smaller_Face = Face_b;
        bigger_Face = Face_a;
        bigger_class = (t8_eclass_t) t8_eclass_face_types[tree_class_a][Face_a->face_number];
    } else
    if (compare < 0) {
        /* tree_class_a is smaller than tree_class_b */
        smaller_Face = Face_a;
        bigger_Face = Face_b;
        bigger_class = (t8_eclass_t) t8_eclass_face_types[tree_class_b][Face_b->face_number];
    } else {
        /* both classes are the same: face with the smaller face id is the smaller one */
        if (Face_a->face_number < Face_b->face_number) {
            smaller_Face = Face_a;
            bigger_Face = Face_b;
            bigger_class = (t8_eclass_t) t8_eclass_face_types[tree_class_b][Face_b->face_number];
        } else {
            smaller_Face = Face_b;
            bigger_Face = Face_a;
            bigger_class = (t8_eclass_t) t8_eclass_face_types[tree_class_a][Face_a->face_number];
        }
    }

    /* number of the first vertex of the smaller face */
    long vertex_zero = smaller_Face->vertices[0];

    /* Find which point in the bigger face is vertex_zero */
    for (iv = 0; iv < t8_eclass_num_vertices[bigger_class]; iv++) {
        if (vertex_zero == bigger_Face->vertices[iv]) {
            /* found corresponding vertex */
            orientation = iv;

            /* set condition to break the loop */
            iv = t8_eclass_num_vertices[bigger_class];
        }
    }
    T8_ASSERT(orientation >= 0);
    return orientation;
}

/* Given the number of vertices and for each element a list of its
 * vertices, find the neighborship relations of each element.
 * This routine only finds neighbors between local trees.
 * Use with care if cmesh is partitioned.
 */
static void
t8wyo_cmesh_mcell_find_neighbors(t8_cmesh_t cmesh,
                                 sc_array_t *vertex_indices){
    sc_hash_t *faces;
    t8_mcell_face_t *Face, **pNeighbor, *Neighbor;
    sc_mempool_t *face_mempool;
    t8_gloidx_t gtree_it;
    t8_gloidx_t gtree_id,gtree_neighbor;
    t8_eclass_t eclass,face_class,neighbor_tclass;
    int num_face_vertices,face_it,vertex_it;
    int retval, orientation;
    long *tree_vertices;
    t8_stash_class_struct_t *class_entry;

    face_mempool = sc_mempool_new(sizeof(t8_mcell_face_t));
    faces = sc_hash_new(t8_mcell_face_hash,
                        t8_mcell_face_equal,
                        cmesh,NULL);

    /* TODO: Does currently not work with partitioned cmesh */
    T8_ASSERT(!cmesh->set_partition);

    /* The cmesh is not allowed to be committed yet */
    T8_ASSERT(t8_cmesh_is_initialized(cmesh));
    t8_debugf("Starting to find tree neighbors\n");

    /* iterate over all local trees */
    for (gtree_it = 0;
         gtree_it < (t8_gloidx_t) cmesh->stash->classes.elem_count;
         gtree_it++) {

        /* We get the class of the current tree.
         * Since we know that the trees were put into the stash in order
         * of their tree id's, we can just read the corresponding entry from the stash.
         *
         * !WARNING: This does not work in general to find the class of a tree
         *    since the order in which the trees are added to the stash is arbitrary.
         *    Use t8_stash_class_bsearch in tat case.
         */
        class_entry = (t8_stash_class_struct_t *) t8_sc_array_index_locidx(&cmesh->stash->classes, gtree_it);
        T8_ASSERT (class_entry->id == gtree_it);
        eclass = class_entry->eclass;

        /* Get the vertices of that tree */
        tree_vertices = *(long **) t8_sc_array_index_locidx(vertex_indices,gtree_it);

        /* loop over all faces of the tree */
        for (face_it = 0; face_it < t8_eclass_num_faces[eclass]; face_it++) {

            Face = (t8_mcell_face_t *) sc_mempool_alloc(face_mempool);

            /* get its eclass and the number of vertices */
            face_class = (t8_eclass_t) t8_eclass_face_types[eclass][face_it];
            num_face_vertices = t8_eclass_num_vertices[face_class];
            Face->vertices = T8_ALLOC(long,num_face_vertices);
            Face->num_vertices = num_face_vertices;
            Face->ltree_id = gtree_it;
            Face->face_number = face_it;

            /* Copy the vertices of the face to the face struct */
            for (vertex_it = 0; vertex_it < num_face_vertices; vertex_it++) {
                Face->vertices[vertex_it] = tree_vertices[t8_face_vertex_to_tree_vertex[eclass][face_it][vertex_it]];
            }

            /* try to insert the face into the hash */
            retval = sc_hash_insert_unique(faces,Face,(void ***) &pNeighbor);
            if (!retval) {
                /* face was already in the hash */
                Neighbor = *pNeighbor;
                T8_ASSERT(Neighbor->ltree_id != gtree_it);

                /* current tree is a neighbor to the tree Neighbor->ltree_id;
                 * need to identify the face number and the orientation
                 */
                gtree_id = gtree_it;
                gtree_neighbor = Neighbor->ltree_id;

                /* compute the orientation of the face connection */
                /* get the element class of the neighbor tree */
                class_entry = (t8_stash_class_struct_t *) t8_sc_array_index_locidx(&cmesh->stash->classes,Neighbor->ltree_id);
                T8_ASSERT(class_entry->id == Neighbor->ltree_id);
                neighbor_tclass = class_entry->eclass;

                /* Calculate the orientation */
                orientation = t8_mcell_face_orientation(Face,Neighbor,eclass,neighbor_tclass);

                /* set the face connection */
                t8_cmesh_set_join(cmesh,
                                  gtree_id,gtree_neighbor,
                                  face_it,Neighbor->face_number,
                                  orientation);

                /* cleanup memory */
                T8_FREE(Face->vertices);
                sc_mempool_free(face_mempool,Face);

                /* free the neighbor here as well */
                sc_hash_remove(faces,Neighbor,NULL);
                T8_FREE(Neighbor->vertices);
                sc_mempool_free(face_mempool,Neighbor);
            }
        }
    }

    /* remaining faces are domain boundaries */
    sc_hash_foreach(faces,t8_mcell_face_set_boundary);

    /* free the faces and the hash */
    sc_hash_foreach(faces,t8_mcell_face_free);
    sc_mempool_destroy(face_mempool);
    sc_hash_destroy(faces);
    printf("Successfully constructed Tree Neighbors.\n");
}

/* This is a helper function to properly register the
 * geometries for the cmesh created in t8_cmesh_from_mcell.
 * It should be called by all processes of the cmesh.
 * Returns 1 on success, 0 on OCC usage error: use_occ_geometry true, but occ not linked.
 * The linear_geometry pointer will point to the newly created linear geometry.
 * The occ_geometry pointer will point to the newly created occ geometry,
 * or to NULL if no occ geometry is used.
 */
static int
t8_cmesh_from_mcell_register_geometries(t8_cmesh_t cmesh,
                                             const int use_occ_geometry,
                                             const int dim,
                                             const char *fileprefix,
                                             const t8_geometry_c **linear_geometry,
                                             const t8_geometry_c **occ_geometry){

    /* register linear geometry */
    const t8_geometry_c *linear_geom = new t8_geometry_linear(dim);
    t8_cmesh_register_geometry(cmesh,linear_geom);
    *linear_geometry = linear_geom;

    if (use_occ_geometry) {
        (void) fileprefix;
#if T8_WITH_OCC
        const t8_geometry_c *occ_geom = t8_geometry_occ_new(dim,fileprefix,"brep_geometry");
        t8_cmesh_register_geometry(cmesh,occ_geom);
        *occ_geometry = occ_geom;
#else /* !T8_WITH_OCC */
        *occ_geometry = NULL;
        return 0;
#endif
    }
    return 1;
}

t8_cmesh_t
t8_cmesh_from_mcell(mcell_t *mcell,
                    int do_bcast,MPI_Comm comm,
                    int dim,int use_occ_geometry){
    t8_cmesh_t cmesh;
    sc_hash_t *vertices = NULL;
    sc_mempool_t *node_mempool = NULL;
    sc_array_t *vertex_indices = NULL;
    long *indices_entry;
    int mpirank,mpisize,mpiret;

    const t8_geometry_c *occ_geometry = NULL;
    const t8_geometry_c *linear_geometry = NULL;

    /* MCELL: 3D only */
    T8_ASSERT(dim==3);

    mpiret = MPI_Comm_size(comm,&mpisize);
    mpiret = MPI_Comm_rank(comm,&mpirank);

    if (!do_bcast || mpirank == 0) {
        /* initialize cmesh structure */
        t8_cmesh_init(&cmesh);

        /* set dimension */
        t8_cmesh_set_dimension(cmesh,dim);

        /* set vertices from mcell data */
        vertices = t8wyo_mcell_nodes(mcell,&node_mempool);

        /* set vertices */
        t8_cmesh_mcell_elements(mcell,cmesh,vertices,&vertex_indices,
                                linear_geometry,use_occ_geometry,occ_geometry);

        /* set faces */
        t8wyo_cmesh_mcell_find_neighbors(cmesh,vertex_indices);

        if(vertices != NULL) sc_hash_destroy(vertices);
        sc_mempool_destroy(node_mempool);

        while (vertex_indices->elem_count > 0) {
            indices_entry = *(long **) sc_array_pop(vertex_indices);
            T8_FREE(indices_entry);
        }
        sc_array_destroy(vertex_indices);
    }

    if (do_bcast) {
        if(mpirank != 0) cmesh = NULL;
        cmesh = t8_cmesh_bcast(cmesh,0,comm);
    }

    /* register the geometries for the cmesh */
    const int registered_geom_success =
        t8_cmesh_from_mcell_register_geometries(cmesh,use_occ_geometry,dim,
                                                mcell->prefix,&linear_geometry,
                                               &occ_geometry);

    if (!registered_geom_success) {
        /* registering failed */
        t8_errorf("OCC is not linked. Cannot use OCC geometry.\n");
        t8_cmesh_destroy(&cmesh);
        return NULL;
    }

    /* commit the cmesh */
    if(cmesh != NULL) t8_cmesh_commit(cmesh,comm);

    /* deallocate fortran data */
    if(mpirank == 0) printf("Successfully constructed cmesh.\n");
    return cmesh;
}