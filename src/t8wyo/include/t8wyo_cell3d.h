/**
 * File:   t8wyo_cell3d.h
 * Author: akirby
 *
 * Created on June 25, 2023, 3:20 PM
 */

#ifndef CELL3D_H
#define CELL3D_H

/* system header files */
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ndf3 boundary face indices */
#define NDF3_N1 0 // node 1
#define NDF3_N2 1 // node 2
#define NDF3_N3 2 // node 3
#define NDF3_E1 3 // element 1
#define NDF3_E2 4 // element 2 (e2<0: boundary face nbf3 index)
#define NDF3_SZ 5 // number of entries per tri face

/* ndf4 boundary face indices */
#define NDF4_N1 0 // node 1
#define NDF4_N2 1 // node 2
#define NDF4_N3 2 // node 3
#define NDF4_N4 3 // node 4
#define NDF4_E1 4 // element 1
#define NDF4_E2 5 // element 2 (e2<0: boundary face nbf4 index)
#define NDF4_SZ 6 // number of entries per quad face

/* data structures */
typedef struct {
    int ntetra,npyr,nprizm,nhex;
    int ntetra_ng,npyr_ng,nprizm_ng,nhex_ng;
    int nnode,ngpt;
    int nbface3,nbface4,nface3,nface4;
    int *ndc4,*ndc5,*ndc6,*ndc8;
    int *ifpat3,*ifpat4;
    int *nbf3,*nbf4;
    int *ndf3,*ndf4;
    double *xgeom;
    char *prefix;
}
mcell_t;

/* external functions */
void cell3d_initialize_();
void cell3d_read_mesh_interface_(int *ntetra,int *npyr,int *nprizm,int *nhex,
                                 int *ntetra_ng,int *npyr_ng,int *nprizm_ng,int *nhex_ng,
                                 int *nnode,int *ngpt,
                                 int *nbface3,int *nbface4,
                                 int *nface3,int *nface4,
                                 int **ndc4,int **ndc5,int **ndc6,int **ndc8,
                                 int **nbf3,int **nbf4,
                                 int **ifpat3,int **ifpat4,
                                 double **xgeom,
                                 int **ndf3,int **ndf4);
void cell3d_deallocate_mesh_();

/* function headers */
void cell3d_initialize();
void cell3d_read_mesh(mcell_t &mcell);
void cell3d_deallocate_mesh();

#ifdef __cplusplus
}
#endif
#endif /* CELL3D_H */