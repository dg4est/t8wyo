/**
 * \file    t8wyo_cell3D.cxx
 * \author  akirby
 *
 * \brief   CELL3D interface.
 */

/* header files */
#include "t8wyo_cell3d.h"

void cell3d_initialize(){
    cell3d_initialize_();
}

void cell3d_read_mesh(mcell_t &mcell){
    cell3d_read_mesh_interface_(&mcell.ntetra,&mcell.npyr,&mcell.nprizm,&mcell.nhex,
                                &mcell.ntetra_ng,&mcell.npyr_ng,&mcell.nprizm_ng,&mcell.nhex_ng,
                                &mcell.nnode,&mcell.ngpt,
                                &mcell.nbface3,&mcell.nbface4,
                                &mcell.nface3,&mcell.nface4,
                                &mcell.ndc4,&mcell.ndc5,&mcell.ndc6,&mcell.ndc8,
                                &mcell.nbf3,&mcell.nbf4,
                                &mcell.ifpat3,&mcell.ifpat4,
                                &mcell.xgeom,
                                &mcell.ndf3,&mcell.ndf4);
}

void cell3d_deallocate_mesh(){
    cell3d_deallocate_mesh_();
}