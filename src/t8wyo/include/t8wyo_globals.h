/**
 * File:   t8wyo_globals.h
 * Author: akirby
 *
 * Created on July 6, 2023, 2:15 PM
 */

#ifndef T8WYO_GLOBALS_H
#define T8WYO_GLOBALS_H

/* header files */
#include "t8wyo_solver.hxx"

/* 3PL header files */
#include <sc.h>
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>

/* global defines */
#define T8WYO_EXTERN_C_BEGIN() SC_EXTERN_C_BEGIN
#define T8WYO_EXTERN_C_END()   SC_EXTERN_C_END

#ifdef _2D_
#  define DIM 2
#else
#  define DIM 3
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* global t8wyo structure for interface code */
t8_forest_t t8wyo_forest;
t8_cmesh_t t8wyo_cmesh;
t8wyo_t t8wyo;

wyo::memory<int> face2cell;
wyo::memory<int> ifacetype;
wyo::memory<int> elem_info;
wyo::memory<Real> elem_vol;
wyo::memory<Real> face_norm;

#ifdef __cplusplus
}
#endif
#endif /* T8WYO_GLOBALS_H */