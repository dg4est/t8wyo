/**
 * \file   amr_test.h
 * \author akirby, mbrazell
 */

#ifndef AMR_TEST_H
#define AMR_TEST_H

/* header files */
#include "t8wyo_solver.hxx"
#include "memory.hxx"

/* 3PL header files */
#include <sc.h>
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

#ifdef __cplusplus
extern "C" {
#endif

/** MPI and p4est initialization function wrapper
 *
 * @param [in]    argc      number of command line arguments
 * @param [in]    argv      command line arguments
 * @param [inout] ctx       context data
 */
int test(int a);

#ifdef __cplusplus
}
#endif
#endif /* AMR_INITIALIZE_H */
