enable_language(Fortran)

set(CMAKE_Fortran_FLAGS
 "${CMAKE_Fortran_FLAGS} -convert big_endian -assume buffered_io")

# source files
set(CELL_MOD_SRC
    my_kinddefs.F90
    mpi_schedule_typedef.F90
    my_typedefs.F90
    mpi_schedules.F90
    io_params.F90
    mp_stuff.F90
    my_mpi_barrier.F90
    params.F90
    get_filename1.F90
    global_arrays_dynamic.F90
    local_arrays_dynamic.F90
    my_allocate.F90
    mpi_schedule_read.F90
    my_deallocate.F90
    all_reduce_max.F90
    all_reduce_min.F90
    all_reduce_sum.F90
    mpi_to_gpt.F90
    triarea.F90
)

set(CELL3D_SRC
    cell3d_initialize.F90
    cell3d_init.F90
    program_mpi_init.F90
    stop_all.F90
    abort_all.F90
    get_time_mp.F90
    readarg.F90
    init_io.F90
    readpm.F90
    char_adjust.F90
    wout.F90
    cell3d_read_mesh.F90
    cell3d_flux_test.F90
    cell3d_test_connectivity.F90
)

# build shared library
add_library(cell_solver_2d SHARED ${CELL_MOD_SRC} ${CELL3D_SRC})
add_library(cell_solver_3d SHARED ${CELL_MOD_SRC} ${CELL3D_SRC})

target_compile_definitions(cell_solver_2d PUBLIC -DMPI_ON)
target_compile_definitions(cell_solver_3d PUBLIC -DMPI_ON)


# install math_kernels
install(TARGETS cell_solver_2d cell_solver_3d
        RUNTIME DESTINATION bin/
        LIBRARY DESTINATION lib/
        ARCHIVE DESTINATION lib/
)