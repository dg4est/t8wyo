# include 3PL header files
include_directories(${t8code_dir}/include
                    ${t8code_dir}/
                    ${metis_dir}/include)

# find 3PL libraries
#            <our_var_name>       <lib(name).so>       <library location>
find_library(t8code_lb      NAMES t8             PATHS ${t8code_dir}/lib)
find_library(sc_lb          NAMES sc             PATHS ${t8code_dir}/lib)
find_library(metis_lb       NAMES metis          PATHS ${metis_dir}/lib)

# git commit file markings
set(GIT_FILE "gitversion.c")
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${GIT_FILE}
           ${CMAKE_CURRENT_BINARY_DIR}/_${GIT_FILE}
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/../../../cmake/version.cmake
    COMMENT "Created ${GIT_FILE}"
)

# header files
file(GLOB SRC_H *.h *.hxx)

# src files
set(SRC
    t8wyo_cell3d.cxx
    t8wyo_cmesh_from_mcell.cxx
    amr_initialize.cxx
)

# Fortran math kernels
add_subdirectory(CELL3D)

# build shared libraries
add_library(t8wyo2d SHARED ${SRC} ${SRC_H} ${GIT_FILE})
add_library(t8wyo3d SHARED ${SRC} ${SRC_H} ${GIT_FILE})
target_compile_definitions(t8wyo2d PUBLIC -D_2D_)

target_compile_definitions(t8wyo2d PUBLIC -DMPI_ON)
target_compile_definitions(t8wyo3d PUBLIC -DMPI_ON)

# link 3PL libraries to shared library
target_link_libraries(t8wyo2d cell_solver_2d ${sc_lb} ${t8code_lb} ${metis_lb} m)
target_link_libraries(t8wyo3d cell_solver_3d ${sc_lb} ${t8code_lb} ${metis_lb} m)

# build executable
add_executable(t8wyo2d.mpi main.cxx)
add_executable(t8wyo3d.mpi main.cxx)

# link shared library to executable
target_link_libraries(t8wyo2d.mpi t8wyo2d)
target_link_libraries(t8wyo3d.mpi t8wyo3d)

# install executables and libraries
install(TARGETS t8wyo2d t8wyo3d t8wyo2d.mpi t8wyo3d.mpi
        RUNTIME DESTINATION bin/
        LIBRARY DESTINATION lib/
        ARCHIVE DESTINATION lib/
)

# install header files
install(FILES ${SRC_H} DESTINATION include)
