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
    t8wyo_cmesh_from_mcell.cxx
    t8wyo_build_lists_ext.cxx
    t8wyo_initialize.cxx
    t8wyo_interface.cxx
)

# build shared libraries
add_library(t8wyo3d SHARED ${SRC} ${SRC_H} ${GIT_FILE})
target_compile_definitions(t8wyo3d PUBLIC -DMPI_ON)
target_link_libraries(t8wyo3d ${sc_lb} ${t8code_lb} ${metis_lb} m)

#add_library(t8wyo2d SHARED ${SRC} ${SRC_H} ${GIT_FILE})
#target_compile_definitions(t8wyo2d PUBLIC -D_2D_)
#target_compile_definitions(t8wyo2d PUBLIC -DMPI_ON)
#target_link_libraries(t8wyo2d ${sc_lb} ${t8code_lb} ${metis_lb} m)

# install executables and libraries
install(TARGETS t8wyo3d
        RUNTIME DESTINATION bin/
        LIBRARY DESTINATION lib/
        ARCHIVE DESTINATION lib/
)

# install header files
install(FILES ${SRC_H} DESTINATION include)
