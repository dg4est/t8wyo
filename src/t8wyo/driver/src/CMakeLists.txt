set(SRC
    my_kinddefs.F90
)

# build driver library
add_library(t8wyo-driver SHARED ${SRC})
target_compile_definitions(t8wyo-driver PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
target_include_directories(t8wyo-driver PUBLIC ${MPI_Fortran_INCLUDE_PATH})

# link 3PL libraries to shared library
target_link_libraries(t8wyo-driver t8wyo3d ${MPI_Fortran_LIBRARIES})

# build executable
add_executable(t8wyo3d.mpi main.F90)

# link shared library to executable
target_link_libraries(t8wyo3d.mpi t8wyo-driver ${MPI_Fortran_LIBRARIES})

# install executables and libraries
install(TARGETS t8wyo-driver t8wyo3d.mpi
        RUNTIME DESTINATION bin/
        LIBRARY DESTINATION lib/
        ARCHIVE DESTINATION lib/
)