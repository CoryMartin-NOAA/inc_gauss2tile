if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model strict -convert big_endian")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fconvert=big-endian -ffree-line-length-0 -ffast-math -fno-second-underscore -frecord-marker=4")
endif()

add_executable(interp_inc.x driver.f90)
set_target_properties(interp_inc.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
target_link_libraries(interp_inc.x PRIVATE NetCDF::NetCDF_Fortran)
target_link_libraries(interp_inc.x PRIVATE MPI::MPI_Fortran)

# Install executable targets
install(TARGETS interp_inc.x RUNTIME DESTINATION bin)
