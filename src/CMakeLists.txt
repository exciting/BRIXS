# Shared modules, compiled once and reused by all three executables. 
# Order doesn't matter: CMake's Fortran scanner infers the real compile
# order from each file's `use` statements.
set(BRIXS_CORE_SOURCES
  m_config.f90
  mod_constants.f90
  mod_mpi.F90
  mod_phdf5.F90
  mod_matmul.F90
  mod_io.F90
  mod_blocks.F90
  mod_blocks_k.F90
)

add_library(brixs_core OBJECT ${BRIXS_CORE_SOURCES})
target_include_directories(brixs_core PUBLIC ${HDF5_Fortran_INCLUDE_DIRS})
target_link_libraries(brixs_core PUBLIC ${HDF5_Fortran_LIBRARIES} ${LAPACK_LIBRARIES})

if(BRIXS_ENABLE_MPI)
  target_compile_definitions(brixs_core PUBLIC MPI)
  target_include_directories(brixs_core PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
  target_link_libraries(brixs_core PUBLIC MPI::MPI_Fortran)
endif()

foreach(brixs_exe rixs_pathway rixs_oscstr rixs_coherence)
  add_executable(${brixs_exe} ${brixs_exe}.F90)
  target_link_libraries(${brixs_exe} PRIVATE brixs_core)
endforeach()

install(TARGETS rixs_pathway rixs_oscstr rixs_coherence RUNTIME DESTINATION bin)
