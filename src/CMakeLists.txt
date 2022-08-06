# Find MPI library
find_package(MPI REQUIRED COMPONENTS Fortran)

set(SOURCE_EXE  mo_kind.f90
                mo_parameters.f90
                mo_mpi.f90
                mo_timer.f90
                mo_utils.f90
                mo_domdcomp.f90
                mo_numerics.f90
                mo_grid.f90
                mo_io.f90
                mo_gridgen.f90
                mo_sponge.f90
                mo_physics.f90
                mo_gcbc.f90
                canard.f90)

add_executable(canard ${SOURCE_EXE})
if(ENABLE_VISCOSITY)
    target_compile_definitions(canard PRIVATE VISCOUS)
endif()

target_compile_options(canard PUBLIC
    $<$<CONFIG:DEBUG>:-cpp -O0>
    $<$<CONFIG:RELEASE>:-cpp -O3>
  )

# Install executables.
install (TARGETS
  canard # executables
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)

