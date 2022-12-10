# Find MPI library
find_package(MPI REQUIRED COMPONENTS Fortran)

option(ENABLE_OPENACC "compile with OpenACC support" OFF)

set(SOURCE_EXE  mo_kind.f90
                mo_parameters.f90
                mo_mpi.f90
                mo_model_ioserver_exchange.f90
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
                mo_io_server.f90
                mo_canard_driver.f90
                mo_aio_driver.f90
                driver.f90)

add_executable(canard ${SOURCE_EXE})
if(ENABLE_VISCOSITY)
    target_compile_definitions(canard PRIVATE VISCOUS)
endif()

# Define compiler options for all supported compilers
if( ${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
  target_compile_options(canard PUBLIC
      $<$<CONFIG:DEBUG>:-cpp>
      $<$<CONFIG:RELEASE>:-cpp>
  )
elseif( ${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )

elseif( ${CMAKE_Fortran_COMPILER_ID} STREQUAL NVHPC )
  target_compile_options(canard PUBLIC
      $<$<CONFIG:DEBUG>:-cpp -Mnofma>
      $<$<CONFIG:RELEASE>:-cpp>
  )
  if(${ENABLE_OPENACC})
    target_compile_options(canard PUBLIC
        $<$<CONFIG:DEBUG>:-cpp -acc -gpu=math_uniform,cuda11.7,cc80 -Minfo=accel,inline>
        $<$<CONFIG:RELEASE>:-cpp -acc -gpu=cc80 -Minfo=accel,inline>
    )
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-acc -gpu=math_uniform,cuda11.7,cc80")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-acc -gpu=cc80")
  endif()
endif()

# Install executables.
install (TARGETS
  canard # executables
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)

