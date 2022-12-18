# Find MPI library
find_package(MPI REQUIRED COMPONENTS Fortran)

option(ENABLE_OPENACC "compile with OpenACC support" OFF)

set(SOURCE_EXE  mo_kind.f90
                mo_parameters.f90
                mo_mpi.f90
                mo_utils.f90
                mo_domdcomp.f90
                mo_io.f90
                dcanard.f90)

add_executable(dcanard ${SOURCE_EXE})

# Define compiler options for all supported compilers
if( ${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
  target_compile_options(dcanard PUBLIC
      $<$<CONFIG:DEBUG>:-cpp>
      $<$<CONFIG:RELEASE>:-cpp>
  )
elseif( ${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )
  target_compile_options(dcanard PUBLIC
      $<$<CONFIG:DEBUG>:-cpp>
      $<$<CONFIG:RELEASE>:-cpp>
  )
elseif( ${CMAKE_Fortran_COMPILER_ID} STREQUAL NVHPC )
  target_compile_options(dcanard PUBLIC
      $<$<CONFIG:DEBUG>:-cpp -Mnofma>
      $<$<CONFIG:RELEASE>:-cpp>
  )
  if(${ENABLE_OPENACC})
    target_compile_options(dcanard PUBLIC
        $<$<CONFIG:DEBUG>:-cpp -acc -gpu=math_uniform,cuda11.7,cc80 -Minfo=accel,inline>
        $<$<CONFIG:RELEASE>:-cpp -acc -gpu=cc80 -Minfo=accel,inline>
    )
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-acc -gpu=math_uniform,cuda11.7,cc80")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-acc -gpu=cc80")
  endif()
endif()

# Install executables.
install (TARGETS
  dcanard # executables
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)

