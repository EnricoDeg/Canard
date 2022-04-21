PGM	= canard.x
FC	= mpifort
OBJ	= mo_kind.f90 mo_parameters.f90 mo_vars.f90 mo_mpi.f90 mo_utils.f90 mo_numerics.f90 \
          mo_gridgen.f90 mo_diagnostics.f90 mo_io.f90 mo_domdcomp.f90 mo_grid.f90 \
          mo_sponge.f90 mo_gcbc.f90 mo_physics.f90 canard.f90
$(PGM): $(OBJ)
	$(FC) -g -O0 -DVISCOUS -cpp -o $(PGM) $(OBJ)
	rm *.mod
