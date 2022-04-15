PGM	= run.exe
FC	= mpifort
OBJ	= MainVar3D.f90 mo_mpi.f90 mo_utils.f90 mo_numerics.f90 \
      GridAerofoil.f90 FanWake.f90 mo_diagnostics.f90 mo_io.f90 mo_domdcomp.f90 mo_grid.f90 \
      mo_sponge.f90 mo_gcbc.f90 Main3D.f90
$(PGM): $(OBJ)
	$(FC) -g -O0 -cpp -o $(PGM) $(OBJ)
	rm *.mod
