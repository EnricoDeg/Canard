PGM	= run.exe
FC	= mpifort
OBJ	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridAerofoil.f90 FanWake.f90 Main3D.f90
$(PGM): $(OBJ)
	$(FC) -g -O0 -o $(PGM) $(OBJ)
	rm *.mod
