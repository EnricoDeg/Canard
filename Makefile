PGM	= run.exe
FC	= mpiifort
OBJ	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridAerofoil.f90 FanWake.f90 Main3D.f90
$(PGM): $(OBJ)
	$(FC) -O3 -o $(PGM) $(OBJ)
	rm *.mod
