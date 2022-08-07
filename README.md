# Canard

The *Compressible Aerodynamics & Aeroacoustics Research Code* is an open source, parallel, high-fidelity CFD code.

Main characteristics:

- multi-block structured grid with body fitted curvilinear mesh
- numerical solution of incompressible / compressible Navier-Stokes equations
- Implicit Large Eddy Simulations (ILES)
- High order compact finite difference schemes
- High order compact filters
- Parallelization with `MPI` library (domain decomposition)
- Generalized Characteristic Boundary conditions (GCBC)
- Non reflective treatments (sponge layer)

#### Download and Install

```bash
git clone https://github.com/EnricoDeg/Canard.git
```

To get a list of all installation options:

```bash
./install -h
```

To get a list of supported machine for installation:

```bash
./install -l
```

If you can not find your machine in the supported ones, you can add a file in the `machine` folder and it will be automatically added to the available machines.

Once the model is installed, the `canard` executable can be found in the `bin` folder. You can then setup an experiment using a script in the `run` folder or writing your own one.

#### Input parameters

The input files are located in `input/parameters` folder.

##### input.main

- `mbk`: domain blocks `(0:mbk)`
- `nts`: initial conditions (**0**: initialize physics with upstream flow -- **1**: initialize physics with restart files)
- `nscrn`: LOG file output frequency (number of time step -- simulation time)
- `ndata`: output files `(0:ndata)`
- `ndatafl`: unsteady fluctuations in output files (**0**: OFF -- **1**: ON)
- `ndataav`: data average in output data (**0**: OFF -- **1**: ON)
- `nkrk`: number of Runge-Kutta stages (time integration)
- `nrestart`: write restart files (**0**: OFF -- **1**: ON)
- `cfl`: CFL number
- `tmax`: end time of the simulation
- `tsam`: start time for output files
- `dto`: time step value (if `dto>0`)
- `nbody`: presence of body in computational domain (**0**: NO -- **1**: YES)
- `ltimer`: enable profiling (`.true.`/`.false.`)

##### input.domdcomp

- `npbxi`: number of `MPI` processes for each domain block in `x` direction (1D array `[0:mbk]`)
- `npbet`: number of `MPI` processes for each domain block in `y` direction (1D array `[0:mbk]`)
- `npbze`: number of `MPI` processes for each domain block in `z` direction (1D array `[0:mbk]`)
- `lximb`: number of cells for each domain block in `x` direction (1D array `[0:mbk]`)
- `letmb`: number of cells for each domain block in `y` direction (1D array `[0:mbk]`)
- `lzemb`: number of cells for each domain block in `z` direction (1D array `[0:mbk]`)

##### input.gridgen

- `lxi0`: 
- `lxi1`: 
- `lxi2`: 
- `let0`: 
- `lze0`: 
- `nthick`: 
- `smg`: 
- `smgvr`: 
- `doml0`: 
- `doml1`: 
- `domh`: 
- `span`: 
- `wlew`: 
- `wlea`: 
- `szth0`: 
- `szth1`: 
- `skew`: 
- `spx`: 
- `gridtype`: 

##### input.numerics

- `fltk`:
- `fltrbc`:
- `nnf1-2-3`:

##### input.physics

- `reoo`: upstream flow Reynolds number
- `tempoo`: upstream flow temperature
- `amach1`: upstream flow Mach number (`x` direction)
- `amach2`: upstream flow Mach number (`y` direction)
- `amach3`: upstream flow Mach number (`z` direction)
- `timf`: moving frame time
- `nsmf`: enable moving frame (**0**: OFF -- **1**: ON)

##### input.sponge

- `szco`: 

##### input.gcbc

- `wtemp`: wall temperature coefficient
- `nextgcic`: 

