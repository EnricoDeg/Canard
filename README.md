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

To get a list of supported machines for installation:

```bash
./install -l
```

If you can not find your machine in the supported ones, you can add a file in the `machine` folder and it will be automatically added to the available machines.

To get a list of compilers and `MPI` implementations combinations for a supported machine:

```bash
./install -m ${MACHINE_NAME} -c help
```

Once the model is installed, the `canard` executable can be found in the `bin` folder. You can then setup an experiment using a script in the `run` folder or writing your own one.

#### Input parameters

The input parameters are in `input.canard` file.

##### nml_driver

- `nio`: number of IO servers
- `mbk`: domain blocks `(0:mbk)`
- `ndata`: output files `(0:ndata)`

##### nml_canard

- `nts`: initial conditions (**0**: initialize physics with upstream flow -- **1**: initialize physics with restart files)
- `nrestart`: write restart files (**0**: OFF -- **1**: ON)
- `cfl`: CFL number
- `tmax`: end time of the simulation
- `ltimer`: enable profiling (`.true.`/`.false.`)

##### nml_domdcomp

- `nbpc`: number of `MPI` processes for each domain block (2D array `[0:mbk,1:3]`)

- `lximb`: number of cells for each domain block in `x` direction (1D array `[0:mbk]`)
- `letmb`: number of cells for each domain block in `y` direction (1D array `[0:mbk]`)
- `lzemb`: number of cells for each domain block in `z` direction (1D array `[0:mbk]`)
- `nbbc`: blocks boundary conditions (3D array `[0:mbk,1:3,1:2]`) 
- `mbcd`: block communication (3D array `[0:mbk,1:3,1:2]`)

##### nml_aio

- `nbpc`: number of `MPI` processes for each domain block (2D array `[0:mbk,1:3]`)

- `lximb`: number of cells for each domain block in `x` direction (1D array `[0:mbk]`)
- `letmb`: number of cells for each domain block in `y` direction (1D array `[0:mbk]`)
- `lzemb`: number of cells for each domain block in `z` direction (1D array `[0:mbk]`)
- `nbbc`: blocks boundary conditions (3D array `[0:mbk,1:3,1:2]`) 
- `mbcd`: block communication (3D array `[0:mbk,1:3,1:2]`)

##### nml_numerics

- `fltk`: filter coefficient
- `fltrbc`: boundary filter coefficient

##### nml_physics

- `reoo`: upstream flow Reynolds number
- `tempoo`: upstream flow temperature
- `amach1`: upstream flow Mach number (`x` direction)
- `amach2`: upstream flow Mach number (`y` direction)
- `amach3`: upstream flow Mach number (`z` direction)
- `timf`: moving frame time
- `nsmf`: enable moving frame (**0**: OFF -- **1**: ON)

##### nml_sponge

- `szco`: magnitude of non reflective sponge
