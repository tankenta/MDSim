# MDSim
MDSim is a simple molecular dynamics simulator written in C++11.

## Features
* Lennard-Jones potential
* Velocity Verlet algorithm
* Temperature control by velocity scaling
* Rectangular cell (periodic or free boundary condition)

## Testing environments
### Ubuntu 14.04
* g++ 4.8.4
* Eigen 3.3.4

## Setup
Install Eigen.
```
./setup.sh
```

## Usage
```
cd src/
make
./mdsim dt total_time temp_cont_time phase bc_mode
```
