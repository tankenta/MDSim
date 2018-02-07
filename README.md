# MDSim
MDSim is a simple molecular dynamics simulator written in C++11.

## Features
* Lennard-Jones potential
* Velocity Verlet algorithm
* Temperature control by velocity scaling
* Rectangular cell (periodic or free boundary condition)
* CSV output

## Testing environments
### Ubuntu 14.04
* g++ 4.8.4
* Eigen 3.3.4

## Build
1. Download Eigen
```
./setup.sh
```

2. Build
```
cd src/
make
```

## Usage
```
./mdsim dt total_time temp_cont_time phase bc_mode csv_dir
```
