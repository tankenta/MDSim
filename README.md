# MDSim
MDSim is a simple molecular dynamics simulator written in C++11.

## Features
* Lennard-Jones potential
* Velocity Verlet algorithm
* Temperature control by velocity scaling
* Rectangular cell (periodic or free boundary condition)
* CSV output
* Python plot script

## Testing environments
### Ubuntu 14.04
* g++ 4.8.4
* Eigen 3.3.4
* Python 2.7.6
* matplotlib 2.0.2

## Build
1. Setup library (Eigen)
```
cd MDSim/
./setup.sh
```

2. Build
```
cd src/
make
```

## Usage
### Simulation
```
mkdir your_csv_dir/
./mdsim dt total_time temp_cont_time phase bc_mode csv_dir
```

### Plot
```
cd MDSim/plot_tools/
./plot_results.py path/to/your_csv_dir/
```

## TODO
- [x] Organize codes into classes
- [ ] Flexible parameter settings by command-line arguments
- [ ] Add GUI
