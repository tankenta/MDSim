#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#define EIGEN_NO_DEBUG
#define EIGEN_MPL2_ONLY
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "file_utils.hpp"
#include "md_utils.hpp"
#include "vec_utils.hpp"

int main(int argc, char const* argv[])
{
    if (argc < 7) {
        std::cout << 
            "Usage: ./mdsim dt total_time temp_cont_time phase bc_mode csv_dir"
            << std::endl;
        exit(1);
    }

    std::cout << "initializing..." << std::endl;
    const double dt = std::stod(std::string(argv[1]));
    const double total_time = std::stod(std::string(argv[2]));
    const double temp_cont_time = std::stod(std::string(argv[3]));
    std::string phase = std::string(argv[4]);
    std::string bc_mode = std::string(argv[5]);
    std::string csv_dir = std::string(argv[6]);

    double number_density;
    if (phase == "solid") {
        number_density = 1.2;
    } else if (phase == "liquid") {
        number_density = 0.8;
    } else if (phase == "gas") {
        number_density = 0.01;
    } else {
        throw std::invalid_argument("phase must be 'solid' or 'liquid' or 'gas'.");
    }

    const double target_temp = 1.0;
    const int num_particles = 256;
    const double ptcl_mass = 1.0;
    const double particle_radius = 0.2;
    const int RDF_hist_size = 200;
    Eigen::Vector3i cube_size = Eigen::Vector3i::Constant(4);

    MDSim md(
            dt, total_time, temp_cont_time, number_density, bc_mode,
            target_temp, num_particles, ptcl_mass, RDF_hist_size, cube_size);

    std::cout << "running the simulation..." << std::endl;
    ProgressBar prog(md.steps.front(), md.steps.back());
    for (const auto& step : md.steps) {
        prog.printProgressBar(step);
        md.next(step);
    }
    prog.finish();

    md.calcTotalEnergy();
    exportPlotData2CSV(
            csv_dir + "total_energy.csv", md.times, md.total_energy_arr,
            "time", "total energy");

    exportPlotData2CSV(
            csv_dir + "potential_energy.csv", md.times, md.potential_arr,
            "time", "potential energy");

    exportPlotData2CSV(
            csv_dir + "kinetic_energy.csv", md.times, md.kinetic_energy_arr,
            "time", "kinetic energy");

    exportPlotData2CSV(
            csv_dir + "temperature.csv", md.times, md.current_temp_arr,
            "time", "temperature");

    auto dr_pair = md.calcRadialDistributionFunction();
    std::vector<double> distance = dr_pair.first;
    std::vector<double> RDF = dr_pair.second;
    exportPlotData2CSV(
            csv_dir + "RDF.csv", distance, RDF,
            "r", "g(r)");

    auto tm_pair = md.calcMeanSquareDisplacement();
    std::vector<double> MSD_time = tm_pair.first;
    std::vector<double> MSD = tm_pair.second;
    exportPlotData2CSV(
            csv_dir + "MSD.csv", MSD_time, MSD,
            "time", "MSD");

    return 0;
}
