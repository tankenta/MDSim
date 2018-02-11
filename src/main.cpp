#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#define EIGEN_NO_DEBUG
#define EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
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

    if (temp_cont_time > total_time) {
        throw std::invalid_argument("total_time must be longer than temp_cont_time.");
    }
    bool control_temp = temp_cont_time > 0;
    std::vector<double> times = generateRange(
            -temp_cont_time+dt, dt, total_time-temp_cont_time);
    const int num_steps = times.size();
    std::vector<int> steps = generateRange(0, 1, num_steps-1);
    std::vector<bool> times_mask = makeVecMask(times, LESR, 0.);
    const int nonneg_step_offset = maskedVec(times, times_mask).size();
    const int num_nonneg_steps = num_steps - nonneg_step_offset;

    const double target_temp = 1.0;
    const int num_particles = 256;
    const double ptcl_mass = 1.0;
    const double particle_radius = 0.2;
    const int RDF_hist_size = 200;

    Eigen::Vector3i cube_size = Eigen::Vector3i::Constant(4);
    const double volume_scalar = num_particles / number_density;
    Eigen::Vector3d volume = Eigen::Vector3d::Constant(std::pow(volume_scalar, 1/3.));

    const double lattice_const = std::pow(volume_scalar, 1/3.)/8.*std::sqrt(2.);
    std::vector<Eigen::Vector3d> ptcls_pos = 
        arrangeParticlesInFCCL(lattice_const, cube_size);
    std::vector<std::vector<Eigen::Vector3d>> ptcls_fpos_allst(
            num_nonneg_steps,
            std::vector<Eigen::Vector3d>(num_particles, Eigen::Vector3d::Zero()));
    std::vector<Eigen::Vector3d> ptcls_velocity = initVelocity(
            num_particles, ptcl_mass, target_temp);
    auto pf_pair = calcLJPotentialAndForce(ptcls_pos, volume, bc_mode);
    std::vector<Eigen::Vector3d> prev_force = pf_pair.second;
    std::vector<Eigen::Vector3d> bc_count(num_particles, Eigen::Vector3d::Zero());
    std::vector<Eigen::Vector3d> next_force(num_particles, Eigen::Vector3d::Zero());
    std::vector<Eigen::Vector3d> bc_count_sum(num_particles, Eigen::Vector3d::Zero());

    std::vector<double> potential(num_steps, 0.);
    std::vector<double> kinetic_energy(num_steps, 0.);
    std::vector<double> current_temp(num_steps, 0.);

    std::cout << "running the simulation..." << std::endl;
    ProgressBar prog(steps.front(), steps.back());
    for (const auto& step : steps) {
        prog.printProgressBar(step);
        for (const auto& ptcl_vel : ptcls_velocity) {
            size_t idx = &ptcl_vel - &ptcls_velocity[0];
            ptcls_pos[idx] += ptcl_vel*dt + prev_force[idx]/(2.*ptcl_mass)*dt*dt;
        }
        auto pc_pair = manageBoundaryCollision(ptcls_pos, volume, bc_mode);
        ptcls_pos = pc_pair.first;
        bc_count = pc_pair.second;
        if (times[step] >= 0.) {
            for (const auto& bcc : bc_count) {
                size_t idx = &bcc - &bc_count[0];
                bc_count_sum[idx] += bcc;
            }
            size_t step_idx = step - nonneg_step_offset;
            for (const auto& ptcl_pos : ptcls_pos) {
                size_t ptcl_idx = &ptcl_pos - &ptcls_pos[0];
                ptcls_fpos_allst[step_idx][ptcl_idx] = ptcl_pos
                    + (bc_count_sum[ptcl_idx].array() * volume.array()).matrix();
            }
        }

        pf_pair = calcLJPotentialAndForce(ptcls_pos, volume, bc_mode);
        potential[step] = pf_pair.first;
        next_force = pf_pair.second;

        kinetic_energy[step] = calcWholeKineticEnergy(ptcls_velocity, ptcl_mass);
        current_temp[step] = 2/3. * kinetic_energy[step]/num_particles;

        for (auto&& ptcl_vel : ptcls_velocity) {
            size_t idx = &ptcl_vel - &ptcls_velocity[0];
            ptcl_vel += dt/(2.*ptcl_mass)*(next_force[idx] + prev_force[idx]);
        }

        if (control_temp) {
            ptcls_velocity = controlTempByScalingVel(
                    ptcls_velocity, ptcl_mass, target_temp);
            if (times[step] > 0.) {
                control_temp = false;
                // TODO: output equilibrium_ptcl_pos to plot
            }
        }

        prev_force = next_force;
    }
    prog.finish();

    std::vector<double> total_energy(num_steps);
    for (const auto& kinetic : kinetic_energy) {
        size_t idx = &kinetic - &kinetic_energy[0];
        total_energy[idx] = potential[idx] + kinetic;
    }
    exportPlotData2CSV(
            csv_dir + "total_energy.csv", times, total_energy,
            "time", "total energy");

    exportPlotData2CSV(
            csv_dir + "potential_energy.csv", times, potential,
            "time", "potential energy");

    exportPlotData2CSV(
            csv_dir + "kinetic_energy.csv", times, kinetic_energy,
            "time", "kinetic energy");

    exportPlotData2CSV(
            csv_dir + "temperature.csv", times, current_temp,
            "time", "temperature");

    auto dr_pair = calcRadialDistributionFunction(
            ptcls_pos, volume, number_density, RDF_hist_size, bc_mode);
    std::vector<double> distance = dr_pair.first;
    std::vector<double> RDF = dr_pair.second;
    exportPlotData2CSV(
            csv_dir + "RDF.csv", distance, RDF,
            "r", "g(r)");

    auto tm_pair = calcMeanSquareDisplacement(ptcls_fpos_allst, dt);
    std::vector<double> MSD_time = tm_pair.first;
    std::vector<double> MSD = tm_pair.second;
    exportPlotData2CSV(
            csv_dir + "MSD.csv", MSD_time, MSD,
            "time", "MSD");

    return 0;
}
