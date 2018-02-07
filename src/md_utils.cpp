#include <cmath>
#include <vector>
#include <utility>
// DEBUG
#include <iostream>
#include <iomanip>

#define EIGEN_NO_DEBUG
#define EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_MPL2_ONLY
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "md_utils.hpp"

// Eigen::Matrix<double, 3, 256> 
std::vector<Eigen::Vector3d> arrangeParticlesInFCCL(
        const double lattice_const, Eigen::Vector3i cube_size) {
    const double side_length = lattice_const * std::sqrt(2);
    Eigen::Matrix<double, 3, 4> pos_criteria;
    pos_criteria << 0., 0., side_length/2., side_length/2.,
        0., side_length/2., side_length/2., 0.,
        0., side_length/2., 0., side_length/2.;
    const int cols_criteria = pos_criteria.cols();
    const int num_particles = cols_criteria
        * cube_size.x() * cube_size.y() * cube_size.z();
    std::vector<Eigen::Vector3d> ptcls_pos(num_particles);

    int n = 0;
    for (int ix = 0; ix < cube_size.x(); ix++) {
        for (int iy = 0; iy < cube_size.y(); iy++) {
            for (int iz = 0; iz < cube_size.z(); iz++) {
                for (int col = 0; col < cols_criteria; col++) {
                    ptcls_pos[n].x() = pos_criteria(1, col) + (ix-1)*side_length;
                    ptcls_pos[n].y() = pos_criteria(2, col) + (iy-1)*side_length;
                    ptcls_pos[n].z() = pos_criteria(3, col) + (iz-1)*side_length;
                    n++;
                }
            }
        }
    }
    return ptcls_pos;
}

double calcWholeKineticEnergy(
        std::vector<Eigen::Vector3d> ptcls_velocity, const double ptcl_mass) {
    const int num_particles = ptcls_velocity.size();
    double sq_v_sum = 0.;
    for (const auto& el : ptcls_velocity) {
        sq_v_sum += el.x()*el.x() + el.y()*el.y() + el.z()*el.z();
    }
    return 1/2. * ptcl_mass * sq_v_sum;
}

std::vector<Eigen::Vector3d> initVelocity(
        const int num_particles, const double ptcl_mass, const double target_temp) {
    std::vector<Eigen::Vector3d> ptcls_velocity(num_particles);
    Eigen::Vector3d ini_v_sum = Eigen::Vector3d::Zero();
    for (auto& el : ptcls_velocity) {
        el = Eigen::Vector3d::Random()*2;
        el = el.array() - 1.;
        ini_v_sum += el;
    }

    Eigen::Vector3d v_average = ini_v_sum / num_particles;
    for (auto& el : ptcls_velocity) {
        el -= v_average;
    }

    ptcls_velocity = controlTempByScalingVel(ptcls_velocity, ptcl_mass, target_temp);
    // DEBUG
    const double kinetic_energy = calcWholeKineticEnergy(
            ptcls_velocity, ptcl_mass) / num_particles;
    const double current_temp = 2/3. * kinetic_energy;
    std::cout << "current temp: " << 
            std::setprecision(3) << current_temp << std::endl;
    return ptcls_velocity;
}

std::vector<Eigen::Vector3d> controlTempByScalingVel(
        std::vector<Eigen::Vector3d> ptcls_velocity,
        const double ptcl_mass, const double target_temp) {
    const int num_particles = ptcls_velocity.size();
    const double kinetic_energy = calcWholeKineticEnergy(
            ptcls_velocity, ptcl_mass) / num_particles;
    const double current_temp = 2/3. * kinetic_energy;
    const double scale_coeff = std::sqrt(target_temp / current_temp);
    for (auto& ptcl_vel : ptcls_velocity) {
        ptcl_vel *= scale_coeff;
    }
    return ptcls_velocity;
}

std::pair<double, std::vector<Eigen::Vector3d>> calcLJPotentialAndForce(
        std::vector<Eigen::Vector3d> ptcls_pos,
        Eigen::Vector3d volume, std::string bc_mode) {
    const int num_particles = ptcls_pos.size();
    const double epsilon = 1.;
    const double sigma = 1.;
    double potential = 0.;
    std::vector<Eigen::Vector3d> ptcls_force(num_particles);

    Eigen::Vector3d pos_diff, trunc, force, tmp;
    for (int j = 1; j < num_particles; j++) {
        for (int i = 0; i < j-1; i++) {
            pos_diff = ptcls_pos[i] - ptcls_pos[j];
            if (bc_mode == "periodic") {
                trunc = 2. * pos_diff.array() / volume.array();
                trunc.x() = std::trunc(trunc.x());
                trunc.y() = std::trunc(trunc.y());
                trunc.z() = std::trunc(trunc.z());
                tmp = trunc.array() * volume.array();
                pos_diff -= tmp;
            } else if (bc_mode == "free") {
                continue;
            } else {
                std::cout << 
                    "Error: bc_mode must be 'periodic' or 'free'." << std::endl;
                exit(1);
            }
            double ptcls_distance = pos_diff.norm();

            potential += 4. * epsilon * (
                    std::pow(sigma/ptcls_distance, 12)
                    - std::pow(sigma/ptcls_distance, 6)
                    );

            double force_coeff = 4. * epsilon * (
                    12 * std::pow(sigma/ptcls_distance, 12)
                    - 6 * std::pow(sigma/ptcls_distance, 6)
                    ) / ptcls_distance;
            force = force_coeff * pos_diff / ptcls_distance;
            ptcls_force[i] += force;
            ptcls_force[j] -= force;
        }
    }
    return std::make_pair(potential, ptcls_force);
}

std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> 
manageBoundaryCollision(
        std::vector<Eigen::Vector3d> ptcls_pos,
        Eigen::Vector3d volume, std::string bc_mode) {
    const int num_particles = ptcls_pos.size();
    std::vector<Eigen::Vector3d> bc_count(num_particles);
    if (bc_mode == "periodic") {
        Eigen::Vector3d tmp;
        for (auto& ptcl_pos : ptcls_pos) {
            size_t idx = &ptcl_pos - &ptcls_pos[0];
            bc_count[idx] = ptcl_pos.array() / volume.array();
            bc_count[idx] = bc_count[idx].array().floor();
            tmp = bc_count[idx].array() * volume.array();
            ptcl_pos -= tmp;
        }
    } else if (bc_mode == "free") {
        for (auto& el : bc_count) {
            el = Eigen::Vector3d::Zero();
        }
    } else {
        std::cout << 
            "Error: bc_mode must be 'periodic' or 'free'." << std::endl;
        exit(1);
    }
    return std::make_pair(ptcls_pos, bc_count);
}
