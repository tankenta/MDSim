#pragma once

#include <vector>
#include <utility>
#include <string>

#define EIGEN_NO_DEBUG
#define EIGEN_MPL2_ONLY
#include <Eigen/Core>
#include <Eigen/Geometry>

class MDSim {
public:
    MDSim(
            double dt, double total_time, double temp_cont_time,
            double number_density, std::string bc_mode, double target_temp,
            int num_particles, double ptcl_mass, int RDF_hist_size,
            const Eigen::Vector3i& cube_size);
    void manageBoundaryCollision();
    double calcLJPotentialAndForce(Eigen::MatrixXd& ptcls_force);
    void next(int step);
    std::pair<std::vector<double>, std::vector<double>> calcRadialDistributionFunction();
    std::pair<std::vector<double>, std::vector<double>> calcMeanSquareDisplacement();

    std::vector<double> times;
    std::vector<int> steps;
    Eigen::MatrixXd ptcls_pos;
    std::vector<double> potential_arr, kinetic_energy_arr,
            total_energy_arr, current_temp_arr;

private:
    Eigen::MatrixXd arrangeParticlesInFCCL(
            double lattice_const, const Eigen::Vector3i& cube_size);
    Eigen::MatrixXd initVelocity();
    std::vector<int> calcHistogram(
            const std::vector<double>& src_arr, int hist_size, 
            const std::pair<double, double>& hist_range);

    const double dt, total_time, temp_cont_time, number_density, target_temp, ptcl_mass;
    const int num_particles, RDF_hist_size;
    int num_steps, nonneg_step_offset; 
    std::string bc_mode;
    bool control_temp;
    Eigen::Vector3d volume;
    Eigen::MatrixXd volume_vecs;
    std::vector<Eigen::MatrixXd> ptcls_fpos_allst;
    Eigen::MatrixXd ptcls_velocity;
    Eigen::MatrixXd prev_force, next_force, bc_count, bc_count_sum;

// inline
public:
    void renewPtclsPos() {
        ptcls_pos += ptcls_velocity*dt + prev_force/(2.*ptcl_mass)*dt*dt;
        return;
    }

    void controlTempByScalingVel() {
        const double kinetic_energy = calcWholeKineticEnergy() / num_particles;
        const double current_temp = 2/3. * kinetic_energy;
        const double scale_coeff = std::sqrt(target_temp / current_temp);
        ptcls_velocity *= scale_coeff;
        return;
    }

    void renewPtclsFreePos(int step) {
        bc_count_sum += bc_count;
        size_t step_idx = step - nonneg_step_offset;
        ptcls_fpos_allst[step_idx] = ptcls_pos
            + (bc_count_sum.array() * volume_vecs.array()).matrix();
        return;
    }

    double calcWholeKineticEnergy() {
        double sq_v_sum = 0.;
        for (int i = 0; i < num_particles; i++) {
            sq_v_sum += ptcls_velocity.col(i).squaredNorm();
        }
        return 1/2. * ptcl_mass * sq_v_sum;
    }

    void renewPtclsVel() {
        ptcls_velocity += dt/(2.*ptcl_mass)*(next_force + prev_force);
        return;
    }

    void calcTotalEnergy() {
        for (const auto& kinetic : kinetic_energy_arr) {
            size_t idx = &kinetic - &kinetic_energy_arr[0];
            total_energy_arr[idx] = potential_arr[idx] + kinetic;
        }
        return;
    }

private:
    std::vector<double> vecZeros(int size) {
        std::vector<double> v(size, 0.);
        return v;
    }
};

