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
    void renewPtclsPos();
    void manageBoundaryCollision();
    void renewPtclsFreePos(int step);
    double calcLJPotentialAndForce(Eigen::MatrixXd& ptcls_force);
    double calcWholeKineticEnergy();
    void renewPtclsVel();
    void next(int step);
    void calcTotalEnergy();
    std::pair<std::vector<double>, std::vector<double>> calcRadialDistributionFunction();
    std::pair<std::vector<double>, std::vector<double>> calcMeanSquareDisplacement();

    std::vector<double> times;
    std::vector<int> steps;
    Eigen::MatrixXd ptcls_pos;
    std::vector<double> potential_arr, kinetic_energy_arr,
            total_energy_arr, current_temp_arr;

private:
    std::vector<double> vecZeros(int size);
    Eigen::MatrixXd arrangeParticlesInFCCL(
            double lattice_const, const Eigen::Vector3i& cube_size);
    Eigen::MatrixXd initVelocity();
    void controlTempByScalingVel();
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
};

