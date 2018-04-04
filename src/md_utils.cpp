#include <cmath>
#include <iostream>
#include <random>
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

MDSim::MDSim(
        double dt, double total_time, double temp_cont_time,
        double number_density, std::string bc_mode_str, double target_temp,
        int num_particles, double ptcl_mass, int RDF_hist_size,
        const Eigen::Vector3i& cube_size)
: dt(dt), total_time(total_time), temp_cont_time(temp_cont_time),
        number_density(number_density), bc_mode(bc_mode),
        target_temp(target_temp), num_particles(num_particles), ptcl_mass(ptcl_mass),
        RDF_hist_size(RDF_hist_size) {
    if (temp_cont_time > total_time) {
        throw std::invalid_argument("total_time must be longer than temp_cont_time.");
    }
    control_temp = temp_cont_time > 0;
    times = VecUtils<double>::generateRange(-temp_cont_time+dt, dt, total_time-temp_cont_time);
    num_steps = times.size();
    steps = VecUtils<int>::generateRange(0, 1, num_steps-1);
    std::vector<bool> times_mask = VecUtils<double>::makeVecMask(times, CompOp::LESS, 0.);
    nonneg_step_offset = VecUtils<double>::maskedVec(times, times_mask).size();
    const int num_nonneg_steps = num_steps - nonneg_step_offset;
    if (bc_mode_str == "periodic") {
        bc_mode = BCMode::Periodic;
    } else if (bc_mode_str == "free") {
        bc_mode = BCMode::Free;
    } else {
        throw std::invalid_argument("bc_mode_str must be 'periodic' or 'free'.");
    }

    const double volume_scalar = num_particles / number_density;
    volume = Eigen::Vector3d::Constant(std::pow(volume_scalar, 1/3.));
    volume_vecs = volume.replicate(1, num_particles);

    const double lattice_const = std::pow(volume_scalar, 1/3.)/8.*std::sqrt(2.);
    ptcls_pos = arrangeParticlesInFCCL(lattice_const, cube_size);
    std::vector<Eigen::MatrixXd> tmp(
            num_nonneg_steps, Eigen::MatrixXd::Zero(3, num_particles));
    ptcls_fpos_allst = tmp;
    ptcls_velocity = initVelocity();
    prev_force = Eigen::MatrixXd::Zero(3, num_particles);
    calcLJPotentialAndForce(prev_force);
    bc_count = Eigen::MatrixXd::Zero(3, num_particles);
    next_force = Eigen::MatrixXd::Zero(3, num_particles);
    bc_count_sum = Eigen::MatrixXd::Zero(3, num_particles);

    potential_arr = vecZeros(num_steps);
    kinetic_energy_arr = vecZeros(num_steps);
    total_energy_arr = vecZeros(num_steps);
    current_temp_arr = vecZeros(num_steps);
}

Eigen::MatrixXd MDSim::arrangeParticlesInFCCL(
        double lattice_const, const Eigen::Vector3i& cube_size) {
    const double side_length = lattice_const * std::sqrt(2);
    Eigen::Matrix<double, 3, 4> pos_criteria;
    pos_criteria << 0., 0., side_length/2., side_length/2.,
        0., side_length/2., side_length/2., 0.,
        0., side_length/2., 0., side_length/2.;
    const int cols_criteria = pos_criteria.cols();
    const int num_particles = cols_criteria
        * cube_size.x() * cube_size.y() * cube_size.z();
    Eigen::MatrixXd ini_ptcls_pos = Eigen::MatrixXd::Zero(3, num_particles);

    int n = 0;
    for (int ix = 0; ix < cube_size.x(); ix++) {
        for (int iy = 0; iy < cube_size.y(); iy++) {
            for (int iz = 0; iz < cube_size.z(); iz++) {
                for (int col = 0; col < cols_criteria; col++) {
                    ini_ptcls_pos(0, n) = pos_criteria(0, col) + ix*side_length;
                    ini_ptcls_pos(1, n) = pos_criteria(1, col) + iy*side_length;
                    ini_ptcls_pos(2, n) = pos_criteria(2, col) + iz*side_length;
                    n++;
                }
            }
        }
    }
    return ini_ptcls_pos;
}

Eigen::MatrixXd MDSim::initVelocity() {
    std::srand((unsigned int) time(0));
    ptcls_velocity = Eigen::MatrixXd::Random(3, num_particles);

    Eigen::Vector3d v_average = ptcls_velocity.rowwise().mean();
    ptcls_velocity -= v_average.replicate(1, num_particles);
    controlTempByScalingVel();
    return ptcls_velocity;
}

void MDSim::manageBoundaryCollision() {
    switch (bc_mode) {
    case BCMode::Periodic:
        bc_count = (ptcls_pos.array() / volume_vecs.array()).floor();
        ptcls_pos -= (bc_count.array() * volume_vecs.array()).matrix();
        break;
    case BCMode::Free:
        bc_count = Eigen::MatrixXd::Zero(3, num_particles);
        break;
    }
    return;
}

double MDSim::calcLJPotentialAndForce(Eigen::MatrixXd& ptcls_force) {
    const double epsilon = 1.;
    const double sigma = 1.;
    double potential = 0.;

    ptcls_force = Eigen::MatrixXd::Zero(3, num_particles);
    Eigen::Vector3d pos_diff = Eigen::Vector3d::Zero();
    Eigen::Vector3d tmp_trunc = Eigen::Vector3d::Zero();
    Eigen::Vector3d force = Eigen::Vector3d::Zero();
    for (int j = 1; j < num_particles; j++) {
        for (int i = 0; i < j; i++) {
            pos_diff = ptcls_pos.col(i) - ptcls_pos.col(j);
            switch (bc_mode) {
            case BCMode::Periodic:
                tmp_trunc = 2. * pos_diff.array() / volume.array();
                tmp_trunc.x() = std::trunc(tmp_trunc.x());
                tmp_trunc.y() = std::trunc(tmp_trunc.y());
                tmp_trunc.z() = std::trunc(tmp_trunc.z());
                pos_diff -= (tmp_trunc.array() * volume.array()).matrix();
                break;
            case BCMode::Free:
                break;
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
            Eigen::Vector3d force = force_coeff * pos_diff / ptcls_distance;
            ptcls_force.col(i) += force;
            ptcls_force.col(j) -= force;
        }
    }
    return potential;
}

void MDSim::next(int step) {
    renewPtclsPos();
    manageBoundaryCollision();
    if (times[step] >= 0.) {
        renewPtclsFreePos(step);
    }

    potential_arr[step] = calcLJPotentialAndForce(next_force);
    kinetic_energy_arr[step] = calcWholeKineticEnergy();
    current_temp_arr[step] = 2/3. * kinetic_energy_arr[step]/num_particles;

    renewPtclsVel();
    if (control_temp) {
        controlTempByScalingVel();
        if (times[step] > 0.) {
            control_temp = false;
            // TODO: output equilibrium_ptcl_pos to plot
        }
    }

    prev_force = next_force;
    return;
}

std::vector<int> MDSim::calcHistogram(
        const std::vector<double>& src_arr, int hist_size, 
        const std::pair<double, double>& hist_range) {
    std::vector<int> hist(hist_size, 0);
    const double bin_width = (hist_range.second - hist_range.first) / hist_size;
    std::vector<int> bin_idxs(src_arr.size(), 0);
    for (const auto& src_el : src_arr) {
        size_t arr_idx = &src_el - &src_arr[0];
        bin_idxs[arr_idx] = std::floor(src_el/bin_width) + 1;
    }

    for (auto&& bin_idx : bin_idxs) {
        bin_idx = std::min(bin_idx, hist_size);
        hist[bin_idx]++;
    }
    return hist;
}

std::pair<std::vector<double>, std::vector<double>> MDSim::calcRadialDistributionFunction() {
    const double pi = 3.14159265358979;
    int ptcls_dist_dim = 0;
    for (int i = 1; i < num_particles; i++) {
        ptcls_dist_dim += i;
    }
    std::vector<double> ptcls_distance(ptcls_dist_dim, 0.);

    Eigen::Vector3d pos_diff, tmp_trunc;
    int dist_idx = 0;
    for (int j = 1; j < num_particles; j++) {
        for (int i = 0; i < j; i++) {
            pos_diff = ptcls_pos.col(i) - ptcls_pos.col(j);
            switch (bc_mode) {
            case BCMode::Periodic:
                tmp_trunc = 2 * pos_diff.array() / volume.array();
                tmp_trunc.x() = std::trunc(tmp_trunc.x());
                tmp_trunc.y() = std::trunc(tmp_trunc.y());
                tmp_trunc.z() = std::trunc(tmp_trunc.z());
                pos_diff -= (tmp_trunc.array() * volume.array()).matrix();
                break;
            case BCMode::Free:
                break;
            }
            ptcls_distance[dist_idx] = pos_diff.norm();
            dist_idx++;
        }
    }

    double max_dist = *std::max_element(ptcls_distance.begin(), ptcls_distance.end());
    std::pair<double, double> hist_range = std::make_pair(0., 1.1*max_dist);
    std::vector<int> hist = calcHistogram(ptcls_distance, RDF_hist_size, hist_range);
    std::vector<double> dist_distn_mean(hist.size(), 0.);
    for (const auto& hist_el : hist) {
        size_t idx = &hist_el - &hist[0];
        dist_distn_mean[idx] = 2. * hist_el / num_particles;
    }
    double bin_width = (hist_range.second - hist_range.first) / RDF_hist_size;
    std::vector<double> dist(dist_distn_mean.size(), 0.);
    std::vector<double> RDF(dist_distn_mean.size(), 0.);
    for (const auto& mean_el : dist_distn_mean) {
        size_t idx = &mean_el - &dist_distn_mean[0];
        dist[idx] = (hist_range.first + idx) * bin_width;
        RDF[idx] = mean_el / (number_density*4.*pi*dist[idx]*dist[idx]*bin_width);
    }
    return std::make_pair(dist, RDF);
}

std::pair<std::vector<double>, std::vector<double>> MDSim::calcMeanSquareDisplacement() {
    const int num_nonneg_steps = ptcls_fpos_allst.size();
    std::vector<double> nonneg_time = VecUtils<double>::generateRange(
            1., 1., num_nonneg_steps-1.);
    for (auto&& el : nonneg_time) {
        el *= dt;
    }

    Eigen::MatrixXd tmp_MSD(num_nonneg_steps-1, num_particles);
    tmp_MSD = Eigen::MatrixXd::Zero(num_nonneg_steps-1, num_particles);
    Eigen::Vector3d displacement;
    std::cout << "calculating MSD..." << std::endl;
    ProgressBar prog(1, num_nonneg_steps-1);
    for (int t = 1; t < num_nonneg_steps; t++) {
        prog.printProgressBar(t);
        for (int t0 = 0; t0 < t; t0++) {
            for (int p = 0; p < num_particles; p++) {
                displacement = ptcls_fpos_allst[t].col(p) - ptcls_fpos_allst[t0].col(p);
                tmp_MSD(t-t0-1, p) += displacement.squaredNorm();
            }
        }
    }

    Eigen::VectorXd MSD = tmp_MSD.rowwise().mean();
    for (int i = 0; i < num_nonneg_steps-1; i++) {
        MSD(i) /= num_nonneg_steps-1. - i;
    }
    prog.finish();
    std::vector<double> ret_MSD(num_nonneg_steps-1, 0.);
    Eigen::Map<Eigen::VectorXd>(&ret_MSD[0], ret_MSD.size()) = MSD;
    return std::make_pair(nonneg_time, ret_MSD);
}
