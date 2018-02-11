#include <cmath>
#include <stdexcept>
#include <vector>
#include <utility>
#include <random>
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
#include "vec_utils.hpp"

std::vector<Eigen::Vector3d> arrangeParticlesInFCCL(
        double lattice_const, const Eigen::Vector3i& cube_size) {
    const double side_length = lattice_const * std::sqrt(2);
    Eigen::Matrix<double, 3, 4> pos_criteria;
    pos_criteria << 0., 0., side_length/2., side_length/2.,
        0., side_length/2., side_length/2., 0.,
        0., side_length/2., 0., side_length/2.;
    const int cols_criteria = pos_criteria.cols();
    const int num_particles = cols_criteria
        * cube_size.x() * cube_size.y() * cube_size.z();
    std::vector<Eigen::Vector3d> ptcls_pos(num_particles, Eigen::Vector3d::Zero());

    int n = 0;
    for (int ix = 0; ix < cube_size.x(); ix++) {
        for (int iy = 0; iy < cube_size.y(); iy++) {
            for (int iz = 0; iz < cube_size.z(); iz++) {
                for (int col = 0; col < cols_criteria; col++) {
                    ptcls_pos[n].x() = pos_criteria(0, col) + ix*side_length;
                    ptcls_pos[n].y() = pos_criteria(1, col) + iy*side_length;
                    ptcls_pos[n].z() = pos_criteria(2, col) + iz*side_length;
                    n++;
                }
            }
        }
    }
    return ptcls_pos;
}

double calcWholeKineticEnergy(
        const std::vector<Eigen::Vector3d>& ptcls_velocity, double ptcl_mass) {
    const int num_particles = ptcls_velocity.size();
    double sq_v_sum = 0.;
    for (const auto& el : ptcls_velocity) {
        sq_v_sum += el.squaredNorm();
    }
    return 1/2. * ptcl_mass * sq_v_sum;
}

std::vector<Eigen::Vector3d> initVelocity(
        int num_particles, double ptcl_mass, double target_temp) {
    std::vector<Eigen::Vector3d> ptcls_velocity(num_particles);
    Eigen::Vector3d ini_v_sum = Eigen::Vector3d::Zero();
    std::srand((unsigned int) time(0));
    for (auto&& el : ptcls_velocity) {
        el = Eigen::Vector3d::Random();
        ini_v_sum += el;
    }

    Eigen::Vector3d v_average = ini_v_sum / num_particles;
    for (auto&& el : ptcls_velocity) {
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
        std::vector<Eigen::Vector3d>& ptcls_velocity,
        double ptcl_mass, double target_temp) {
    const int num_particles = ptcls_velocity.size();
    const double kinetic_energy = calcWholeKineticEnergy(
            ptcls_velocity, ptcl_mass) / num_particles;
    const double current_temp = 2/3. * kinetic_energy;
    const double scale_coeff = std::sqrt(target_temp / current_temp);
    for (auto&& ptcl_vel : ptcls_velocity) {
        ptcl_vel *= scale_coeff;
    }
    return ptcls_velocity;
}

std::pair<double, std::vector<Eigen::Vector3d>> calcLJPotentialAndForce(
        const std::vector<Eigen::Vector3d>& ptcls_pos,
        const Eigen::Vector3d& volume, const std::string& bc_mode) {
    const int num_particles = ptcls_pos.size();
    const double epsilon = 1.;
    const double sigma = 1.;
    double potential = 0.;
    std::vector<Eigen::Vector3d> ptcls_force(num_particles, Eigen::Vector3d::Zero());

    Eigen::Vector3d pos_diff = Eigen::Vector3d::Zero();
    Eigen::Vector3d tmp_trunc = Eigen::Vector3d::Zero();
    for (int j = 1; j < num_particles; j++) {
        for (int i = 0; i < j; i++) {
            pos_diff = ptcls_pos[i] - ptcls_pos[j];
            if (bc_mode == "periodic") {
                tmp_trunc = 2. * pos_diff.array() / volume.array();
                tmp_trunc.x() = std::trunc(tmp_trunc.x());
                tmp_trunc.y() = std::trunc(tmp_trunc.y());
                tmp_trunc.z() = std::trunc(tmp_trunc.z());
                pos_diff -= (tmp_trunc.array() * volume.array()).matrix();
            } else if (bc_mode == "free") {
                continue;
            } else {
                throw std::invalid_argument("bc_mode must be 'periodic' or 'free'.");
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
            ptcls_force[i] += force;
            ptcls_force[j] -= force;
        }
    }
    return std::make_pair(potential, ptcls_force);
}

std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> 
manageBoundaryCollision(
        std::vector<Eigen::Vector3d>& ptcls_pos,
        const Eigen::Vector3d& volume, const std::string& bc_mode) {
    const int num_particles = ptcls_pos.size();
    std::vector<Eigen::Vector3d> bc_count(num_particles, Eigen::Vector3d::Zero());
    if (bc_mode == "periodic") {
        for (auto&& ptcl_pos : ptcls_pos) {
            size_t idx = &ptcl_pos - &ptcls_pos[0];
            bc_count[idx] = (ptcl_pos.array() / volume.array()).floor();
            ptcl_pos -= (bc_count[idx].array() * volume.array()).matrix();
        }
    } else if (bc_mode == "free") {
        for (auto&& el : bc_count) {
            el = Eigen::Vector3d::Zero();
        }
    } else {
        throw std::invalid_argument("bc_mode must be 'periodic' or 'free'.");
    }
    return std::make_pair(ptcls_pos, bc_count);
}

std::vector<int> calcHistogram(
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

std::pair<std::vector<double>, std::vector<double>> calcRadialDistributionFunction(
        const std::vector<Eigen::Vector3d>& ptcls_pos,
        const Eigen::Vector3d& volume,
        double number_density, int hist_size, const std::string bc_mode) {
    const double pi = 3.14159265358979;
    const int num_particles = ptcls_pos.size();
    int ptcls_dist_dim = 0;
    for (int i = 1; i < num_particles; i++) {
        ptcls_dist_dim += i;
    }
    std::vector<double> ptcls_distance(ptcls_dist_dim, 0.);

    Eigen::Vector3d pos_diff, tmp_trunc;
    int dist_idx = 0;
    for (int j = 1; j < num_particles; j++) {
        for (int i = 0; i < j; i++) {
            pos_diff = ptcls_pos[i] - ptcls_pos[j];
            if (bc_mode == "periodic") {
                tmp_trunc = 2 * pos_diff.array() / volume.array();
                tmp_trunc.x() = std::trunc(tmp_trunc.x());
                tmp_trunc.y() = std::trunc(tmp_trunc.y());
                tmp_trunc.z() = std::trunc(tmp_trunc.z());
                pos_diff -= (tmp_trunc.array() * volume.array()).matrix();
            } else if (bc_mode == "free") {
                continue;
            } else {
                throw std::invalid_argument("bc_mode must be 'periodic' or 'free'.");
            }
            ptcls_distance[dist_idx] = pos_diff.norm();
            dist_idx++;
        }
    }

    double max_dist = *std::max_element(ptcls_distance.begin(), ptcls_distance.end());
    std::pair<double, double> hist_range = std::make_pair(0., 1.1*max_dist);
    std::vector<int> hist = calcHistogram(ptcls_distance, hist_size, hist_range);
    std::vector<double> dist_distn_mean(hist.size(), 0.);
    for (const auto& hist_el : hist) {
        size_t idx = &hist_el - &hist[0];
        dist_distn_mean[idx] = 2. * hist_el / num_particles;
    }
    double bin_width = (hist_range.second - hist_range.first) / hist_size;
    std::vector<double> dist(dist_distn_mean.size(), 0.);
    std::vector<double> RDF(dist_distn_mean.size(), 0.);
    for (const auto& mean_el : dist_distn_mean) {
        size_t idx = &mean_el - &dist_distn_mean[0];
        dist[idx] = (hist_range.first + idx) * bin_width;
        RDF[idx] = mean_el / (number_density*4.*pi*dist[idx]*dist[idx]*bin_width);
    }
    return std::make_pair(dist, RDF);
}

std::pair<std::vector<double>, std::vector<double>>
calcMeanSquareDisplacement(
        const std::vector<std::vector<Eigen::Vector3d>>& ptcls_fpos_allst, double dt) {
    const int num_nonneg_steps = ptcls_fpos_allst.size();
    const int num_particles = ptcls_fpos_allst.front().size();
    std::vector<double> time = generateRange(1., 1., num_nonneg_steps-1.);
    for (auto&& el : time) {
        el *= dt;
    }
    Eigen::MatrixXd tmp_MSD(num_nonneg_steps-1, num_particles);
    tmp_MSD = Eigen::MatrixXd::Zero(num_nonneg_steps-1, num_particles);
    Eigen::Vector3d displacement;
    for (int t = 1; t < num_nonneg_steps; t++) {
        for (int t0 = 0; t0 < t; t0++) {
            for (int p = 0; p < num_particles; p++) {
                displacement = ptcls_fpos_allst[t][p] - ptcls_fpos_allst[t0][p];
                tmp_MSD(t-t0-1, p) += displacement.squaredNorm();
            }
        }
    }
    Eigen::VectorXd MSD = tmp_MSD.rowwise().mean();
    for (int i = 0; i < num_nonneg_steps-1; i++) {
        MSD(i) /= num_nonneg_steps-1. - i;
    }
    std::vector<double> ret_MSD(num_nonneg_steps-1, 0.);
    Eigen::Map<Eigen::VectorXd>(&ret_MSD[0], ret_MSD.size()) = MSD;
    return std::make_pair(time, ret_MSD);
}
