#include <cmath>
#include <iostream>
#include <random>
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

MDSim::MDSim(
        double dt, double total_time, double temp_cont_time,
        double number_density, std::string bc_mode, double target_temp,
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
    times = generateRange(-temp_cont_time+dt, dt, total_time-temp_cont_time);
    num_steps = times.size();
    steps = generateRange(0, 1, num_steps-1);
    std::vector<bool> times_mask = makeVecMask(times, LESR, 0.);
    nonneg_step_offset = maskedVec(times, times_mask).size();
    const int num_nonneg_steps = num_steps - nonneg_step_offset;

    const double volume_scalar = num_particles / number_density;
    volume = Eigen::Vector3d::Constant(std::pow(volume_scalar, 1/3.));

    const double lattice_const = std::pow(volume_scalar, 1/3.)/8.*std::sqrt(2.);
    ptcls_pos = arrangeParticlesInFCCL(lattice_const, cube_size);
    std::vector<std::vector<Eigen::Vector3d>> tmp(
            num_nonneg_steps,
            std::vector<Eigen::Vector3d>(num_particles, Eigen::Vector3d::Zero()));
    ptcls_fpos_allst = tmp;
    ptcls_velocity = fillVecWithZeros(num_particles);
    initVelocity();
    prev_force = fillVecWithZeros(num_particles);
    calcLJPotentialAndForce(prev_force);
    bc_count = fillVecWithZeros(num_particles);
    next_force = fillVecWithZeros(num_particles);
    bc_count_sum = fillVecWithZeros(num_particles);

    potential_arr = fillVecWithZero(num_steps);
    kinetic_energy_arr = fillVecWithZero(num_steps);
    total_energy_arr = fillVecWithZero(num_steps);
    current_temp_arr = fillVecWithZero(num_steps);
    // printIniVal();   // DEBUG
}

void MDSim::printIniVal() {
    std::cout << "times len  : " << times.size() << std::endl;
    std::cout << "times front: " << times.front() << std::endl;
    std::cout << "times back : " << times.back() << std::endl;

    std::cout << "steps len  : " << steps.size() << std::endl;
    std::cout << "steps front: " << steps.front() << std::endl;
    std::cout << "steps back : " << steps.back() << std::endl;

    std::cout << "p_pos len  : " << ptcls_pos.size() << std::endl;
    std::cout << "p_pos front: " << ptcls_pos.front() << std::endl;
    std::cout << "p_pos back : " << ptcls_pos.back() << std::endl;

    std::cout << "poten len  : " << potential_arr.size() << std::endl;
    std::cout << "poten front: " << potential_arr.front() << std::endl;
    std::cout << "poten back : " << potential_arr.back() << std::endl;

    std::cout << "kinet len  : " << kinetic_energy_arr.size() << std::endl;
    std::cout << "kinet front: " << kinetic_energy_arr.front() << std::endl;
    std::cout << "kinet back : " << kinetic_energy_arr.back() << std::endl;

    std::cout << "total len  : " << total_energy_arr.size() << std::endl;
    std::cout << "total front: " << total_energy_arr.front() << std::endl;
    std::cout << "total back : " << total_energy_arr.back() << std::endl;

    std::cout << "c_temp len  : " << current_temp_arr.size() << std::endl;
    std::cout << "c_temp front: " << current_temp_arr.front() << std::endl;
    std::cout << "c_temp back : " << current_temp_arr.back() << std::endl;

    std::cout << "p_force len  : " << prev_force.size() << std::endl;
    std::cout << "p_force front: " << prev_force.front() << std::endl;
    std::cout << "p_force back : " << prev_force.back() << std::endl;

    std::cout << "n_force len  : " << next_force.size() << std::endl;
    std::cout << "n_force front: " << next_force.front() << std::endl;
    std::cout << "n_force back : " << next_force.back() << std::endl;

    std::cout << "bcc len  : " << bc_count.size() << std::endl;
    std::cout << "bcc front: " << bc_count.front() << std::endl;
    std::cout << "bcc back : " << bc_count.back() << std::endl;

    std::cout << "bcs len  : " << bc_count_sum.size() << std::endl;
    std::cout << "bcs front: " << bc_count_sum.front() << std::endl;
    std::cout << "bcs back : " << bc_count_sum.back() << std::endl;

    std::cout << "p_vel len  : " << ptcls_velocity.size() << std::endl;
    std::cout << "p_vel front: " << ptcls_velocity.front() << std::endl;
    std::cout << "p_vel back : " << ptcls_velocity.back() << std::endl;

    std::cout << "f_pos len    : " << ptcls_fpos_allst.size() << std::endl;
    std::cout << "f_pos fr len : " << ptcls_fpos_allst.front().size() << std::endl;
    std::cout << "f_pos fr fr  : " << ptcls_fpos_allst.front().front() << std::endl;
    std::cout << "f_pos bc bc  : " << ptcls_fpos_allst.back().back() << std::endl;

    std::cout << "dt : " << dt << std::endl;
    std::cout << "total_time : " << total_time << std::endl;
    std::cout << "temp_cont_time : " << temp_cont_time << std::endl;
    std::cout << "number_density : " << number_density << std::endl;
    std::cout << "target_temp : " << target_temp << std::endl;
    std::cout << "ptcl_mass : " << ptcl_mass << std::endl;
    std::cout << "num_particles : " << num_particles << std::endl;
    std::cout << "RDF_hist_size : " << RDF_hist_size << std::endl;
    std::cout << "num_steps : " << num_steps << std::endl;
    std::cout << "nonneg_step_offset : " << nonneg_step_offset << std::endl;
    std::cout << "bc_mode : " << bc_mode << std::endl;
    std::cout << "control_temp : " << control_temp << std::endl;
    std::cout << "volume : " << volume << std::endl;
    return;
}

std::vector<Eigen::Vector3d> MDSim::fillVecWithZeros(int size) {
    std::vector<Eigen::Vector3d> v(size, Eigen::Vector3d::Zero());
    return v;
}

std::vector<double> MDSim::fillVecWithZero(int size) {
    std::vector<double> v(size, 0.);
    return v;
}

std::vector<Eigen::Vector3d> MDSim::arrangeParticlesInFCCL(
        double lattice_const, const Eigen::Vector3i& cube_size) {
    const double side_length = lattice_const * std::sqrt(2);
    Eigen::Matrix<double, 3, 4> pos_criteria;
    pos_criteria << 0., 0., side_length/2., side_length/2.,
        0., side_length/2., side_length/2., 0.,
        0., side_length/2., 0., side_length/2.;
    const int cols_criteria = pos_criteria.cols();
    const int num_particles = cols_criteria
        * cube_size.x() * cube_size.y() * cube_size.z();
    std::vector<Eigen::Vector3d> ini_ptcls_pos(num_particles, Eigen::Vector3d::Zero());

    int n = 0;
    for (int ix = 0; ix < cube_size.x(); ix++) {
        for (int iy = 0; iy < cube_size.y(); iy++) {
            for (int iz = 0; iz < cube_size.z(); iz++) {
                for (int col = 0; col < cols_criteria; col++) {
                    ini_ptcls_pos[n].x() = pos_criteria(0, col) + ix*side_length;
                    ini_ptcls_pos[n].y() = pos_criteria(1, col) + iy*side_length;
                    ini_ptcls_pos[n].z() = pos_criteria(2, col) + iz*side_length;
                    n++;
                }
            }
        }
    }
    return ini_ptcls_pos;
}

void MDSim::initVelocity() {
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

    controlTempByScalingVel();
    return;
}

void MDSim::renewPtclsPos() {
    for (const auto& ptcl_vel : ptcls_velocity) {
        size_t idx = &ptcl_vel - &ptcls_velocity[0];
        ptcls_pos[idx] += ptcl_vel*dt + prev_force[idx]/(2.*ptcl_mass)*dt*dt;
    }
    return;
}

void MDSim::controlTempByScalingVel() {
    const double kinetic_energy = calcWholeKineticEnergy() / num_particles;
    const double current_temp = 2/3. * kinetic_energy;
    const double scale_coeff = std::sqrt(target_temp / current_temp);
    for (auto&& ptcl_vel : ptcls_velocity) {
        ptcl_vel *= scale_coeff;
    }
    return;
}

void MDSim::manageBoundaryCollision() {
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
    return;
}

void MDSim::renewPtclsFreePos(int step) {
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
    return;
}

double MDSim::calcLJPotentialAndForce(std::vector<Eigen::Vector3d>& ptcls_force) {
    const double epsilon = 1.;
    const double sigma = 1.;
    double potential = 0.;

    ptcls_force = fillVecWithZeros(num_particles);
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
    return potential;
}

double MDSim::calcWholeKineticEnergy() {
    double sq_v_sum = 0.;
    for (const auto& el : ptcls_velocity) {
        sq_v_sum += el.squaredNorm();
    }
    return 1/2. * ptcl_mass * sq_v_sum;
}

void MDSim::renewPtclsVel() {
    for (auto&& ptcl_vel : ptcls_velocity) {
        size_t idx = &ptcl_vel - &ptcls_velocity[0];
        ptcl_vel += dt/(2.*ptcl_mass)*(next_force[idx] + prev_force[idx]);
    }
    return;
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

void MDSim::calcTotalEnergy() {
    for (const auto& kinetic : kinetic_energy_arr) {
        size_t idx = &kinetic - &kinetic_energy_arr[0];
        total_energy_arr[idx] = potential_arr[idx] + kinetic;
    }
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
    std::vector<double> time = generateRange(1., 1., num_nonneg_steps-1.);
    for (auto&& el : time) {
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
                displacement = ptcls_fpos_allst[t][p] - ptcls_fpos_allst[t0][p];
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
    return std::make_pair(time, ret_MSD);
}
