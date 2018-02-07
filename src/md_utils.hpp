#ifndef MD_UTILS_H_
#define MD_UTILS_H_

#include <vector>
#include <utility>
#include <string>

#define EIGEN_NO_DEBUG
#define EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_MPL2_ONLY
#include <Eigen/Core>
#include <Eigen/Geometry>

std::vector<Eigen::Vector3d> arrangeParticlesInFCCL(
        double lattice_const, const Eigen::Vector3i& cube_size);

double calcWholeKineticEnergy(
        const std::vector<Eigen::Vector3d>& ptcls_velocity, double ptcl_mass);

std::vector<Eigen::Vector3d> initVelocity(
        int num_particles, double ptcl_mass, double target_temp);

std::vector<Eigen::Vector3d> controlTempByScalingVel(
        std::vector<Eigen::Vector3d>& ptcls_velocity,
        double ptcl_mass, double target_temp);

std::pair<double, std::vector<Eigen::Vector3d>> calcLJPotentialAndForce(
        const std::vector<Eigen::Vector3d>& ptcls_pos,
        const Eigen::Vector3d& volume, const std::string& bc_mode);

std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> 
manageBoundaryCollision(
        std::vector<Eigen::Vector3d>& ptcls_pos,
        const Eigen::Vector3d& volume, const std::string& bc_mode);

std::vector<int> calcHistogram(
        const std::vector<double>& src_arr, int hist_size, 
        const std::pair<double, double>& hist_range);

std::pair<std::vector<double>, std::vector<double>> calcRadialDistributionFunction(
        const std::vector<Eigen::Vector3d>& ptcls_pos,
        const Eigen::Vector3d& volume,
        double number_density, int RDF_hist_size, const std::string bc_mode);

std::pair<std::vector<double>, std::vector<double>>
calcMeanSquareDisplacement(
        const std::vector<std::vector<Eigen::Vector3d>>& ptcls_fpos_allst, double dt);

#endif
