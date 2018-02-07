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
        const double lattice_const, Eigen::Vector3i cube_size);

double calcWholeKineticEnergy(
        std::vector<Eigen::Vector3d> ptcls_velocity, const double ptcl_mass);

std::vector<Eigen::Vector3d> initVelocity(
        const int num_particles, const double ptcl_mass, const double target_temp);

std::vector<Eigen::Vector3d> controlTempByScalingVel(
        std::vector<Eigen::Vector3d> ptcls_velocity,
        const double ptcl_mass, const double target_temp);

std::pair<double, std::vector<Eigen::Vector3d>> calcLJPotentialAndForce(
        std::vector<Eigen::Vector3d> ptcls_pos,
        Eigen::Vector3d volume, std::string bc_mode);

std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> 
manageBoundaryCollision(
        std::vector<Eigen::Vector3d> ptcls_pos,
        Eigen::Vector3d volume, std::string bc_mode);

#endif
