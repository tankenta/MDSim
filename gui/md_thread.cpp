#include <thread>
#include <atomic>

#include "md_utils.hpp"

MDSimThread::MDSimThread(
        double dt, double total_time, double temp_cont_time,
        double number_density, std::string bc_mode_str, double target_temp,
        int num_particles, double ptcl_mass, int RDF_hist_size,
        const Eigen::Vector3i& cube_size)
: MDSim(
        dt, total_time, temp_cont_time,
        number_density, bc_mode_str, target_temp,
        num_particles, ptcl_mass, RDF_hist_size,
        cube_size) {
    curr_step.store(0, std::memory_order_release);
    th = std::thread(runThread);
}

MDSimThread::~MDSimThread() {
    th.join();
}

MDSimThread::runThread() {
    for (const auto& step : steps) {
        curr_step.store(step, std::memory_order_release);
        next();
    }
}
