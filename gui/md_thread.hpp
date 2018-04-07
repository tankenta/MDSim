#pragma once
#include <thread>
#include <atomic>

#include "md_utils.hpp"

class MDSimThread : public MDSim {
private:
    std::thread th;

    void runThread();

public:
    std::atomic<int> curr_step;

    MDSimThread();
    ~MDSimThread();
};
