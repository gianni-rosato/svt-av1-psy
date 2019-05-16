/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file PerformanceCollect.cc
 *
 * @brief Impelmentation of performance tool for timing collection
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#include <thread>
#include "PerformanceCollect.h"
#include "gtest/gtest.h"
#include "../util.h"

uint64_t PerformanceCollect::get_time_tick() {
    using namespace std::chrono;
    steady_clock::time_point tp = steady_clock::now();
    steady_clock::duration dtn = tp.time_since_epoch();
    return dtn.count() / 1000000;
}
