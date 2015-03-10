/**
 * @file		test_copyable_atomics.cpp
 * @ingroup
 * @author	Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include <vector>
#include <unordered_map>
#include "concurrent/copyable_atomic.hpp"

#include "utils/logging.h"

using namespace bliss::concurrent;

int main(int argc, char** argv) {

        std::vector<std::atomic<int> > atoms(12);   // can't add, but can preallocate.

        for (int i = 0; i < 12; ++i) {
                atoms[i].store(i);
        }
        INFOF("Atomics preallocated in vector\n");
        for (int i = 0; i < 12; ++i) {
                auto x = atoms[i].fetch_add(1);
                auto y = atoms[i].load();
                INFOF("%d: %d->%d\n", i, x, y);
        }


        std::unordered_map<int, copyable_atomic<int> > atoms2;

        for (int i = 0; i < 12; ++i) {
                atoms2.emplace(std::make_pair(int(i*2), std::move(copyable_atomic<int>(i))));
                atoms2.emplace(i*2+1, copyable_atomic<int>(-i));
                atoms2[i+24] = i+24;
        }

        INFOF("Copyable atomics in unordered_map\n");
        for (int i = 0; i < 12; ++i) {
                auto x= atoms2.at(i*2).fetch_add(1);
                auto y = atoms2.at(i*2).load();
                INFOF("%d: %d->%d\n", i*2, x, y);
                INFOF("%d: %d\n", i*2+1, atoms2.at(i*2+1).load());
                INFOF("%d: %d\n", i+24, atoms2.at(i+24).load());
        }






}
