#include "CubeCG.hpp"
#include "LCGShuffle.hpp"
#include "QuickBBS.hpp"
#include "util.hpp"
#include <iostream>
#include <ctime>
#include <vector>
#include <cstdint>
#include <numeric>
#include <chrono>

#ifdef _WIN32
#include <windows.h>
#endif
#include <random>

int main() {
    #ifdef _WIN32
    SetConsoleCP(65001);
    SetConsoleOutputCP(65001);
    #endif


    uint32_t seed = std::time(nullptr);
    std::vector<std::vector<uint32_t>> gen1_results(25, std::vector<uint32_t>(2000));
    std::vector<std::vector<uint32_t>> gen2_results(25, std::vector<uint32_t>(2000));
    std::vector<std::vector<uint32_t>> gen3_results(25, std::vector<uint32_t>(2000));

    CubeCG generator1 = CubeCG(seed);
    LCGShuffle generator2 = LCGShuffle(seed);
    QuickBBS generator3 = QuickBBS(seed);

    for (unsigned int i = 0; i < 25; ++i) {
        for (unsigned int j = 0; j < 2000; ++j) {
            gen1_results[i][j] = generator1.generate();
            gen2_results[i][j] = generator2.generate();
            gen3_results[i][j] = generator3.generate();
        }
    }

    for (unsigned int i = 0; i < 25; ++i) {
        util::printStats(gen1_results[i], 1, i + 1);
        util::printStats(gen2_results[i], 2, i + 1);
        util::printStats(gen3_results[i], 3, i + 1);
        std::cout << "==========================\n";
    }

    seed = 54321;
    std::vector<int> sizes = {1000, 5000, 10000, 50000, 100000, 250000, 500000, 1000000, 100000000};
    std::vector<double> time1, time2, time3, time_std;

    for (int n : sizes) {
        CubeCG gen1(seed);
        LCGShuffle gen2(seed);
        QuickBBS gen3(seed);
        std::mt19937 gen_std(seed);

        volatile uint32_t dump1 = 0;
        volatile uint32_t dump2 = 0;
        volatile uint32_t dump3 = 0;
        volatile uint32_t dump_std = 0;

        // 1. CubeCG
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n; ++i) dump1 = gen1.generate();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed1 = end - start;
        time1.push_back(elapsed1.count());

        // 2. LCGShuffle
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n; ++i) dump2 = gen2.generate();
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed2 = end - start;
        time2.push_back(elapsed2.count());

        // 3. QuickBBS
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n; ++i) dump3 = gen3.generate();
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed3 = end - start;
        time3.push_back(elapsed3.count());

        // 4. MT19937
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n; ++i) dump_std = gen_std();
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed_std = end - start;
        time_std.push_back(elapsed_std.count());

        std::cout << "Создано: " << n << " элементов\n";
    }

    util::write_csv("../data/out_times.csv", sizes, time1, time2, time3, time_std);
    return 0;
}