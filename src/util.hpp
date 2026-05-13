#pragma once
#include <iostream>
#include <ctime>
#include <vector>
#include <cstdint>
#include <numeric>
#include <cmath>
#include <fstream>
#include <string>
#include "tests.hpp"

namespace util {
    /**
     * @brief Функция для проверки на равномерность
     */
    void chiSquareUniform(const std::vector<uint32_t>& sample) {
        int n = sample.size();
        const int num_bins = 10;
        
        std::vector<int> observed(num_bins, 0);
        
        double expected = static_cast<double>(n) / num_bins;

        for (uint32_t val : sample) {
            double norm_val = val / 4294967296.0; 
            int bin = static_cast<int>(norm_val * num_bins);
            if (bin >= num_bins) bin = num_bins - 1; 
            observed[bin]++;
        }

        double chi2_stat = 0.0;
        for (int i = 0; i < num_bins; ++i) {
            double diff = observed[i] - expected;
            chi2_stat += (diff * diff) / expected;
        }

        double chi2_crit_min = 3.325;
        double chi2_crit_max = 16.919;

        std::cout << " | Хи-квадрат U: " << chi2_stat;
        if (chi2_stat > chi2_crit_max) {
            std::cout << " (НЕ равномерно: разброс слишком велик)";
        } else if (chi2_stat < chi2_crit_min) {
            std::cout << " (НЕ равномерно: разброс слишком мал)";
        } else {
            std::cout << " (Равномерно)";
        }
    }

    /**
     * @brief Функция для проверки на случайность
     */
    void chiSquareRandomness(const std::vector<uint32_t>& sample) {
        int n = sample.size();
        const int k = 5;
        const int num_cells = k * k; 
        
        std::vector<int> observed(num_cells, 0);

        int pairs_count = n / 2;
        double expected = static_cast<double>(pairs_count) / num_cells;

        for (int i = 0; i < n - 1; i += 2) {
            double norm1 = sample[i] / 4294967296.0;
            double norm2 = sample[i+1] / 4294967296.0;

            int bin1 = static_cast<int>(norm1 * k);
            int bin2 = static_cast<int>(norm2 * k);

            if (bin1 >= k) bin1 = k - 1;
            if (bin2 >= k) bin2 = k - 1;

            int cell_index = bin1 * k + bin2;
            observed[cell_index]++;
        }

        double chi2_stat = 0.0;
        for (int i = 0; i < num_cells; ++i) {
            double diff = observed[i] - expected;
            chi2_stat += (diff * diff) / expected;
        }

        double chi2_crit_min = 13.848;
        double chi2_crit_max = 36.415;

        std::cout << " | Хи-квадрат R: " << chi2_stat;
        if (chi2_stat > chi2_crit_max) {
            std::cout << " (НЕ случайно: разброс слишком велик)\n";
        } else if (chi2_stat < chi2_crit_min) {
            std::cout << " (НЕ случайно: разброс слишком мал)\n";
        } else {
            std::cout << " (Случайно)\n";
        }
    }

    /**
     * Функция, выводящая на экран информацию о выборке: среднее, дисперсию, коэффициент ковариации; результаты проверки
     * на равномерность и случайность; результаты проверки тестов NIST
     */
    void printStats(const std::vector<uint32_t>& sample, int gen_id, int sample_id) {
        uint64_t sum = std::accumulate(sample.begin(), sample.end(), 0ULL);
        double mean = static_cast<double>(sum) / sample.size();

        double variance_sum = 0.0;
        for (uint32_t val : sample) {
            double diff = val - mean;
            variance_sum += diff * diff;
        }
        double stddev = std::sqrt(variance_sum / sample.size());

        double cv = stddev / mean;

        std::cout << "Генератор: " << gen_id << ", Выборка: " << sample_id 
                << " | Среднее: " << mean 
                << " | Отклонение: " << stddev 
                << " | Коэф. Ков.: " << cv;

        util::chiSquareUniform(sample); // TODO(zhabka_zhaba): проверить
        util::chiSquareRandomness(sample); // TODO(zhabka_zhaba): проверить
        tests::runAllNistTests(sample); 
    }

    /**
     * Функция для записи времен замеров в файл
     */
    void write_csv(const std::string& filename, 
                   const std::vector<int>& sizes, 
                   const std::vector<double>& times_gen1,
                   const std::vector<double>& times_gen2,
                   const std::vector<double>& times_gen3,
                   const std::vector<double>& times_std) {
        
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Ошибка: не удалось открыть файл " << filename << " для записи.\n";
            return;
        }

        file << "Size,CubeCG,LCGShuffle,QuickBBS,MT19937\n";

        for (size_t i = 0; i < sizes.size(); ++i) {
            file << sizes[i] << ","
                 << times_gen1[i] << ","
                 << times_gen2[i] << ","
                 << times_gen3[i] << ","
                 << times_std[i] << "\n";
        }

        file.close();
        std::cout << "Данные успешно сохранены в " << filename << "\n";
    }
}