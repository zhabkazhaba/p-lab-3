#pragma once
#include <cmath>
#include <cstdint>

/**
 * Вариация LCG с кубической формулой вместо линейной
 */
class CubeCG {
    private:
    uint32_t value;
    const uint32_t A = 16;
    const uint32_t B = 36;
    const uint32_t C = 11298535;
    const uint32_t D = 11111;

    public:
    CubeCG(uint32_t seed) : value(seed) {};

    uint32_t generate() {
        uint32_t value2 = value * value;
        uint32_t value3 = value2 * value;

        value = A * value3 + B * value2 + C * value + D;

        return value;
    }
};