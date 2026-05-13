#pragma once
#include <cmath>
#include <cstdint>

/**
 * Вариация LCG с битовым смещением и операцией XOR
 */
class LCGShuffle {
    private:
    uint32_t value;
    const uint32_t A = 874235315;
    const uint32_t C = 1321313;
    uint32_t shuffle1 = 9;
    uint32_t shuffle2 = 13;

    public:
    LCGShuffle(uint32_t seed) : value(seed) {};

    uint32_t generate() {
        value = A * value + C;

        uint32_t output = value ^ (value >> shuffle1);
        output = output ^ (output << shuffle2);

        return output;
    }
};