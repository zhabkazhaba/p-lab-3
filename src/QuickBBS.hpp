#pragma once
#include <cstdint>

/**
 * Быстрая (но не устойчивая криптографически) вариация BBS 
 */
class QuickBBS {
private:
    uint64_t value;
    
    const uint64_t M = 564995933; // 22739 * 24847

public:
    QuickBBS(uint32_t seed) {
        value = seed % M;
        if (value == 0 || value == 1) {
            value = 111; 
        }
    }

    uint32_t generate() {
        value = (value * value) % M;

        uint32_t output = static_cast<uint32_t>(value);

        output = output ^ (output >> 13);
        output = output ^ (output << 9);

        return output;
    }
};