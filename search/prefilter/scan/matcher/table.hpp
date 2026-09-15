#pragma once

// One table byte per code, one load per code. Scalar, any target, any cover.
// The oracle the vector kernels are tested against.

#include <array>

#include "../scan.hpp"

namespace fsst::search::prefilter::scan::matcher {

class Table {
   public:
    explicit Table(const ProbeCover& cover) : admits_{} {
        for (uint8_t t : cover.points) admits_[t] = 0xFF;
        for (CodeRange r : cover.ranges)
            for (unsigned c = r.begin; c <= r.last; ++c) admits_[c] = 0xFF;
    }

    bool check(const uint8_t* codes, Mask& bits) const {
        uint64_t any = 0;
        for (size_t word = 0; word < bits.size(); ++word) {
            const uint8_t* at = codes + word * 64;
            uint64_t packed = 0;
            for (size_t group = 0; group < 8; ++group) {
                uint64_t byte = 0;
                for (size_t bit = 0; bit < 8; ++bit)
                    byte |= uint64_t{admits_[at[group * 8 + bit]] & 1u} << bit;
                packed |= byte << (8 * group);
            }
            bits[word] = packed;
            any |= packed;
        }
        return any != 0;
    }

   private:
    std::array<uint8_t, 256> admits_;
};

}  // namespace fsst::search::prefilter::scan::matcher
