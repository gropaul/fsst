#pragma once

// Load, OR, any and movemask around each kernel's compare. NEON, u8 lanes:
// 64 codes are four vectors, hits are 0xFF/0x00 bytes.

#ifndef __ARM_NEON
#error "the vector matchers are NEON only for now"
#endif

#include <arm_neon.h>

#include <array>

#include "../scan.hpp"

namespace fsst::search::prefilter::scan::matcher {

using Vectors = std::array<uint8x16_t, 4>;
using Hits = std::array<uint8x16_t, 4>;

inline Vectors load(const uint8_t* at) {
    return {vld1q_u8(at), vld1q_u8(at + 16), vld1q_u8(at + 32), vld1q_u8(at + 48)};
}

inline Hits or_hits(Hits a, Hits b) {
    return {vorrq_u8(a[0], b[0]), vorrq_u8(a[1], b[1]), vorrq_u8(a[2], b[2]),
            vorrq_u8(a[3], b[3])};
}

inline bool any(Hits a, Hits b) {
    uint8x16_t x = vorrq_u8(vorrq_u8(a[0], a[1]), vorrq_u8(a[2], a[3]));
    uint8x16_t y = vorrq_u8(vorrq_u8(b[0], b[1]), vorrq_u8(b[2], b[3]));
    return vmaxvq_u8(vorrq_u8(x, y)) != 0;
}

// One bit per lane, then three pairwise adds; both words leave in one store.
inline void movemask(uint64_t* bits, Hits a, Hits b) {
    alignas(16) static const uint8_t LANE_BIT[16] = {1, 2, 4, 8, 16, 32, 64, 128,
                                                     1, 2, 4, 8, 16, 32, 64, 128};
    const uint8x16_t lane_bit = vld1q_u8(LANE_BIT);
    auto half = [lane_bit](Hits h) {
        uint8x16_t lo = vpaddq_u8(vandq_u8(h[0], lane_bit), vandq_u8(h[1], lane_bit));
        uint8x16_t hi = vpaddq_u8(vandq_u8(h[2], lane_bit), vandq_u8(h[3], lane_bit));
        return vpaddq_u8(lo, hi);
    };
    vst1q_u64(bits, vreinterpretq_u64_u8(vpaddq_u8(half(a), half(b))));
}

// Per pair of mask words: load 128 codes, hits on each 64, movemask both.
// Returns whether any pair was written; a skipped pair is zeroed.
template <bool SKIP_MOVEMASK_IF_NO_MATCH, typename HitsOf>
bool words(const uint8_t* codes, Mask& bits, HitsOf hits) {
    bool written = false;
    for (size_t pair = 0; pair < BLOCK / 128; ++pair) {
        const uint8_t* at = codes + pair * 128;
        Hits a = hits(load(at));
        Hits b = hits(load(at + 64));
        if (SKIP_MOVEMASK_IF_NO_MATCH && !any(a, b)) {
            bits[2 * pair] = 0;
            bits[2 * pair + 1] = 0;
            continue;
        }
        written = true;
        movemask(&bits[2 * pair], a, b);
    }
    return written;
}

}  // namespace fsst::search::prefilter::scan::matcher
