#pragma once

// The universal 16x16 bitmap of the cover, addressed by nibble. Any K, any
// R, flat cost: ranges are set into the tables and the kernel never sees them.
// `low` holds the rows for high nibbles 0..7, `high` those for 8..15.

#include "shared.hpp"

namespace fsst::search::prefilter::scan::matcher {

template <bool SKIP_MOVEMASK_IF_NO_MATCH>
class Nibble {
   public:
    explicit Nibble(const ProbeCover& cover) {
        alignas(16) uint8_t low[16] = {};
        alignas(16) uint8_t high[16] = {};
        auto set = [&](uint8_t code) {
            uint8_t row = code & 0x0f;
            uint8_t bit = code >> 4;
            if (bit < 8)
                low[row] |= static_cast<uint8_t>(1u << bit);
            else
                high[row] |= static_cast<uint8_t>(1u << (bit - 8));
        };
        for (uint8_t t : cover.points) set(t);
        for (CodeRange r : cover.ranges)
            for (unsigned c = r.begin; c <= r.last; ++c) set(static_cast<uint8_t>(c));
#if defined(__ARM_NEON)
        low_ = vld1q_u8(low);
        high_ = vld1q_u8(high);
#elif defined(__AVX512BW__)
        low_ = broadcast_table(low);
        high_ = broadcast_table(high);
#endif
    }

    // 1 << (high nibble % 8), indexed by the whole high nibble.
    static constexpr uint8_t BIT[16] = {1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128};

#if defined(__ARM_NEON)
    bool check(const uint8_t* codes, Mask& bits) const {
        const uint8x16_t bit = vld1q_u8(BIT);
        const uint8x16_t row = vdupq_n_u8(0x8f);
        const uint8x16_t half = vdupq_n_u8(0x80);
        return words<SKIP_MOVEMASK_IF_NO_MATCH>(codes, bits, [&](Vectors v, const uint8_t*) {
            Hits out;
            for (size_t i = 0; i < 4; ++i) {
                // bit 7 kept in the index: vqtbl1q returns 0 for the wrong half.
                uint8x16_t index = vandq_u8(v[i], row);
                uint8x16_t rows = vorrq_u8(vqtbl1q_u8(low_, index),
                                           vqtbl1q_u8(high_, veorq_u8(index, half)));
                out[i] = vtstq_u8(rows, vqtbl1q_u8(bit, vshrq_n_u8(v[i], 4)));
            }
            return out;
        });
    }

   private:
    uint8x16_t low_;
    uint8x16_t high_;
#elif defined(__AVX512BW__)
    bool check(const uint8_t* codes, Mask& bits) const {
        const __m512i bit = broadcast_table(BIT);
        const __m512i row = _mm512_set1_epi8(static_cast<char>(0x8f));
        const __m512i half = _mm512_set1_epi8(static_cast<char>(0x80));
        const __m512i nibble = _mm512_set1_epi8(0x0f);
        return words<SKIP_MOVEMASK_IF_NO_MATCH>(codes, bits, [&](Vectors v, const uint8_t*) {
            // bit 7 kept in the index: vpshufb returns 0 for the wrong half.
            __m512i index = _mm512_and_si512(v, row);
            __m512i rows = _mm512_or_si512(_mm512_shuffle_epi8(low_, index),
                                           _mm512_shuffle_epi8(high_, _mm512_xor_si512(index, half)));
            __m512i n1 = _mm512_and_si512(_mm512_srli_epi16(v, 4), nibble);
            return _mm512_test_epi8_mask(rows, _mm512_shuffle_epi8(bit, n1));
        });
    }

   private:
    __m512i low_;
    __m512i high_;
#endif
};

}  // namespace fsst::search::prefilter::scan::matcher
