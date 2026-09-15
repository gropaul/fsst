#pragma once

// Two nibble-table shuffles per code, ANDed: a token bit that survives both
// matched the whole code. Eight tokens per batch, batches ORed, ranges beside.
// At u8 width one batch is all the planner asks for; the bitmap in nibble.hpp
// takes over from nine tokens.

#include <array>
#include <vector>

#include "limits.hpp"
#include "range.hpp"
#include "shared.hpp"

namespace fsst::search::prefilter::scan::matcher {


// Token k owns bit 1 << k in the row each of its nibbles indexes.
inline void batch_tables(const uint8_t* tokens, size_t count, uint8_t* low, uint8_t* high) {
    for (size_t k = 0; k < count; ++k) {
        low[tokens[k] & 0x0f] |= static_cast<uint8_t>(1u << k);
        high[tokens[k] >> 4] |= static_cast<uint8_t>(1u << k);
    }
}

#if defined(__ARM_NEON)

struct Batch {
    uint8x16_t low;
    uint8x16_t high;
};

inline Batch batch(const uint8_t* tokens, size_t count) {
    alignas(16) uint8_t low[16] = {};
    alignas(16) uint8_t high[16] = {};
    batch_tables(tokens, count, low, high);
    return {vld1q_u8(low), vld1q_u8(high)};
}

template <size_t BATCHES>
inline Hits probe(const std::array<Batch, BATCHES>& batches, Vectors codes) {
    const uint8x16_t nibble = vdupq_n_u8(0x0f);
    Hits out;
    for (size_t i = 0; i < 4; ++i) {
        uint8x16_t n0 = vandq_u8(codes[i], nibble);
        uint8x16_t n1 = vshrq_n_u8(codes[i], 4);
        uint8x16_t low = vqtbl1q_u8(batches[0].low, n0);
        uint8x16_t high = vqtbl1q_u8(batches[0].high, n1);
        if constexpr (BATCHES == 1) {
            out[i] = vtstq_u8(low, high);
        } else {
            uint8x16_t hit = vandq_u8(low, high);
            for (size_t b = 1; b < BATCHES; ++b)
                hit = vorrq_u8(hit, vandq_u8(vqtbl1q_u8(batches[b].low, n0),
                                             vqtbl1q_u8(batches[b].high, n1)));
            out[i] = vtstq_u8(hit, hit);
        }
    }
    return out;
}

#elif defined(__AVX512BW__)

struct Batch {
    __m512i low;
    __m512i high;
};

inline Batch batch(const uint8_t* tokens, size_t count) {
    alignas(16) uint8_t low[16] = {};
    alignas(16) uint8_t high[16] = {};
    batch_tables(tokens, count, low, high);
    return {broadcast_table(low), broadcast_table(high)};
}

// No per-byte shift on x86: the high nibble is a word shift masked back.
template <size_t BATCHES>
inline Hits probe(const std::array<Batch, BATCHES>& batches, Vectors codes) {
    const __m512i nibble = _mm512_set1_epi8(0x0f);
    __m512i n0 = _mm512_and_si512(codes, nibble);
    __m512i n1 = _mm512_and_si512(_mm512_srli_epi16(codes, 4), nibble);
    __m512i low = _mm512_shuffle_epi8(batches[0].low, n0);
    __m512i high = _mm512_shuffle_epi8(batches[0].high, n1);
    if constexpr (BATCHES == 1) return _mm512_test_epi8_mask(low, high);
    __m512i hit = _mm512_and_si512(low, high);
    for (size_t b = 1; b < BATCHES; ++b)
        hit = _mm512_or_si512(hit, _mm512_and_si512(_mm512_shuffle_epi8(batches[b].low, n0),
                                                    _mm512_shuffle_epi8(batches[b].high, n1)));
    return _mm512_test_epi8_mask(hit, hit);
}

#elif defined(__AVX2__)

struct Batch {
    __m256i low;
    __m256i high;
};

inline Batch batch(const uint8_t* tokens, size_t count) {
    alignas(16) uint8_t low[16] = {};
    alignas(16) uint8_t high[16] = {};
    batch_tables(tokens, count, low, high);
    return {broadcast_table(low), broadcast_table(high)};
}

template <size_t BATCHES>
inline Hits probe(const std::array<Batch, BATCHES>& batches, Vectors codes) {
    const __m256i nibble = _mm256_set1_epi8(0x0f);
    Hits out;
    for (size_t i = 0; i < 2; ++i) {
        __m256i n0 = _mm256_and_si256(codes[i], nibble);
        __m256i n1 = _mm256_and_si256(_mm256_srli_epi16(codes[i], 4), nibble);
        __m256i hit = _mm256_and_si256(_mm256_shuffle_epi8(batches[0].low, n0),
                                       _mm256_shuffle_epi8(batches[0].high, n1));
        for (size_t b = 1; b < BATCHES; ++b)
            hit = _mm256_or_si256(hit, _mm256_and_si256(_mm256_shuffle_epi8(batches[b].low, n0),
                                                        _mm256_shuffle_epi8(batches[b].high, n1)));
        out[i] = nonzero(hit);
    }
    return out;
}

#endif

template <size_t BATCHES, bool SKIP_MOVEMASK_IF_NO_MATCH>
class NibbleN8 {
   public:
    explicit NibbleN8(const ProbeCover& cover) {
        const size_t n = cover.points.size();
        for (size_t b = 0; b < BATCHES; ++b) {
            size_t from = std::min(b * PER_BATCH, n);
            size_t to = std::min(from + PER_BATCH, n);
            batches_[b] = batch(cover.points.data() + from, to - from);
        }
        for (CodeRange r : cover.ranges) ranges_.push_back(hold(r));
    }

    bool check(const uint8_t* codes, Mask& bits) const {
        return words<SKIP_MOVEMASK_IF_NO_MATCH>(codes, bits, [&](Vectors v, const uint8_t*) {
            return check_ranges(probe(batches_, v), ranges_, v);
        });
    }

   private:
    std::array<Batch, BATCHES> batches_;
    std::vector<Held> ranges_;
};

}  // namespace fsst::search::prefilter::scan::matcher
