#pragma once

// One compare per token per vector, ORed, plus the ranges, plus the pairs:
// a pair is the compare of the first code ANDed with the compare of the
// second on the codes one lane later. K + P >= 1, any R.

#include <vector>

#include "range.hpp"
#include "shared.hpp"

namespace fsst::search::prefilter::scan::matcher {

template <bool SKIP_MOVEMASK_IF_NO_MATCH>
class EqOr {
   public:
    explicit EqOr(const ProbeCover& cover) {
        for (uint8_t t : cover.points) tokens_.push_back(broadcast(t));
        for (CodeRange r : cover.ranges) ranges_.push_back(hold(r));
        for (CodePair p : cover.pairs) pairs_.push_back({broadcast(p.first), broadcast(p.second)});
    }

    // The pair path is chosen once per block, so a cover without pairs runs
    // the token loop alone.
    bool check(const uint8_t* codes, Mask& bits) const {
        if (pairs_.empty())
            return words<SKIP_MOVEMASK_IF_NO_MATCH>(codes, bits, [&](Vectors v, const uint8_t*) {
                return check_ranges(token_hits(v), ranges_, v);
            });
        return words<SKIP_MOVEMASK_IF_NO_MATCH>(codes, bits, [&](Vectors v, const uint8_t* at) {
            Hits hit = tokens_.empty() ? no_hits() : token_hits(v);
            Vectors next = load(at + 1);
            for (const Pair& p : pairs_) hit = or_hits(hit, pair_hits(v, next, p));
            return check_ranges(hit, ranges_, v);
        });
    }

   private:
#if defined(__ARM_NEON)
    using Token = uint8x16_t;

    static Token broadcast(uint8_t code) { return vdupq_n_u8(code); }

    static Hits no_hits() { return {vdupq_n_u8(0), vdupq_n_u8(0), vdupq_n_u8(0), vdupq_n_u8(0)}; }

    // K >= 1 here.
    Hits token_hits(Vectors codes) const {
        Hits hit;
        for (size_t i = 0; i < 4; ++i) hit[i] = vceqq_u8(codes[i], tokens_[0]);
        for (size_t k = 1; k < tokens_.size(); ++k)
            for (size_t i = 0; i < 4; ++i) hit[i] = vorrq_u8(hit[i], vceqq_u8(codes[i], tokens_[k]));
        return hit;
    }
#elif defined(__AVX512BW__)
    using Token = __m512i;

    static Token broadcast(uint8_t code) { return _mm512_set1_epi8(static_cast<char>(code)); }

    static Hits no_hits() { return 0; }

    Hits token_hits(Vectors codes) const {
        Hits hit = _mm512_cmpeq_epi8_mask(codes, tokens_[0]);
        for (size_t k = 1; k < tokens_.size(); ++k) hit |= _mm512_cmpeq_epi8_mask(codes, tokens_[k]);
        return hit;
    }
#endif

    struct Pair {
        Token first, second;
    };

#if defined(__ARM_NEON)
    static Hits pair_hits(Vectors v, Vectors next, const Pair& p) {
        Hits hit;
        for (size_t i = 0; i < 4; ++i) hit[i] = vandq_u8(vceqq_u8(v[i], p.first), vceqq_u8(next[i], p.second));
        return hit;
    }
#elif defined(__AVX512BW__)
    static Hits pair_hits(Vectors v, Vectors next, const Pair& p) {
        return _mm512_cmpeq_epi8_mask(v, p.first) & _mm512_cmpeq_epi8_mask(next, p.second);
    }
#endif

    std::vector<Token> tokens_;
    std::vector<Held> ranges_;
    std::vector<Pair> pairs_;
};

}  // namespace fsst::search::prefilter::scan::matcher
