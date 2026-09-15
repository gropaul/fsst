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
        for (uint8_t t : cover.points) tokens_.push_back(vdupq_n_u8(t));
        for (CodeRange r : cover.ranges) ranges_.push_back(hold(r));
        for (CodePair p : cover.pairs) pairs_.push_back({vdupq_n_u8(p.first), vdupq_n_u8(p.second)});
    }

    // The pair path is chosen once per block, so a cover without pairs runs
    // the token loop alone.
    bool check(const uint8_t* codes, Mask& bits) const {
        if (pairs_.empty())
            return words<SKIP_MOVEMASK_IF_NO_MATCH>(codes, bits, [&](Vectors v, const uint8_t*) {
                return check_ranges(token_hits(v), ranges_, v);
            });
        return words<SKIP_MOVEMASK_IF_NO_MATCH>(codes, bits, [&](Vectors v, const uint8_t* at) {
            Hits hit = tokens_.empty() ? Hits{vdupq_n_u8(0), vdupq_n_u8(0), vdupq_n_u8(0), vdupq_n_u8(0)}
                                       : token_hits(v);
            Vectors next = load(at + 1);
            for (const Pair& p : pairs_)
                for (size_t i = 0; i < 4; ++i)
                    hit[i] = vorrq_u8(hit[i], vandq_u8(vceqq_u8(v[i], p.first), vceqq_u8(next[i], p.second)));
            return check_ranges(hit, ranges_, v);
        });
    }

   private:
    struct Pair {
        uint8x16_t first, second;
    };

    // K >= 1 here.
    Hits token_hits(Vectors codes) const {
        Hits hit;
        for (size_t i = 0; i < 4; ++i) hit[i] = vceqq_u8(codes[i], tokens_[0]);
        for (size_t k = 1; k < tokens_.size(); ++k)
            for (size_t i = 0; i < 4; ++i) hit[i] = vorrq_u8(hit[i], vceqq_u8(codes[i], tokens_[k]));
        return hit;
    }

    std::vector<uint8x16_t> tokens_;
    std::vector<Held> ranges_;
    std::vector<Pair> pairs_;
};

}  // namespace fsst::search::prefilter::scan::matcher
