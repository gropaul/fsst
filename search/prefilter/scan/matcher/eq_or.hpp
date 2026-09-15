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

    bool check(const uint8_t* codes, Mask& bits) const {
        return words<SKIP_MOVEMASK_IF_NO_MATCH>(codes, bits, [&](Vectors v, const uint8_t* at) {
            return check_ranges(hits(v, at), ranges_, v);
        });
    }

   private:
    struct Pair {
        uint8x16_t first, second;
    };

    Hits hits(Vectors codes, const uint8_t* at) const {
        Hits hit;
        for (size_t i = 0; i < 4; ++i) hit[i] = vdupq_n_u8(0);
        for (const uint8x16_t& t : tokens_)
            for (size_t i = 0; i < 4; ++i) hit[i] = vorrq_u8(hit[i], vceqq_u8(codes[i], t));
        if (!pairs_.empty()) {
            Vectors next = load(at + 1);
            for (const Pair& p : pairs_)
                for (size_t i = 0; i < 4; ++i)
                    hit[i] = vorrq_u8(hit[i], vandq_u8(vceqq_u8(codes[i], p.first), vceqq_u8(next[i], p.second)));
        }
        return hit;
    }

    std::vector<uint8x16_t> tokens_;
    std::vector<Held> ranges_;
    std::vector<Pair> pairs_;
};

}  // namespace fsst::search::prefilter::scan::matcher
