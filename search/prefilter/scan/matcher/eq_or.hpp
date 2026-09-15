#pragma once

// One compare per token per vector, ORed, plus the ranges. K >= 1, any R.

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
    }

    bool check(const uint8_t* codes, Mask& bits) const {
        return words<SKIP_MOVEMASK_IF_NO_MATCH>(codes, bits, [&](Vectors v) {
            return check_ranges(hits(v), ranges_, v);
        });
    }

   private:
    Hits hits(Vectors codes) const {
        Hits hit;
        for (size_t i = 0; i < 4; ++i) hit[i] = vceqq_u8(codes[i], tokens_[0]);
        for (size_t k = 1; k < tokens_.size(); ++k)
            for (size_t i = 0; i < 4; ++i)
                hit[i] = vorrq_u8(hit[i], vceqq_u8(codes[i], tokens_[k]));
        return hit;
    }

    std::vector<uint8x16_t> tokens_;
    std::vector<Held> ranges_;
};

}  // namespace fsst::search::prefilter::scan::matcher
