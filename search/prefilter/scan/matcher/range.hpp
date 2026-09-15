#pragma once

// `code - begin <= last - begin` unsigned, ORed over R ranges; folded into
// every kernel and run alone here. K = 0, R >= 1.

#include <vector>

#include "shared.hpp"

namespace fsst::search::prefilter::scan::matcher {

struct Held {
    uint8x16_t lo;
    uint8x16_t width;
};

inline Held hold(CodeRange r) {
    return {vdupq_n_u8(r.begin), vdupq_n_u8(static_cast<uint8_t>(r.last - r.begin))};
}

inline Hits inside(Held h, Vectors codes) {
    Hits out;
    for (size_t i = 0; i < 4; ++i) out[i] = vcleq_u8(vsubq_u8(codes[i], h.lo), h.width);
    return out;
}

inline Hits check_ranges(Hits hit, const std::vector<Held>& held, Vectors codes) {
    for (const Held& h : held) hit = or_hits(hit, inside(h, codes));
    return hit;
}

template <bool SKIP_MOVEMASK_IF_NO_MATCH>
class Range {
   public:
    explicit Range(const ProbeCover& cover) : first_(hold(cover.ranges[0])) {
        for (size_t i = 1; i < cover.ranges.size(); ++i) rest_.push_back(hold(cover.ranges[i]));
    }

    bool check(const uint8_t* codes, Mask& bits) const {
        return words<SKIP_MOVEMASK_IF_NO_MATCH>(codes, bits, [&](Vectors v) {
            return check_ranges(inside(first_, v), rest_, v);
        });
    }

   private:
    Held first_;
    std::vector<Held> rest_;
};

}  // namespace fsst::search::prefilter::scan::matcher
