#pragma once

// How often each code stands in the stream, one counting pass, with the
// cumulative sum so a code range is priced in one subtraction. Literal bytes
// behind escapes count as the code they look like, which is what the
// matcher sees. Advisory weights only; soundness never reads them.

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "cover.hpp"

namespace fsst::search::prefilter {

struct Frequency {
    std::array<uint32_t, 256> count{};
    std::array<uint32_t, 257> cum{};  // cum[c] = sum of count[0..c)
    uint32_t total = 0;

    static Frequency of(const uint8_t* codes, size_t len) {
        Frequency f;
        for (size_t i = 0; i < len; ++i) ++f.count[codes[i]];
        for (size_t c = 0; c < 256; ++c) f.cum[c + 1] = f.cum[c] + f.count[c];
        f.total = static_cast<uint32_t>(len);
        return f;
    }

    uint32_t of_range(CodeRange r) const { return cum[static_cast<size_t>(r.last) + 1] - cum[r.begin]; }

    // `codes` distinct, so the sum stays within total.
    uint32_t of_codes(const std::vector<uint8_t>& codes) const {
        uint32_t sum = 0;
        for (uint8_t c : codes) sum += count[c];
        return sum;
    }

    uint32_t of_cover(const ProbeCover& cover) const {
        uint32_t sum = of_codes(cover.points);
        for (const CodeRange& r : cover.ranges) sum += of_range(r);
        return sum;
    }
};

}  // namespace fsst::search::prefilter
