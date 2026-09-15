#pragma once

#include <cstdint>
#include <vector>

namespace fsst::search::prefilter {

// Inclusive run of codes.
struct CodeRange {
    uint8_t begin;
    uint8_t last;

    bool contains(uint8_t code) const { return begin <= code && code <= last; }
};

// The codes the scan probes for: single codes and inclusive runs, disjoint.
struct ProbeCover {
    std::vector<uint8_t> points;
    std::vector<CodeRange> ranges;

    bool empty() const { return points.empty() && ranges.empty(); }

    bool contains(uint8_t code) const {
        for (uint8_t p : points)
            if (p == code) return true;
        for (const CodeRange& r : ranges)
            if (r.contains(code)) return true;
        return false;
    }
};

}  // namespace fsst::search::prefilter
