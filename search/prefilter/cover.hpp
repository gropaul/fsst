#pragma once

#include <algorithm>
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

    // Merge runs that overlap or abut, then file single-code runs as points
    // and the rest as ranges. Input in any order.
    static ProbeCover from_runs(std::vector<CodeRange> runs) {
        std::sort(runs.begin(), runs.end(), [](CodeRange a, CodeRange b) { return a.begin < b.begin; });
        std::vector<CodeRange> merged;
        for (CodeRange run : runs) {
            if (!merged.empty() && run.begin <= static_cast<size_t>(merged.back().last) + 1)
                merged.back().last = std::max(merged.back().last, run.last);
            else
                merged.push_back(run);
        }
        ProbeCover cover;
        for (CodeRange run : merged) {
            if (run.begin == run.last)
                cover.points.push_back(run.begin);
            else
                cover.ranges.push_back(run);
        }
        return cover;
    }
};

}  // namespace fsst::search::prefilter
