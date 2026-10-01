#pragma once

// A LIKE pattern as the literal runs between its wildcards. One run is a substring, or a prefix when it is
// anchored, and both are exact. Several runs are answered as a superset: every run is required to occur in
// the row, but not to occur in that order, so the caller confirms the survivors with the exact predicate.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "prefilter.hpp"

namespace fsst::search::prefilter {

struct PatternPlan {
    std::vector<Analysis> runs;
    size_t driver = 0;    // the run whose cover drives the scan: the anchored one, else the rarest
    bool anchored = false;  // runs[0] must start at the row's first byte
    bool exact = false;     // one run: the rows are the answer, not a superset
    bool usable = false;
};

// One analysis per run. Fails when a run is empty or the symbol table cannot spell it.
inline PatternPlan plan_pattern(const Dictionary& dict, const Frequency& freq, const std::vector<std::string>& runs,
                                bool anchor_start, bool anchor_end, size_t row_count) {
    PatternPlan plan;
    if (runs.empty()) return plan;
    for (size_t i = 0; i < runs.size(); ++i) {
        const Anchor anchor = (i == 0 && anchor_start) ? Anchor::Start : Anchor::None;
        Analysis a = analyze(reinterpret_cast<const uint8_t*>(runs[i].data()), runs[i].size(), dict, freq, row_count,
                             anchor);
        if (a.matches_all) return plan;
        if (anchor == Anchor::None && a.cover.empty()) return plan;
        plan.runs.push_back(std::move(a));
    }
    plan.anchored = anchor_start;
    if (!anchor_start)
        for (size_t i = 1; i < plan.runs.size(); ++i)
            if (plan.runs[i].covered_frequency < plan.runs[plan.driver].covered_frequency) plan.driver = i;
    plan.exact = plan.runs.size() == 1 && !anchor_end;
    plan.usable = true;
    return plan;
}

// Whether the run occurs anywhere in codes[start..end): its cover names the positions worth testing and its
// walk decides each one.
inline bool row_has(const Dictionary& dict, const uint8_t* codes, size_t start, size_t end, const Analysis& run) {
    for (size_t i = start; i < end; ++i)
        if (run.cover.matches(codes, i, end) && run.walk.check(dict, codes, start, end, i)) return true;
    return false;
}

// Stage two for a multi-run pattern: the driver's hit must verify, and every other run must occur in the row.
struct PatternCheck {
    const PatternPlan& plan;
    const Dictionary& dict;
    const uint8_t* codes;

    bool passes(size_t hit, size_t row_start, size_t row_end) const {
        if (!plan.runs[plan.driver].walk.check(dict, codes, row_start, row_end, hit)) return false;
        for (size_t i = 0; i < plan.runs.size(); ++i)
            if (i != plan.driver && !row_has(dict, codes, row_start, row_end, plan.runs[i])) return false;
        return true;
    }
};

// Ascending rows the pattern can match.
template <typename O>
void pattern_rows(const uint8_t* codes, size_t len, const O* row_offsets, size_t offsets_len, const Dictionary& dict,
                  const PatternPlan& plan, std::vector<size_t>& out) {
    if (!plan.usable) return;
    if (plan.anchored) {
        for (size_t row = 0; row + 1 < offsets_len; ++row) {
            const size_t start = row_offsets[row], end = row_offsets[row + 1];
            if (start == end) continue;
            if (!plan.runs[0].walk.check(dict, codes, start, end, start)) continue;
            bool all = true;
            for (size_t i = 1; i < plan.runs.size() && all; ++i) all = row_has(dict, codes, start, end, plan.runs[i]);
            if (all) out.push_back(row);
        }
        return;
    }
    const Analysis& driver = plan.runs[plan.driver];
    scan::execute(driver.cover, codes, len, row_offsets, offsets_len, driver.covered_frequency,
                  driver.total_frequency, PatternCheck{plan, dict, codes}, out);
}

}  // namespace fsst::search::prefilter
