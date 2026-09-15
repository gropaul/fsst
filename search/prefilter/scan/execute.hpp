#pragma once

// The scan as the planner enters it: the region's facts, the plan, the run.
// The rest of what onpair keeps in scan/mod.rs beside the seam in scan.hpp.

#include <optional>
#include <vector>

#include "dispatch.hpp"
#include "policy/policy.hpp"

namespace fsst::search::prefilter::scan {

// Covered codes expected in this region: exact for the indexed population,
// a projection of it for a subset. Planning only, never correctness.
inline size_t expected_hits(size_t covered_frequency, size_t total_frequency, size_t code_count) {
    if (total_frequency == code_count) return covered_frequency;
    if (total_frequency == 0) return 0;
    unsigned __int128 projected = static_cast<unsigned __int128>(covered_frequency) * code_count / total_frequency;
    return projected > static_cast<unsigned __int128>(SIZE_MAX) ? SIZE_MAX : static_cast<size_t>(projected);
}

// None where the region admits nothing whatever runs over it: an empty
// cover proves no row matches, and a region with no rows has nothing to emit.
inline std::optional<policy::Facts> facts(const ProbeCover& cover, size_t code_count, size_t offsets_len,
                                          size_t covered_frequency, size_t total_frequency) {
    size_t row_count = offsets_len == 0 ? 0 : offsets_len - 1;
    if (cover.empty() || row_count == 0) return std::nullopt;
    return policy::Facts{expected_hits(covered_frequency, total_frequency, code_count), code_count, row_count};
}

// The planned scan over one region under whatever stage two asks of a hit.
// `covered_frequency` and `total_frequency` are the frequency index's counts.
template <typename O, typename Check>
void execute(const ProbeCover& cover, const uint8_t* codes, size_t len, const O* row_offsets, size_t offsets_len,
             size_t covered_frequency, size_t total_frequency, const Check& check, std::vector<size_t>& out) {
    auto f = facts(cover, len, offsets_len, covered_frequency, total_frequency);
    if (!f) return;
    dispatch::run<O>(policy::select(cover, *f), cover, codes, len, row_offsets, offsets_len, check, out);
}

}  // namespace fsst::search::prefilter::scan
