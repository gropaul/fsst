#pragma once

// The prefilter's entry points: analyze a pattern into a sound probe cover
// and its walk, then run the cover over a code stream, exact through the
// walk or as the superset of rows holding a covered code.

#include <cstdint>
#include <vector>

#include "plan.hpp"
#include "scan/execute.hpp"

namespace fsst::search::prefilter {

// One node per needle offset plus the sink, inside the walk's u16 ids.
constexpr size_t MAX_PATTERN_LEN = 65535;

struct Analysis {
    ProbeCover cover;
    uint32_t covered_frequency;
    uint32_t total_frequency;
    double scan_ns;
    scan::Walk walk;
    bool matches_all;  // the empty pattern: every row, no scan
};

inline Analysis analyze(const uint8_t* pattern, size_t n, const Dictionary& dict, const Frequency& freq,
                        size_t row_count) {
    if (n == 0) return Analysis{ProbeCover{}, 0, freq.total, 0.0, {}, true};
    assert(n <= MAX_PATTERN_LEN);
    Planned p = plan(dict, pattern, n, freq, row_count);
    return Analysis{std::move(p.cover), p.covered, freq.total, p.scan_ns, std::move(p.walk), false};
}

inline void all_rows(size_t offsets_len, std::vector<size_t>& out) {
    for (size_t row = 0; row + 1 < offsets_len; ++row) out.push_back(row);
}

// Ascending rows containing the pattern: every hit goes through the walk.
// `dict` is the table the analysis was built over. A non-empty pattern with
// an empty cover appends nothing; the empty pattern appends every row.
template <typename O>
void candidate_rows(const uint8_t* codes, size_t len, const O* row_offsets, size_t offsets_len, const Dictionary& dict,
                    const Analysis& analysis, std::vector<size_t>& out) {
    if (analysis.matches_all) return all_rows(offsets_len, out);
    scan::execute(analysis.cover, codes, len, row_offsets, offsets_len, analysis.covered_frequency,
                  analysis.total_frequency, scan::WalkCheck{analysis.walk, dict, codes}, out);
}

// Ascending rows holding a covered code, without the walk: a sound superset
// for a caller that decodes the survivors.
template <typename O>
void superset_rows(const uint8_t* codes, size_t len, const O* row_offsets, size_t offsets_len,
                   const Analysis& analysis, std::vector<size_t>& out) {
    if (analysis.matches_all) return all_rows(offsets_len, out);
    scan::execute(analysis.cover, codes, len, row_offsets, offsets_len, analysis.covered_frequency,
                  analysis.total_frequency, scan::Superset{}, out);
}

}  // namespace fsst::search::prefilter
