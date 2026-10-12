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
    ProbeCover cover;  // empty when anchored: the row loop needs no cover
    uint32_t covered_frequency;
    uint32_t total_frequency;
    double scan_ns;
    scan::Walk walk;
    bool matches_all;  // the empty pattern: every row, no scan
};

inline Analysis analyze(const uint8_t* pattern, size_t n, const Dictionary& dict, const Frequency& freq,
                        size_t row_count, Anchor anchor = Anchor::None) {
    if (n == 0) return Analysis{ProbeCover{}, 0, freq.total, 0.0, {}, true};
    assert(n <= MAX_PATTERN_LEN);
    if (anchor != Anchor::None) {
        // one layout, so nothing to cut and no cover to price: the walk alone decides, per row
        AlignmentGraph graph = build_alignment_graph(dict, pattern, n, freq, anchor);
        return Analysis{ProbeCover{}, 0, freq.total, 0.0, scan::Walk(graph, pattern), false};
    }
    Planned p = plan(dict, pattern, n, freq, row_count);
    return Analysis{std::move(p.cover), p.covered, freq.total, p.scan_ns, std::move(p.walk), false};
}

inline void all_rows(size_t offsets_len, std::vector<size_t>& out) {
    for (size_t row = 0; row + 1 < offsets_len; ++row) out.push_back(row);
}

// Ascending rows containing the pattern: every hit goes through the walk.
// `dict` is the table the analysis was built over. A non-empty pattern with
// an empty cover appends nothing; the empty pattern appends every row.
inline void candidate_rows(const uint8_t* codes, size_t len, const uint32_t* row_offsets, size_t offsets_len,
                           const Dictionary& dict, const Analysis& analysis, std::vector<size_t>& out) {
    if (analysis.matches_all) return all_rows(offsets_len, out);
    scan::execute(analysis.cover, codes, len, row_offsets, offsets_len, analysis.covered_frequency,
                  analysis.total_frequency, scan::WalkCheck{analysis.walk, dict, codes}, out);
}

namespace detail {

template <Anchor A>
size_t anchored_rows(const uint8_t* codes, const int32_t* row_starts, const uint32_t* row_lengths,
                     const uint32_t* sel, size_t count, const Dictionary& dict, const scan::Walk& walk, uint32_t* out) {
    size_t out_count = 0;
    for (size_t index = 0; index < count; ++index) {
        const uint32_t row = sel ? sel[index] : static_cast<uint32_t>(index);
        const size_t start = static_cast<size_t>(row_starts[row]);
        const size_t end = start + static_cast<size_t>(row_lengths[row]);
        const bool match = walk.check_anchored<A>(dict, codes, start, end);
        out[out_count] = row;
        out_count += match;
    }
    return out_count;
}

}  // namespace detail

// The rows that start (Anchor::Start) or end (Anchor::End) with the walk's pattern, written to
// out[0 .. returned count) in the order checked: exact, and without a byte scan, since the anchor names the only
// position the pattern can sit at. The rows checked are sel[0 .. count), or 0 .. count when sel is null; out must
// hold count entries and may be sel itself. Each row is read on its own,
// codes[row_starts[row] .. row_starts[row] + row_lengths[row]), so the rows need not be adjacent or in order;
// an empty row matches nothing.
inline size_t anchored_rows(const uint8_t* codes, const int32_t* row_starts, const uint32_t* row_lengths,
                            const uint32_t* sel, size_t count, const Dictionary& dict, const scan::Walk& walk,
                            uint32_t* out) {
    if (walk.anchor() == Anchor::End) {
        return detail::anchored_rows<Anchor::End>(codes, row_starts, row_lengths, sel, count, dict, walk, out);
    }
    return detail::anchored_rows<Anchor::Start>(codes, row_starts, row_lengths, sel, count, dict, walk, out);
}

}  // namespace fsst::search::prefilter
