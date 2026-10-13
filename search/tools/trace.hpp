#pragma once

// The planner's alignment graph and chosen cut, for logging and figures.
// Paper tooling, removable as a whole: this directory, the plan_log target in
// search/CMakeLists.txt, and the two hooks it reads, CutGraph::origin and
// Planned::cut.

#include <cassert>
#include <cstdint>
#include <vector>

#include "../prefilter/prefilter.hpp"

namespace fsst::search::prefilter {

// The alignment graph, which edges the cut selects outright, and for each
// edge the codes before it the cut pairs it with instead.
struct Trace {
    AlignmentGraph graph;
    std::vector<bool> cut;
    std::vector<std::vector<uint8_t>> paired;
    Analysis analysis;  // what analyze returns for the same pattern
};

// Plans a non-empty pattern as analyze does, keeping the graph and the cut.
inline Trace trace(const uint8_t* pattern, size_t n, const Dictionary& dict, const Frequency& freq, size_t row_count) {
    assert(n > 0 && n <= MAX_PATTERN_LEN);
    Trace t{build_alignment_graph(dict, pattern, n, freq), {}, {}, {}};
    const CutGraph cg = CutGraph::with_pairs(t.graph, freq);
    Planned p = cheapest_cover(cg, freq, scan::policy::Region{freq.total, row_count});
    t.cut.assign(t.graph.edges.size(), false);
    t.paired.assign(t.graph.edges.size(), {});
    for (uint32_t at : p.cut) {
        const Edge& e = cg.edges[at];
        if (e.probe == Probe::Pair)
            t.paired[cg.origin[at]].push_back(e.pair.first);
        else
            t.cut[cg.origin[at]] = true;
    }
    t.analysis = Analysis{std::move(p.cover), p.covered, freq.total, p.scan_ns, scan::Walk(t.graph, pattern), false};
    return t;
}

// Every hit passes: the rows holding a covered code, a sound superset.
struct AnyHit {
    bool passes(size_t, size_t, size_t) const { return true; }
};

// Ascending rows holding a covered code of a non-empty pattern's analysis.
template <typename O>
void superset_rows(const uint8_t* codes, size_t len, const O* row_offsets, size_t offsets_len,
                   const Analysis& analysis, std::vector<size_t>& out) {
    scan::execute(analysis.cover, codes, len, row_offsets, offsets_len, analysis.covered_frequency,
                  analysis.total_frequency, AnyHit{}, out);
}

}  // namespace fsst::search::prefilter
