#pragma once

// Pattern to probe cover: build the alignment graph, unroll it for pairs,
// cut it under frequency + lambda * (points + 2 * ranges + 2 * pairs) for a
// sweep of lambda, price every distinct cut with scan_ns and keep the
// cheapest. See docs/prefilter.md section 5.

#include <algorithm>
#include <cstdint>
#include <vector>

#include "graph.hpp"
#include "mincut.hpp"
#include "scan/policy/policy.hpp"
#include "scan/walk/walk.hpp"

namespace fsst::search::prefilter {

struct Planned {
    ProbeCover cover;
    uint32_t covered;  // codes the cover matches in the indexed stream
    double scan_ns;
    scan::Walk walk;
};

inline uint64_t weight_at(const Edge& e, uint64_t lambda) {
    return e.frequency + lambda * (e.points + 2 * e.ranges + 2 * e.pairs);
}

// The narrowest cut first, since no lambda passes it, then the ladder from
// zero until it arrives there.
inline Planned cheapest_cover(const CutGraph& graph, const Frequency& freq, scan::policy::Region region) {
    const uint64_t ceiling = freq.total;
    MinCut solver(graph);
    Planned best{};
    bool priced = false;
    auto price = [&](const std::vector<uint32_t>& cut) {
        std::vector<const Edge*> edges;
        for (uint32_t at : cut) edges.push_back(&graph.edges[at]);
        ProbeCover cover = from_edge_cut(edges);
        uint32_t covered = freq.of_cover(cover);
        double ns = scan::policy::scan_ns(cover, covered, region);
        if (!priced || ns < best.scan_ns) {
            best = Planned{std::move(cover), covered, ns, {}};
            priced = true;
        }
    };
    std::vector<uint32_t> narrowest = solver.solve(graph, [&](const Edge& e) { return weight_at(e, ceiling + 1); });
    price(narrowest);
    std::vector<uint32_t> last;
    uint64_t lambda = 0;
    while (lambda <= ceiling) {
        const std::vector<uint32_t>& cut = solver.solve(graph, [&](const Edge& e) { return weight_at(e, lambda); });
        if (cut == narrowest) break;
        if (cut != last) {
            last = cut;
            price(last);
        }
        lambda = std::max<uint64_t>(lambda * 4, 1);
    }
    return best;
}

inline Planned plan(const Dictionary& dict, const uint8_t* pattern, size_t n, const Frequency& freq,
                    size_t row_count) {
    AlignmentGraph graph = build_alignment_graph(dict, pattern, n, freq);
    Planned p = cheapest_cover(CutGraph::with_pairs(graph, freq), freq, scan::policy::Region{freq.total, row_count});
    p.walk = scan::Walk(graph, pattern);
    return p;
}

}  // namespace fsst::search::prefilter
