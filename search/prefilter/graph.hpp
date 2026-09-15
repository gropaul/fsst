#pragma once

// The alignment graph: every layout of the needle across symbol boundaries
// as one DAG, whose source-to-sink cuts are the sound covers. Node o is the
// parse position at needle[o..], node 0 the source, node n the sink. See
// docs/prefilter.md sections 4 and 9.

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <vector>

#include "cover.hpp"
#include "dictionary.hpp"
#include "frequency.hpp"

namespace fsst::search::prefilter {

constexpr size_t PROBE_SET_SIZE_LIMIT = 512;
constexpr size_t PROBE_SET_SIZE_LIMIT_K1 = 16;

// Point: the greedy symbol at an interior offset. Escape: the marker 255
// with the literal needle byte behind it. Range: the symbols a needle suffix
// is a prefix of, into the sink. Set: an alignment's first symbols out of
// the source, or the symbols holding the whole needle past their start.
// SetTooBig: an alignment whose set was not enumerated; a cut may not
// select it.
enum class Probe : uint8_t { Point, Escape, Range, Set, SetTooBig };

struct Edge {
    uint32_t from;
    uint32_t to;
    Probe probe;
    uint8_t byte;                // Escape only: the literal after the marker
    std::vector<uint8_t> codes;  // ascending; Point and Escape hold one
    uint32_t frequency;          // codes the probe matches in the stream
    uint32_t points;             // runs the codes merge to: the lambda weight's terms
    uint32_t ranges;

    bool cuttable() const { return probe != Probe::SetTooBig; }
};

struct AlignmentGraph {
    std::vector<Edge> edges;
    size_t needle_len;

    size_t node_count() const { return needle_len + 1; }
    uint32_t source() const { return 0; }
    uint32_t sink() const { return static_cast<uint32_t>(needle_len); }
};

// The cover a cut's probes form: every code a run of one, merged by from_runs.
inline ProbeCover from_edge_cut(const std::vector<const Edge*>& cut) {
    std::vector<CodeRange> runs;
    for (const Edge* e : cut) {
        assert(e->cuttable());
        for (uint8_t c : e->codes) runs.push_back({c, c});
    }
    return ProbeCover::from_runs(std::move(runs));
}

// For each alignment k >= 1, the symbols whose last k bytes are needle[..k]
// and are longer than k, and the symbols holding the whole needle at an
// offset past their start. Counts saturate one past the limit.
struct Candidates {
    std::array<size_t, MAX_TOKEN_SIZE> count{};
    std::array<std::vector<uint8_t>, MAX_TOKEN_SIZE> first;
    std::vector<uint8_t> contained;

    void record(size_t k, uint8_t code) {
        size_t limit = k == 1 ? PROBE_SET_SIZE_LIMIT_K1 : PROBE_SET_SIZE_LIMIT;
        if (count[k] < limit) {
            first[k].push_back(code);
            ++count[k];
        } else {
            count[k] = PROBE_SET_SIZE_LIMIT + 1;
        }
    }
};

inline Candidates alignment_candidates(const Dictionary& dict, const uint8_t* needle, size_t n) {
    Candidates c;
    size_t kmax = std::min(n, MAX_TOKEN_SIZE);
    for (size_t code = 0; code < dict.count; ++code) {
        const uint8_t* s = dict.symbol(code);
        size_t len = dict.length(code);
        for (size_t k = 1; k < kmax && k < len; ++k)
            if (std::memcmp(s + len - k, needle, k) == 0) c.record(k, static_cast<uint8_t>(code));
        for (size_t p = 1; p + n <= len; ++p)
            if (std::memcmp(s + p, needle, n) == 0) {
                c.contained.push_back(static_cast<uint8_t>(code));
                break;
            }
    }
    return c;
}

namespace detail {

struct Builder {
    const Dictionary& dict;
    const uint8_t* needle;
    size_t n;
    const Frequency& freq;
    std::vector<Edge> edges;
    std::vector<bool> built;

    void add(uint32_t from, uint32_t to, Probe probe, std::vector<uint8_t> codes, uint8_t byte = 0) {
        std::vector<CodeRange> runs;
        for (uint8_t c : codes) runs.push_back({c, c});
        ProbeCover shape = ProbeCover::from_runs(runs);
        Edge e{from, to, probe, byte, std::move(codes), 0, static_cast<uint32_t>(shape.points.size()),
               static_cast<uint32_t>(shape.ranges.size())};
        e.frequency = freq.of_codes(e.codes);
        edges.push_back(std::move(e));
    }

    // The symbols needle[o..] is a prefix of, the exact one included.
    std::vector<uint8_t> terminal_set(size_t o) const {
        std::vector<uint8_t> set;
        size_t m = n - o;
        if (m > MAX_TOKEN_SIZE) return set;
        for (size_t code = 0; code < dict.count; ++code)
            if (dict.length(code) >= m && std::memcmp(dict.symbol(code), needle + o, m) == 0)
                set.push_back(static_cast<uint8_t>(code));
        return set;
    }

    // The encoder's longest match restricted to the needle: the longest
    // symbol that is a prefix of needle[o..]. Symbols of three or more bytes
    // have distinct three-byte prefixes, so at most one can match and the
    // encoder takes it before any shorter one.
    bool greedy(size_t o, uint8_t& code, size_t& len) const {
        size_t m = n - o;
        len = 0;
        for (size_t c = 0; c < dict.count; ++c) {
            size_t l = dict.length(c);
            if (l > len && l <= m && std::memcmp(dict.symbol(c), needle + o, l) == 0) {
                code = static_cast<uint8_t>(c);
                len = l;
            }
        }
        return len > 0;
    }

    // The steps out of offset o; returns where the parse lands. No symbol
    // matching means the encoder escapes needle[o] and continues at o + 1.
    size_t build_state(size_t o) {
        uint32_t state = static_cast<uint32_t>(o);
        std::vector<uint8_t> terminal = terminal_set(o);
        if (!terminal.empty()) add(state, sink(), Probe::Range, terminal);
        uint8_t code;
        size_t len;
        if (greedy(o, code, len)) {
            size_t next = o + len;
            if (next < n)
                add(state, static_cast<uint32_t>(next), Probe::Point, {code});
            else
                assert(std::find(terminal.begin(), terminal.end(), code) != terminal.end());
            return next;
        }
        add(state, static_cast<uint32_t>(o + 1), Probe::Escape, {ESCAPE}, needle[o]);
        return o + 1;
    }

    void ensure_chain(size_t start) {
        size_t o = start;
        while (!built[o]) {
            built[o] = true;
            size_t next = build_state(o);
            if (next >= n) break;
            o = next;
        }
    }

    uint32_t sink() const { return static_cast<uint32_t>(n); }
};

}  // namespace detail

inline AlignmentGraph build_alignment_graph(const Dictionary& dict, const uint8_t* needle, size_t n,
                                            const Frequency& freq) {
    assert(n > 0);
    detail::Builder b{dict, needle, n, freq, {}, std::vector<bool>(n, false)};
    Candidates cand = alignment_candidates(dict, needle, n);
    size_t kmax = std::min(n, MAX_TOKEN_SIZE);
    for (size_t k = 0; k < kmax; ++k) {
        if (k != 0 && cand.count[k] == 0) continue;
        b.ensure_chain(k);
        if (k != 0) {
            if (cand.count[k] <= PROBE_SET_SIZE_LIMIT)
                b.add(0, static_cast<uint32_t>(k), Probe::Set, cand.first[k]);
            else
                b.add(0, static_cast<uint32_t>(k), Probe::SetTooBig, {});
        }
    }
    if (!cand.contained.empty()) b.add(0, b.sink(), Probe::Set, cand.contained);
    assert(b.edges.size() <= 2 * n + 16);
    return AlignmentGraph{std::move(b.edges), n};
}

}  // namespace fsst::search::prefilter
