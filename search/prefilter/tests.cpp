// The graph, the cut and the planned cover against byte containment on
// FSST-compressed rows, plus the min cut on synthetic graphs.

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

#include "fsst.h"
#include "prefilter.hpp"
#include "scan/tests.hpp"

using namespace fsst::search::prefilter;

struct Column {
    std::vector<uint8_t> codes;
    std::vector<uint32_t> offsets;
    Dictionary dict;
    Frequency freq;
    std::vector<std::string> rows;

    size_t row_count() const { return offsets.size() - 1; }
};

static Column compress_rows(const std::vector<std::string>& rows) {
    const size_t n = rows.size();
    std::vector<size_t> len(n);
    std::vector<const unsigned char*> str(n);
    size_t raw_bytes = 0;
    for (size_t i = 0; i < n; ++i) {
        len[i] = rows[i].size();
        str[i] = reinterpret_cast<const unsigned char*>(rows[i].data());
        raw_bytes += len[i];
    }
    fsst_encoder_t* encoder = fsst_create(n, len.data(), str.data(), 0);
    fsst_sort_codes(encoder);
    std::vector<unsigned char> out(2 * raw_bytes + 8 * n + 8);
    std::vector<size_t> len_out(n);
    std::vector<unsigned char*> str_out(n);
    size_t done = 0;
    unsigned char* cursor = out.data();
    while (done < n) {
        size_t count = fsst_compress(encoder, n - done, len.data() + done, str.data() + done,
                                     static_cast<size_t>(out.data() + out.size() - cursor), cursor,
                                     len_out.data() + done, str_out.data() + done);
        CHECK(count > 0);
        if (count == 0) break;
        done += count;
        cursor = str_out[done - 1] + len_out[done - 1];
    }
    Column col;
    col.rows = rows;
    col.codes.assign(out.data(), cursor);
    for (size_t i = 0; i < n; ++i) col.offsets.push_back(static_cast<uint32_t>(str_out[i] - out.data()));
    col.offsets.push_back(static_cast<uint32_t>(col.codes.size()));

    unsigned char header[FSST_MAXHEADER];
    fsst_export(encoder, header);
    fsst_decoder_t decoder;
    fsst_import(&decoder, header);
    col.dict = Dictionary::from(decoder.symbol, decoder.len, symbol_count(header));
    CHECK(col.dict.sorted());
    fsst_destroy(encoder);
    col.freq = Frequency::of(col.codes.data(), col.codes.size());

    std::vector<unsigned char> plain(raw_bytes + 8);
    for (size_t i = 0; i < n; ++i) {
        size_t got = fsst_decompress(&decoder, len_out[i], str_out[i], plain.size(), plain.data());
        CHECK_MSG(got == len[i] && std::memcmp(plain.data(), str[i], got) == 0, "row %zu does not round trip", i);
    }
    return col;
}

static bool byte_contains(const std::string& hay, const std::string& needle) {
    return hay.find(needle) != std::string::npos;
}

static std::vector<size_t> rows_containing(const Column& col, const std::string& pat) {
    std::vector<size_t> want;
    for (size_t r = 0; r < col.row_count(); ++r)
        if (byte_contains(col.rows[r], pat)) want.push_back(r);
    return want;
}

static const uint8_t* bytes(const std::string& s) { return reinterpret_cast<const uint8_t*>(s.data()); }

template <typename Blocked>
static bool sink_reachable_avoiding(const AlignmentGraph& g, Blocked blocked) {
    std::vector<std::vector<const Edge*>> adjacency(g.node_count());
    for (const Edge& e : g.edges) adjacency[e.from].push_back(&e);
    std::vector<bool> seen(g.node_count(), false);
    std::vector<size_t> stack{g.source()};
    seen[g.source()] = true;
    while (!stack.empty()) {
        size_t node = stack.back();
        stack.pop_back();
        if (node == g.sink()) return true;
        for (const Edge* e : adjacency[node])
            if (!seen[e->to] && !blocked(*e)) {
                seen[e->to] = true;
                stack.push_back(e->to);
            }
    }
    return false;
}

static bool selected(const std::vector<const Edge*>& selection, const Edge& e) {
    return std::find(selection.begin(), selection.end(), &e) != selection.end();
}

static bool covers_every_match(const Column& col, const std::vector<const Edge*>& selection,
                               const std::vector<size_t>& want) {
    ProbeCover cover = from_edge_cut(selection);
    for (size_t row : want) {
        bool holds = false;
        for (size_t i = col.offsets[row]; i < col.offsets[row + 1]; ++i) holds |= cover.contains(col.codes[i]);
        if (!holds) return false;
    }
    return true;
}

static uint64_t by_frequency(const Edge& e) { return e.frequency; }

static size_t g_escape_edges = 0;

static void check_graph(const Column& col, const std::string& pat, const std::vector<size_t>& want) {
    const size_t n = pat.size();
    AlignmentGraph g = build_alignment_graph(col.dict, bytes(pat), n, col.freq);
    std::vector<const Edge*> probes;
    for (const Edge& e : g.edges) {
        CHECK_MSG(e.from < e.to && e.to <= n, "%s: edge %u -> %u", pat.c_str(), e.from, e.to);
        if (e.probe == Probe::Escape) {
            ++g_escape_edges;
            CHECK(e.to == e.from + 1 && e.codes.size() == 1 && e.codes[0] == ESCAPE && e.byte == bytes(pat)[e.from]);
        }
        if (e.probe == Probe::Range) {
            CHECK(e.codes.empty() && e.points == 0 && e.ranges == 1 && e.range.begin <= e.range.last);
            // Exactly the symbols the suffix is a prefix of, contiguous because the codes are sorted.
            size_t m = n - e.from;
            for (size_t c = 0; c < col.dict.count; ++c) {
                bool prefixed = col.dict.length(c) >= m && std::memcmp(col.dict.symbol(c), bytes(pat) + e.from, m) == 0;
                CHECK_MSG(prefixed == e.range.contains(static_cast<uint8_t>(c)), "%s: range %u-%u at code %zu",
                          pat.c_str(), e.range.begin, e.range.last, c);
            }
        } else {
            ProbeCover shape = ProbeCover::from_runs([&] {
                std::vector<CodeRange> runs;
                for (uint8_t c : e.codes) runs.push_back({c, c});
                return runs;
            }());
            CHECK(e.points == shape.points.size() && e.ranges == shape.ranges.size());
        }
        if (e.cuttable()) probes.push_back(&e);
    }
    CHECK_MSG(!sink_reachable_avoiding(g, [](const Edge& e) { return e.cuttable(); }),
              "%s: a source-to-sink path carries no probe", pat.c_str());
    for (const Edge* e : probes) {
        ProbeCover cover = from_edge_cut({e});
        size_t matched = 0;
        for (uint8_t c : col.codes) matched += cover.contains(c);
        CHECK_MSG(e->frequency == matched, "%s: probe %u -> %u reports %u, matches %zu", pat.c_str(), e->from, e->to,
                  e->frequency, matched);
    }
    CHECK_MSG(covers_every_match(col, probes, want), "%s: a matching row holds no probe code at all", pat.c_str());
    CHECK_MSG(g.edges.size() <= 2 * n + 16, "%s: %zu edges", pat.c_str(), g.edges.size());

    std::vector<const Edge*> cut = min_cut(g, by_frequency);
    CHECK_MSG(!sink_reachable_avoiding(g, [&](const Edge& e) { return selected(cut, e); }),
              "%s: the minimum cut leaves a path open", pat.c_str());
    CHECK_MSG(covers_every_match(col, cut, want), "%s: the minimum cut misses a matching row", pat.c_str());

    if (probes.size() <= 12) {
        uint64_t best = 0;
        for (const Edge* e : cut) best += e->frequency;
        for (uint32_t mask = 0; mask < (1u << probes.size()); ++mask) {
            std::vector<const Edge*> subset;
            uint64_t weight = 0;
            for (size_t bit = 0; bit < probes.size(); ++bit)
                if ((mask >> bit) & 1) {
                    subset.push_back(probes[bit]);
                    weight += probes[bit]->frequency;
                }
            if (weight < best)
                CHECK_MSG(sink_reachable_avoiding(g, [&](const Edge& e) { return selected(subset, e); }),
                          "%s: a cover of weight %llu beats the cut's %llu", pat.c_str(),
                          static_cast<unsigned long long>(weight), static_cast<unsigned long long>(best));
        }
    }
}

// Every pattern against byte containment: the graph's invariants, then the
// planned cover run as a superset.
static void check(const Column& col, const std::vector<std::string>& patterns) {
    for (const std::string& pat : patterns) {
        std::vector<size_t> want = rows_containing(col, pat);
        if (!pat.empty()) check_graph(col, pat, want);
        Analysis a = analyze(bytes(pat), pat.size(), col.dict, col.freq, col.row_count());
        std::vector<size_t> got;
        superset_rows(col.codes.data(), col.codes.size(), col.offsets.data(), col.offsets.size(), a, got);
        if (pat.empty()) {
            CHECK(a.matches_all && got.size() == col.row_count());
            for (size_t r = 0; r < got.size(); ++r) CHECK(got[r] == r);
            continue;
        }
        CHECK(a.covered_frequency == col.freq.of_cover(a.cover));
        CHECK(a.total_frequency == col.codes.size());
        for (size_t i = 1; i < got.size(); ++i) CHECK(got[i - 1] < got[i]);
        CHECK_MSG(std::includes(got.begin(), got.end(), want.begin(), want.end()), "%s: the superset misses a row",
                  pat.c_str());
        for (size_t row : got) {
            bool holds = false;
            for (size_t i = col.offsets[row]; i < col.offsets[row + 1]; ++i) holds |= a.cover.contains(col.codes[i]);
            CHECK_MSG(holds, "%s: row %zu emitted without a covered code", pat.c_str(), row);
        }
        if (!want.empty()) CHECK_MSG(!a.cover.empty(), "%s: empty cover with matching rows", pat.c_str());
        std::vector<size_t> exact;
        candidate_rows(col.codes.data(), col.codes.size(), col.offsets.data(), col.offsets.size(), col.dict, a, exact);
        CHECK_MSG(exact == want, "%s: the walk's rows differ from byte containment (%zu against %zu)", pat.c_str(),
                  exact.size(), want.size());
    }
}

static std::vector<std::string> url_rows(size_t count, Rng& rng) {
    static const char* hosts[] = {"www.example.com", "shop.example.org", "cdn.static.net", "api.service.io"};
    static const char* paths[] = {"page", "user", "item", "search", "img", "docs"};
    std::vector<std::string> rows;
    for (size_t i = 0; i < count; ++i) {
        std::string r = rng.below(4) == 0 ? "http://" : "https://";
        r += hosts[rng.below(4)];
        size_t depth = rng.below(4);
        for (size_t d = 0; d < depth; ++d) r += "/" + std::string(paths[rng.below(6)]) + std::to_string(rng.below(1000));
        if (rng.below(2)) r += "?user=" + std::to_string(rng.below(50)) + "&q=" + std::string(paths[rng.below(6)]);
        rows.push_back(r);
    }
    return rows;
}

static void sound_on_url_text() {
    Rng rng{3};
    std::vector<std::string> rows = url_rows(400, rng);
    rows.push_back("");
    rows.push_back(std::string(300, 'x') + "https://www.example.com/page1" + std::string(20, 'y'));
    Column col = compress_rows(rows);
    std::vector<std::string> patterns = {"",
                                         "e",
                                         "/",
                                         "://",
                                         "http",
                                         "https://",
                                         ".com/page",
                                         "https://www.example.com",
                                         "user=1",
                                         "q=docs",
                                         "zzz absent",
                                         "example.org/item",
                                         rows[5],
                                         rows.back(),
                                         rows.back().substr(1)};
    check(col, patterns);
}

static void sound_on_repetitive_text() {
    std::vector<std::string> rows = {std::string(16, 'a'),   "aaaaaaaa aaaaaaaa", "xaaaaay", "aaa",  "aa",
                                     "zzaazz",               "abababababab",      "ba",     "aab",  "baa",
                                     std::string(600, 'a'),  "appappapple",       "apple",  "app",  "papp",
                                     "abcabcabcabc",         "cab",               "bca",    "abcab"};
    Column col = compress_rows(rows);
    check(col, {"a", "aa", "aaa", "aaaa", "aaaaaaaaa", "ab", "ba", "aba", "bab", "zz", "za", "app", "apple", "pp",
                "ppa", "abc", "cabc", "bcab", std::string(64, 'a'), std::string(300, 'a'), std::string(601, 'a')});
}

static void sound_on_random_bytes() {
    Rng rng{7};
    std::vector<std::string> rows;
    for (size_t i = 0; i < 300; ++i) {
        std::string r(1 + rng.below(40), '\0');
        for (char& c : r) c = static_cast<char>(rng.below(256));
        rows.push_back(r);
    }
    rows.push_back(std::string("\xff\xff\xff\xff", 4));
    rows.push_back(std::string("\x00\xff\x00\xff", 4));
    Column col = compress_rows(rows);
    CHECK_MSG(col.freq.count[ESCAPE] > 0, "random bytes should escape");
    std::vector<std::string> patterns = {std::string(1, '\xff'), std::string("\xff\xff", 2), std::string(1, '\0'),
                                         std::string("\x00\xff", 2)};
    for (size_t i = 0; i < 40; ++i) {
        const std::string& row = rows[rng.below(300)];
        size_t len = 1 + rng.below(std::min<size_t>(row.size(), 12));
        size_t at = rng.below(row.size() - len + 1);
        patterns.push_back(row.substr(at, len));
    }
    for (size_t i = 0; i < 8; ++i) patterns.push_back(std::string(1, static_cast<char>(rng.below(256))));
    check(col, patterns);
}

// Text with a handful of bytes the table has no symbol for, so occurrences
// begin, end and continue through escapes.
static void sound_with_escapes_inside_the_needle() {
    Rng rng{11};
    std::vector<std::string> rows = url_rows(300, rng);
    std::string rare;
    for (int b = 0x80; b < 0x100; ++b) rare.push_back(static_cast<char>(b));
    rows.push_back(rare);
    rows.push_back("page\x80\x81page");
    rows.push_back("\x90user=1\x91");
    rows.push_back("x\xff\xffy");
    rows.push_back("\x80");
    rows.push_back("ab\x80");
    rows.push_back("\x80\x80\x80");
    Column col = compress_rows(rows);
    CHECK_MSG(col.freq.count[ESCAPE] > 0, "rare bytes should escape");
    size_t before = g_escape_edges;
    check(col, {std::string(1, '\x80'), "page\x80", "\x80\x81page", "e\x80\x81p", "\x90user", "user=1\x91",
                std::string("\xff\xff", 2), "x\xff", "\xffy", "\x80\x80", "b\x80", rare.substr(3, 9)});
    CHECK_MSG(g_escape_edges > before, "no pattern built an escape edge");
}

static void empty_probe_cover_appends_nothing() {
    Column col = compress_rows({"alpha", "beta", "gamma"});
    Analysis a = analyze(bytes("\x01\x02\x03"), 3, col.dict, col.freq, col.row_count());
    std::vector<size_t> got;
    superset_rows(col.codes.data(), col.codes.size(), col.offsets.data(), col.offsets.size(), a, got);
    CHECK(rows_containing(col, "\x01\x02\x03").empty());
    CHECK(std::includes(got.begin(), got.end(), got.begin(), got.end()));
}

static void false_zero_frequencies_cannot_hide_a_true_match() {
    std::vector<std::string> rows = {"alpha beta", "gamma", "alphabet soup", "delta"};
    Column col = compress_rows(rows);
    Frequency zeros;
    zeros.total = col.freq.total;
    std::string pat = "alpha";
    Analysis a = analyze(bytes(pat), pat.size(), col.dict, zeros, col.row_count());
    CHECK(!a.cover.empty());
    CHECK(a.covered_frequency == 0);
    std::vector<size_t> got;
    superset_rows(col.codes.data(), col.codes.size(), col.offsets.data(), col.offsets.size(), a, got);
    std::vector<size_t> want = rows_containing(col, pat);
    CHECK(std::includes(got.begin(), got.end(), want.begin(), want.end()));
    std::vector<size_t> exact;
    candidate_rows(col.codes.data(), col.codes.size(), col.offsets.data(), col.offsets.size(), col.dict, a, exact);
    CHECK(exact == want);
}

static void cover_probes_the_maximal_runs_of_its_membership() {
    ProbeCover pf = ProbeCover::from_runs({{7, 8}, {0, 0}, {3, 3}, {6, 7}, {1, 1}, {8, 8}});
    CHECK(pf.points == std::vector<uint8_t>{3});
    CHECK(pf.ranges.size() == 2 && pf.ranges[0].begin == 0 && pf.ranges[0].last == 1 && pf.ranges[1].begin == 6 &&
          pf.ranges[1].last == 8);
    bool members[] = {true, true, false, true, false, false, true, true, true};
    for (uint8_t code = 0; code < 9; ++code) CHECK_MSG(pf.contains(code) == members[code], "membership at %u", code);
    CHECK(ProbeCover::from_runs({}).empty());
    ProbeCover top = ProbeCover::from_runs({{254, 254}, {255, 255}});
    CHECK(top.points.empty() && top.ranges.size() == 1 && top.ranges[0].last == 255);
}

static void sweep_prices_at_or_below_the_frequency_cut() {
    Rng rng{5};
    Column col = compress_rows(url_rows(500, rng));
    scan::policy::Region region{col.codes.size(), col.row_count()};
    for (std::string pat : {"e", "://", ".com/page", "https://www.example.com", "user=1"}) {
        AlignmentGraph g = build_alignment_graph(col.dict, bytes(pat), pat.size(), col.freq);
        ProbeCover baseline = from_edge_cut(min_cut(g, by_frequency));
        double baseline_ns = scan::policy::scan_ns(baseline, col.freq.of_cover(baseline), region);
        Planned p = cheapest_cover(g, col.freq, region);
        CHECK(p.covered == col.freq.of_cover(p.cover));
        CHECK(p.scan_ns == scan::policy::scan_ns(p.cover, p.covered, region));
        CHECK_MSG(p.scan_ns <= baseline_ns, "%s: sweep %f against %f", pat.c_str(), p.scan_ns, baseline_ns);
    }
}

// Synthetic graphs for the cut alone. Uncuttable steps are SetTooBig.
static Edge synthetic(uint32_t from, uint32_t to, int64_t frequency) {
    if (frequency < 0) return Edge{from, to, Probe::SetTooBig, 0, {}, CodeRange{0, 0}, 0, 0, 0};
    return Edge{from, to, Probe::Point, 0, {0}, CodeRange{0, 0}, static_cast<uint32_t>(frequency), 1, 0};
}

using Steps = std::vector<std::pair<uint32_t, uint32_t>>;

static Steps steps(const std::vector<const Edge*>& cut) {
    Steps s;
    for (const Edge* e : cut) s.emplace_back(e->from, e->to);
    return s;
}

static void shared_suffix_beats_two_local_choices() {
    AlignmentGraph g{{synthetic(0, 1, -1), synthetic(0, 2, -1), synthetic(1, 3, 4), synthetic(2, 3, 4),
                      synthetic(3, 4, 6)},
                     4};
    CHECK(steps(min_cut(g, by_frequency)) == Steps({{3, 4}}));
}

static void disjoint_paths_are_cut_separately() {
    AlignmentGraph g{{synthetic(0, 1, -1), synthetic(1, 4, 5), synthetic(0, 2, -1), synthetic(2, 3, 9),
                      synthetic(3, 4, 0)},
                     4};
    CHECK(steps(min_cut(g, by_frequency)) == Steps({{1, 4}, {3, 4}}));
}

static void deep_chain_does_not_exhaust_the_stack() {
    const uint32_t len = 100000;
    AlignmentGraph g{{}, len - 1};
    for (uint32_t v = 0; v + 1 < len; ++v) g.edges.push_back(synthetic(v, v + 1, 7));
    g.edges[len / 2] = synthetic(len / 2, len / 2 + 1, 3);
    CHECK(steps(min_cut(g, by_frequency)) == Steps({{len / 2, len / 2 + 1}}));
}

int main() {
    shared_suffix_beats_two_local_choices();
    disjoint_paths_are_cut_separately();
    deep_chain_does_not_exhaust_the_stack();
    cover_probes_the_maximal_runs_of_its_membership();
    sound_on_url_text();
    sound_on_repetitive_text();
    sound_on_random_bytes();
    sound_with_escapes_inside_the_needle();
    empty_probe_cover_appends_nothing();
    false_zero_frequencies_cannot_hide_a_true_match();
    sweep_prices_at_or_below_the_frequency_cut();
    return finish("prefilter tests");
}
