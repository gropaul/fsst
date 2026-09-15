// Stage one: every kernel against every probe shape, timed over a prefix of
// the stream, one CSV row per pair. Mirrors onpair's matcher_fit sweep so the
// paper's plot_matcher.py reads the file unchanged.
//
//   prefilter_matcher_sweep
//
// A probe is K codes from the catalog, R ranges spread over the code space,
// or both. Only R is swept for ranges; the codes per range are the control.

#include <algorithm>
#include <cstdio>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "../matcher/eq_or.hpp"
#include "../matcher/nibble.hpp"
#include "../matcher/nibble_n8.hpp"
#include "../matcher/range.hpp"
#include "../matcher/table.hpp"
#include "../resolver/linear_seek.hpp"
#include "../scan.hpp"
#include "loader.hpp"
#include "utils.hpp"

using namespace fsst::search::prefilter;
using namespace fsst::search::prefilter::scan;
using namespace fsst::search::prefilter::scan::bench;

constexpr size_t RANGE_COUNTS[] = {1, 2, 3, 4, 6, 8, 12, 16, 24, 32};
constexpr size_t RANGE_WIDTHS[] = {1, 16, 256};
constexpr size_t RANGED_TOKEN_COUNTS[] = {1, 8, 16};
constexpr double RANGED_TARGET = 0.01;
constexpr size_t RANGED_WIDTH = 16;
constexpr size_t CODE_SPACE = 256;
constexpr size_t MAX_MEASURED_BATCHES = 32;

struct Probe {
    const NeedleSet* set;
    std::vector<CodeRange> ranges;
    size_t width;

    ProbeCover cover() const {
        return ProbeCover{set ? set->needles : std::vector<uint8_t>{}, ranges};
    }
};

// R disjoint ranges of `width` codes spread over the code space, or none
// where they do not fit.
static std::optional<std::vector<CodeRange>> spread_ranges(size_t count, size_t width) {
    size_t stride = CODE_SPACE / count;
    if (width > stride) return std::nullopt;
    std::vector<CodeRange> ranges;
    for (size_t at = 0; at < count; ++at) {
        size_t lo = at * stride;
        ranges.push_back({static_cast<uint8_t>(lo), static_cast<uint8_t>(lo + width - 1)});
    }
    return ranges;
}

static std::vector<Probe> probes(const std::vector<NeedleSet>& sets) {
    std::vector<Probe> out;
    for (const NeedleSet& set : sets) out.push_back({&set, {}, 0});
    for (size_t count : RANGE_COUNTS)
        for (size_t width : RANGE_WIDTHS)
            if (auto ranges = spread_ranges(count, width)) out.push_back({nullptr, *ranges, width});
    for (const NeedleSet& set : sets) {
        bool beside = std::find(std::begin(RANGED_TOKEN_COUNTS), std::end(RANGED_TOKEN_COUNTS), set.count) !=
                          std::end(RANGED_TOKEN_COUNTS) &&
                      std::abs(set.target - RANGED_TARGET) < 1e-9;
        if (!beside) continue;
        for (size_t count : RANGE_COUNTS)
            if (auto ranges = spread_ranges(count, RANGED_WIDTH)) out.push_back({&set, *ranges, RANGED_WIDTH});
    }
    return out;
}

struct Stream {
    const uint8_t* codes;
    size_t len;
    const uint8_t* checked;
    size_t checked_len;
    const std::vector<uint32_t>* check_rows;
};

struct Timed {
    double seconds;
    std::vector<size_t> found;
};

template <typename M>
static Timed time_kernel(const Probe& probe, const Stream& s) {
    ProbeCover cover = probe.cover();
    M matcher(cover);
    Mask bits{};
    double seconds = best([&] {
        blocks(s.codes, s.len, [&](const uint8_t* block, size_t, size_t) {
            matcher.check(block, bits);
            // The mask is thrown away; keep the compiler from throwing away the kernel with it.
            asm volatile("" : : "r"(bits.data()) : "memory");
        });
    });
    std::vector<size_t> found;
    both_stages<M, resolver::LinearSeek<uint32_t>>(cover, s.checked, s.checked_len, s.check_rows->data(),
                                                   s.check_rows->size(), Superset{}, found);
    return {seconds, std::move(found)};
}

using Run = std::function<std::optional<Timed>(const Probe&, const Stream&)>;

struct Kernel {
    const char* name;
    Run run;
};

template <typename M>
static Run only_if(std::function<bool(size_t, size_t)> takes) {
    return [takes](const Probe& probe, const Stream& s) -> std::optional<Timed> {
        ProbeCover cover = probe.cover();
        if (!takes(cover.points.size(), cover.ranges.size())) return std::nullopt;
        return time_kernel<M>(probe, s);
    };
}

// NibbleN8 at ceil(K / 8) batches, measured past the one batch the planner
// picks at u8 so the figure shows where the batches stop paying.
template <bool SKIP, size_t... B>
static std::optional<Timed> n8k_at(size_t batches, const Probe& probe, const Stream& s,
                                   std::index_sequence<B...>) {
    std::optional<Timed> out;
    (void)((batches == B + 1 ? (out = time_kernel<matcher::NibbleN8<B + 1, SKIP>>(probe, s), true) : false) ||
           ...);
    return out;
}

template <bool SKIP>
static Run nibble_n8k() {
    return [](const Probe& probe, const Stream& s) -> std::optional<Timed> {
        size_t k = probe.cover().points.size();
        if (k == 0) return std::nullopt;
        size_t batches = (k + matcher::PER_BATCH - 1) / matcher::PER_BATCH;
        return n8k_at<SKIP>(batches, probe, s, std::make_index_sequence<MAX_MEASURED_BATCHES>{});
    };
}

static std::vector<Kernel> kernels() {
    auto any = [](size_t, size_t) { return true; };
    auto tokens = [](size_t k, size_t) { return k > 0; };
    auto ranges_alone = [](size_t k, size_t r) { return k == 0 && r > 0; };
    return {
        {"table", only_if<matcher::Table>(any)},
        {"eq_or", only_if<matcher::EqOr<false>>(tokens)},
        {"eq_or_skip", only_if<matcher::EqOr<true>>(tokens)},
        {"range", only_if<matcher::Range<false>>(ranges_alone)},
        {"range_skip", only_if<matcher::Range<true>>(ranges_alone)},
        {"nibble_n8k", nibble_n8k<false>()},
        {"nibble_n8k_skip", nibble_n8k<true>()},
        {"nibble", only_if<matcher::Nibble<false>>(any)},
        {"nibble_skip", only_if<matcher::Nibble<true>>(any)},
    };
}

static std::string number(double v) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.15g", v);
    return buf;
}

static void measure(const std::string& stream, const std::vector<Kernel>& kernels, const std::string& machine,
                    std::vector<std::string>& out) {
    auto [corpus_path, needles_path] = paths(stream, ENCODING);
    Corpus corpus = load_corpus(corpus_path);
    size_t len = std::min(CODES, corpus.codes.size());
    std::vector<uint32_t> row_offsets = row_layer(corpus.row_offsets, len);
    std::vector<NeedleSet> sets;
    for (NeedleSet& set : load_needles(needles_path))
        if (set.sample == SAMPLE) sets.push_back(std::move(set));
    std::vector<Probe> all = probes(sets);
    std::printf("\n%s %s: %zu codes, %zu rows, %zu probes over %zu needle sets\n", stream.c_str(), ENCODING,
                len, row_offsets.size() - 1, all.size(), sets.size());

    size_t checked_len = std::min(CHECK_CODES, len);
    std::vector<uint32_t> check_rows = row_layer(corpus.row_offsets, checked_len);
    Stream s{corpus.codes.data(), len, corpus.codes.data(), checked_len, &check_rows};

    for (const Probe& probe : all) {
        ProbeCover cover = probe.cover();
        std::vector<size_t> expected;
        for (size_t row = 0; row + 1 < check_rows.size(); ++row)
            for (size_t i = check_rows[row]; i < check_rows[row + 1]; ++i)
                if (cover.contains(corpus.codes[i])) {
                    expected.push_back(row);
                    break;
                }
        for (const Kernel& kernel : kernels) {
            std::optional<Timed> timed = kernel.run(probe, s);
            if (!timed) continue;
            if (timed->found != expected) {
                std::fprintf(stderr, "%s differs from the cover on K=%zu R=%zu\n", kernel.name, cover.points.size(),
                             cover.ranges.size());
                std::exit(1);
            }
            size_t length = probe.set ? 1 : 0;
            size_t count = probe.set ? probe.set->count : 0;
            double target = probe.set ? probe.set->target : 0.0;
            double achieved = probe.set ? probe.set->achieved : 0.0;
            double gbs = static_cast<double>(len) / timed->seconds / 1e9;
            std::string row = stream + "," + ENCODING + "," + std::to_string(len) + "," +
                              std::to_string(row_offsets.size() - 1) + "," + kernel.name + "," +
                              std::to_string(length) + "," + std::to_string(count) + "," +
                              std::to_string(probe.ranges.size()) + "," + std::to_string(probe.width) + "," +
                              number(target) + "," + number(achieved) + "," +
                              number(static_cast<double>(timed->found.size()) / (check_rows.size() - 1)) + "," +
                              number(gbs) + "," + number(gbs) + "," + machine;
            out.push_back(std::move(row));
        }
    }
}

int main() {
    std::string name = machine();
    std::vector<std::string> rows;
    std::vector<Kernel> all = kernels();
    for (const char* stream : STREAMS) measure(stream, all, name, rows);
    auto path = write_csv("novel_mask",
                          "stream,encoding,codes,rows,matcher,length,count,ranges,width,target,achieved,"
                          "prefix_selectivity,gbs,gcodes,machine",
                          rows);
    std::printf("%zu rows -> %s\n", rows.size(), path.c_str());
    return 0;
}
