// One pattern over one raw corpus, end to end: compress with sorted codes,
// analyze, run the exact and the superset scan, time them, check the exact
// rows against memmem on the raw rows.
//
//   prefilter_search <raw.csv> <pattern>

#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "../../prefilter.hpp"
#include "fsst.h"
#include "utils.hpp"

using namespace fsst::search::prefilter;
using fsst::search::prefilter::scan::bench::best;

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: %s <raw.csv> <pattern>\n", argv[0]);
        return 2;
    }
    std::ifstream in(argv[1], std::ios::binary);
    if (!in) {
        std::fprintf(stderr, "cannot read %s\n", argv[1]);
        return 1;
    }
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::string pattern = argv[2];

    std::vector<size_t> len;
    std::vector<const unsigned char*> str;
    for (size_t at = 0; at < text.size();) {
        size_t end = text.find('\n', at);
        if (end == std::string::npos) end = text.size();
        str.push_back(reinterpret_cast<const unsigned char*>(text.data() + at));
        len.push_back(end - at);
        at = end + 1;
    }
    const size_t n = str.size();
    size_t raw_bytes = 0;
    for (size_t l : len) raw_bytes += l;

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
        if (count == 0) {
            std::fprintf(stderr, "fsst_compress made no progress at row %zu\n", done);
            return 1;
        }
        done += count;
        cursor = str_out[done - 1] + len_out[done - 1];
    }
    const uint8_t* codes = out.data();
    const size_t code_count = static_cast<size_t>(cursor - out.data());
    std::vector<uint32_t> offsets(n + 1);
    for (size_t i = 0; i < n; ++i) offsets[i] = static_cast<uint32_t>(str_out[i] - out.data());
    offsets[n] = static_cast<uint32_t>(code_count);

    unsigned char header[FSST_MAXHEADER];
    fsst_export(encoder, header);
    fsst_decoder_t decoder;
    fsst_import(&decoder, header);
    Dictionary dict = Dictionary::from(decoder.symbol, decoder.len, symbol_count(header));
    fsst_destroy(encoder);
    Frequency freq = Frequency::of(codes, code_count);

    const uint8_t* pat = reinterpret_cast<const uint8_t*>(pattern.data());
    Analysis a = analyze(pat, pattern.size(), dict, freq, n);

    std::vector<size_t> exact, superset;
    double exact_s = best([&] {
        exact.clear();
        candidate_rows(codes, code_count, offsets.data(), offsets.size(), dict, a, exact);
    });
    double superset_s = best([&] {
        superset.clear();
        superset_rows(codes, code_count, offsets.data(), offsets.size(), a, superset);
    });

    std::vector<size_t> want;
    double memmem_s = best([&] {
        want.clear();
        for (size_t r = 0; r < n; ++r)
            if (memmem(str[r], len[r], pattern.data(), pattern.size()) != nullptr) want.push_back(r);
    });

    std::printf("%s: %zu rows, %zu bytes, %zu codes, %zu symbols\n", argv[1], n, raw_bytes, code_count, dict.count);
    std::printf("pattern '%s': cover %zu points, %zu ranges, %zu pairs, covers %u codes (%.4f%%), planned %.0f us\n",
                pattern.c_str(), a.cover.points.size(), a.cover.ranges.size(), a.cover.pairs.size(), a.covered_frequency,
                100.0 * a.covered_frequency / static_cast<double>(a.total_frequency), a.scan_ns / 1e3);
    for (uint8_t c : a.cover.points)
        std::printf("  point %3u '%.*s' (%u codes)\n", c, static_cast<int>(dict.length(c)),
                    reinterpret_cast<const char*>(dict.symbol(c)), freq.count[c]);
    for (CodeRange r : a.cover.ranges)
        std::printf("  range %3u-%3u '%.*s' .. '%.*s' (%u codes)\n", r.begin, r.last, static_cast<int>(dict.length(r.begin)),
                    reinterpret_cast<const char*>(dict.symbol(r.begin)), static_cast<int>(dict.length(r.last)),
                    reinterpret_cast<const char*>(dict.symbol(r.last)), freq.of_range(r));
    for (CodePair p : a.cover.pairs)
        std::printf("  pair  %3u %3u '%.*s' '%.*s' (%u codes estimated)\n", p.first, p.second,
                    static_cast<int>(dict.length(p.first)), reinterpret_cast<const char*>(dict.symbol(p.first)),
                    static_cast<int>(dict.length(p.second)), reinterpret_cast<const char*>(dict.symbol(p.second)),
                    freq.pair_estimate(p));
    std::printf("exact    %8.0f us  %6.2f GB/s raw  %6.2f GB/s codes  %zu rows  %s\n", exact_s * 1e6,
                raw_bytes / exact_s / 1e9, code_count / exact_s / 1e9, exact.size(), exact == want ? "== memmem" : "!= memmem");
    std::printf("superset %8.0f us  %6.2f GB/s raw  %6.2f GB/s codes  %zu rows\n", superset_s * 1e6,
                raw_bytes / superset_s / 1e9, code_count / superset_s / 1e9, superset.size());
    std::printf("memmem   %8.0f us  %6.2f GB/s raw  %zu rows (uncompressed rows, the reference)\n", memmem_s * 1e6,
                raw_bytes / memmem_s / 1e9, want.size());
    return exact == want ? 0 : 1;
}
