// Encodes a raw corpus, one value per line, into an FSST code stream in the
// layout the sweeps and onpair's bench_needles read: one row per line, the
// codes as comma-separated decimals, escape markers and their literal bytes
// both written as codes. The symbol table goes beside it as <out>.symtab.
//
//   prefilter_corpus <raw.csv> <out.fsst.csv>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "fsst.h"

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: %s <raw.csv> <out.fsst.csv>\n", argv[0]);
        return 2;
    }
    std::ifstream in(argv[1], std::ios::binary);
    if (!in) {
        std::fprintf(stderr, "cannot read %s\n", argv[1]);
        return 1;
    }
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    std::vector<size_t> len;
    std::vector<const unsigned char*> str;
    size_t at = 0;
    while (at < text.size()) {
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

    std::string csv;
    csv.reserve(raw_bytes * 2);
    size_t codes = 0, escapes = 0;
    char digits[8];
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < len_out[i]; ++j) {
            unsigned char c = str_out[i][j];
            if (c == FSST_ESC) ++escapes;
            int w = std::snprintf(digits, sizeof digits, "%u", static_cast<unsigned>(c));
            if (j) csv.push_back(',');
            csv.append(digits, static_cast<size_t>(w));
        }
        csv.push_back('\n');
        codes += len_out[i];
    }
    std::ofstream(argv[2], std::ios::binary) << csv;

    unsigned char table[FSST_MAXHEADER];
    unsigned table_len = fsst_export(encoder, table);
    std::ofstream(std::string(argv[2]) + ".symtab", std::ios::binary)
        .write(reinterpret_cast<const char*>(table), table_len);
    fsst_destroy(encoder);

    std::printf("%zu rows, %zu bytes -> %zu codes (%zu escapes, %.1f%% of codes), ratio %.2f\n", n,
                raw_bytes, codes, escapes, 100.0 * escapes / codes, double(raw_bytes) / codes);
    return 0;
}
