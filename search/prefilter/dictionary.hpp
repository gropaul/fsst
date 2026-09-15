#pragma once

// The FSST symbol table as the planner reads it: code to bytes and length.
// See docs/prefilter.md section 3.

#include <array>
#include <cstddef>
#include <cstdint>

namespace fsst::search::prefilter {

constexpr size_t MAX_TOKEN_SIZE = 8;
constexpr uint8_t ESCAPE = 255;

struct Dictionary {
    size_t count = 0;  // codes 0..count-1 name symbols
    std::array<uint8_t, 255> len{};
    std::array<std::array<uint8_t, 8>, 255> bytes{};

    // From fsst_decoder_t::symbol (little-endian words) and ::len. fsst_import
    // keeps no symbol count; symbol_count reads it off the fsst_export header.
    static Dictionary from(const uint64_t* symbol, const uint8_t* length, size_t count) {
        Dictionary d;
        d.count = count;
        for (size_t c = 0; c < count; ++c) {
            d.len[c] = length[c];
            for (size_t i = 0; i < 8; ++i) d.bytes[c][i] = static_cast<uint8_t>(symbol[c] >> (8 * i));
        }
        return d;
    }

    const uint8_t* symbol(size_t code) const { return bytes[code].data(); }
    size_t length(size_t code) const { return len[code]; }
};

inline size_t symbol_count(const uint8_t* export_header) { return export_header[1]; }

}  // namespace fsst::search::prefilter
