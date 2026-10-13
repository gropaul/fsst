#pragma once

// The FSST symbol table as the planner reads it: code to bytes and length.
// See docs/prefilter.md section 3.

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <cstring>

#include "cover.hpp"

namespace fsst::search::prefilter {

constexpr size_t MAX_TOKEN_SIZE = 8;
constexpr uint8_t ESCAPE = 255;

// The first min(m, 8) bytes at `p` as a little-endian word, the rest zero.
inline uint64_t load_word(const uint8_t* p, size_t m) {
    uint64_t w = 0;
    std::memcpy(&w, p, m < 8 ? m : 8);
    return w;
}

// The low l bytes of a word.
inline uint64_t low_bytes(size_t l) { return l >= 8 ? ~uint64_t{0} : (uint64_t{1} << (8 * l)) - 1; }

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

    // The symbol as one little-endian word, byte i at bits 8i; a symbol is at
    // most 8 bytes, so every compare against it is a masked word compare.
    uint64_t word(size_t code) const {
        uint64_t w;
        std::memcpy(&w, bytes[code].data(), 8);
        return w;
    }

    // Whether the codes are in byte order of their symbols, a symbol before
    // its extensions: what fsst_sort_codes gives and the planner relies on.
    bool sorted() const {
        for (size_t c = 1; c < count; ++c) {
            size_t m = std::min(length(c - 1), length(c));
            int cmp = std::memcmp(symbol(c - 1), symbol(c), m);
            if (cmp > 0 || (cmp == 0 && length(c - 1) >= length(c))) return false;
        }
        return true;
    }

    // The symbols prefix[0..m) is a prefix of, the exact one included:
    // contiguous in a sorted table, found by binary search.
    CodeRange prefix_range(const uint8_t* prefix, size_t m) const {
        assert(m <= MAX_TOKEN_SIZE);
        const uint64_t key = order_word(load_word(prefix, m));
        const uint64_t top = ~order_word(~low_bytes(m));  // the high m bytes of an order word
        size_t first = lower_bound(key, m);
        size_t end = first;
        for (size_t lo = first, hi = count; lo < hi;) {  // the first code past the prefix
            size_t mid = (lo + hi) / 2;
            if ((order(mid) & top) <= key) lo = end = mid + 1; else hi = mid;
        }
        if (first >= end) return CodeRange::none();
        return CodeRange{static_cast<uint8_t>(first), static_cast<uint8_t>(end - 1)};
    }

    // The code of the symbol s[0..l), or `count` if the table has none.
    size_t find(const uint8_t* s, size_t l) const {
        const uint64_t key = order_word(load_word(s, l));
        size_t c = lower_bound(key, l);
        return c < count && order(c) == key && length(c) == l ? c : count;
    }

  private:
    // Byte order as integer order: the symbol's bytes big-endian, unused bytes zero.
    static uint64_t order_word(uint64_t little) { return __builtin_bswap64(little); }
    uint64_t order(size_t code) const { return order_word(word(code) & low_bytes(length(code))); }

    // The first code at or after (key, l) in table order: bytes, then a shorter symbol first.
    size_t lower_bound(uint64_t key, size_t l) const {
        size_t lo = 0, hi = count;
        while (lo < hi) {
            size_t mid = (lo + hi) / 2;
            uint64_t o = order(mid);
            if (o < key || (o == key && length(mid) < l)) lo = mid + 1; else hi = mid;
        }
        return lo;
    }
};

inline size_t symbol_count(const uint8_t* export_header) { return export_header[1]; }

}  // namespace fsst::search::prefilter
