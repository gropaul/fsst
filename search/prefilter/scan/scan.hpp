#pragma once

// Codes to bit mask (matcher/), bit mask to rows (resolver/), and the block
// driver joining them. See docs/prefilter.md section 6.

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include "../cover.hpp"

namespace fsst::search::prefilter::scan {

constexpr size_t BLOCK = 4096;
using Block = std::array<uint8_t, BLOCK>;
using Mask = std::array<uint64_t, BLOCK / 64>;

// Stage two asks this of every hit before its row goes out.
struct Superset {
    bool passes(size_t, size_t, size_t) const { return true; }
};

inline void clear_from(Mask& bits, size_t valid) {
    if (valid >= BLOCK) return;
    bits[valid / 64] &= (uint64_t{1} << (valid % 64)) - 1;
    for (size_t w = valid / 64 + 1; w < bits.size(); ++w) bits[w] = 0;
}

// step(block, at, valid): BLOCK codes starting at stream index `at`, of which
// the first `valid` are the stream's and the rest zero padding.
template <typename Step>
void blocks(const uint8_t* codes, size_t len, Step step) {
    size_t at = 0;
    for (; at + BLOCK <= len; at += BLOCK) step(codes + at, at, BLOCK);
    if (at < len) {
        Block tail{};
        std::memcpy(tail.data(), codes + at, len - at);
        step(tail.data(), at, len - at);
    }
}

// row_offsets has offsets_len entries, the last one equal to len.
template <typename M, typename R, typename Check>
void both_stages(const ProbeCover& cover, const uint8_t* codes, size_t len,
                 const typename R::Offset* row_offsets, size_t offsets_len,
                 const Check& check, std::vector<size_t>& out) {
    M matcher(cover);
    R resolver(row_offsets, offsets_len);
    Mask bits{};
    blocks(codes, len, [&](const uint8_t* block, size_t at, size_t valid) {
        if (matcher.check(block, bits)) {
            clear_from(bits, valid);
            resolver.rows(bits, at, check, out);
        }
    });
}

}  // namespace fsst::search::prefilter::scan
