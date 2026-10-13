#pragma once

// Codes to bit mask (matcher/), bit mask to rows (resolver/), and the block
// driver joining them. See docs/prefilter.md section 6.

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <vector>

#include "../cover.hpp"

namespace fsst::search::prefilter::scan {

constexpr size_t BLOCK = 512;
// A block buffer carries 16 codes of padding: the pair compare reads one
// code past the block.
using Block = std::array<uint8_t, BLOCK + 16>;
// One mask word per 64 codes of a block.
template <size_t W>
using Bits = std::array<uint64_t, W>;
using Mask = Bits<BLOCK / 64>;
// Mask words per block of a streamed selection, whose rows are short spans.
constexpr size_t SELECTED_WORDS = 1;

// A check over an interleaved selection layer (OddRows): it is told the
// row, and `streamed` places the blocks at the selected rows.
template <typename C, typename = void>
struct on_layer : std::false_type {};
template <typename C>
struct on_layer<C, std::void_t<decltype(C::streamed)>> : std::true_type {};

template <size_t W>
void clear_from(Bits<W>& bits, size_t valid) {
    if (valid >= 64 * W) return;
    bits[valid / 64] &= (uint64_t{1} << (valid % 64)) - 1;
    for (size_t w = valid / 64 + 1; w < bits.size(); ++w) bits[w] = 0;
}

// step(block, at, valid): 64 * W codes starting at stream index `at`, of
// which the first `valid` are the ones to scan. The code after a block is
// readable: the stream's next code, or the padding of the copied last block.
template <size_t W, typename Step>
void block_at(const uint8_t* codes, size_t len, size_t at, size_t valid, Step step) {
    if (at + 64 * W < len) return step(codes + at, at, valid);
    std::array<uint8_t, 64 * W + 16> tail{};
    std::memcpy(tail.data(), codes + at, len - at);
    step(tail.data(), at, valid);
}

template <size_t W, typename Step>
void blocks(const uint8_t* codes, size_t len, Step step) {
    for (size_t at = 0; at < len; at += 64 * W) block_at<W>(codes, len, at, std::min(64 * W, len - at), step);
}

// blocks over the odd rows of interleaved offsets only. A block runs from the
// first selected row not yet scanned to the end of the last selected row
// starting inside it, and `valid` stops there.
template <size_t W, typename O, typename Step>
void selected_blocks(const uint8_t* codes, size_t len, const O* row_offsets, size_t offsets_len, Step step) {
    const size_t rows = offsets_len / 2 - 1;
    auto start = [&](size_t k) { return static_cast<size_t>(row_offsets[2 * k + 1]); };
    auto end = [&](size_t k) { return static_cast<size_t>(row_offsets[2 * k + 2]); };
    size_t k = 0, scanned = 0;
    while (k < rows) {
        const size_t from = std::max(start(k), scanned);
        if (end(k) <= from) {
            ++k;
            continue;
        }
        while (k + 1 < rows && start(k + 1) < from + 64 * W) ++k;
        const size_t to = std::min(end(k), from + 64 * W);
        block_at<W>(codes, len, from, to - from, step);
        scanned = to;
    }
}

// row_offsets has offsets_len entries, the last one equal to len.
template <typename M, typename R, typename Check>
void both_stages(const ProbeCover& cover, const uint8_t* codes, size_t len,
                 const typename R::Offset* row_offsets, size_t offsets_len,
                 const Check& check, std::vector<size_t>& out) {
    M matcher(cover);
    R resolver(row_offsets, offsets_len);
    const bool pairs = !cover.pairs.empty();
    auto step = [&](auto& bits) {
        return [&](const uint8_t* block, size_t at, size_t valid) {
            if (matcher.check(block, bits)) {
                clear_from(bits, valid);
                // The stream's last code pairs with padding; only a covered code counts there.
                if (pairs && valid > 0 && at + valid == len) {
                    size_t last = valid - 1;
                    uint64_t bit = uint64_t{1} << (last % 64);
                    if (cover.contains(block[last]))
                        bits[last / 64] |= bit;
                    else
                        bits[last / 64] &= ~bit;
                }
                resolver.rows(bits, at, check, out);
            }
        };
    };
    if constexpr (on_layer<Check>::value) {
        if (check.streamed) {
            Bits<SELECTED_WORDS> bits{};
            selected_blocks<SELECTED_WORDS>(codes, len, row_offsets, offsets_len, step(bits));
            return;
        }
    }
    Mask bits{};
    blocks<BLOCK / 64>(codes, len, step(bits));
}

}  // namespace fsst::search::prefilter::scan
