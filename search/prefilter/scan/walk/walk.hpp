#pragma once

// Exact verification of a hit in the compressed domain. The graph lists the
// steps the hit code can be; forward takes the greedy step out of each node
// until a terminal set, backward the one step into each node until the
// source or an entry the planner did not enumerate. Escapes are units of
// two codes, marker then literal; a position is a literal iff the run of
// 255s before it has odd length. See docs/prefilter.md sections 7 and 9.

#include <algorithm>
#include <array>
#include <cstring>
#include <vector>

#include "../../dictionary.hpp"
#include "../../graph.hpp"

namespace fsst::search::prefilter::scan {

class Walk {
  public:
    Walk() = default;

    Walk(const AlignmentGraph& graph, const uint8_t* needle) : nodes_(graph.node_count()) {
        std::memcpy(head_bytes_.data(), needle, std::min(graph.needle_len, MAX_TOKEN_SIZE));
        std::vector<Step> steps;
        for (const Edge& e : graph.edges) {
            uint16_t from = static_cast<uint16_t>(e.from), to = static_cast<uint16_t>(e.to);
            switch (e.probe) {
                case Probe::Point:
                    nodes_[from].has_step = true;
                    nodes_[from].step = e.codes[0];
                    nodes_[from].next = to;
                    steps.push_back({e.codes[0], false, 0, from, to});
                    break;
                case Probe::Escape:
                    nodes_[from].has_step = true;
                    nodes_[from].escape = true;
                    nodes_[from].step = e.byte;
                    nodes_[from].next = to;
                    steps.push_back({ESCAPE, true, e.byte, from, to});
                    break;
                case Probe::Range:
                    nodes_[from].has_terminal = true;
                    nodes_[from].terminal = e.range;
                    for (size_t c = e.range.begin; c <= e.range.last; ++c)
                        steps.push_back({static_cast<uint8_t>(c), false, 0, from, to});
                    break;
                case Probe::Set:
                    for (uint8_t c : e.codes) steps.push_back({c, false, 0, from, to});
                    break;
                case Probe::SetTooBig: nodes_[to].entry_len = static_cast<uint8_t>(to); break;
                case Probe::Pair: break;  // the cut graph's; the walk reads the alignment graph
            }
        }
        std::sort(steps.begin(), steps.end(), [](const Step& a, const Step& b) {
            return std::tie(a.code, a.from, a.to, a.byte) < std::tie(b.code, b.from, b.to, b.byte);
        });
        steps.erase(std::unique(steps.begin(), steps.end(),
                                [](const Step& a, const Step& b) {
                                    return a.code == b.code && a.from == b.from && a.to == b.to && a.byte == b.byte;
                                }),
                    steps.end());
        for (const Step& s : steps) ++head_[s.code + 1];
        for (size_t c = 0; c < 256; ++c) head_[c + 1] += head_[c];
        steps_ = std::move(steps);
    }

    // Whether an occurrence inside codes[row_start..row_end) has the unit at
    // `hit` as one of its parse steps.
    bool check(const Dictionary& dict, const uint8_t* codes, size_t row_start, size_t row_end, size_t hit) const {
        if (escapes_before(codes, row_start, hit) & 1) return false;
        uint8_t c = codes[hit];
        if (c == ESCAPE) {
            if (hit + 1 >= row_end) return false;
            uint8_t byte = codes[hit + 1];
            for (uint32_t s = head_[ESCAPE]; s < head_[ESCAPE + 1]; ++s) {
                const Step& step = steps_[s];
                if (step.byte == byte && (step.to == sink() || forward(step.to, hit + 2, codes, row_end)) &&
                    backward(step.from, hit, codes, row_start, dict))
                    return true;
            }
            return false;
        }
        for (uint32_t s = head_[c]; s < head_[c + 1]; ++s) {
            const Step& step = steps_[s];
            if ((step.to == sink() || forward(step.to, hit + 1, codes, row_end)) &&
                backward(step.from, hit, codes, row_start, dict))
                return true;
        }
        return false;
    }

  private:
    struct Node {
        bool has_step = false;
        bool escape = false;  // the step is the marker with literal `step`
        uint8_t step = 0;     // Point: the code; Escape: the literal byte
        uint16_t next = 0;
        uint8_t entry_len = 0;  // > 0 where the set into this node was too big
        bool has_terminal = false;
        CodeRange terminal{0, 0};

        bool terminal_has(uint8_t c) const { return has_terminal && terminal.contains(c); }
    };

    struct Step {
        uint8_t code;
        bool escape;
        uint8_t byte;
        uint16_t from;
        uint16_t to;
    };

    uint16_t sink() const { return static_cast<uint16_t>(nodes_.size() - 1); }

    // 255s in codes[row_start..i), counted back from i - 1.
    static size_t escapes_before(const uint8_t* codes, size_t row_start, size_t i) {
        size_t run = 0;
        while (i > row_start && codes[i - 1] == ESCAPE) {
            ++run;
            --i;
        }
        return run;
    }

    // The greedy path from `node` over codes[i..row_end) reaches the sink,
    // through a terminal set or, when the needle ends in a literal, an
    // escape step landing on it.
    bool forward(uint16_t node, size_t i, const uint8_t* codes, size_t row_end) const {
        while (i < row_end) {
            uint8_t c = codes[i];
            const Node& at = nodes_[node];
            if (at.terminal_has(c)) return true;
            if (!at.has_step) return false;
            if (at.escape) {
                if (c != ESCAPE || i + 1 >= row_end || codes[i + 1] != at.step) return false;
                i += 2;
            } else {
                if (c != at.step) return false;
                i += 1;
            }
            node = at.next;
            if (node == sink()) return true;
        }
        return false;
    }

    // Some path from the source reaches `node` with the unit ending at
    // end - 1 as its last step. At most one step of a code enters a node.
    bool backward(uint16_t node, size_t end, const uint8_t* codes, size_t row_start, const Dictionary& dict) const {
        while (node != 0) {
            if (end <= row_start) return false;
            uint8_t c = codes[end - 1];
            if (escapes_before(codes, row_start, end - 1) & 1) {
                const Step* step = step_into(ESCAPE, node, c);
                if (step == nullptr) return false;
                node = step->from;
                end -= 2;
                continue;
            }
            if (c == ESCAPE) return false;
            const Step* step = step_into(c, node, 0);
            if (step == nullptr) {
                const Node& at = nodes_[node];
                return at.entry_len > 0 && dict.length(c) >= at.entry_len &&
                       std::memcmp(dict.symbol(c) + dict.length(c) - at.entry_len, head_bytes_.data(), at.entry_len) == 0;
            }
            node = step->from;
            end -= 1;
        }
        return true;
    }

    const Step* step_into(uint8_t code, uint16_t node, uint8_t byte) const {
        for (uint32_t s = head_[code]; s < head_[code + 1]; ++s)
            if (steps_[s].to == node && steps_[s].byte == byte) return &steps_[s];
        return nullptr;
    }

    std::vector<Node> nodes_;
    std::vector<Step> steps_;
    std::array<uint32_t, 257> head_{};
    std::array<uint8_t, MAX_TOKEN_SIZE> head_bytes_{};  // needle[..8], what an entry test compares
};

// The walk as stage two's check over one stream.
struct WalkCheck {
    const Walk& walk;
    const Dictionary& dict;
    const uint8_t* codes;

    bool passes(size_t hit, size_t row_start, size_t row_end) const {
        return walk.check(dict, codes, row_start, row_end, hit);
    }
};

}  // namespace fsst::search::prefilter::scan
