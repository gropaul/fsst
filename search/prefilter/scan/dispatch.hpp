#pragma once

// A plan as the template parameters both_stages wants.

#include <vector>

#include "matcher/eq_or.hpp"
#include "matcher/nibble.hpp"
#include "matcher/nibble_n8.hpp"
#include "matcher/range.hpp"
#include "policy/policy.hpp"
#include "resolver/gallop_seek.hpp"
#include "resolver/linear_seek.hpp"
#include "scan.hpp"

namespace fsst::search::prefilter::scan::dispatch {

template <typename O, typename M, typename Check>
void with_resolver(policy::Plan plan, const ProbeCover& cover, const uint8_t* codes, size_t len, const O* row_offsets,
                   size_t offsets_len, const Check& check, std::vector<size_t>& out) {
    if (plan.resolver == policy::Resolve::LinearSeek)
        both_stages<M, resolver::LinearSeek<O>>(cover, codes, len, row_offsets, offsets_len, check, out);
    else
        both_stages<M, resolver::GallopSeek<O>>(cover, codes, len, row_offsets, offsets_len, check, out);
}

// The same kernel compiled with the movemask skipped, S, and without, P.
template <typename O, typename S, typename P, typename Check>
void with_skip(policy::Plan plan, const ProbeCover& cover, const uint8_t* codes, size_t len, const O* row_offsets,
               size_t offsets_len, const Check& check, std::vector<size_t>& out) {
    if (plan.skip)
        with_resolver<O, S>(plan, cover, codes, len, row_offsets, offsets_len, check, out);
    else
        with_resolver<O, P>(plan, cover, codes, len, row_offsets, offsets_len, check, out);
}

template <typename O, typename Check>
void run(policy::Plan plan, const ProbeCover& cover, const uint8_t* codes, size_t len, const O* row_offsets,
         size_t offsets_len, const Check& check, std::vector<size_t>& out) {
    using namespace matcher;
    switch (plan.matcher) {
        case policy::Match::EqOr:
            return with_skip<O, EqOr<true>, EqOr<false>>(plan, cover, codes, len, row_offsets, offsets_len, check, out);
        case policy::Match::Range:
            return with_skip<O, Range<true>, Range<false>>(plan, cover, codes, len, row_offsets, offsets_len, check,
                                                           out);
        case policy::Match::NibbleN8K:
            return with_skip<O, NibbleN8<1, true>, NibbleN8<1, false>>(plan, cover, codes, len, row_offsets,
                                                                       offsets_len, check, out);
        case policy::Match::Nibble:
            return with_skip<O, Nibble<true>, Nibble<false>>(plan, cover, codes, len, row_offsets, offsets_len, check,
                                                             out);
    }
}

}  // namespace fsst::search::prefilter::scan::dispatch
