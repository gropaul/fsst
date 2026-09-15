#pragma once

// One step per row crossed: what a dense mask over long rows wants.

#include "shared.hpp"

namespace fsst::search::prefilter::scan::resolver {

template <typename O>
class LinearSeek {
   public:
    using Offset = O;

    LinearSeek(const O* row_offsets, size_t offsets_len) : cursor_(row_offsets, offsets_len) {}

    template <typename Check>
    void rows(const Mask& bits, size_t first_code, const Check& check, std::vector<size_t>& out) {
        cursor_.rows(bits, first_code, check, out, walk);
    }

   private:
    static size_t walk(const O* row_offsets, size_t, size_t from, size_t code) {
        size_t row = from;
        while (static_cast<size_t>(row_offsets[row + 1]) <= code) ++row;
        return row;
    }

    Cursor<O> cursor_;
};

}  // namespace fsst::search::prefilter::scan::resolver
