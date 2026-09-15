// Both resolvers against the rows a mask names, worked out from the row
// layer directly, over every density and row shape that decides between them.

#include <algorithm>
#include <string>
#include <vector>

#include "../tests.hpp"
#include "gallop_seek.hpp"
#include "linear_seek.hpp"

using namespace fsst::search::prefilter::scan;
using namespace fsst::search::prefilter::scan::resolver;

// `mask` covers several blocks back to back; `check` filters the bits.
template <typename Offset, typename Check>
static std::vector<size_t> expected(const std::vector<uint64_t>& mask, const std::vector<Offset>& offsets,
                                    const Check& check) {
    std::vector<size_t> rows;
    for (size_t word = 0; word < mask.size(); ++word) {
        uint64_t set = mask[word];
        while (set) {
            size_t code = word * 64 + static_cast<size_t>(__builtin_ctzll(set));
            set &= set - 1;
            size_t row = static_cast<size_t>(
                std::upper_bound(offsets.begin(), offsets.end(), code,
                                 [](size_t c, Offset o) { return c < static_cast<size_t>(o); }) -
                offsets.begin()) - 1;
            if (!check.passes(code, offsets[row], offsets[row + 1])) continue;
            if (rows.empty() || rows.back() != row) rows.push_back(row);
        }
    }
    return rows;
}

template <typename R, typename Check>
static std::vector<size_t> resolve(const std::vector<uint64_t>& mask,
                                   const std::vector<typename R::Offset>& offsets, const Check& check) {
    R resolver(offsets.data(), offsets.size());
    std::vector<size_t> out;
    for (size_t block = 0; block * (BLOCK / 64) < mask.size(); ++block) {
        Mask bits{};
        for (size_t w = 0; w < bits.size(); ++w) {
            size_t at = block * (BLOCK / 64) + w;
            bits[w] = at < mask.size() ? mask[at] : 0;
        }
        resolver.rows(bits, block * BLOCK, check, out);
    }
    return out;
}

struct RejectEven {
    bool passes(size_t code, size_t, size_t) const { return code % 2 == 1; }
};

template <typename Offset>
static void every_resolver(const std::string& name, const std::vector<uint64_t>& mask,
                           const std::vector<Offset>& offsets) {
    auto want = expected(mask, offsets, Superset{});
    CHECK_MSG(resolve<LinearSeek<Offset>>(mask, offsets, Superset{}) == want, "linear: %s", name.c_str());
    CHECK_MSG(resolve<GallopSeek<Offset>>(mask, offsets, Superset{}) == want, "gallop: %s", name.c_str());
    auto odd = expected(mask, offsets, RejectEven{});
    CHECK_MSG(resolve<LinearSeek<Offset>>(mask, offsets, RejectEven{}) == odd, "linear odd: %s",
              name.c_str());
    CHECK_MSG(resolve<GallopSeek<Offset>>(mask, offsets, RejectEven{}) == odd, "gallop odd: %s",
              name.c_str());
}

static const size_t CODES = 3 * BLOCK + 500;

static std::vector<uint32_t> uniform(size_t per_row) {
    std::vector<uint32_t> offsets{0};
    while (offsets.back() < CODES)
        offsets.push_back(static_cast<uint32_t>(std::min(CODES, offsets.back() + per_row)));
    return offsets;
}

static std::vector<uint32_t> ragged(uint64_t seed) {
    Rng rng{seed};
    std::vector<uint32_t> offsets{0};
    while (offsets.back() < CODES) {
        size_t len = rng.below(8) == 0 ? 0 : rng.below(40) + 1;
        offsets.push_back(static_cast<uint32_t>(std::min(CODES, offsets.back() + len)));
    }
    return offsets;
}

static std::vector<uint64_t> strided(size_t stride) {
    std::vector<uint64_t> mask(CODES / 64 + 1, 0);
    for (size_t code = 0; code < CODES; code += stride) mask[code / 64] |= uint64_t{1} << (code % 64);
    return mask;
}

static std::vector<uint64_t> random_mask(uint64_t seed, unsigned one_in) {
    Rng rng{seed};
    std::vector<uint64_t> mask(CODES / 64 + 1, 0);
    for (size_t code = 0; code < CODES; ++code)
        if (rng.below(one_in) == 0) mask[code / 64] |= uint64_t{1} << (code % 64);
    return mask;
}

static std::vector<std::pair<std::string, std::vector<uint32_t>>> layers() {
    return {{"uniform 1", uniform(1)},     {"uniform 3", uniform(3)},   {"uniform 12", uniform(12)},
            {"uniform 257", uniform(257)}, {"uniform 5000", uniform(5000)}, {"ragged a", ragged(1)},
            {"ragged b", ragged(2)},       {"one row", {0, static_cast<uint32_t>(CODES)}}};
}

static void agree_on_every_density() {
    for (auto& [name, offsets] : layers()) {
        every_resolver(name + " empty", std::vector<uint64_t>(CODES / 64 + 1, 0), offsets);
        for (size_t stride : {1, 2, 3, 7, 64, 65, 1000, 4095, 4096, 4097})
            every_resolver(name + " stride " + std::to_string(stride), strided(stride), offsets);
        for (unsigned one_in : {1u, 2u, 10u, 100u, 3000u})
            every_resolver(name + " random 1/" + std::to_string(one_in), random_mask(one_in, one_in),
                           offsets);
    }
}

static void a_row_hit_in_two_blocks() {
    std::vector<uint32_t> offsets{0, 100, BLOCK - 10, BLOCK + 10, BLOCK + 11, static_cast<uint32_t>(CODES)};
    std::vector<uint64_t> mask(CODES / 64 + 1, 0);
    auto set = [&](size_t code) { mask[code / 64] |= uint64_t{1} << (code % 64); };
    set(BLOCK - 5);
    set(BLOCK - 1);
    set(BLOCK);
    set(BLOCK + 9);
    set(BLOCK + 10);
    every_resolver("two blocks", mask, offsets);
    CHECK((resolve<LinearSeek<uint32_t>>(mask, offsets, Superset{}) == std::vector<size_t>{2, 3}));
}

static void the_ends_of_the_stream() {
    for (auto& [name, offsets] : layers()) {
        std::vector<uint64_t> mask(CODES / 64 + 1, 0);
        mask[0] |= 1;
        mask[(CODES - 1) / 64] |= uint64_t{1} << ((CODES - 1) % 64);
        every_resolver(name + " ends", mask, offsets);
    }
}

static void one_hit_after_many_rows() {
    std::vector<uint32_t> offsets = uniform(1);
    std::vector<uint64_t> mask(CODES / 64 + 1, 0);
    mask[(CODES - 1) / 64] |= uint64_t{1} << ((CODES - 1) % 64);
    every_resolver("one hit at the end of 12k rows", mask, offsets);
    CHECK((resolve<GallopSeek<uint32_t>>(mask, offsets, Superset{}) == std::vector<size_t>{CODES - 1}));
}

static void wide_offsets() {
    for (auto& [name, offsets32] : layers()) {
        std::vector<uint64_t> offsets(offsets32.begin(), offsets32.end());
        every_resolver(name + " u64", strided(13), offsets);
        every_resolver(name + " u64 dense", strided(1), offsets);
    }
}

int main() {
    agree_on_every_density();
    a_row_hit_in_two_blocks();
    the_ends_of_the_stream();
    one_hit_after_many_rows();
    wide_offsets();
    return finish("prefilter resolver tests");
}
