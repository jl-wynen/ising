#ifndef ISLE_TEST_UTIL_HPP
#define ISLE_TEST_UTIL_HPP

#include <vector>
#include <algorithm>

/// Create and return a vector of n elements produced using generator.
template <typename Gen>
auto produce(Gen generator, size_t const n)
{
    using std::begin, std::end;

    std::vector<decltype(generator())> values(n);
    std::generate(begin(values), end(values), generator);
    return values;
}

#endif  // ndef ISLE_TEST_UTIL_HPP
