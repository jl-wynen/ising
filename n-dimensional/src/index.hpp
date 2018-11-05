#ifndef ISING_INDEX_HPP
#define ISING_INDEX_HPP

#include <cstddef>

#include "arithmetic.hpp"

/// Type for indexes into configurations.
using Index = ArithmeticType<size_t, struct IndexTag>;

constexpr Index operator""_i(unsigned long long const x) noexcept
{
    return Index{x};
}

#endif  // ndef ISING_INDEX_HPP
