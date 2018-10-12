#ifndef ISING_LATTICE_HPP
#define ISING_LATTICE_HPP

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "arithmetic.hpp"
#include "ndebug.hpp"

/// Type for indexes into configurations.
using Index = ArithmeticType<size_t, struct IndexTag>;

constexpr Index operator""_i(unsigned long long const x) noexcept
{
    return Index{x};
}

namespace {
    // overly complicated but awesome way to compute the prouct of all elements in an array
    template <typename Array, size_t ...Idxs>
    constexpr Index prod_impl(Array const &array, std::index_sequence<Idxs...>)
    {
        return (array[Idxs]*...);
    }
    template <typename T, size_t N>
    constexpr Index prod(std::array<T, N> const &array)
    {
        return prod_impl(array, std::make_index_sequence<N>());
    }
}


constexpr Index NDIM {2};
constexpr std::array<Index, NDIM.get()> LATSHAPE{Index{4ul},
                                                 Index{3ul}};
constexpr Index LATSIZE = prod(LATSHAPE);


namespace {
    template <size_t dim, typename ...Idxs>
    constexpr Index totalIdx_impl(Index const head, Idxs const ...tail) noexcept(ndebug)
    {
        if constexpr (not ndebug) {
            if (head >= LATSHAPE[dim])  // Index is unsigned => this checks also for 'head < 0'.
                throw std::out_of_range("Index is out of range.");
        }

        if constexpr (dim < NDIM.get()-1) {
             return head + LATSHAPE[dim]*(totalIdx_impl<dim+1>(tail...));
        }
        else {
            return head;
        }
    }
}

/// Compute the total lattice index from individual indices for each dimension.
/**
 * The layout is 'column-major' in the sense that the first index is fastest running.
 */
template <typename ...Idxs>
constexpr Index totalIdx(Idxs const ...idxs) noexcept(ndebug)
{
    static_assert(sizeof...(Idxs) == NDIM.get(),
                  "Need to pass as many indices as there are dimensions.");
    static_assert((std::is_same_v<Idxs, Index>&&...),
                  "Function only accepts Index as parameter.");

    return totalIdx_impl<0>(idxs...);
}

/// Apply periodic boundary conditions to idx.
constexpr Index applyPeriodicBC(Index const idx, Index const n) noexcept
{
    if (idx == n)
        return 0_i;
    else if (idx == static_cast<Index>(-1))
        return n-1_i;
    return idx;
}

#endif  // ndef ISING_LATTICE_HPP
