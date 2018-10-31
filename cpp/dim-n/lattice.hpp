#ifndef ISING_LATTICE_HPP
#define ISING_LATTICE_HPP

#include <array>
#include <vector>

#include "index.hpp"

/// Number of dimensions.
constexpr Index NDIM{2};
/// Shape of the lattice, NDIM numbers.
constexpr std::array<Index, NDIM.get()> LATSHAPE{Index{3ul},
                                                 Index{3ul}};

/// Total lattice size.
constexpr Index LATSIZE = [](auto &shape) {
                              Index size{1};
                              for (Index d : shape) {
                                  size = size*d;
                              }
                              return size;
                          }(LATSHAPE);

/// Compute the flat index from a set of NDIM indices.
inline Index totalIndex(std::array<Index, NDIM.get()> const &index)
{
    Index total{index[0]};
    for (Index i = 1_i; i < NDIM; ++i) {
        total = total*LATSHAPE[i.get()] + index[i.get()];
    }
    return total;
}

namespace detail_ {
    /// Build the neighbour index list.
    std::vector<Index> makeNeighbourList();
}

/// Return the list of neares neighbour indices.
inline std::vector<Index> const &neighbourList()
{
    // TODO global variable to avoid static initialisation overhead?
    static std::vector<Index> neighbours{detail_::makeNeighbourList()};
    return neighbours;
}

#endif  // ndef ISING_LATTICE_HPP
