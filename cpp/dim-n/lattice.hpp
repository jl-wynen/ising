#ifndef ISING_LATTICE_HPP
#define ISING_LATTICE_HPP

#include <array>
#include <vector>

#include "index.hpp"

/// Number of dimensions.
constexpr Index NDIM{2};
/// Shape of the lattice, NDIM numbers.
constexpr std::array<Index, NDIM.get()> LATSHAPE{Index{4ul},
                                                 Index{3ul}};

/// Total lattice size.
constexpr Index LATSIZE = [](auto &shape) {
                              Index size{1};
                              for (Index d : shape) {
                                  size = size*d;
                              }
                              return size;
                          }(LATSHAPE);

namespace detail_ {
    /// Build the neighbour index list.
    std::vector<Index> makeNeighbourList();
}

/// Return the list of neares neighbour indices.
inline std::vector<Index> const &neighbourList()
{
    static std::vector<Index> neighbours{detail_::makeNeighbourList()};
    return neighbours;
}

#endif  // ndef ISING_LATTICE_HPP
