#include "lattice.hpp"

#include <cmath>

#include "ndebug.hpp"

namespace
{
    Index latticeSize(MultiIndex const &shape)
    {
        Index size{1};
        for (Index d : shape) {
            size = size*d;
        }
        return size;
    }

    /// Increment multi dimensional index.
    void increment(MultiIndex &index,
                   MultiIndex const &shape)
    {
        Index const ndim = Index{std::size(index)};

        for (Index i = 0_i; i < ndim; ++i) {
            // go through in inverse order to get row-major layout
            Index const j = ndim-1_i-i;

            // try to increment at position i
            if (++index[j.get()] == shape[j.get()]) {
                // back off, increment was too much
                index[j.get()] = 0_i;
            }
            else {
                // it is fine, keep the index
                return;
            }
        }
        // wrap around and start from 0 again
        index[0] = 0_i;
    }

    /// Increment the index at a given position taking PBCs into account.
    MultiIndex incrementedAt(MultiIndex const &index,
                             Index const pos,
                             MultiIndex const &shape)
    {
        if constexpr (not ndebug) {
            if (pos.get() > std::size(shape)) {
                throw std::out_of_range("Index position to increment out of range");
            }
        }

        MultiIndex incremented{index};
        if (index[pos.get()] == shape[pos.get()]-1_i) {
            incremented[pos.get()] = 0_i;
        }
        else {
            incremented[pos.get()] = index[pos.get()]+1_i;
        }
        return incremented;
    }

    /// Decrement the index at a given position taking PBCs into account.
    MultiIndex decrementedAt(MultiIndex const &index,
                             Index const pos,
                             MultiIndex const &shape)
    {
        if constexpr (not ndebug) {
            if (pos.get() > std::size(shape)) {
                throw std::out_of_range("Index position to decrement out of range");
            }
        }

        MultiIndex decremented{index};
        if (index[pos.get()] == 0_i) {
            decremented[pos.get()] = shape[pos.get()]-1_i;
        }
        else {
            decremented[pos.get()] = index[pos.get()]-1_i;
        }
        return decremented;
    }

    /// Build list of nearest neighbours.
    std::vector<Index> makeNeighbourList(MultiIndex const &shape)
    {
        Index const latsize = latticeSize(shape);
        Index const ndim = Index{std::size(shape)};

        std::vector<Index> neighbours((2_i*ndim*latsize).get());
        MultiIndex index(ndim.get(), 0_i);

        for (Index i = 0_i; i < latsize; ++i) {
            for (Index d = 0_i; d < ndim; ++d) {
                neighbours[(2_i*ndim*i + 2_i*d).get()] = totalIndex(incrementedAt(index, d, shape), shape);
                neighbours[(2_i*ndim*i + 2_i*d+1_i).get()] = totalIndex(decrementedAt(index, d, shape), shape);
            }

            increment(index, shape);
        }

        return neighbours;
    }

    /// Return the minimum distance between two sites in a specific dimension.
    /// Parameters are passed as indices / number of sites for this one dimension.
    constexpr int mindist1d(Index const x0, Index const x1, Index const nx) noexcept
    {
        int const forward = static_cast<int>(((nx + x1 - x0) % nx).get());
        int const backward = static_cast<int>(nx.get()) - forward;
        return std::min(forward, backward);
    }

    int sqMindist(MultiIndex const &site0, MultiIndex const &site1,
                  MultiIndex const &shape,
                  int (*sqdist)(std::vector<int> const &))
    {
        std::vector<int> individual;
        individual.reserve(std::size(shape));
        for (size_t dim = 0; dim < std::size(shape); ++dim) {
            individual.emplace_back(mindist1d(site0[dim], site1[dim], shape[dim]));
        }
        return sqdist(individual);
    }

    /// Return square of Euclidean distance given individual differences per dimension.
    int sqEuclidean(std::vector<int> const &individual) noexcept
    {
        int sqDistance = 0;
        for (int const d : individual) {
            sqDistance += d*d;
        }
        return sqDistance;
    }

    /// Return square of Manhattan distance given individual differences per dimension.
    int sqManhattan(std::vector<int> const &individual) noexcept
    {
        int distance = 0;
        for (int const d : individual) {
            distance += d;
        }
        return distance * distance;
    }

    /// Cosntruct a map of all distances on a lattice to pairs of indices of sites with that distance.
    DistMap buildDistMap(MultiIndex const &shape,
                         std::optional<double> const maxDist,
                         int (*sqdist)(std::vector<int> const &))
    {
        size_t const ndim = std::size(shape);
        Index const latsize = latticeSize(shape);

        DistMap distmap;
        MultiIndex site0(ndim, 0_i);  // first index

        for (Index i = 0_i; i < latsize-1_i; ++i) {  // move site0 through the lattice
            MultiIndex site1 = site0;  // compute distances between this and site0
            for (Index j = i; j < latsize; ++j) {  // move site1 through the lattice
                int const dist = sqMindist(site0, site1, shape, sqdist);

                if (maxDist and std::sqrt(static_cast<double>(dist)) < maxDist.value()) {
                    if (distmap.find(dist) != std::end(distmap)) {
                        // had this distance before => append new pair
                        distmap.at(dist).emplace_back(i, j);
                    }
                    else {
                        // new distance => insert new vector
                        distmap.emplace(dist, std::vector(1, std::make_pair(i, j)));
                    }
                }

                increment(site1, shape);
            }

            increment(site0, shape);
        }

        return distmap;
    }
}

Lattice::Lattice(std::vector<Index> const &shape,
                 std::optional<double> const maxDist,
                 DistanceFn const distfn)
    : neighbourList_{makeNeighbourList(shape)},
      shape_{shape},
      size_{latticeSize(shape)},
      distMap_{buildDistMap(shape, maxDist,
                            distfn == DistanceFn::EUCLIDEAN ? sqEuclidean : sqManhattan)}
{ }
