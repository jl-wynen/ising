#ifndef ISING_LATTICE_HPP
#define ISING_LATTICE_HPP

#include <vector>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <optional>
#include <algorithm>

#include "index.hpp"
#include "ndebug.hpp"


/// N-dimensional index of a lattice site.
using MultiIndex = std::vector<Index>;

/// Map of squared distances to vectors of pairs of lattice sites with that separation.
using DistMap = std::unordered_map<int, std::vector<std::pair<Index, Index>>>;

/// Represent an n-dimensional lattice with an arbitrary (but hyperrectangular) shape.
struct Lattice
{
    /// Identify a function to calculate distances on a lattice.
    enum class DistanceFn { EUCLIDEAN, MANHATTAN };

    /// Construct from a shape and configure distance map.
    /**
     * \param shape N-dimensional shape of the lattice to construct.
     * \param maxDist Construct distance map only up to this maximum (non squared) distance.
     * \param distfn Function to use to compute distances on the lattice.
     */
    explicit Lattice(MultiIndex const &shape,
                     std::optional<double> maxDist = std::optional<double>{},
                     DistanceFn distfn = DistanceFn::EUCLIDEAN);

    /// Construct from a shape.
    template <typename... Shape>
    explicit Lattice(Shape ...shape)
        : Lattice{{shape...}}
    {
        static_assert((... and std::is_same_v<Shape, Index>),
                      "Lattice shape must be given as variables of type Index");
    }

    /// Return total lattice size.
    Index size() const noexcept
    {
        return size_;
    }

    /// Return the lattice extend in dimention dim.
    Index extend(Index const dim) const noexcept(ndebug)
    {
        if constexpr (not ndebug) {
            if (dim.get() >= shape_.size()) {
                throw std::out_of_range("Dimension for Lattice::extend is out of range.");
            }
        }

        return shape_[dim.get()];
    }

    /// Return the full lattice shape.
    auto const &shape() const noexcept
    {
        return shape_;
    }

    /// Return the number of dimensions.
    Index ndim() const noexcept
    {
        return Index{std::size(shape_)};
    }

    /// Return the list of nearest neighbour indices.
    auto const &neighbourList() const noexcept
    {
        return neighbourList_;
    }

    /// Return the index of neighbour number `neigh` of site `site`.
    Index neighbour(Index const site, Index const neigh) const noexcept(ndebug)
    {
        if constexpr (not ndebug) {
            if (site >= size_) {
                throw std::out_of_range("Index site out of range in Lattice::neighbour().");
            }
            if (neigh >= 2_i*ndim()) {
                throw std::out_of_range("Neihbour number out of range in Lattice::neighbour().");
            }
        }

        return neighbourList_[(2_i*ndim()*site+neigh).get()];
    }

    /// Return a tuple of iterators to the beginning and one past end of neighbours of a site.
    auto neighbours(Index const site) const
    {
        // type to add onto interators
        using diff = decltype(std::cbegin(neighbourList_))::difference_type;

        Index const NDIM = ndim();
        auto begin = std::cbegin(neighbourList_)+static_cast<diff>((2_i*NDIM*site).get());
        auto end = begin + static_cast<diff>((2_i*NDIM).get());
        return std::make_tuple(std::move(begin), std::move(end));
    }

    /// Return sorted vector of all squared distances stored in the map.
    auto sqDistances() const
    {
        std::vector<int> distances;
        for (auto const [key, _] : distMap_) {
            distances.emplace_back(key);
        }
        std::sort(std::begin(distances), std::end(distances));
        return distances;
    }

    /// Return vector of pairs of all sites with given squared distance.
    auto const &pairsWithSqDistance(int const sqDistance) const
    {
        return distMap_.at(sqDistance);
    }

private:
    /// Indices of nearest neighbours.
    std::vector<Index> const neighbourList_;
    /// Shape of the lattice.
    MultiIndex const shape_;
    /// Total size of the lattice.
    Index const size_;
    /// Squared distances with all pairs of lattice sites witht hat separation up to some max distance.
    DistMap const distMap_;
};

/// Return the total size of a lattice.
inline Index size(Lattice const &lat) noexcept
{
    return lat.size();
}

/// Compute the flat index from a set of Ndim indices.
inline Index totalIndex(std::vector<Index> const &index,
                        std::vector<Index> const &shape) noexcept(ndebug)
{
    if constexpr (not ndebug) {
        if (std::size(index) != std::size(shape)) {
            throw std::runtime_error("Invalid number of indices");
        }
    }

    Index total{index[0]};
    for (Index i = 1_i; i < Index{std::size(shape)}; ++i) {
        total = total*shape[i.get()] + index[i.get()];
    }
    return total;
}

/// Compute the flat index from a set of Ndim indices.
inline Index totalIndex(std::vector<Index> const &index,
                        Lattice const &lat) noexcept(ndebug)
{
    return totalIndex(index, lat.shape());
}

#endif  // ndef ISING_LATTICE_HPP
