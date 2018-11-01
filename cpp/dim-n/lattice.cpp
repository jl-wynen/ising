#include "lattice.hpp"

#include "ndebug.hpp"


namespace
{
    // TODO is there a std lib function for this?
    Index latticeSize(std::vector<Index> const &shape)
    {
        Index size{1};
        for (Index d : shape) {
            size = size*d;
        }
        return size;
    }

    /// Increment multi dimensional index.
    void increment(std::vector<Index> &index,
                   std::vector<Index> const &shape)
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
    std::vector<Index> incrementedAt(std::vector<Index> const &index,
                                     Index const pos,
                                     std::vector<Index> const &shape)
    {
        if constexpr (not ndebug) {
            if (pos.get() > std::size(shape)) {
                throw std::out_of_range("Index position to increment out of range");
            }
        }

        std::vector<Index> incremented{index};
        if (index[pos.get()] == shape[pos.get()]-1_i) {
            incremented[pos.get()] = 0_i;
        }
        else {
            incremented[pos.get()] = index[pos.get()]+1_i;
        }
        return incremented;
    }

    /// Decrement the index at a given position taking PBCs into account.
    std::vector<Index> decrementedAt(std::vector<Index> const &index,
                                     Index const pos,
                                     std::vector<Index> const &shape)
    {
        if constexpr (not ndebug) {
            if (pos.get() > std::size(shape)) {
                throw std::out_of_range("Index position to decrement out of range");
            }
        }

        std::vector<Index> decremented{index};
        if (index[pos.get()] == 0_i) {
            decremented[pos.get()] = shape[pos.get()]-1_i;
        }
        else {
            decremented[pos.get()] = index[pos.get()]-1_i;
        }
        return decremented;
    }
}

/// Build list of nearest neighbours.
std::vector<Index> makeNeighbourList(std::vector<Index> const &shape)
{
    Index const latsize = latticeSize(shape);
    Index const ndim = Index{std::size(shape)};

    std::vector<Index> neighbours((2_i*ndim*latsize).get());
    std::vector<Index> index(ndim.get(), 0_i);

    for (Index i = 0_i; i < latsize; ++i) {
        for (Index d = 0_i; d < ndim; ++d) {
            neighbours[(2_i*ndim*i + 2_i*d).get()] = totalIndex(incrementedAt(index, d, shape), shape);
            neighbours[(2_i*ndim*i + 2_i*d+1_i).get()] = totalIndex(decrementedAt(index, d, shape), shape);
        }

        increment(index, shape);
    }

    return neighbours;
}


Lattice::Lattice(std::vector<Index> const &shape)
    : neighbourList_{makeNeighbourList(shape)},
      shape_{shape},
      size_{latticeSize(shape)}
{ }

