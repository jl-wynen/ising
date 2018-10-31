#include "lattice.hpp"

#include <array>

#include "ndebug.hpp"

namespace detail_
{
    namespace
    {
        /// Increment multi dimensional index.
        void increment(std::array<Index, NDIM.get()> &index)
        {
            for (Index i = 0_i; i < NDIM; ++i) {
                // try to increment at position i
                if (++index[i.get()] == LATSHAPE[i.get()]) {
                    // back off, increment was too much
                    index[i.get()] = 0_i;
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
        std::array<Index, NDIM.get()> incrementedAt(
            std::array<Index, NDIM.get()> const &index,
            Index const pos)
        {
            if constexpr (not ndebug) {
                if (pos > NDIM) {
                    throw std::out_of_range("Index position to increment out of range");
                }
            }

            std::array<Index, NDIM.get()> incremented{index};
            if (index[pos.get()] == LATSHAPE[pos.get()]-1_i) {
                incremented[pos.get()] = 0_i;
            }
            else {
                incremented[pos.get()] = index[pos.get()]+1_i;
            }
            return incremented;
        }

        /// Decrement the index at a given position taking PBCs into account.
        std::array<Index, NDIM.get()> decrementedAt(
            std::array<Index, NDIM.get()> const &index,
            Index const pos)
        {
            if constexpr (not ndebug) {
                if (pos > NDIM) {
                    throw std::out_of_range("Index position to decrement out of range");
                }
            }

            std::array<Index, NDIM.get()> decremented{index};
            if (index[pos.get()] == 0_i) {
                decremented[pos.get()] = LATSHAPE[pos.get()]-1_i;
            }
            else {
                decremented[pos.get()] = index[pos.get()]-1_i;
            }
            return decremented;
        }
    }

    std::vector<Index> makeNeighbourList()
    {
        std::vector<Index> neighbours((2_i*NDIM*LATSIZE).get());
        std::array<Index, NDIM.get()> index;
        std::fill(std::begin(index), std::end(index), 0_i);

        for (Index i = 0_i; i < LATSIZE; ++i) {
            for (Index d = 0_i; d < NDIM; ++d) {
                neighbours[(2_i*NDIM*i + d).get()] = totalIndex(incrementedAt(index, d));
                neighbours[(2_i*NDIM*i + d+1_i).get()] = totalIndex(decrementedAt(index, d));
            }

            increment(index);
        }

        return neighbours;
    }
}
