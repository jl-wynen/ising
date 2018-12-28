#include "lattice.hpp"

#include <cassert>

#include "catch.hpp"

using IVec = std::vector<Index>;

TEST_CASE("Construction of Lattice", "[Lattice]")
{
    SECTION("Constructing just stores the shape")
    {
        IVec shape;

        for (size_t i = 0; i < 5; ++i) {
            shape.emplace_back((i-2)*(i-2)+6);

            Lattice const lat{shape, 0.0};
            REQUIRE(lat.shape() == shape);
            REQUIRE(lat.ndim().get() == std::size(shape));

            Index product = 1_i;
            for (Index i : shape)
                product = product * i;
            REQUIRE(size(lat) == product);
        }
    }
}

namespace {
    /// Return a sorted vector of all neighbour indices of site.
    IVec sortedNeighbours(Lattice const &lat, Index const site)
    {
        auto const [begin, end] = lat.neighbours(site);
        IVec neighbours{begin, end};
        std::sort(std::begin(neighbours), std::end(neighbours));
        return neighbours;
    }
}

TEST_CASE("Nearest Neighbour Indices", "[Lattice]")
{
    using Catch::Matchers::VectorContains;

    SECTION("Manual simple 1D case")
    {
        Lattice const lat{{5_i}, 0.0};
        REQUIRE(sortedNeighbours(lat, 0_i) == IVec{1_i, 4_i});
        REQUIRE(sortedNeighbours(lat, 1_i) == IVec{0_i, 2_i});
        REQUIRE(sortedNeighbours(lat, 2_i) == IVec{1_i, 3_i});
        REQUIRE(sortedNeighbours(lat, 3_i) == IVec{2_i, 4_i});
        REQUIRE(sortedNeighbours(lat, 4_i) == IVec{0_i, 3_i});
    }

    SECTION("Manual simple 2D case")
    {
        Lattice const lat{{3_i, 3_i}, 0.0};
        REQUIRE(sortedNeighbours(lat, 0_i) == IVec{1_i, 2_i, 3_i, 6_i});
        REQUIRE(sortedNeighbours(lat, 1_i) == IVec{0_i, 2_i, 4_i, 7_i});
        REQUIRE(sortedNeighbours(lat, 2_i) == IVec{0_i, 1_i, 5_i, 8_i});
        REQUIRE(sortedNeighbours(lat, 3_i) == IVec{0_i, 4_i, 5_i, 6_i});
        REQUIRE(sortedNeighbours(lat, 4_i) == IVec{1_i, 3_i, 5_i, 7_i});
        REQUIRE(sortedNeighbours(lat, 5_i) == IVec{2_i, 3_i, 4_i, 8_i});
        REQUIRE(sortedNeighbours(lat, 6_i) == IVec{0_i, 3_i, 7_i, 8_i});
        REQUIRE(sortedNeighbours(lat, 7_i) == IVec{1_i, 4_i, 6_i, 8_i});
        REQUIRE(sortedNeighbours(lat, 8_i) == IVec{2_i, 5_i, 6_i, 7_i});
    }


    std::vector<IVec> const shapes{
        {8_i},
        {32_i, 16_i},
        {16_i, 16_i, 8_i, 24_i}
    };

    SECTION("All ways of getting neighbours are equivalent")
    {
        for (auto const &shape : shapes) {
            Lattice const lat{shape, 0.0};

            for (Index i = 0_i; i < lat.size(); ++i) {
                auto const [begin, end] = lat.neighbours(i);
                REQUIRE(end-begin == 2*lat.ndim().get());

                for (Index n = 0_i; n < 2_i*lat.ndim(); ++n) {
                    REQUIRE(lat.neighbour(i, n) == lat.neighbourList()[(2_i*lat.ndim()*i+n).get()]);
                    REQUIRE(lat.neighbour(i, n) == *(begin+n.get()));
                }
            }
        }
    }

    SECTION("Neighbour indices are symmetric")
    {
        for (auto const &shape : shapes) {
            Lattice const lat{shape, 0.0};

            for (Index i = 0_i; i < lat.size(); ++i) {

                for (Index n = 0_i; n < 2_i*lat.ndim(); ++n) {
                    Index const neighbour = lat.neighbour(i, n);
                    auto const [begin, end] = lat.neighbours(neighbour);
                    IVec const neighboursOfNeighbour{begin, end};

                    REQUIRE_THAT(neighboursOfNeighbour, VectorContains(i));
                }
            }
        }
    }
}

TEST_CASE("Lattice layout", "[Lattice]")
{
    SECTION("1D lattice layout is a flat vector")
    {
        IVec const shape{4_i};
        IVec const index{2_i};
        IVec const incremented{3_i};
        IVec const decremented{1_i};
        REQUIRE(totalIndex(incremented, shape) == totalIndex(index, shape)+1_i);
        REQUIRE(totalIndex(decremented, shape) == totalIndex(index, shape)-1_i);
    }

    SECTION("2D lattice layout is row-major")
    {
        Lattice const lat{{4_i, 7_i}, 0.0};
        IVec const index{{2_i, 4_i}};

        // in dimension 0 (row)
        IVec const incremented0{3_i, 4_i};
        IVec const decremented0{1_i, 4_i};
        REQUIRE(totalIndex(incremented0, lat) == totalIndex(index, lat)+1_i*lat.extend(1_i));
        REQUIRE(totalIndex(decremented0, lat) == totalIndex(index, lat)-1_i*lat.extend(1_i));

        // in dimension 1(column)
        IVec const incremented1{2_i, 5_i};
        IVec const decremented1{2_i, 3_i};
        REQUIRE(totalIndex(incremented1, lat) == totalIndex(index, lat)+1_i);
        REQUIRE(totalIndex(decremented1, lat) == totalIndex(index, lat)-1_i);
    }

    SECTION("ND lattice layout is 'generalised row-major'")
    {
        // test for those shapes
        std::vector<IVec> const shapes{
            {16_i, 16_i, 8_i},
            {32_i, 3_i, 4_i, 5_i},
            {8_i, 4_i, 8_i, 16_i, 32_i, 5_i}
        };
        // indices to centre tests around
        std::vector<IVec> const indices{
            {5_i, 2_i, 6_i},
            {17_i, 1_i, 2_i, 4_i},
            {1_i, 2_i, 4_i, 9_i, 24_i, 2_i}
        };
        assert(std::size(shapes) == std::size(indices));

        for (size_t i = 0; i < std::size(shapes); ++i) {
            auto const &shape = shapes.at(i);
            auto const &index = indices.at(i);
            assert(std::size(shape) == std::size(index));

            Index stride = 1_i;  // The distance index and increment are apart from each other.
                                 // Should just be prod of all shape[d'] for d' > d.
            for (size_t d = std::size(shape)-1; d != static_cast<size_t>(-1); --d) {
                IVec incremented{index};
                incremented.at(d) = incremented.at(d) + 1_i;
                IVec decremented{index};
                decremented.at(d) = decremented.at(d) - 1_i;

                REQUIRE(totalIndex(incremented, shape) == totalIndex(index, shape)+stride);
                REQUIRE(totalIndex(decremented, shape) == totalIndex(index, shape)-stride);

                stride = stride * shape.at(d);
            }
        }
    }
}
