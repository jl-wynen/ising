#include "rng.hpp"

#include "catch.hpp"
#include "util.hpp"

using IVec = std::vector<Index>;

namespace {

    /// Produce n values using generator and check that all satify predicate.
    template <typename Gen, typename Pred>
    bool produceAndCheck(Gen generator, Pred predicate, size_t const n)
    {
        using std::begin, std::end;

        auto values = produce(generator, n);
        return std::all_of(begin(values), end(values), predicate);
    }

    /// Check that and index is in the proper range.
    struct IndexPredicate
    {
        Index const latsize;

        bool operator()(Index const idx) const
        {
            // > 0 not needed for unsigned, it is there just to be sure.
            return idx >= 0_i and idx < latsize;
        }
    };
}

TEST_CASE("Rng produces numbers in the correct range", "Rng")
{
    constexpr Index latsize = 143_i;
    constexpr size_t ncheck = 500;

    Rng rng(latsize, 538);

    SECTION("Indices")
    {
        REQUIRE(produceAndCheck([&rng](){ return rng.genIndex(); },
                                IndexPredicate{latsize},
                                ncheck));
    }

    SECTION("Real")
    {
        REQUIRE(produceAndCheck([&rng](){ return rng.genReal(); },
                                [](double const real) { return real >= 0.0 and real < 1.0; },
                                ncheck));
    }

    SECTION("Spin")
    {
        REQUIRE(produceAndCheck([&rng](){ return rng.genSpin(); },
                                [](Spin const spin) { return spin == Spin{-1} or spin == Spin{+1}; },
                                ncheck));
    }
}

SCENARIO("Rng latsize can be reset", "[Rng]")
{
    constexpr Index firstLatsize = 143_i;
    constexpr Index greaterLatsize = 187_i;
    constexpr Index smallerLatsize = 187_i;
    constexpr size_t ncheck = 500;

    GIVEN("An rng with set lattice size")
    {
        Rng rng(firstLatsize, 538);
        auto generator = [&rng](){ return rng.genIndex(); };
        auto const firstSet = produce(generator, ncheck);

        WHEN("The lattice size is increased")
        {
            rng.setLatsize(greaterLatsize);

            THEN("genIndex produces values in a greater range")
            {
                REQUIRE(produceAndCheck(generator,
                                        IndexPredicate{greaterLatsize},
                                        ncheck));
            }

            // make sure the internal rng state does not get reset
            THEN("genIndex produces different values")
            {
                auto const greaterSet = produce(generator, ncheck);
                REQUIRE(firstSet != greaterSet);
            }
        }

        WHEN("The lattice size is decreased")
        {
            rng.setLatsize(smallerLatsize);

            THEN("genIndex produces values in a smaller range")
            {
                REQUIRE(produceAndCheck(generator,
                                        IndexPredicate{smallerLatsize},
                                        ncheck));
            }

            // make sure the internal rng state does not get reset
            THEN("genIndex produces different values")
            {
                auto const smallerSet = produce(generator, ncheck);
                REQUIRE(firstSet != smallerSet);
            }
        }
    }
}
