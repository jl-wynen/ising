#include "ising.hpp"

#include "catch.hpp"

#include "rng.hpp"

TEST_CASE("Parameters can be constructed properly", "[Parameters]")
{
    REQUIRE(Parameters{1.0, 0.0}.JT == 1.0);
    REQUIRE(Parameters{1.0, 0.0}.hT == 0.0);

    REQUIRE(Parameters{-3.1, 2.6}.JT == -3.1);
    REQUIRE(Parameters{-3.1, 2.6}.hT == 2.6);

    REQUIRE(Parameters{0.0, -1.23}.JT == 0.0);
    REQUIRE(Parameters{0.0, -1.23}.hT == -1.23);
}

TEST_CASE("Hamiltonian", "[Ising]")
{
    std::vector<std::vector<Index>> const shapes{
        {8_i},
        {32_i, 16_i},
        {16_i, 16_i, 16_i, 16_i},
        {8_i, 4_i, 8_i, 16_i, 5_i}
    };
    constexpr size_t nsamples = 10;

    SECTION("For J=0, hamiltonian is like magnetisation")
    {
        Rng rng(1_i, 6274);
        for (auto const &shape : shapes) {
            Lattice const lat{shape};
            rng.setLatsize(size(lat));

            for (size_t sample = 0; sample < nsamples; ++sample) {
                Parameters params{0.0, -0.7+sample*0.13};
                auto cfg = randomCfg(size(lat), rng);
                REQUIRE(hamiltonian(cfg, params, lat)
                        == Approx(-params.hT*magnetisation(cfg)*size(lat).get()));
            }
        }
    }

    SECTION("All spins aligned gives -(ndim*J+h)*sum s")
    {
        for (auto const &shape : shapes) {
            Lattice const lat{shape};

            for (size_t sample = 0; sample < nsamples; ++sample) {
                Parameters params{1.1-2*sample, -0.7+sample*0.13};
                Configuration cfg{size(lat), Spin{+1}};

                REQUIRE(hamiltonian(cfg, params, lat)
                        == Approx(-(lat.ndim().get()*params.JT + params.hT)
                                  * std::accumulate(begin(cfg), end(cfg), Spin{0}).get()));
            }
        }
    }
}
