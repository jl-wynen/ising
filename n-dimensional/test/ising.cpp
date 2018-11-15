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
    Rng rng(1_i, 6274);
    constexpr size_t nsamples = 10;

    SECTION("Manual 1 - Checkerboard")
    {
        // - + -
        // + - +
        // - + -
        Lattice const lat{{3_i, 3_i}};
        Configuration cfg(size(lat), Spin{+1});
        for (Index i = 0_i; i < size(cfg); i=i+2_i) {
            cfg[i] = Spin{-1};
        }

        for (size_t sample = 0; sample < nsamples; ++sample) {
            Parameters const params{rng.genReal()*4.0-2.0, rng.genReal()-0.5};
            REQUIRE(hamiltonian(cfg, params, lat) == Approx(params.JT*6.0 + params.hT));
        }
    }

    SECTION("Manual 2 - Cluster")
    {
        // + + - -
        // + + - -
        // - - - -
        // - - - -
        Lattice const lat{{4_i, 4_i}};
        Configuration cfg(size(lat), Spin{-1});
        cfg[0_i] = Spin{+1};
        cfg[1_i] = Spin{+1};
        cfg[4_i] = Spin{+1};
        cfg[5_i] = Spin{+1};

        for (size_t sample = 0; sample < nsamples; ++sample) {
            Parameters const params{rng.genReal()*3.0-1.5, rng.genReal()*2.0-1.0};
            INFO("Sample " << sample << " with JT=" << params.JT << " hT=" << params.hT);
            REQUIRE(hamiltonian(cfg, params, lat) == Approx(-params.JT*16.0 + params.hT*8.0));
        }
    }

    SECTION("Manual 3 - 3D stripes")
    {
        // + - +  |  + + +  |  - - -
        // + - +  |  + + +  |  + + +
        // + - +  |  + + +  |  - - -
        Lattice const lat{{3_i, 3_i, 3_i}};
        Configuration cfg(size(lat), Spin{+1});

        cfg[totalIndex({0_i, 1_i, 0_i}, lat.shape())] = Spin{-1};
        cfg[totalIndex({1_i, 1_i, 0_i}, lat.shape())] = Spin{-1};
        cfg[totalIndex({2_i, 1_i, 0_i}, lat.shape())] = Spin{-1};

        cfg[totalIndex({0_i, 0_i, 2_i}, lat.shape())] = Spin{-1};
        cfg[totalIndex({0_i, 1_i, 2_i}, lat.shape())] = Spin{-1};
        cfg[totalIndex({0_i, 2_i, 2_i}, lat.shape())] = Spin{-1};

        cfg[totalIndex({2_i, 0_i, 2_i}, lat.shape())] = Spin{-1};
        cfg[totalIndex({2_i, 1_i, 2_i}, lat.shape())] = Spin{-1};
        cfg[totalIndex({2_i, 2_i, 2_i}, lat.shape())] = Spin{-1};

        for (size_t sample = 0; sample < nsamples; ++sample) {
            Parameters const params{rng.genReal()*2.0-1.0, rng.genReal()*4.2-2.1};
            INFO("Sample " << sample << " with JT=" << params.JT << " hT=" << params.hT);
            REQUIRE(hamiltonian(cfg, params, lat) == Approx(-params.JT*29.0 - params.hT*9.0));
        }
    }

    std::vector<std::vector<Index>> const shapes{
        {8_i},
        {32_i, 16_i},
        {16_i, 16_i, 16_i, 16_i},
        {8_i, 4_i, 8_i, 16_i, 5_i}
    };


    SECTION("For J=0, hamiltonian is like magnetisation")
    {
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

TEST_CASE("Delta E", "[Ising]")
{
    std::vector<std::vector<Index>> const shapes{
        {3_i, 3_i},
        {32_i, 16_i},
        {5_i, 5_i, 5_i},
        {8_i, 4_i, 8_i, 5_i}
    };
    Rng rng(1_i, 6274);
    constexpr size_t nsamples = 10;

    SECTION("Delta E gives the same result as difference of hamiltonian")
    {
        for (auto const &shape : shapes) {
            Lattice const lat{shape};
            rng.setLatsize(size(lat));

            for (size_t sample = 0; sample < nsamples; ++sample) {
                Parameters const params{rng.genReal()*2.0-1.0, rng.genReal()*3.2-1.6};
                INFO("Sample " << sample << " with JT=" << params.JT << " hT=" << params.hT);
                Configuration cfg = randomCfg(size(lat), rng);
                double const energy = hamiltonian(cfg, params, lat);

                for (Index site = 0_i; site < size(cfg); ++site) {
                    cfg.flip(site);
                    double const flippedEnergy = hamiltonian(cfg, params, lat);
                    cfg.flip(site);
                    REQUIRE(deltaE(cfg, site, params, lat) == Approx(flippedEnergy-energy));
                }
            }
        }
    }
}
