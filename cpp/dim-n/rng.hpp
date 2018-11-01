#ifndef ISING_RNG_HPP
#define ISING_RNG_HPP

#include <random>

#include "lattice.hpp"
#include "configuration.hpp"

/// Helper class to handle a random number generator.
struct Rng
{
    /// Seed the rng and set up distributions.
    explicit Rng(Index const latsize,
                 unsigned long const seed)
        : rng{seed}, indexDist{0, latsize.get()-1}, realDist{0, 1}, spinDist{0, 1}
    { }

    /// Generate a random index into a configuration.
    Index genIndex()
    {
        return Index{indexDist(rng)};
    }

    /// Generate a random double in [0, 1).
    double genReal()
    {
        return realDist(rng);
    }

    /// Generate a random spin, one of {-1, +1}.
    Spin genSpin()
    {
        return spinDist(rng)==0 ? Spin{-1} : Spin{1};
    }

private:
    /// The generator.
    std::mt19937 rng;

    /// Distribution to generator lattice indices.
    std::uniform_int_distribution<typename Index::Underlying> indexDist;

    /// Distribution to generate floating point numbers in [0, 1).
    std::uniform_real_distribution<double> realDist;

    /// Distribution to generate spins, i.e. values 0 or 1.
    std::uniform_int_distribution<typename Spin::Underlying> spinDist;
};


/// Generate a random spin configuration.
inline Configuration randomCfg(Index const latsize, Rng &rng)
{
    Configuration cfg(latsize);
    for (Spin &s : cfg) {
        s = rng.genSpin();
    }
    return cfg;
}

#endif  // ndef ISING_RNG_HPP
