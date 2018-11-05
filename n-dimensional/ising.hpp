#ifndef ISING_ISING_HPP
#define ISING_ISING_HPP

#include <numeric>

#include "configuration.hpp"
#include "lattice.hpp"


/*
 * The energy is computed using each link only once.
 * That mean that if <x,y> is a pair of nearest neighbours,
 * the link x->y is counted but not y->x.
 * This can be absorbed into parameter J.
 */


/// Physical dimensionless parameters of the model.
struct Parameters
{
    double JT;  // J / (k_B T)
    double hT;  // h / (k_B T)
};

/// Sum spins of all neighbours of a given site.
inline Spin sumOfNeighbours(Configuration const &cfg,
                            Index const site,
                            Lattice const &lat) noexcept(ndebug)
{
    Spin neighbourSum{0};
    for (auto [it, end] = lat.neighbours(site);
         it != end;
         ++it)
    {
        neighbourSum = neighbourSum + cfg[*it];
    }
    return neighbourSum;
}

/// Evaluate the Hamiltonian on a configuration.
inline double hamiltonian(Configuration const &cfg,
                          Parameters const &params,
                          Lattice const &lat) noexcept(ndebug)
{
    Spin coupling{0};  // the nearest-neighbour coupling
    Spin magn{0};  // the sum over all sites

    for (Index i = 0_i; i < size(lat); ++i) {
        coupling = coupling + cfg[i]*sumOfNeighbours(cfg, i, lat);
        magn = magn + cfg[i];
    }

    // divide by 2 so each link is counted only once
    return -params.JT*coupling.get()/2.0 - params.hT*magn.get();
}

/// Compute the change in energy if the spin at site idx were flipped.
inline double deltaE(Configuration const &cfg, Index const site,
                     Parameters const &params,
                     Lattice const &lat) noexcept(ndebug)
{
    return 2.0*cfg[site].get()*(params.JT*sumOfNeighbours(cfg, site, lat).get()
                                + params.hT);
}

/// Compute the magnetisation on a configuration.
inline double magnetisation(Configuration const &cfg) noexcept(ndebug)
{
    return std::accumulate(begin(cfg), end(cfg), Spin{0}).get() /
        static_cast<double>(size(cfg).get());
}

#endif  // ndef ISING_ISING_HPP
