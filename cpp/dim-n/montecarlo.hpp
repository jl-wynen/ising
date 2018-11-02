#ifndef ISING_MONTECARLO_HPP
#define ISING_MONTECARLO_HPP

#include <vector>
#include <tuple>

#include "configuration.hpp"
#include "ising.hpp"
#include "rng.hpp"

/// Store Monte-Carlo history of observables.
struct Observables
{
    std::vector<double> energy;
    std::vector<double> magnetisation;
};


std::tuple<Configuration, double, size_t>
evolve(Configuration cfg, double energy, Parameters const& params,
       Lattice const &lat, Rng &rng, size_t const nsweep, Observables * const obs);

#endif  // ndef ISING_MONTECARLO_HPP
