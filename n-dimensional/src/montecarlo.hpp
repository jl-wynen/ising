#ifndef ISING_MONTECARLO_HPP
#define ISING_MONTECARLO_HPP

#include <vector>
#include <tuple>
#include <functional>

#include "configuration.hpp"
#include "ising.hpp"
#include "rng.hpp"

using Measurement = std::function<void(Configuration const&, double)>;

/// Store Monte-Carlo history of observables.
struct Observables
{
    std::vector<double> energy;
    std::vector<double> magnetisation;

    struct Correlator
    {
        std::vector<int> const sqDistances;
        std::vector<std::vector<double>> correlator;

        explicit Correlator(std::vector<int> &&sqd);
    } corr;

    explicit Observables(Lattice const &lat);
};

/// Evolve a configuration in Monte-Carlo time.
/**
 * \param cfg Starting configuration.
 * \param energy Starting energy.
 * \param params Physical parameters of the ensemble.
 * \param lat Lattice to run on, must be consistent with cfg.
 * \param rng Random number generator to use. Its internal
 *            state is advanced by this function.
 * \param nsweed Number of sweeps to perform. A sweep is a set of size(lat) updates.
 * \param obs Storage for measuring observables.
 *            Can be nullptr in which case no measurements are performed.
 * \param extraMeas Additional measurements to perform.
 *                  Each vector element is called after every sweep.
 *
 * \returns Tuple of
 *   - final configuration
 *   - final energy
 *   - acceptance rate.
 */
std::tuple<Configuration, double, double>
evolve(Configuration cfg, double energy, Parameters const& params,
       Lattice const &lat, Rng &rng, size_t const nsweep,
       Observables *obs, std::vector<Measurement> const & extraMeas={});

#endif  // ndef ISING_MONTECARLO_HPP
