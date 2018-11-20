#include "montecarlo.hpp"

#include <cmath>

using std::exp;

Observables::Observables(Lattice const &lat)
    : energy(), magnetisation(), corr(lat.sqDistances())
{ }

Observables::Correlator::Correlator(std::vector<int> &&sqd)
    : sqDistances(std::move(sqd)), correlator()
{
    correlator.resize(std::size(sqDistances));
}


namespace {
    void measureCorrelator(Observables::Correlator &corr, Lattice const &lat, Configuration const &cfg)
    {
        for (size_t sqdi = 0; sqdi < std::size(corr.sqDistances); ++sqdi) {
            Spin aux = Spin{0};
            size_t n = 0;
            for (auto const [i, j] : lat.pairsWithSqDistance(corr.sqDistances[sqdi])) {
                aux = aux + cfg[i]*cfg[j];
                ++n;
            }
            corr.correlator[sqdi].emplace_back(static_cast<double>(aux)/static_cast<double>(n));
        }
    }

    /// Measure observables if an instance of Observables is given.
    void measure(Observables * const obs, Lattice const &lat,
                 Configuration const &cfg, double const energy)
    {
        if (obs) {
            obs->energy.emplace_back(energy);
            obs->magnetisation.emplace_back(magnetisation(cfg));
            measureCorrelator(obs->corr, lat, cfg);
        }
    }
}


std::tuple<Configuration, double, double>
evolve(Configuration cfg, double energy, Parameters const& params,
       Lattice const &lat, Rng &rng, size_t const nsweep,
       Observables * const obs, std::vector<Measurement> const & extraMeas)
{
    size_t naccept = 0;  // running number of accepted spin flips

    for (size_t sweep = 0; sweep < nsweep; ++sweep) {
        for (size_t step = 0; step < size(lat).get(); ++step) {
            Index const site = rng.genIndex();  // flip spin at this site

            double const delta = deltaE(cfg, site, params, lat);  // proposed change in energy

            // Metropolis-Hastings accept-reject
            // The first check is not necessary for this to be correct but avoids
            // evaluating the costly exponential and RNG.
            if (delta <= 0 or exp(-delta) > rng.genReal()) {
                // accept change
                cfg.flip(site);
                energy += delta;
                ++naccept;
            }
            // else: discard
        }

        measure(obs, lat, cfg, energy);

        // perform extra measurements
        for (auto const &meas : extraMeas) {
            meas(cfg, energy);
        }
    }

    return std::make_tuple(std::move(cfg), energy,
                           static_cast<double>(naccept)
                           / static_cast<double>(nsweep)
                           / static_cast<double>(size(lat).get()));
}
