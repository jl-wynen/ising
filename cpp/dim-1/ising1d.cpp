/**
 * C++ implementation of the 1D Ising Model simulation.
 */

#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>


/// Type for indices into configurations.
using Index = size_t;

//--------------------------
// Set run parameters here.

constexpr size_t NTHERM_INIT = 1000;  // number of thermalisation sweeps in the beginning
constexpr size_t NTHERM = 1000;  // number of thermalisation sweeps per temperature
constexpr size_t NPROD = 10000;  // number of production sweeps (with measurements) per temperature

constexpr Index N = 16;  // number of lattice sites

// seed for the random number generator
constexpr int SEED = 538;

// directory for output data
std::string const datadir = "data";

// Return an array of temperatures to run the simulation with.
auto listTemperatures() noexcept {
    // linearly interpolate between 0.4 and 6.0
    constexpr size_t N = 10;
    std::vector<double> temps(N);
    for (size_t i = 0; i < N; ++i) {
        temps[i] = 6. - (6.-0.4)/N*i;
    }
    return temps;

    // just one element
    // return std::array{2.};
}

// End of run parameters.
//------------------------

/// Helper class to handle a random number generator.
struct Rng
{
    std::ranlux24 rng;  ///< The generator.
    /// Distribution to generator lattice indices.
    std::uniform_int_distribution<Index> indexDist;
    /// Distribution to generate floating point numbers in [0, 1).
    std::uniform_real_distribution<double> realDist;
    /// Distribution to generate spins, i.e. values 0 or 1.
    std::uniform_int_distribution<int> spinDist;

    /// Seed the rng and set up distributions.
    explicit Rng(int const seed) : rng{seed}, indexDist{0, N-1},
                                   realDist{0, 1}, spinDist{0, 1}
    { }

    /// Generate random indices into a configuration.
    Index genIndex()
    {
        return indexDist(rng);
    }

    /// Generate a random double in [0, 1).
    double genReal()
    {
        return realDist(rng);
    }

    /// Generate a random spin, one of {-1, +1}.
    int genSpin()
    {
        return spinDist(rng)==0 ? -1 : 1;
    }
};


/// Apply anti-periodic boundary conditions to idx.
constexpr Index applyPeriodicBC(Index const idx, Index const n) noexcept
{
    if (idx == n)
        return 0;
    else if (idx == static_cast<Index>(-1))
        return n-1;
    return idx;
}

/// Hold a spin configuration on the lattice.
struct Configuration
{
    /// Number of sites.
    Index const n;

    /// The actual configuration.
    std::vector<int> cfg;

    /// Initialise with given lattice size.
    /**
     * Sets up the configuration to spin +1 at every site.
     */
    explicit Configuration(Index const n)
        : n(n), cfg(n, 1) { }

    /// Access site.
    int &operator()(Index const x) {
        return cfg.at(applyPeriodicBC(x, n));
    }

    /// Access site.
    int const &operator()(Index const x) const {
        return cfg.at(applyPeriodicBC(x, n));
    }
};


/// Store Monte-Carlo history of observables.
struct Observables
{
    std::vector<double> energy;
    std::vector<double> magnetisation;
};

void writeTemperatures(std::string const &fname, std::vector<double> const &temperatures)
{
    // write temperatures and indices
    std::ofstream ofs;
    ofs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    ofs.open(fname);
    for (size_t i = 0; i < temperatures.size(); ++i) {
        ofs << i << ": " << temperatures[i] << '\n';
    }
}

/// Write observables to a data file.
void writeObservables(std::string const &fname, Observables const &obs)
{
    std::ofstream ofs;
    ofs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    ofs.open(fname);

    for (auto const &energy : obs.energy)
        ofs << static_cast<double>(energy) << ' ';
    ofs << '\n';

    for (auto const &magn : obs.magnetisation)
        ofs << magn << ' ';
    ofs << '\n';
}

/// Generate a random spin configuration.
Configuration randomCfg(Rng &rng, Index const n)
{
    Configuration cfg(n);
    for (Index idx = 0; idx < cfg.n; ++idx) {
        cfg(idx) = rng.genSpin();
    }
    return cfg;
}

/// Evaluate the Hamiltonian on a configuration.
int hamiltonian(Configuration const &cfg)
{
    int energy = 0;
    for (Index idx = 0; idx < cfg.n; ++idx) {
        energy -= cfg(idx) * (cfg(idx+1) + cfg(idx-1));
    }
    return energy;
}

/// Compute the magnetisation on a configuration.
double magnetisation(Configuration const &cfg)
{
    int magn = 0;
    for (Index idx = 0; idx < cfg.n; ++idx) {
        magn += cfg(idx);
    }
    return magn / static_cast<double>(cfg.n);
}

/// Compute the change in energy if the spin at site idx were flipped.
int deltaE(Configuration const &cfg, Index const idx)
{
    return 2*cfg(idx) * (cfg(idx+1) + cfg(idx-1));
}

/// Perform N local MCMC updates on a configuration.
/**
 * Flips spins at random sites N times and accepting or
 * rejecting the change using the Metropolis-Hastings algorithm.
 *
 * \param cfg In: starting configuration;
 *            Out: final configuration.
 * \param energy In: starting energy;
 *               Out: final energy.
 * \param beta Inverse temperature to simulate with.
 * \param rng Random number generator to use for picking sites and
 *            accept-reject.
 *
 * \returns Number of accepted changes.
 */
size_t localUpdateSweep(Configuration &cfg, double &energy, double const beta, Rng &rng)
{
    size_t naccept = 0;  // running number of accepted spin flips

    for (size_t step = 0; step < cfg.n; ++step) {
        // flip spin at this site
        Index const idx = rng.genIndex();

        // proposed change in energy
        int const delta = deltaE(cfg, idx);

        // Metropolis-Hastings accept-reject
        if (std::exp(-beta*delta) > rng.genReal()) {
            // accept change
            cfg(idx) *= -1;
            energy += delta;
            ++naccept;
        }
        // else: discard
    }

    return naccept;
}

/// Evolve a configuration in Monte-Carlo time and measure observables.
/**
 * Performs a given number of sweeps of local updates and measures
 * given observables after every sweep.
 *
 * \param cfg In: starting configuration;
 *            Out: final configuration.
 * \param energy In: starting energy;
 *               Out: final energy.
 * \param beta Inverse temperature to simulate with.
 * \param rng Random number generator to use for picking sites and
 *            accept-reject.
 * \param nsweep Number of sweeps of the lattice to perform.
 * \param obs Instance of Observables to store results of measurements.
 *            May be nullptr in which case no measurements are performed.
 *
 * \returns Number of accepted changes.
 */
size_t integrate(Configuration &cfg, double &energy, double const beta,
                 Rng &rng, size_t const nsweep, Observables * const obs)
{
    size_t naccept = 0;  // running number of accepted spin flips

    for (size_t sweep = 0; sweep < nsweep; ++sweep) {
        naccept += localUpdateSweep(cfg, energy, beta, rng);

        // measure observables if an instance of Observables is given
        if (obs) {
            obs->energy.emplace_back(energy);
            obs->magnetisation.emplace_back(magnetisation(cfg));
        }
    }

    return naccept;
}

int main()
{
    Rng rng{SEED};  // one rng to rule them all

    auto const temperatures = listTemperatures();
    writeTemperatures(datadir+"/temperatures.dat", temperatures);

    // initial condition (hot start)
    auto cfg = randomCfg(rng, N);
    double energy = hamiltonian(cfg);

    // start measuring time now
    auto startTime = std::chrono::high_resolution_clock::now();

    // initial thermalisation
    auto naccept = integrate(cfg, energy, 1./temperatures[0], rng, NTHERM_INIT, nullptr);
    std::cout << "Initial thermalisation acceptance rate: "
              << static_cast<double>(naccept)/NTHERM_INIT/N << '\n';

    // iterate over temperatures
    for (size_t itemp = 0; itemp < temperatures.size(); ++itemp) {
        std::cout << "Running T = " << temperatures[itemp] << '\n';
        auto const beta = 1./temperatures[itemp];

        // thermalise
        naccept = integrate(cfg, energy, beta, rng, NTHERM, nullptr);
        std::cout << "  Thermalisation acceptance rate: "
                  << static_cast<double>(naccept)/NTHERM/N << '\n';

        // measure
        Observables obs;
        naccept = integrate(cfg, energy, beta, rng, NPROD, &obs);
        std::cout << "  Production acceptance rate: "
                  << static_cast<double>(naccept)/NPROD/N << '\n';

        // save result
        writeObservables(datadir+"/"+std::to_string(itemp)+".dat", obs);
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << "Duration in wall clock time: " << std::setprecision(4)
              << std::chrono::duration<float>(endTime-startTime).count() << "s\n";
}
