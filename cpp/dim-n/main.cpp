/**
 * C++ implementation of the Ising Model simulation.
 */

// #include <algorithm>
// #include <array>
// #include <chrono>
// #include <cmath>
// #include <filesystem>
// #include <fstream>
// #include <iomanip>
#include <iostream>
// #include <numeric>
// #include <random>
// #include <tuple>
// #include <vector>
// #include <exception>
#include <cassert>

// namespace fs = std::filesystem;

#include "lattice.hpp"
#include "configuration.hpp"
#include "rng.hpp"
#include "ising.hpp"

std::ostream &operator<<(std::ostream &os, Spin const s)
{
    assert(s == Spin{1} or s == Spin{-1});
    os << (s==Spin{+1} ? '+' : '-');
    return os;
}

std::ostream &operator<<(std::ostream &os, Configuration const &cfg)
{
    os << '(';
    for (Spin s : cfg) {
        os << s << ' ';
    }
    os << ')';

    return os;
}

std::ostream &operator<<(std::ostream &os, std::vector<Index> const &vec)
{
    os << '[';
    for (Index i: vec) {
        os << i.get() << ' ';
    }
    os << ']';

    return os;
}


int main()
{
    Lattice lat{3_i, 4_i};


    // Rng rng{541};
    // Configuration cfg = randomCfg(rng);
    // Parameters params{1.1, 3.5};

    // Index flipSite{4};

    // std::cout << cfg << '\n';

    // Configuration cfg2{cfg};
    // cfg2.flip(flipSite);
    // std::cout << cfg2 << '\n';

    // // std::cout << magnetisation(cfg) << '\n';
    // // std::cout << hamiltonian(cfg, params) << '\n';
    // std::cout << hamiltonian(cfg, params) << ", " << hamiltonian(cfg2, params) << '\n';



    // std::cout << hamiltonian(cfg2, params)-hamiltonian(cfg, params) << '\n';

    // std::cout << deltaE(cfg, flipSite, params) << '\n';
}


//--------------------------
// Set run parameters here.

// constexpr size_t NTHERM_INIT = 1000;  // number of thermalisation sweeps in the beginning
// constexpr size_t NTHERM = 1000;  // number of thermalisation sweeps per temperature
// constexpr size_t NPROD = 10000;  // number of production sweeps (with measurements) per temperature

// // seed for the random number generator
// constexpr unsigned long SEED = 538;

// // Return an array of temperatures to run the simulation with.
// constexpr auto listTemperatures() noexcept {
//     // linearly interpolate between 0.4 and 6.0
//     // constexpr size_t N = 10;
//     // std::array<double, N> temps{};
//     // for (size_t i = 0; i < N; ++i) {
//     //     temps[i] = 6. - (6.-0.4)/N*i;
//     // }
//     // return temps;

//     // just one element
//     return std::array{2.};
// }

// End of run parameters.
//------------------------



// /// Store Monte-Carlo history of observables.
// struct Observables
// {
//     std::vector<double> energy;
//     std::vector<double> magnetisation;
// };

// /// Create the output data directory and write the temperature file.
// /**
//  * Deletes the directory and all its contents if it exists.
//  */
// template <size_t NTEMP>
// void prepareDatadir(fs::path const &datadir, std::array<double, NTEMP> const &temperatures)
// {
//     // delete if exists and (re-)create
//     if (fs::exists(datadir)) {
//         std::cerr << "Data directory " << datadir << " exists, deleting!\n";
//         fs::remove_all(datadir);
//     }
//     fs::create_directory(datadir);

//     // write temperatures and indices
//     std::ofstream ofs;
//     ofs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
//     ofs.open(datadir/"temperatures.dat");
//     for (size_t i = 0; i < temperatures.size(); ++i) {
//         ofs << i << ": " << temperatures[i] << '\n';
//     }
// }

// /// Write observables to a data file.
// void write(fs::path const &fname, Observables const &obs)
// {
//     std::ofstream ofs;
//     ofs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
//     ofs.open(fname);

//     for (auto const &energy : obs.energy)
//         ofs << static_cast<double>(energy) << ' ';
//     ofs << '\n';

//     for (auto const &magn : obs.magnetisation)
//         ofs << magn << ' ';
//     ofs << '\n';
// }

// /// Exponential function for acceptance probability.
// struct Exp
// {
//     /// Precompute possible values.
//     /**
//      * \param beta Coefficient in fron to energy difference (J/(k_B T)).
//      */
//     explicit Exp(double const beta) noexcept
//         : exp4{std::exp(-beta*4)}, exp8{std::exp(-beta*8)} { }

//     /// Evaluate function.
//     /**
//      * \param delta Energy difference.
//      * \returns exp(-beta*delta).
//      * \attention Only allows values delta = 4,8. Checks in debug build.
//      */
//     double operator()(int const delta) const noexcept(ndebug)
//     {
//         if constexpr (not ndebug) {
//             if (delta != 4 and delta != 8)
//                 throw std::invalid_argument("Unsupported parameter in Exp::operator().");
//         }
//         return (delta==4) ? exp4 : exp8;
//     }

// private:
//     double const exp4, exp8;
// };

// /// Evolve a configuration in Monte-Carlo time.
// /**
//  * Flips spins at random sites nsweep*NX*NY times and accepting or
//  * rejecting the change using the Metropolis-Hastings algorithm.
//  * Measures observables every NX*NY steps, i.e. once per sweep.
//  *
//  * \param cfg Starting configuration.
//  * \param energy Starting energy.
//  * \param beta Inverse temperature to simulate with.
//  * \param rng Random number generator to use for picking sites and
//  *            accept-reject.
//  * \param nsweep Number of sweeps of the lattice to perform.
//  * \param obs Instance of Observables to store results of measurements.
//  *            May be nullptr in which case no measurements are performed.
//  *
//  * \returns Tuple of
//  *          - Final configuration
//  *          - Final energy
//  *          - Number of accepted changes.
//  */
// auto evolve(Configuration cfg, double energy, double const beta,
//             Rng &rng, size_t const nsweep, Observables * const obs)
// {
//     size_t naccept = 0;  // running number of accepted spin flips
//     Exp exp(beta);  // fast way to compute exponentials

//     for (size_t sweep = 0; sweep < nsweep; ++sweep) {
//         for (size_t step = 0; step < NX*NY; ++step) {
//             Index const idx = rng.genIndex();  // flip spin at this site

//             int const delta = deltaE(cfg, idx);  // proposed change in energy

//             // Metropolis-Hastings accept-reject
//             // The first check is not necessary for this to be correct but avoids
//             // evaluating the costly exponential and RNG.
//             if (delta <= 0 or exp(delta) > rng.genReal()) {
//                 // accept change
//                 cfg[idx] *= -1;
//                 energy += delta;
//                 ++naccept;
//             }
//             // else: discard
//         }

//         // measure observables if an instance of Observables is given
//         if (obs) {
//             obs->energy.emplace_back(energy);
//             obs->magnetisation.emplace_back(magnetisation(cfg));
//         }
//     }

//     return std::make_tuple(std::move(cfg), energy, naccept);
// }

// int main(int const argc, char const * const argv[])
// {
//     Rng rng{SEED};  // one rng for all

//     constexpr auto temperatures = listTemperatures();

//     // make output directory and store temperatures
//     fs::path datadir = fs::weakly_canonical(argc == 2 ? argv[1] : "data");
//     prepareDatadir(datadir, temperatures);

//     // initial condition (hot start)
//     auto cfg = randomCfg(rng);
//     auto energy = hamiltonian(cfg);
//     size_t naccept;

//     // start measuring time now
//     auto const startTime = std::chrono::high_resolution_clock::now();

//     // initial thermalisation
//     std::tie(cfg, energy, naccept) = evolve(cfg, energy, 1./temperatures[0], rng, NTHERM_INIT, nullptr);
//     std::cout << "Initial thermalisation acceptance rate: "
//               << static_cast<double>(naccept)/NTHERM_INIT/NX/NY << '\n';

//     // iterate over temperatures
//     for (size_t itemp = 0; itemp < temperatures.size(); ++itemp) {
//         std::cout << "Running T = " << temperatures[itemp] << '\n';
//         auto const beta = 1./temperatures[itemp];

//         // re-thermalise
//         std::tie(cfg, energy, naccept) = evolve(cfg, energy, beta, rng, NTHERM, nullptr);
//         std::cout << "  Thermalisation acceptance rate: "
//                   << static_cast<double>(naccept)/NTHERM/NX/NY << '\n';

//         // measure
//         Observables obs;
//         std::tie(cfg, energy, naccept) = evolve(cfg, energy, beta, rng, NPROD, &obs);
//         std::cout << "  Production acceptance rate: "
//                   << static_cast<double>(naccept)/NPROD/NX/NY << '\n';

//         // save result
//         write(datadir/(std::to_string(itemp)+".dat"), obs);
//     }

//     auto const endTime = std::chrono::high_resolution_clock::now();
//     std::cout << "Duration in wall clock time: " << std::setprecision(4)
//               << std::chrono::duration<float>(endTime-startTime).count() << "s\n";
// }
