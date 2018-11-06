/**
 * Main program of N dimensional Ising model simulation-
 */

#include <chrono>
#include <iomanip>
#include <iostream>

#include <yaml-cpp/yaml.h>

#include "lattice.hpp"
#include "configuration.hpp"
#include "rng.hpp"
#include "ising.hpp"
#include "montecarlo.hpp"
#include "fileio.hpp"

using Clock = std::chrono::steady_clock;
using Milliseconds = std::chrono::milliseconds;


/// Parse command line arguments.
auto parseArgs(int const argc, char const * const argv[])
{
    if (argc != 3) {
        throw std::runtime_error("Need two parameters, in order: input file, output directory!");
    }

    return std::make_tuple(fs::path{argv[1]}, fs::path{argv[2]});
}


int main(int const argc, char const * const argv[])
{
    // load / prepare files
    auto const [infile, outdir] = parseArgs(argc, argv);
    auto const input = YAML::LoadFile(infile).as<ProgConfig>();
    prepareOutdir(outdir);

    Lattice const lat{input.latticeShape};
    Rng rng{size(lat), input.rngSeed};

    // initial state (hot)
    Configuration cfg = randomCfg(size(lat), rng);
    double energy = hamiltonian(cfg, input.params.at(0), lat);
    double accRate;

    // initial thermalisation
    auto startTime = Clock::now();
    std::tie(cfg, energy, accRate) = evolve(cfg, energy, input.params.at(0), lat,
                                            rng, input.nthermInit, nullptr);
    auto endTime = Clock::now();
    std::cout << "Initial thermalisation acceptance rate: " << std::setprecision(4)
              << accRate << '\n'
              << "Run time: " << std::chrono::duration_cast<Milliseconds>(endTime-startTime).count()
              << "ms\n";

    for (size_t i = 0; i < std::size(input.params); ++i) {
        auto const params = input.params.at(i);
        auto const ntherm = input.ntherm.at(i);
        auto const nprod = input.nprod.at(i);

        std::vector<Measurement> meas;
        if (input.writeCfg) {
            meas.emplace_back([&](Configuration const &cfg, double const)
                              {
                                  write(outdir, i, cfg, params, lat);
                              });
        }

        std::cout << "Running with {J/kT = " << params.JT
                  << ", h/kT = " << params.hT << "}\n";

        // (re-)thermalise
        startTime = Clock::now();
        std::tie(cfg, energy, accRate) = evolve(cfg, energy, params,
                                                lat, rng, ntherm, nullptr);
        std::cout << "  Thermalisation acceptance rate: " << std::setprecision(4)
                  << accRate << '\n';

        // measure
        Observables obs;
        std::tie(cfg, energy, accRate) = evolve(cfg, energy, params,
                                                lat, rng, nprod, &obs, meas);
        endTime = Clock::now();
        std::cout << "  Production acceptance rate: " << std::setprecision(4)
                  << accRate << '\n'
                  << "  Run time: " << std::chrono::duration_cast<Milliseconds>(endTime-startTime).count()
                  << "ms\n";

        write(outdir, i, obs, params, lat);
    }
}
