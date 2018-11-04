/**
 * Main program of N dimensional Ising model simulation-
 */

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <cassert>

namespace fs = std::filesystem;

#include <yaml-cpp/yaml.h>

#include "lattice.hpp"
#include "configuration.hpp"
#include "rng.hpp"
#include "ising.hpp"
#include "montecarlo.hpp"
#include "fileio.hpp"

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
    std::tie(cfg, energy, accRate) = evolve(cfg, energy, input.params.at(0), lat,
                                            rng, input.nthermInit, nullptr);
    std::cout << "Initial thermalisation acceptance rate: " << std::setprecision(4)
              << accRate << '\n';

    for (size_t i = 0; i < std::size(input.params); ++i) {
        auto const params = input.params.at(i);
        auto const ntherm = input.ntherm.at(i);
        auto const nprod = input.ntherm.at(i);

        std::cout << "Running with {J/kT = " << params.JT << ", h/kT = " << params.hT << "}\n";

        // (re-)thermalise
        std::tie(cfg, energy, accRate) = evolve(cfg, energy, params, lat, rng, ntherm, nullptr);
        std::cout << "  Thermalisation acceptance rate: " << std::setprecision(4)
                  << accRate << '\n';

        // measure
        Observables obs;
        std::tie(cfg, energy, accRate) = evolve(cfg, energy, params, lat, rng, nprod, &obs);
        std::cout << "  Production acceptance rate: " << std::setprecision(4)
                  << accRate << '\n';

        write(outdir, i, obs, params, lat);
    }
}
