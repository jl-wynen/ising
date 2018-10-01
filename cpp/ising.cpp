#include <array>
#include <random>
#include <fstream>
#include <numeric>
#include <cmath>
#include <tuple>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <filesystem>

namespace fs = std::filesystem;

using index = int;

constexpr size_t NTHERM_INIT = 1000;
constexpr size_t NTHERM = 1000;
constexpr size_t NPROD = 10000;

constexpr index NX = 5;
constexpr index NY = 5;

constexpr int SEED = 538;

constexpr auto listTemperatures() noexcept {
    // return std::array{0.5, 1., 1.5, 2.};
    return std::array{2.};
}

struct Rng
{
    std::ranlux24 rng;
    std::uniform_int_distribution<index> indexDist;
    std::uniform_real_distribution<double> realDist;
    std::uniform_int_distribution<int> spinDist;

    explicit Rng(int const seed) : rng{seed}, indexDist{0, NX*NY-1},
                                   realDist{0, 1}, spinDist{0, 1}
    { }

    index genIndex()
    {
        return indexDist(rng);
    }

    double genReal()
    {
        return realDist(rng);
    }

    int genSpin()
    {
        return spinDist(rng)==0 ? -1 : 1;
    }
};


struct Configuration
{
    std::array<int, NX*NY> cfg;
    std::array<int, 4*NX*NY> neighbours;

    Configuration()
    {
        for (index y = 0; y < NY; ++y) {
            for (index x = 0; x < NX; ++x) {
                neighbours[(y*NX+x)*4 + 0] = x == NX-1 ? (y*NX)        : (y*NX + x+1);
                neighbours[(y*NX+x)*4 + 1] = x == 0    ? (y*NX + NX-1) : (y*NX + x-1);
                neighbours[(y*NX+x)*4 + 2] = y == NY-1 ? x             : ((y+1)*NX + x);
                neighbours[(y*NX+x)*4 + 3] = y == 0    ? (NY-1)*NX+x   : ((y-1)*NX + x);
            }
        }
    }

    int &operator[](index const idx) noexcept {
        return cfg[idx];
    }

    int const &operator[](index const idx) const noexcept {
        return cfg[idx];
    }

    int &operator()(index const x, index const y) noexcept {
        return cfg[y*NX+x];
    }

    int const &operator()(index const x, index const y) const noexcept {
        return cfg[y*NX+x];
    }
};

struct Observables
{
    std::vector<double> energy;
    std::vector<double> magnetisation;
};

template <size_t NTEMP>
void prepareDatadir(fs::path const &datadir, std::array<double, NTEMP> const &temperatures)
{
    if (fs::exists(datadir)) {
        std::cerr << "Data directory " << datadir << " exists, deleting!\n";
        fs::remove_all(datadir);
    }
    fs::create_directories(datadir);

    std::ofstream ofs;
    ofs.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
    ofs.open(datadir/"temperatures.dat");

    for (size_t i = 0; i < temperatures.size(); ++i) {
        ofs << i << ": " << temperatures[i] << '\n';
    }
}


void write(fs::path const &fname, Observables const &obs)
{
    std::ofstream ofs;
    ofs.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
    ofs.open(fname);

    for (auto const &energy : obs.energy)
        ofs << static_cast<double>(energy) << ' ';
    ofs << '\n';

    for (auto const &magn : obs.magnetisation)
        ofs << magn << ' ';
    ofs << '\n';
}


Configuration randomCfg(Rng &rng)
{
    Configuration cfg;
    for (index i = 0; i < NX*NY; ++i) {
        cfg[i] = rng.genSpin();
    }
    return cfg;
}

int hamiltonian(Configuration const &cfg) noexcept
{
    int energy = 0;
    for (index idx = 0; idx < NX*NY; ++idx) {
        energy += cfg[idx] * (cfg[cfg.neighbours[4*idx]]
                              + cfg[cfg.neighbours[4*idx+1]]
                              + cfg[cfg.neighbours[4*idx+2]]
                              + cfg[cfg.neighbours[4*idx+3]]);
    }
    return -energy;
}

double magnetisation(Configuration const &cfg) noexcept
{
    return std::accumulate(std::begin(cfg.cfg), std::end(cfg.cfg), 0) /
        static_cast<double>(NX*NY);
}

int deltaE(Configuration const &cfg, index const idx) noexcept
{
    return 2*cfg[idx] * (cfg[cfg.neighbours[4*idx]]
                         + cfg[cfg.neighbours[4*idx+1]]
                         + cfg[cfg.neighbours[4*idx+2]]
                         + cfg[cfg.neighbours[4*idx+3]]);
}

auto evolve(Configuration cfg, double energy, double const beta,
            Rng &rng, size_t const nsweep, Observables * const obs)
{
    size_t naccept = 0;

    for (size_t sweep = 0; sweep < nsweep; ++sweep) {
        for (size_t step = 0; step < NX*NY; ++step) {
            index const idx = rng.genIndex();

            int const delta = deltaE(cfg, idx);

            if (delta <= 0 or std::exp(-beta*delta) > rng.genReal()) {
                cfg[idx] *= -1;
                energy += delta;
                ++naccept;
            }
        }

        if (obs) {
            obs->energy.emplace_back(energy);
            obs->magnetisation.emplace_back(magnetisation(cfg));
        }
    }

    return std::make_tuple(cfg, energy, naccept);
}

int main(int const argc, char const * const argv[])
{
    Rng rng{SEED};

    constexpr auto temperatures = listTemperatures();

    fs::path datadir = argc == 2 ? argv[1] : "data";
    prepareDatadir(datadir, temperatures);

    auto cfg = randomCfg(rng);
    auto energy = hamiltonian(cfg);
    size_t naccept;

    auto startTime = std::chrono::high_resolution_clock::now();

    std::tie(cfg, energy, naccept) = evolve(cfg, energy, 1./temperatures[0], rng, NTHERM_INIT, nullptr);
    std::cout << "Initial thermalisation acceptance rate: "
              << static_cast<double>(naccept)/NTHERM_INIT/NX/NY << '\n';

    for (size_t itemp = 0; itemp < temperatures.size(); ++itemp) {
        std::cout << "Running T = " << temperatures[itemp] << '\n';
        auto const beta = 1./temperatures[itemp];

        // thermalise
        std::tie(cfg, energy, naccept) = evolve(cfg, energy, beta, rng, NTHERM, nullptr);
        std::cout << "  Thermalisation acceptance rate: "
                  << static_cast<double>(naccept)/NTHERM/NX/NY << '\n';

        // measure
        Observables obs;
        std::tie(cfg, energy, naccept) = evolve(cfg, energy, beta, rng, NPROD, &obs);
        std::cout << "  Production acceptance rate: "
                  << static_cast<double>(naccept)/NPROD/NX/NY << '\n';

        write(datadir/(std::to_string(itemp)+".dat"), obs);
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << "Duration in wall clock time: " << std::setprecision(4)
              << std::chrono::duration<float>(endTime-startTime).count() << "s\n";
}
