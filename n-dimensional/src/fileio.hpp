#ifndef ISLE_FILEIO_HPP
#define ISLE_FILEIO_HPP

#include <vector>
#include <filesystem>

#include <yaml-cpp/yaml.h>

#include "ising.hpp"
#include "montecarlo.hpp"
#include "index.hpp"

namespace fs = std::filesystem;


/// Program configuration as read from input file.
struct ProgConfig
{
    std::vector<Index> latticeShape;
    unsigned long rngSeed;
    std::vector<Parameters> params;
    size_t nthermInit;
    std::vector<size_t> ntherm;
    std::vector<size_t> nprod;
    bool writeCfg;
};


namespace YAML {
    /// Allow reading of Indices from YAML.
    template <>
    struct convert<Index>
    {
        static bool decode(Node const &node, Index &idx);
    };

    /// Allow reading of program configurations from YAML.
    template <>
    struct convert<ProgConfig>
    {
        static bool decode(Node const &node, ProgConfig &pc);
    };
}


/// Load a vector from a YAML sequence or scalar node.
template <typename T>
std::vector<T> loadVector(YAML::Node const &node)
{
    std::vector<T> values;

    if (node.IsSequence()) {
        values = node.as<std::vector<T>>();
    }
    else if (node.IsScalar()) {
        values.emplace_back(node.as<T>());
    }
    else {
        throw std::invalid_argument("Invalid YAML node type to read vector");
    }

    return values;
}

/// Load parameters from YAML.
std::vector<Parameters> loadParams(YAML::Node const &node);


/// Create the output data directory.
/**
 * Deletes the directory and all its contents if it exists.
 */
void prepareOutdir(fs::path const &outdir);

/// Write observables to a data file.
void write(fs::path const &outdir, size_t ensemble,
           Observables const &obs, Parameters const &params,
           Lattice const &lat);

/// Write a configuration to a file.
/**
 * Appends the config if the file already exists.
 */
void write(fs::path const &outdir, size_t ensemble,
           Configuration const &cfg,
           Parameters const &params, Lattice const &lat);

#endif  // nde ISLE_FILEIO_HPP
