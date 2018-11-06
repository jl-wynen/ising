#include "fileio.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace {
    /// Load parameters as individual arrays.
    auto loadParamsIndividual(YAML::Node const &node)
    {
        auto J = loadVector<double>(node["J"]);
        auto h = loadVector<double>(node["h"]);

        if (std::size(J) > 1 and std::size(h) > 1) {
            if (std::size(J) != std::size(h)) {
                throw std::invalid_argument("Both J and h are given as sequences but have different lengths");
            }
            return std::make_tuple(std::move(J), std::move(h));
        }
        // at mode one param is a sequence, broadcast the other

        if (std::size(J) > 1) {
            h.resize(std::size(J), h.front());
        }
        else if (std::size(h) > 1) {
            J.resize(std::size(h), J.front());
        }

        return std::make_tuple(std::move(J), std::move(h));
    }

    /// Make sure that size(vec) == 1 or desired, broadcast to desired iff 1.
    template <typename T>
    void checkSizeAndBroadcast(std::vector<T> &vec, std::size_t const desired)
    {
        if (std::size(vec) > 1) {
            if (std::size(vec) != desired) {
                throw std::invalid_argument("Inconsistent lengths in sequences in input");
            }
            return;  // sizes match, nothing to be done
        }

        // broadcast
        vec.resize(desired, vec.front());
    }

    /// Return the output file name for given ensemble number.
    fs::path outFname(size_t const ensemble, char const extension[]=".dat")
    {
        std::ostringstream oss;
        oss << std::setfill('0') << std::setw(4) << ensemble << extension;
        return {oss.str()};
    }

    /// Output a vector of arbitrary elements.
    template <typename T>
    std::ostream &operator<<(std::ostream &os, std::vector<T> const &vec)
    {
        for (size_t i = 0; i < std::size(vec)-1; ++i) {
            os << vec[i] << ", ";
        }
        os << vec.back();
        return os;
    }

    /// Output an Index.
    std::ostream &operator<<(std::ostream &os, Index const idx)
    {
        return (os << idx.get());
    }

    /// Write ensemble metadata to given file.
    std::ofstream writeMetadata(std::ofstream &ofs, Parameters const &params,
                                Lattice const &lat)
    {
        ofs << "# J=" << params.JT << " h=" << params.hT
            << " shape=[" << lat.shape() << "]\n";
        return std::move(ofs);
    }

    /// Write ensemble metadata to a new file.
    std::ofstream writeMetadata(fs::path const &fname, Parameters const &params,
                                Lattice const &lat)
    {
        std::ofstream ofs;
        ofs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        ofs.open(fname, std::ios::trunc);
        return writeMetadata(ofs, params, lat);
    }
}

namespace YAML {
    bool convert<Index>::decode(Node const &node, Index &idx)
    {
        if (not node.IsScalar()) {
            return false;
        }
        idx = Index{node.as<typename Index::Underlying>()};
        return true;
    }

    bool convert<ProgConfig>::decode(Node const &node, ProgConfig &pc)
    {
        if (not node.IsMap()) {
            return false;
        }

        pc.latticeShape = node["Lattice"]["shape"].as<std::vector<Index>>();
        pc.rngSeed = node["RNG"]["seed"].as<unsigned long>();
        pc.params = loadParams(node["Parameters"]);

        auto const &mcNode = node["MC"];
        pc.nthermInit = mcNode["ntherm_init"].as<size_t>();
        pc.ntherm = loadVector<size_t>(mcNode["ntherm"]);
        pc.nprod = loadVector<size_t>(mcNode["nprod"]);

        checkSizeAndBroadcast(pc.ntherm, std::size(pc.params));
        checkSizeAndBroadcast(pc.nprod, std::size(pc.params));

        std::string const startStr = mcNode["start"].as<std::string>();
        if (startStr == "hot") {
            pc.start = ProgConfig::HOT;
        }
        else if (startStr == "cold") {
            pc.start = ProgConfig::COLD;
        }
        else {
            throw std::invalid_argument("Invalid argument to input param 'start'");
        }

        pc.writeCfg = node["write_cfg"].as<bool>();

        return true;
    }
}

std::vector<Parameters> loadParams(YAML::Node const &node)
{
    std::vector<Parameters> params;
    auto const [J, h] = loadParamsIndividual(node);
    for (size_t i = 0; i < std::size(J); ++i) {
        params.emplace_back(Parameters{J[i], h[i]});
    }
    return params;
}

void prepareOutdir(fs::path const &outdir)
{
    if (fs::exists(outdir)) {
        std::cerr << "Output directory " << outdir << " exists, deleting!\n";
        fs::remove_all(outdir);
    }
    fs::create_directory(outdir);
}

void write(fs::path const &outdir, size_t const ensemble,
           Observables const &obs, Parameters const &params,
           Lattice const &lat)
{
    auto ofs = writeMetadata(outdir/outFname(ensemble), params, lat);
    ofs << obs.energy << '\n'
        << obs.magnetisation << '\n';
}

void write(fs::path const &outdir, size_t const ensemble,
           Configuration const &cfg,
           Parameters const &params, Lattice const &lat)
{
    fs::path const outfile = outdir/outFname(ensemble, ".cfg");
    std::ofstream ofs;
    if (not fs::exists(outfile)) {
        ofs = writeMetadata(outfile, params, lat);
    }
    else {
        ofs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        ofs.open(outfile, std::ios::app);
    }

    for (Index i = 0_i; i < size(cfg)-1_i; ++i) {
        ofs << cfg[i].get() << ", ";
    }
    ofs << cfg[size(cfg)-1_i].get() << '\n';
}
