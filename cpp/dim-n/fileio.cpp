#include "fileio.hpp"

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

        return true;
    }
}
