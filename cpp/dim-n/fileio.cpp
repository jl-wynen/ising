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
