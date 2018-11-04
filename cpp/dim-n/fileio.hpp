#ifndef ISLE_FILEIO_HPP
#define ISLE_FILEIO_HPP

#include <vector>

#include <yaml-cpp/yaml.h>

#include "ising.hpp"

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

#endif  // nde ISLE_FILEIO_HPP
