#include "fileio.hpp"

#include <random>
#include <string>

#include "catch.hpp"

namespace {
    template <typename T>
    std::string writeVector(std::vector<T> const &vec)
    {
        std::ostringstream oss;
        oss << '[';
        for (size_t i = 0; i < std::size(vec)-1; ++i) {
            oss << vec[i] << ", ";
        }
        oss << vec.back() << ']';
        return oss.str();
    }
}


TEST_CASE("Loading vectors from YAML", "[YAML]")
{
    std::mt19937 rng{919};
    // use int to avoid floating point issues
    std::uniform_int_distribution<int> numberDist{-50, 50};
    std::uniform_int_distribution<size_t> lengthDist{1, 30};

    constexpr size_t nsamples = 50;

    SECTION("Vectors can be loaded from scalar nodes")
    {
        for (size_t i = 0; i < nsamples; ++i) {
            int const number = numberDist(rng);
            auto const node = YAML::Load(std::to_string(number));
            auto const vec = loadVector<int>(node);

            REQUIRE(std::size(vec) == 1);
            REQUIRE(vec[0] == number);
        }
    }

    SECTION("Vectors can be loaded from sequence nodes")
    {
        for (size_t i = 0; i < nsamples; ++i) {
            std::vector<int> vecOut(lengthDist(rng));
            std::generate(std::begin(vecOut), std::end(vecOut),
                          [&](){ return numberDist(rng); });

            auto const node = YAML::Load(writeVector(vecOut));
            auto const vecIn = loadVector<int>(node);

            REQUIRE(vecIn == vecOut);
        }
    }
}
