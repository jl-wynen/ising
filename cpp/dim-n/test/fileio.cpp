#include "fileio.hpp"

#include <random>
#include <string>

#include "catch.hpp"
#include "util.hpp"

namespace {
    // TODO use YAML to output
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

    std::string writeParams(Parameters const &params)
    {
        std::ostringstream oss;
        oss << "J: " << params.JT << '\n';
        oss << "h: " << params.hT << '\n';
        return oss.str();
    }

    /// Check if two Parameter instances are exactly equal.
    // Not very useful outside of these tests, so don't provide it globally.
    bool operator==(Parameters const &a, Parameters const &b) noexcept
    {
        return a.JT == b.JT and a.hT == b.hT;
    }

    /// Generate random ints and return them as doubles.
    // This should avoid problems when comparing floats with operator==.
    struct NumberRNG
    {
        std::mt19937 &rng;
        std::uniform_int_distribution<int> &dist;

        double operator()()
        {
            return static_cast<double>(dist(rng));
        }
    };
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
        for (size_t sample = 0; sample < nsamples; ++sample) {
            int const number = numberDist(rng);
            auto const node = YAML::Load(std::to_string(number));
            auto const vec = loadVector<int>(node);

            REQUIRE(std::size(vec) == 1);
            REQUIRE(vec[0] == number);
        }
    }

    SECTION("Vectors can be loaded from sequence nodes")
    {
        for (size_t sample = 0; sample < nsamples; ++sample) {
            auto vecOut = produce(NumberRNG{rng, numberDist}, lengthDist(rng));
            auto const node = YAML::Load(writeVector(vecOut));
            auto const vecIn = loadVector<double>(node);

            REQUIRE(vecIn == vecOut);
        }
    }
}

TEST_CASE("Loading Parameters from YAML", "[YAML]")
{
    std::mt19937 rng{1532};
    // generate int to avoid floating point issues
    // but use floats because Parameters needs those
    std::uniform_int_distribution<int> numberDist{-50, 50};
    std::uniform_int_distribution<size_t> lengthDist{1, 30};

    constexpr size_t nsamples = 10;

    SECTION("For only scalar parameters, gives just one object")
    {
        for (size_t sample = 0; sample < nsamples; ++sample) {
            Parameters const paramsOut{static_cast<double>(numberDist(rng)),
                                       static_cast<double>(numberDist(rng))};

            auto const node = YAML::Load(writeParams(paramsOut));
            auto const paramsInVec = loadParams(node);

            REQUIRE(std::size(paramsInVec) == 1);
            REQUIRE(paramsInVec[0] == paramsOut);
        }
    }

    SECTION("If all parameters are sequences, they are combined")
    {
        for (size_t sample = 0; sample < nsamples; ++sample) {
            auto const length = lengthDist(rng);
            auto vecJOut = produce(NumberRNG{rng, numberDist}, length);
            auto vechOut = produce(NumberRNG{rng, numberDist}, length);

            auto const node = YAML::Load("J: "+writeVector(vecJOut)+"\nh: "+writeVector(vechOut));
            auto const paramsInVec = loadParams(node);

            REQUIRE(std::size(paramsInVec) == length);
            for (size_t i = 0; i < length; ++i) {
                REQUIRE(paramsInVec[i].JT == vecJOut[i]);
                REQUIRE(paramsInVec[i].hT == vechOut[i]);
            }
        }
    }

    SECTION("If there are multiple sequences, they must have equal lengths")
    {
        for (size_t sample = 0; sample < nsamples; ++sample) {
            auto const lengthJ = lengthDist(rng)+1;  // +1 because we don't want saclar nodes here
            auto vecJOut = produce(NumberRNG{rng, numberDist}, lengthJ);
            auto lengthh = lengthDist(rng)+1;  // +1 because we don't want saclar nodes here
            if (lengthh == lengthJ) lengthh++;  // make sure the lengths are different
            auto vechOut = produce(NumberRNG{rng, numberDist}, lengthh);

            auto const node = YAML::Load("J: "+writeVector(vecJOut)+"\nh: "+writeVector(vechOut));
            REQUIRE_THROWS(loadParams(node));
        }
    }

    SECTION("If one parameter is a sequence, the others are broadcast")
    {
        SECTION("J is a sequence")
        {
            for (size_t sample = 0; sample < nsamples; ++sample) {
                auto const length = lengthDist(rng);
                auto const vecJOut = produce(NumberRNG{rng, numberDist}, length);
                double const hOut = numberDist(rng);

                auto const node = YAML::Load("J: "+writeVector(vecJOut)+"\nh: "+std::to_string(hOut));
                auto const paramsInVec = loadParams(node);

                REQUIRE(std::size(paramsInVec) == length);
                for (size_t i = 0; i < length; ++i) {
                    REQUIRE(paramsInVec[i].JT == vecJOut[i]);
                    REQUIRE(paramsInVec[i].hT == hOut);
                }
            }
        }

        SECTION("h is a sequence")
        {
            for (size_t sample = 0; sample < nsamples; ++sample) {
                auto const length = lengthDist(rng);
                double const JOut = numberDist(rng);
                auto const vechOut = produce(NumberRNG{rng, numberDist}, length);

                auto const node = YAML::Load("J: "+std::to_string(JOut)+"\nh: "+writeVector(vechOut));
                auto const paramsInVec = loadParams(node);

                REQUIRE(std::size(paramsInVec) == length);
                for (size_t i = 0; i < length; ++i) {
                    REQUIRE(paramsInVec[i].JT == JOut);
                    REQUIRE(paramsInVec[i].hT == vechOut[i]);
                }
            }
        }
    }
}
