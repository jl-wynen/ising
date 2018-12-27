#ifndef ISING_CONFIGURATION_HPP
#define ISING_CONFIGURATION_HPP

#include <vector>

#include "lattice.hpp"
#include "arithmetic.hpp"
#include "ndebug.hpp"

using Spin = ArithmeticType<int, struct SpinTag>;

/// Hold a spin configuration on the lattice.
struct Configuration
{
    /// Initialise with a given initial spin.
    Configuration(Index const size,
                  Spin const &initial=Spin{+1})
        : cfg_(size.get(), initial)
    {
        if constexpr (not ndebug) {
            if (initial != Spin{+1} and initial != Spin{-1})
                throw std::runtime_error("Invlaid initial spin. Must be +1 or -1");
        }
    }

    /// Access site by total index.
    Spin &operator[](Index const idx) noexcept(ndebug)
    {
        if constexpr (not ndebug) {
            if (idx.get() >= std::size(cfg_))  // Index is unsigned => this checks for 'idx < 0' too.
                throw std::out_of_range("Configuration index is out of range.");
        }
        return cfg_[idx.get()];
    }

    /// Access site by total index.
    Spin const &operator[](Index const idx) const noexcept(ndebug)
    {
        if constexpr (not ndebug) {
            if (idx.get() >= std::size(cfg_))  // Index is unsigned => this checks for 'idx < 0' too.
                throw std::out_of_range("Configuration index is out of range.");
        }
        return cfg_[idx.get()];
    }

    /// Return the size of this configuration.
    Index size() const noexcept
    {
        return Index{cfg_.size()};
    }

    /// Return iterator to beginning of config.
    friend auto begin(Configuration &cfg)
    {
        return cfg.cfg_.begin();
    }

    /// Return iterator to beginning of config.
    friend auto begin(Configuration const &cfg)
    {
        return cfg.cfg_.begin();
    }

    /// Return iterator to end of config.
    friend auto end(Configuration &cfg)
    {
        return cfg.cfg_.end();
    }

    /// Return iterator to end of config.
    friend auto end(Configuration const &cfg)
    {
        return cfg.cfg_.end();
    }

    /// Flip the sign at site idx.
    void flip(Index const idx) noexcept(ndebug)
    {
        (*this)[idx] = (*this)[idx] * Spin(-1);
    }

private:
    /// The actual configuration.
    std::vector<Spin> cfg_;
};


/// Return the size of a configuration.
inline Index size(Configuration const &cfg) noexcept
{
    return cfg.size();
}

#endif  // ndef ISING_CONFIGURATION_HPP
