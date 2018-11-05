#ifndef ISING_NDEBUG_HPP
#define ISING_NDEBUG_HPP

/// Indicate whether NDEBUG macro was set.
#ifdef NDEBUG
    constexpr bool ndebug = true;
#else
    constexpr bool ndebug = false;
#endif

#endif  // ndef ISING_NDEBUG_HPP
