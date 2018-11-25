#ifndef ISING_ARITHMETIC_HPP
#define ISING_ARITHMETIC_HPP

#include <type_traits>
#include <utility>

template <typename T, typename Tag>
struct ArithmeticType
{
    using Underlying = T;

    constexpr explicit ArithmeticType()
        noexcept(std::is_nothrow_default_constructible_v<Underlying>)
            : value_{}
    { }

    constexpr explicit ArithmeticType(Underlying const &x)
        noexcept(std::is_nothrow_copy_constructible_v<Underlying>)
        : value_{x}
    { }

    constexpr explicit ArithmeticType(Underlying &&x)
        noexcept(std::is_nothrow_move_constructible_v<Underlying>)
        : value_{std::move(x)}
    { }

    constexpr ArithmeticType(ArithmeticType const &other)
        noexcept(std::is_nothrow_copy_constructible_v<Underlying>) = default;

    constexpr ArithmeticType(ArithmeticType &&other)
        noexcept(std::is_nothrow_move_constructible_v<Underlying>) = default;

    constexpr ArithmeticType &operator=(ArithmeticType const &other)
        noexcept(std::is_nothrow_copy_assignable_v<Underlying>) = default;

    constexpr ArithmeticType &operator=(ArithmeticType &&other)
        noexcept(std::is_nothrow_move_assignable_v<Underlying>) = default;


    constexpr Underlying &get() noexcept
    {
        return value_;
    }

    constexpr Underlying const &get() const noexcept
    {
        return value_;
    }

    template <typename U>
    constexpr explicit operator U const() const noexcept
    {
        return static_cast<U>(value_);
    }

private:
    Underlying value_;
};


template <typename T, typename Tag>
constexpr auto operator+(ArithmeticType<T, Tag> const &a, ArithmeticType<T, Tag> const &b)
{
    return ArithmeticType<T, Tag>{a.get() + b.get()};
}

template <typename T, typename Tag>
constexpr auto operator-(ArithmeticType<T, Tag> const &a, ArithmeticType<T, Tag> const &b)
{
    return ArithmeticType<T, Tag>{a.get() - b.get()};
}

template <typename T, typename Tag>
constexpr auto operator*(ArithmeticType<T, Tag> const &a, ArithmeticType<T, Tag> const &b)
{
    return ArithmeticType<T, Tag>{a.get() * b.get()};
}

template <typename T, typename Tag>
constexpr auto operator/(ArithmeticType<T, Tag> const &a, ArithmeticType<T, Tag> const &b)
{
    return ArithmeticType<T, Tag>{a.get() / b.get()};
}

template <typename T, typename Tag>
constexpr auto operator%(ArithmeticType<T, Tag> const &a, ArithmeticType<T, Tag> const &b)
{
    return ArithmeticType<T, Tag>{a.get() % b.get()};
}

template <typename T, typename Tag>
constexpr auto &operator++(ArithmeticType<T, Tag> &x)
{
    ++x.get();
    return x;
}

template <typename T, typename Tag>
constexpr auto &operator++(ArithmeticType<T, Tag> &x, int)
{
    decltype(x) c{x.get()++};
    return c;
}

template <typename T, typename Tag>
constexpr auto &operator--(ArithmeticType<T, Tag> &x)
{
    --x.get();
    return x;
}

template <typename T, typename Tag>
constexpr auto &operator--(ArithmeticType<T, Tag> &x, int)
{
    decltype(x) c{x.get()--};
    return c;
}


template <typename T, typename Tag>
constexpr bool operator==(ArithmeticType<T, Tag> const &a, ArithmeticType<T, Tag> const &b)
{
    return a.get() == b.get();
}

template <typename T, typename Tag>
constexpr bool operator!=(ArithmeticType<T, Tag> const &a, ArithmeticType<T, Tag> const &b)
{
    return not (a == b);
}

template <typename T, typename Tag>
constexpr bool operator>(ArithmeticType<T, Tag> const &a, ArithmeticType<T, Tag> const &b)
{
    return a.get() > b.get();
}

template <typename T, typename Tag>
constexpr bool operator>=(ArithmeticType<T, Tag> const &a, ArithmeticType<T, Tag> const &b)
{
    return a.get() >= b.get();
}

template <typename T, typename Tag>
constexpr bool operator<(ArithmeticType<T, Tag> const &a, ArithmeticType<T, Tag> const &b)
{
    return a.get() < b.get();
}

template <typename T, typename Tag>
constexpr bool operator<=(ArithmeticType<T, Tag> const &a, ArithmeticType<T, Tag> const &b)
{
    return a.get() <= b.get();
}

#endif  // ndef ISING_ARITHMETIC_HPP
