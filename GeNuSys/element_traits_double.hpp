#include <cmath>
#include <complex>
#include <stdexcept>

namespace GeNuSys
{

    template<>
    struct ElementTypeTraits<double>
    {

        typedef double RealType;

        typedef double RationalType;

        typedef std::complex<double> ComplexType;

        typedef double AbsType;

        typedef double AbsSqrType;

    };

    template<>
    inline
    const double& ElementTraits<double>::zero()
    {
        static const double zero = 0.0;
        return zero;
    }

    template<>
    inline
    const double& ElementTraits<double>::one()
    {
        static const double one = 1.0;
        return one;
    }

    template<>
    inline
    const double& ElementTraits<double>::epsilon()
    {
        static const double epsilon = 0.000001;
        return epsilon;
    }

    template<>
    template<>
    inline
    double ElementTraits<double>::asType<double>(const double& value)
    {
        return value;
    }

    template<>
    template<>
    inline
    std::complex<double> ElementTraits<double>::asType<std::complex<double>>(const double& value)
    {
        return std::complex<double>(value, 0);
    }

    template<>
    template<>
    inline
    int ElementTraits<double>::asTypeUnsafe<int>(const double& value)
    {
        return (int)floor(value + 0.5);
    }

    template<>
    template<>
    inline
    long long ElementTraits<double>::asTypeUnsafe<long long>(const double& value)
    {
        return (long long)floor(value + 0.5);
    }

    template<>
    inline
    double ElementTraits<double>::abs(const double& value)
    {
        return std::abs(value);
    }

    template<>
    inline
    double ElementTraits<double>::absSqr(const double& value)
    {
        return value * value;
    }

    template<>
    inline
    double ElementTraits<double>::sqrt(const double& value)
    {
        return std::sqrt(value);
    }

    template<>
    inline
    double ElementTraits<double>::root(const double& value, int n)
    {
        return std::pow(value, 1.0 / n);
    }

    template<>
    inline
    double ElementTraits<double>::div(const double& a, const double& b)
    {
        return a / b;
    }

    //TODO:

    template<>
    inline
    bool ElementTraits<double>::divisible(const double&, const double&)
    {
        ASSERT_EXCEPTION(false, std::logic_error);
    }

    template<>
    inline
    double ElementTraits<double>::idiv(const double&, const double&)
    {
        ASSERT_EXCEPTION(false, std::logic_error);
    }

    template<>
    inline
    double ElementTraits<double>::mod(const double&, const double&)
    {
        ASSERT_EXCEPTION(false, std::logic_error);
    }

    template<>
    inline
    double ElementTraits<double>::mods(const double&, const double&)
    {
        ASSERT_EXCEPTION(false, std::logic_error);
    }
}
