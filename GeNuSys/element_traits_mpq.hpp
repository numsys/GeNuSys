#include <complex>

#include <gmp.h>
#include <gmpxx.h>

namespace GeNuSys
{

    template<>
    struct ElementTypeTraits<mpq_class>
    {

        typedef mpf_class RealType;

        typedef mpq_class RationalType;

        typedef std::complex<mpq_class> ComplexType;

        typedef mpq_class AbsType;

        typedef mpq_class AbsSqrType;

    };

    template<>
    inline
    const mpq_class& ElementTraits<mpq_class>::zero()
    {
        static const mpq_class zero(0);
        return zero;
    }

    template<>
    inline
    const mpq_class& ElementTraits<mpq_class>::one()
    {
        static const mpq_class one(1);
        return one;
    }
    
    template<>
    inline
    const mpq_class& ElementTraits<mpq_class>::epsilon()
    {
        ASSERT_EXCEPTION(!"INVALID_OPERATION", std::runtime_error);        
    }

    template<>
    template<>
    inline
    mpq_class ElementTraits<mpq_class>::asType<mpq_class>(const mpq_class& value)
    {
        return value;
    }

    template<>
    template<>
    inline
    mpf_class ElementTraits<mpq_class>::asType<mpf_class>(const mpq_class& value)
    {
        return mpf_class(value);
    }

    template<>
    template<>
    inline
    std::complex<mpq_class> ElementTraits<mpq_class>::asType<std::complex<mpq_class>>(const mpq_class& value)
    {
        return std::complex<mpq_class>(value, mpq_class(0));
    }

    template<>
    template<>
    inline
    double ElementTraits<mpq_class>::asType<double>(const mpq_class& value)
    {
        return value.get_d();
    }
    
    template<>
    template<>
    inline
    std::complex<mpf_class> ElementTraits<mpq_class>::asType<std::complex<mpf_class>>(const mpq_class& value)
    {
        return std::complex<mpf_class>(value, mpq_class(0));
    }

    template<>
    template<>
    inline
    mpz_class ElementTraits<mpq_class>::asTypeUnsafe<mpz_class>(const mpq_class& value)
    {
        mpz_class result;
        mpz_fdiv_q(result.get_mpz_t(), value.get_num_mpz_t(), value.get_den_mpz_t());

        return result;
    }
    
    template<>
    template<>
    inline
    int ElementTraits<mpq_class>::asTypeUnsafe<int>(const mpq_class& value)
    {
        mpz_class result;
        mpz_fdiv_q(result.get_mpz_t(), value.get_num_mpz_t(), value.get_den_mpz_t());
        if (result.fits_sint_p())
        {
            return result.get_si();
        }
        else
        {
            throw std::out_of_range{"Conversion to int failed: out of range"};
        }
    }

    template<>
    template<>
    inline
    mpq_class ElementTraits<int>::asType<mpq_class>(const int& value)
    {
        return mpq_class(value);
    }

    template<>
    template<>
    inline
    std::complex<mpq_class> ElementTraits<int>::asType<std::complex<mpq_class>>(const int& value)
    {
        return std::complex<mpq_class>(value, mpq_class(0));
    }

    template<>
    inline
    mpq_class ElementTraits<mpq_class>::abs(const mpq_class& value)
    {
        mpq_class result;
        mpq_abs(result.get_mpq_t(), value.get_mpq_t());

        return result;
    }

    template<>
    inline
    mpq_class ElementTraits<mpq_class>::absSqr(const mpq_class& value)
    {
        return value * value;
    }

    template<>
    inline
    mpq_class ElementTraits<mpq_class>::pow(const mpq_class& value, unsigned int n)
    {
        mpq_class result;
        mpz_pow_ui(mpq_numref(result.get_mpq_t()), mpq_numref(value.get_mpq_t()), n);
        mpz_pow_ui(mpq_denref(result.get_mpq_t()), mpq_denref(value.get_mpq_t()), n);

        return result;
    }

    template<>
    inline
    mpf_class ElementTraits<mpq_class>::sqrt(const mpq_class& value)
    {
        mpf_class sqrt_num, sqrt_den;
        mpf_sqrt(sqrt_num.get_mpf_t(), mpf_class(mpz_class(mpq_numref(value.get_mpq_t()))).get_mpf_t());
        mpf_sqrt(sqrt_den.get_mpf_t(), mpf_class(mpz_class(mpq_denref(value.get_mpq_t()))).get_mpf_t());

        return sqrt_num / sqrt_den;
    }

    template<>
    inline
    mpq_class ElementTraits<mpq_class>::div(const mpq_class& a, const mpq_class& b)
    {
        mpq_class result;
        mpq_div(result.get_mpq_t(), a.get_mpq_t(), b.get_mpq_t());

        return result;
    }

}
