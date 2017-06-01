#include <GeNuSys/element_traits.h>

#include <GeNuSys/vector.h>
#include <GeNuSys/sparse_vector.h>
#include <GeNuSys/matrix.h>
#include <GeNuSys/sparse_matrix.h>

#include <GeNuSys/p_norm.h>
#include <GeNuSys/frobenius_norm.h>
#include <GeNuSys/operator_norm.h>

#include <GeNuSys/linalg_algorithms.h>

#include <GeNuSys/polynomial.h>
#include <GeNuSys/lehmer_schur.h>

#include <GeNuSys/number_system.h>
#include <GeNuSys/radix_properties.h>
#include <GeNuSys/digit_set.h>
#include <GeNuSys/smith_hash.h>
#include <GeNuSys/numsys_traits.h>
#include <GeNuSys/simultaneous.h>

int main()
{
    //A randomly chosen base
    GeNuSys::LinAlg::Matrix<long long> matrix(3 , 3, std::vector<long long>{1, 2, 4, 5, 10, 6, 7, 1, -5});
    //GeNuSys::LinAlg::Matrix<long long> matrix(3 , 3, std::vector<long long>{1, 12, 12, 5, 10, 6, 7, 1, -5});
    //GeNuSys::LinAlg::Matrix<long long> matrix(3 , 3, std::vector<long long>{1, 3, 12, 2, 10, 6, 7, 1, -5});
    std::cout << "Base:" << std::endl;
    std::cout << matrix << std::endl;

    //Compute some properties
    GeNuSys::NumSys::RadixProperties<long long> props(matrix);
    std::cout << "Abs determinant:" << std::endl;
    std::cout << props.getAbsDet() << std::endl;
    std::cout << "Adjoint:" << std::endl;
    std::cout << props.getAdjoint() << std::endl;

    //Compute some digit sets
    std::vector<GeNuSys::LinAlg::SparseVector<long long>> symmetric = GeNuSys::NumSys::DigitSet::getJSymmetric(props, 0);
    std::cout << "Symmetric digit set:" << std::endl;
    for(const auto& d: symmetric)
    {
        std::cout << d << std::endl;
    }
    std::vector<GeNuSys::LinAlg::SparseVector<long long>> canonical = GeNuSys::NumSys::DigitSet::getJCanonical(props, 0);
    std::cout << "Canonical digit set:" << std::endl;
    for(const auto& d: canonical)
    {
        std::cout << d << std::endl;
    }
    std::vector<GeNuSys::LinAlg::Vector<long long>> adjoint = GeNuSys::NumSys::DigitSet::getAdjoint(props);
    std::cout << "Adjoint digit set:" << std::endl;
    for(const auto& d: adjoint)
    {
        std::cout << d << std::endl;
    }
    std::vector<GeNuSys::LinAlg::Vector<long long>> dense = GeNuSys::NumSys::DigitSet::getDense(props);
    std::cout << "Dense digit set:" << std::endl;
    for(const auto& d: dense)
    {
        std::cout << d << std::endl;
    }

    std::cout << "Volume:" << std::endl;
    std::cout << GeNuSys::NumSys::Traits::getVolume(props.getInverse(), symmetric) << std::endl;

    //Construct a number system using the symmetric digit set and the operator norm
    GeNuSys::NumSys::NumberSystem<long long, GeNuSys::LinAlg::SparseVector, GeNuSys::LinAlg::Matrix, GeNuSys::LinAlg::OperatorNorm<typename GeNuSys::ElementTraits<long long>::RationalType>> numSys1(props, symmetric, props.getOperatorNorm());

    //Compute cycles
    auto cycles1 = numSys1.getCycles();

    std::cout << "Cycles:" << std::endl;
    int cc = 0;
    for(const auto& c: cycles1)
    {
        std::cout << "Cycle " << cc << ":" << std::endl;
        ++cc;
        for(const auto& n: c)
        {
            std::cout << n << std::endl;
        }
    }

    //Optimize the base matrix to reduce volume
    GeNuSys::LinAlg::Matrix<long long> T = GeNuSys::NumSys::Traits::findBasisTransformation(props.getInverse(), symmetric, 15, 5, 2);

    //Compute the improved matrix and digit set (apply the previously computed T transformation matrix)
    auto imprM = T * matrix * GeNuSys::LinAlg::Traits::template convertUnsafe<typename GeNuSys::ElementTraits<long long int>::RationalType, long long>(GeNuSys::LinAlg::Algorithms::invert(T));
    std::vector<GeNuSys::LinAlg::SparseVector<long long>> imprDigits;
    for (unsigned int i = 0; i < symmetric.size(); ++i)
    {
        imprDigits.push_back(T * symmetric[i]);
    }
    GeNuSys::NumSys::RadixProperties<long long> imprProps(imprM);

    std::cout << "Reduced volume:" << GeNuSys::NumSys::Traits::getVolume(imprProps.getInverse(), imprDigits) << std::endl;

    //Construct a number system using the improved matrix and digit set
    GeNuSys::NumSys::NumberSystem<long long, GeNuSys::LinAlg::SparseVector, GeNuSys::LinAlg::Matrix, GeNuSys::LinAlg::OperatorNorm<typename GeNuSys::ElementTraits<long long>::RationalType>> numSys2(imprProps, imprDigits, imprProps.getOperatorNorm());

    //Compute cycles
    auto cycles2 = numSys2.getCycles();

    std::cout << "Cycles (improved matrix):" << std::endl;
    cc = 0;
    for(const auto& c: cycles2)
    {
        std::cout << "Cycle " << cc << ":" << std::endl;
        ++cc;
        for(const auto& n: c)
        {
            std::cout << n << std::endl;
        }
    }

    return 0;
}
