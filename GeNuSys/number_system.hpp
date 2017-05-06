#include "numsys_traits.h"

#include "bit_vector.h"
#include "vector_coder.h"

namespace GeNuSys
{
    namespace NumSys
    {

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType,
            typename Norm
            >
        NumberSystem<ElementType, VectorType, MatrixType, Norm>::NumberSystem(const RadixProperties<ElementType>& props,
                                                                              const std::vector<VectorType<ElementType>>& digitSet, const Norm& norm): props(props), digitSet(digitSet), hash(props), hashTable(hash, digitSet), norm(norm)
        {
        }

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType,
            typename Norm
            >
        GeNuSys::LinAlg::Vector<ElementType> NumberSystem<ElementType, VectorType, MatrixType, Norm>::phi(const GeNuSys::LinAlg::Vector<ElementType>& z) const
        {
            GeNuSys::LinAlg::Vector<ElementType> phiZ = props.getAdjoint() * (z - hashTable(z));
            GeNuSys::LinAlg::Operations::vct_idiv(phiZ, props.getDet());

            return phiZ;
        }

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType,
            typename Norm
            >
        const VectorType<ElementType>& NumberSystem<ElementType, VectorType, MatrixType, Norm>::phi(GeNuSys::LinAlg::Vector<ElementType>& z, GeNuSys::LinAlg::Vector<ElementType>& phiZ, GeNuSys::LinAlg::Vector<ElementType>& Uz) const
        {
            const VectorType<ElementType>& digit = hashTable(z, Uz);
            GeNuSys::LinAlg::Operations::vct_sub(z, digit);
            GeNuSys::LinAlg::Operations::mat_mul(props.getAdjoint(), z, phiZ);
            GeNuSys::LinAlg::Operations::vct_idiv(phiZ, props.getDet());

            return digit;
        }

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType,
            typename Norm
            >
        std::vector<VectorType<ElementType>> NumberSystem<ElementType, VectorType, MatrixType, Norm>::getExpansion(const GeNuSys::LinAlg::Vector<ElementType>& z)
        {
            std::vector<VectorType<ElementType>> result;

            GeNuSys::LinAlg::Vector<ElementType> Uz = hash.createCache();
            GeNuSys::LinAlg::Vector<ElementType> act[2] = { z, z };
            int i = 0;
            do
            {
                result.push_back(phi(act[i], act[(i + 1) % 2], Uz));
                i = (i + 1) % 2;
            }
            while (GeNuSys::LinAlg::PNorm<00>::norm(act[i]) != ElementTraits<ElementType>::zero());

            return result;
        }

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType,
            typename Norm
            >
        std::vector<GeNuSys::LinAlg::Vector<ElementType>> NumberSystem<ElementType, VectorType, MatrixType, Norm>::getOrbit(const GeNuSys::LinAlg::Vector<ElementType>& z)
        {
            std::vector<GeNuSys::LinAlg::Vector<ElementType>> result;

            result.push_back(z);

            if (GeNuSys::LinAlg::PNorm<00>::norm(result.back()) == ElementTraits<ElementType>::zero())
            {
                return result;
            }

            for (;;)
            {
                result.push_back(phi(result.back()));
                if (GeNuSys::LinAlg::PNorm<00>::norm(result.back()) == ElementTraits<ElementType>::zero())
                {
                    return result;
                }
                for (size_t i = 0; i < result.size() - 1; ++i)
                {
                    if (result[i] == result.back())
                    {
                        result.pop_back();
                        return result;
                    }
                }
            }
        }

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType,
            typename Norm
            >
        std::vector<std::vector<GeNuSys::LinAlg::Vector<ElementType>>> NumberSystem<ElementType, VectorType, MatrixType, Norm>::getCycles()
        {
            const unsigned int N = props.getBase().getRows();

            std::vector<int> lowerBound, upperBound;
            Traits::getBounds(props.getInverse(), digitSet, lowerBound, upperBound);

            VectorCoder coder(lowerBound, upperBound);

            GeNuSys::LinAlg::Vector<ElementType> Uz = hash.createCache();
            GeNuSys::LinAlg::Vector<ElementType> act[2] = { GeNuSys::LinAlg::Vector<ElementType>(N), GeNuSys::LinAlg::Vector<ElementType>(N) };

            std::vector<std::vector<GeNuSys::LinAlg::Vector<ElementType>>> result;

            BitVector touched(coder.getSize());
            std::vector<unsigned long long> path;

            for (unsigned long long i = 0; i < coder.getSize(); ++i)
            {
                if ((i + 1) % (coder.getSize() / 80) == 0)
                {
                    // std::cout << "#";
                    //std::cout << i << " / " << coder.getSize() <<  std::endl;
                }

                if (touched[i])
                {
                    continue;
                }

                path.clear();

                coder.decode(i, act[0]);

                unsigned long long idx = i;
                path.push_back(idx);

                bool valid = true;

                int actIdx = 0;
                do
                {
                    touched.set(idx);

                    phi(act[actIdx], act[(actIdx + 1) % 2], Uz);
                    actIdx = (actIdx + 1) % 2;

                    idx = coder.encode(act[actIdx], valid);

                    path.push_back(idx);
                }
                while (valid && !touched[idx]);

                if (!valid)
                {
                    continue;
                }

                for (int j = path.size() - 2; j >= 0; --j)
                {
                    if (path[j] == path[path.size() - 1])
                    {
                        std::vector<GeNuSys::LinAlg::Vector<ElementType>> loop;
                        GeNuSys::LinAlg::Vector<ElementType> vct(N);
                        for (unsigned int k = j; k < path.size(); ++k)
                        {
                            coder.decode(path[k], vct);
                            loop.push_back(vct);
                        }

                        result.push_back(loop);

                        break;
                    }
                }
            }

            return result;
        }

    }
}
