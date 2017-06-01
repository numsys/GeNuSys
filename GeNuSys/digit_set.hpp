#include "linalg_traits.h"

namespace GeNuSys
{
    namespace NumSys
    {

        template<typename ElementType>
        std::vector<GeNuSys::LinAlg::SparseVector<ElementType>> DigitSet::getJCanonical(const RadixProperties<ElementType>& props, unsigned int j)
        {
            std::vector<GeNuSys::LinAlg::SparseVector<ElementType>> result;
            GeNuSys::LinAlg::SparseVector<ElementType> ej(props.getSize());
            for (int i = 0; i < props.getAbsDet(); ++i)
            {
                ej.set(j, ElementTraits<ElementType>::one() * i);
                result.push_back(ej);
            }

            return result;
        }

        template<typename ElementType>
        std::vector<GeNuSys::LinAlg::SparseVector<ElementType>> DigitSet::getJSymmetric(const RadixProperties<ElementType>& props, unsigned int j)
        {
            std::vector<GeNuSys::LinAlg::SparseVector<ElementType>> result;
            ElementType mid = props.getAbsDet() / 2;
            GeNuSys::LinAlg::SparseVector<ElementType> ej(props.getSize());
            for (int i = 0; i < props.getAbsDet(); ++i)
            {
                ej.set(j, ElementTraits<ElementType>::one() * (i - mid));
                result.push_back(ej);
            }

            return result;
        }

        template<typename ElementType>
        std::vector<GeNuSys::LinAlg::Vector<ElementType>> DigitSet::getAdjoint(const RadixProperties<ElementType>& props)
        {
            /////
            typedef typename ElementTraits<ElementType>::RationalType RationalType;

            GeNuSys::LinAlg::Matrix<ElementType> invU =
                GeNuSys::LinAlg::Traits::template convertUnsafe<RationalType, ElementType>(GeNuSys::LinAlg::Algorithms::invert(props.getSmithNormalForm().U));

            GeNuSys::LinAlg::Vector<ElementType> smithDiag = GeNuSys::LinAlg::Vector<ElementType>(props.getSize());
            for (unsigned int i = 0; i < props.getSize(); ++i)
            {
                smithDiag.set(i, props.getSmithNormalForm().S(i, i));
            }
            /////

            std::vector<GeNuSys::LinAlg::Vector<ElementType>> digits;
            GeNuSys::LinAlg::Vector<ElementType> v(props.getSize());
            digits.push_back(v);
            bool flag = true;
            while (flag && (int)digits.size() < props.getAbsDet())
            {
                int i;
                for (i = v.getLength() - 1; i >= 0 && v[i] == smithDiag[i] - 1; --i)
                {
                    v.set(i, ElementTraits<ElementType>::zero());
                }
                flag = (i >= 0);

                if (flag)
                {
                    v.set(i, v[i] + 1);
                    digits.push_back(invU * v);
                }
            }

            std::vector<GeNuSys::LinAlg::Vector<ElementType>> result;
            for (unsigned int i = 0; i < digits.size(); ++i)
            {
                GeNuSys::LinAlg::Vector<ElementType> aV = props.getAdjoint() * digits[i];
                GeNuSys::LinAlg::Operations::vct_mods(aV, props.getAbsDet());
                GeNuSys::LinAlg::Vector<ElementType> MaV = props.getBase() * aV;
                GeNuSys::LinAlg::Operations::vct_idiv(MaV, props.getDet());

                result.push_back(MaV);
            }

            return result;
        }

        template<typename ElementType>
        std::vector<GeNuSys::LinAlg::Vector<ElementType>> DigitSet::getDense(const RadixProperties<ElementType>& props)
        {
            std::vector<GeNuSys::LinAlg::Vector<ElementType>> adjoint = getAdjoint(props);
            GeNuSys::LinAlg::OperatorNorm<typename ElementTraits<ElementType>::RationalType> norm = props.getOperatorNorm();
            ElementType t = props.getAbsDet();

            std::vector<GeNuSys::LinAlg::Vector<ElementType>> result;
            for (unsigned int i = 0; i < adjoint.size(); ++i)
            {
                GeNuSys::LinAlg::Vector<ElementType> v = adjoint[i];

                bool flag;
                do
                {
                    flag = false;

                    for (unsigned int j = 0; j < v.getLength(); ++j)
                    {
                        ElementType start = v[j];

                        typename GeNuSys::LinAlg::OperatorNorm<ElementType>::template NormType<ElementType>::Type actNorm;
                        typename GeNuSys::LinAlg::OperatorNorm<ElementType>::template NormType<ElementType>::Type nextNorm = norm.norm(v);

                        do
                        {
                            actNorm = nextNorm;
                            v.set(j, v[j] + t);
                            nextNorm = norm.norm(v);
                        }
                        while (nextNorm < actNorm);
                        v.set(j, v[j] - t);
                        nextNorm = norm.norm(v);

                        do
                        {
                            actNorm = nextNorm;
                            v.set(j, v[j] - t);
                            nextNorm = norm.norm(v);
                        }
                        while (nextNorm < actNorm);
                        v.set(j, v[j] + t);

                        if (v[j] != start)
                        {
                            flag = true;
                        }
                    }
                }
                while (flag);

                result.push_back(v);
            }

            return result;
        }

    }
}
