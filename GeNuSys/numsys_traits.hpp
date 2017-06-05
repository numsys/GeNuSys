/*
GeNuSys - computations with generalized number systems
Copyright (C) 2015-2017  Bence Németh
Copyright (C) 2017  Tamás Krutki

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>

namespace GeNuSys
{
    namespace NumSys
    {

        template <
            typename ElementType,
            template<typename> class VectorType
            >
        void Traits::getBounds(const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM,
                               const std::vector<VectorType<ElementType>>& digitSet,
                               const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& T,
                               std::vector<int>& lowerBound, std::vector<int>& upperBound)
        {

            typedef typename ElementTraits<ElementType>::RationalType RationalType;
            typedef typename GeNuSys::LinAlg::PNorm<00>::template NormType<RationalType>::Type NormType;

            const unsigned int N = invM.getRows();

            std::vector<VectorType<RationalType>> rationalDigits(digitSet.size());
            for (unsigned int i = 0; i < digitSet.size(); ++i)
            {
                rationalDigits[i] = digitSet[i];
                rationalDigits[i] = T * rationalDigits[i];
            }

            GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType> base = T * invM * GeNuSys::LinAlg::Algorithms::invert(T);

            std::vector<RationalType> low(N, 0);
            std::vector<RationalType> high(N, 0);

            GeNuSys::LinAlg::Matrix<RationalType> X[2] = { base, GeNuSys::LinAlg::Matrix<RationalType>::identity(N, N) };
            std::vector<VectorType<RationalType>> multipliedDigits(digitSet.size(), VectorType<RationalType>(N));

            int actXIdx = 0;
            do
            {
                for (unsigned int i = 0; i < rationalDigits.size(); ++i)
                {
                    GeNuSys::LinAlg::Operations::mat_mul(X[actXIdx], rationalDigits[i], multipliedDigits[i]);
                }
                for (unsigned int i = 0; i < N; ++i)
                {
                    RationalType min = multipliedDigits[0][i];
                    RationalType max = multipliedDigits[0][i];
                    for (unsigned int j = 1; j < multipliedDigits.size(); ++j)
                    {
                        if (multipliedDigits[j][i] < min)
                        {
                            min = multipliedDigits[j][i];
                        }
                        if (multipliedDigits[j][i] > max)
                        {
                            max = multipliedDigits[j][i];
                        }
                    }
                    low[i] += max;
                    high[i] += min;
                }
                GeNuSys::LinAlg::Operations::mat_mul(base, X[actXIdx], X[(actXIdx + 1) % 2]);
                actXIdx = (actXIdx + 1) % 2;
            }
            while (GeNuSys::LinAlg::PNorm<00>::norm(X[actXIdx]) > ElementTraits<NormType>::epsilon());

            typename ElementTraits<NormType>::RationalType coef =
                ElementTraits<NormType>::div(ElementTraits<NormType>::one(), ElementTraits<NormType>::one() - GeNuSys::LinAlg::PNorm<00>::norm(X[actXIdx]));

            lowerBound = std::vector<int>(N);
            upperBound = std::vector<int>(N);
            for (unsigned int i = 0; i < N; ++i)
            {
                lowerBound[i] = ElementTraits<RationalType>::template asTypeUnsafe<int>(std::ceil(ElementTraits<RationalType>::template asTypeUnsafe<double>(-low[i] * coef)));
                upperBound[i] = ElementTraits<RationalType>::template asTypeUnsafe<int>(std::floor(ElementTraits<RationalType>::template asTypeUnsafe<double>(-high[i] * coef)));
            }
        }

        template <
            typename ElementType,
            template<typename> class VectorType
            >
        void Traits::getBounds(const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM,
                               const std::vector<VectorType<ElementType>>& digitSet,
                               std::vector<int>& lowerBound, std::vector<int>& upperBound)
        {
            getBounds(invM, digitSet, GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>::identity(invM.getRows(), invM.getCols()), lowerBound, upperBound);
        }

        template <
            typename ElementType,
            template<typename> class VectorType
            >
        unsigned long long Traits::getVolume(const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM,
                                             const std::vector<VectorType<ElementType>>& digitSet,
                                             const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& T)
        {

            std::vector<int> lowerBound;
            std::vector<int> upperBound;
            getBounds(invM, digitSet, T, lowerBound, upperBound);

            unsigned long long size = 1;
            for (unsigned int i = 0; i < lowerBound.size(); ++i)
            {
                size *= upperBound[i] - lowerBound[i] + 1;
            }

            return size;
        }

        template <
            typename ElementType,
            template<typename> class VectorType
            >
        unsigned long long Traits::getVolume(const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM,
                                             const std::vector<VectorType<ElementType>>& digitSet)
        {

            std::vector<int> lowerBound;
            std::vector<int> upperBound;
            getBounds(invM, digitSet, lowerBound, upperBound);

            unsigned long long size = 1;
            for (unsigned int i = 0; i < lowerBound.size(); ++i)
            {
                size *= upperBound[i] - lowerBound[i] + 1;
            }

            return size;
        }

        template<typename ElementType>
        struct Transformation
        {

            GeNuSys::LinAlg::Matrix<ElementType> T;

            unsigned long long vol;

            Transformation() {}

            Transformation(const GeNuSys::LinAlg::Matrix<ElementType>& T, unsigned long long vol): T(T), vol(vol) {}

            bool operator <(const Transformation<ElementType>& tr) const
            {
                return vol < tr.vol;
            }

        };

        template <
            typename ElementType,
            template<typename> class VectorType
            >
        GeNuSys::LinAlg::Matrix<ElementType> Traits::findBasisTransformation(
            const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM, const std::vector<VectorType<ElementType>>& digitSet,
            const unsigned int candNum, const unsigned int mutateNum, const unsigned int noImprLimit)
        {
            typedef typename ElementTraits<ElementType>::RationalType RationalType;

            const unsigned int N = invM.getRows();

            unsigned long long origVol = getVolume(invM, digitSet);

            std::vector<Transformation<RationalType>> candidates(candNum, Transformation<RationalType>(GeNuSys::LinAlg::Matrix<RationalType>::identity(N, N), origVol));

            unsigned int noImpr = 0;
            while (noImpr < noImprLimit)
            {
                std::vector<Transformation<RationalType>> newCandidates;
                for (unsigned int i = 0; i < candidates.size(); ++i)
                {
                    for (unsigned int j = 0; j < mutateNum; ++j)
                    {
                        unsigned int x, y;
                        do
                        {
                            x = rand() % N;
                            y = rand() % N;
                        }
                        while (x == y);

                        if (x > y)
                        {
                            unsigned int t = x;
                            x = y;
                            y = t;
                        }

                        GeNuSys::LinAlg::Matrix<RationalType> U = candidates[i].T;
                        unsigned long long origVolU;
                        unsigned long long volU = getVolume(invM, digitSet, U);
                        do
                        {
                            origVolU = volU;
                            U.set(x, y, U(x, y) + 1);
                            volU = getVolume(invM, digitSet, U);
                        }
                        while (origVolU > volU);
                        U.set(x, y, U(x, y) - 1);

                        newCandidates.push_back(Transformation<RationalType>(U, origVolU));

                        GeNuSys::LinAlg::Matrix<RationalType> V = candidates[i].T;
                        unsigned long long origVolV;
                        unsigned long long volV = getVolume(invM, digitSet, V);
                        do
                        {
                            origVolV = volV;
                            V.set(x, y, V(x, y) - 1);
                            volV = getVolume(invM, digitSet, V);
                        }
                        while (origVolV > volV);
                        V.set(x, y, V(x, y) + 1);

                        newCandidates.push_back(Transformation<RationalType>(V, origVolV));
                    }
                }

                std::sort(newCandidates.begin(), newCandidates.end());

                if (newCandidates[0].vol < candidates[0].vol)
                {
                    noImpr = 0;
                }
                else
                {
                    noImpr += 1;
                }

                for (unsigned int i = 0; i < candidates.size(); ++i)
                {
                    candidates[i] = newCandidates[i];
                }
            }

            //std::cout << "REDUCED : " << origVol << " -> " << candidates[0].vol << std::endl;
            return GeNuSys::LinAlg::Traits::template convertUnsafe<RationalType, ElementType>(candidates[0].T);
        }

    }
}
