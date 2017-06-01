namespace GeNuSys
{
    namespace NumSys
    {

        template<typename ElementType>
        GeNuSys::LinAlg::Matrix<ElementType> Simultaneous::getSimultaneous(ElementType a, ElementType b)
        {
            GeNuSys::LinAlg::Matrix<ElementType> simultaneous(6, 6);

            simultaneous.set(0, 0, a);
            simultaneous.set(0, 1, -b);
            simultaneous.set(1, 0, b);
            simultaneous.set(1, 1, a - b);

            simultaneous.set(2, 2, a + 1);
            simultaneous.set(2, 3, -b);
            simultaneous.set(3, 2, b);
            simultaneous.set(3, 3, a - b + 1);

            simultaneous.set(4, 4, a + 1);
            simultaneous.set(4, 5, -b - 1);
            simultaneous.set(5, 4, b + 1);
            simultaneous.set(5, 5, a - b);

            return simultaneous;
        }

        template<typename ElementType, typename Norm>
        void Simultaneous::tryDigitCandidate(
            const ElementType i, const ElementType j, GeNuSys::LinAlg::Vector<ElementType>& vct,
            std::vector<bool>& digitFound, std::vector<typename Norm::template NormType<ElementType>::Type>& digitNorm, std::vector<GeNuSys::LinAlg::Vector<ElementType>>& digits,
            int& numFound,
            const Norm& norm, typename Norm::template NormType<ElementType>::Type& maxNorm, bool& flag,
            const SmithHash<ElementType, GeNuSys::LinAlg::Matrix>& hash, GeNuSys::LinAlg::Vector<ElementType>& Uz)
        {
            typedef typename Norm::template NormType<ElementType>::Type NormType;

            vct.set(0, i);
            vct.set(2, i);
            vct.set(4, i);
            vct.set(1, j);
            vct.set(3, j);
            vct.set(5, j);

            NormType n = norm(vct);
            if (n <= maxNorm)
            {
                flag = true;
            }

            ElementType h = hash(vct, Uz);
            if (!digitFound[ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(h)])
            {
                digitFound[ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(h)] = true;
                ++numFound;

                digitNorm[ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(h)] = n;
                digits[ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(h)] = vct;
                if (maxNorm < n)
                {
                    maxNorm = n;
                }
            }
            else if (n <= digitNorm[ElementTraits<ElementType>::template asType<long int>(h)])
            {
                digitNorm[ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(h)] = n;
                digits[ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(h)] = vct;
                if (maxNorm < n)
                {
                    maxNorm = n;
                }
            }
        }

        template<typename ElementType, typename Norm>
        std::vector<GeNuSys::LinAlg::Vector<ElementType>> Simultaneous::getDigitSet(const RadixProperties<ElementType>& props, const Norm& norm, bool global)
        {
            typedef typename Norm::template NormType<ElementType>::Type NormType;

            SmithHash<ElementType, GeNuSys::LinAlg::Matrix> hash(props);
            GeNuSys::LinAlg::Vector<ElementType> Uz = hash.createCache();

            int numFound = 0;
            NormType maxNorm = ElementTraits<NormType>::zero();

            std::vector<bool> digitFound(ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(props.getAbsDet()), false);
            std::vector<NormType> digitNorm(ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(props.getAbsDet()));
            std::vector<GeNuSys::LinAlg::Vector<ElementType>> digits(ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(props.getAbsDet()));

            bool flag = false;
            GeNuSys::LinAlg::Vector<ElementType> cand(6);
            tryDigitCandidate(ElementTraits<ElementType>::zero(), ElementTraits<ElementType>::zero(), cand, digitFound, digitNorm, digits, numFound, norm, maxNorm, flag, hash, Uz);
            for (ElementType i = 1; numFound < props.getAbsDet() || (global && flag); ++i)
            {
                flag = false;
                for (ElementType j = -i + 1; j < i; ++j)
                {
                    tryDigitCandidate(ElementType(-i), j, cand, digitFound, digitNorm, digits, numFound, norm, maxNorm, flag, hash, Uz);
                    tryDigitCandidate(i, j, cand, digitFound, digitNorm, digits, numFound, norm, maxNorm, flag, hash, Uz);
                    tryDigitCandidate(j, i, cand, digitFound, digitNorm, digits, numFound, norm, maxNorm, flag, hash, Uz);
                    tryDigitCandidate(j, ElementType(-i), cand, digitFound, digitNorm, digits, numFound, norm, maxNorm, flag, hash, Uz);
                }
                tryDigitCandidate(ElementType(-i), ElementType(-i), cand, digitFound, digitNorm, digits, numFound, norm, maxNorm, flag, hash, Uz);
                tryDigitCandidate(ElementType(-i), i, cand, digitFound, digitNorm, digits, numFound, norm, maxNorm, flag, hash, Uz);
                tryDigitCandidate(i, i, cand, digitFound, digitNorm, digits, numFound, norm, maxNorm, flag, hash, Uz);
                tryDigitCandidate(i, ElementType(-i), cand, digitFound, digitNorm, digits, numFound, norm, maxNorm, flag, hash, Uz);
            }

            return digits;
        }

    }
}
