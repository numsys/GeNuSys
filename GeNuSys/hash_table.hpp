namespace GeNuSys
{
    namespace NumSys
    {

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType
            >
        HashTable<ElementType, VectorType, MatrixType>::HashTable(const SmithHash<ElementType, MatrixType>& hash, const std::vector<VectorType<ElementType>>& digitSet): hash(hash), digitSet(digitSet)
        {
            // TODO: VERIFY CRS PROPERTY!!!
            hashTable = std::vector<VectorType<ElementType>>(digitSet.size());
            GeNuSys::LinAlg::Vector<ElementType> Uz = hash.createCache();
            for (unsigned int i = 0; i < digitSet.size(); ++i)
            {
                hashTable[ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(hash(digitSet[i], Uz))] = digitSet[i];
            }
        }

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType
            >
        const VectorType<ElementType>& HashTable<ElementType, VectorType, MatrixType>::operator()(const GeNuSys::LinAlg::Vector<ElementType>& z) const
        {
            return hashTable[hash(z, hash.createCache())];
        }

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType
            >
        const VectorType<ElementType>& HashTable<ElementType, VectorType, MatrixType>::operator()(const GeNuSys::LinAlg::Vector<ElementType>& z, GeNuSys::LinAlg::Vector<ElementType>& Uz) const
        {
            return hashTable[ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(hash(z, Uz))];
        }

    }
}
