namespace GeNuSys
{
    namespace Tests
    {

        template<typename T, typename S>
        bool TestUtils::equals(const GeNuSys::LinAlg::Vector<T>& vctA, const GeNuSys::LinAlg::Vector<S>& vctB)
        {
            if (vctA.getLength() != vctB.getLength())
            {
                return false;
            }

            bool equals = true;
            for (unsigned int i = 0; i < vctA.getLength() && equals; ++i)
            {
                equals = (vctA[i] == GeNuSys::ElementTraits<S>::template asType<T>(vctB[i]));
            }
            return equals;
        }

        template<typename T, typename S>
        bool TestUtils::equals(const GeNuSys::LinAlg::Vector<T>& vctA, const GeNuSys::LinAlg::SparseVector<S>& vctB)
        {
            if (vctA.getLength() != vctB.getLength())
            {
                return false;
            }

            bool equals = true;
            for (unsigned int i = 0; i < vctA.getLength() && equals; ++i)
            {
                equals = (vctA[i] == GeNuSys::ElementTraits<S>::template asType<T>(vctB[i]));
            }
            return equals;
        }

        template<typename T, typename S>
        bool TestUtils::equals(const GeNuSys::LinAlg::SparseVector<T>& vctA, const GeNuSys::LinAlg::Vector<S>& vctB)
        {
            if (vctA.getLength() != vctB.getLength())
            {
                return false;
            }

            bool equals = true;
            for (unsigned int i = 0; i < vctA.getLength() && equals; ++i)
            {
                equals = (vctA[i] == GeNuSys::ElementTraits<S>::template asType<T>(vctB[i]));
            }
            return equals;
        }

        template<typename T, typename S>
        bool TestUtils::equals(const GeNuSys::LinAlg::SparseVector<T>& vctA, const GeNuSys::LinAlg::SparseVector<S>& vctB)
        {
            if (vctA.getLength() != vctB.getLength())
            {
                return false;
            }

            bool equals = true;
            for (unsigned int i = 0; i < vctA.getLength() && equals; ++i)
            {
                equals = (vctA[i] == GeNuSys::ElementTraits<S>::template asType<T>(vctB[i]));
            }
            return equals;
        }

        template<typename T, typename S>
        bool TestUtils::equals(const GeNuSys::LinAlg::Matrix<T>& matA, const GeNuSys::LinAlg::Matrix<S>& matB)
        {
            if (matA.getRows() != matB.getRows() || matA.getCols() != matB.getCols())
            {
                return false;
            }

            bool equals = true;
            for (unsigned int i = 0; i < matA.getRows() && equals; ++i)
            {
                for (unsigned int j = 0; j < matA.getCols() && equals; ++j)
                {
                    equals = (matA(i, j) == GeNuSys::ElementTraits<S>::template asType<T>(matB(i, j)));
                }
            }
            return equals;
        }

        template<typename T, typename S>
        bool TestUtils::equals(const GeNuSys::LinAlg::Matrix<T>& matA, const GeNuSys::LinAlg::SparseMatrix<S>& matB)
        {
            if (matA.getRows() != matB.getRows() || matA.getCols() != matB.getCols())
            {
                return false;
            }

            bool equals = true;
            for (unsigned int i = 0; i < matA.getRows() && equals; ++i)
            {
                for (unsigned int j = 0; j < matA.getCols() && equals; ++j)
                {
                    equals = (matA(i, j) == GeNuSys::ElementTraits<S>::template asType<T>(matB(i, j)));
                }
            }
            return equals;
        }

        template<typename T, typename S>
        bool TestUtils::equals(const GeNuSys::LinAlg::SparseMatrix<T>& matA, const GeNuSys::LinAlg::Matrix<S>& matB)
        {
            if (matA.getRows() != matB.getRows() || matA.getCols() != matB.getCols())
            {
                return false;
            }

            bool equals = true;
            for (unsigned int i = 0; i < matA.getRows() && equals; ++i)
            {
                for (unsigned int j = 0; j < matA.getCols() && equals; ++j)
                {
                    equals = (matA(i, j) == GeNuSys::ElementTraits<S>::template asType<T>(matB(i, j)));
                }
            }
            return equals;
        }

        template<typename T, typename S>
        bool TestUtils::equals(const GeNuSys::LinAlg::SparseMatrix<T>& matA, const GeNuSys::LinAlg::SparseMatrix<S>& matB)
        {
            if (matA.getRows() != matB.getRows() || matA.getCols() != matB.getCols())
            {
                return false;
            }

            bool equals = true;
            for (unsigned int i = 0; i < matA.getRows() && equals; ++i)
            {
                for (unsigned int j = 0; j < matA.getCols() && equals; ++j)
                {
                    equals = (matA(i, j) == ElementTraits<S>::template asType<T>(matB(i, j)));
                }
            }
            return equals;
        }

    }
}
