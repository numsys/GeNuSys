namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename ElementType>
        typename FrobeniusNorm::template NormType<ElementType>::Type FrobeniusNorm::norm(const Matrix<ElementType>& mat)
        {
            typedef typename ElementTraits<ElementType>::AbsSqrType AbsSqrType;

            AbsSqrType sum = ElementTraits<AbsSqrType>::zero();
            for (unsigned int idx = 0; idx < mat.size(); ++idx)
            {
                sum += ElementTraits<ElementType>::absSqr(mat.elem[idx]);
            }

            return ElementTraits<AbsSqrType>::sqrt(sum);
        }

        template<typename ElementType>
        typename FrobeniusNorm::template NormType<ElementType>::Type FrobeniusNorm::norm(const SparseMatrix<ElementType>& mat)
        {
            typedef typename ElementTraits<ElementType>::AbsSqrType AbsSqrType;

            AbsSqrType sum = ElementTraits<AbsSqrType>::zero();
            for (unsigned int idx = 0; idx < mat.size(); ++idx)
            {
                sum += ElementTraits<AbsSqrType>::absSqr(mat.elem[idx].value);
            }

            return ElementTraits<AbsSqrType>::sqrt(sum);
        }

        template<typename ElementType>
        typename FrobeniusNorm::template NormType<ElementType>::Type FrobeniusNorm::operator()(const Matrix<ElementType>& mat) const
        {
            return norm(mat);
        }

        template<typename ElementType>
        typename FrobeniusNorm::template NormType<ElementType>::Type FrobeniusNorm::operator()(const SparseMatrix<ElementType>& mat) const
        {
            return norm(mat);
        }

    }
}
