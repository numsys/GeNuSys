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
