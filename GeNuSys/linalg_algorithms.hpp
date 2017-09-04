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

#include <algorithm>

#include "linalg_traits.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename ElementType>
        LU<ElementType>::LU(
            const Matrix<typename ElementTraits<ElementType>::RationalType>& L,
            const Matrix<typename ElementTraits<ElementType>::RationalType>& U): L(L), U(U)
        {
        }

        template<typename ElementType>
        QR<ElementType>::QR(
            const Matrix<typename ElementTraits<ElementType>::RealType>& Q,
            const Matrix<typename ElementTraits<ElementType>::RealType>& R): Q(Q), R(R)
        {
        }

        template<typename ElementType>
        SmithNormalForm<ElementType>::SmithNormalForm(const Matrix<ElementType>& S,
                                                      const Matrix<ElementType>& U,
                                                      const Matrix<ElementType>& V): S(S), U(U), V(V)
        {
        }

        template<typename ElementType>
        HessenbergForm<ElementType>::HessenbergForm(
            const Matrix<typename ElementTraits<ElementType>::RealType>& Q,
            const Matrix<typename ElementTraits<ElementType>::RealType>& A,
            const Matrix<typename ElementTraits<ElementType>::RealType>& QT): Q(Q), A(A), QT(QT)
        {
        }

        template<typename ElementType>
        SchurForm<ElementType>::SchurForm(
            const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType>& Q,
            const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType>& U,
            const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType>& QT): Q(Q), U(U), QT(QT)
        {
        }

        template<typename ElementType>
        JordanForm<ElementType>::JordanForm(
            const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType>& P,
            const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType>& J,
            const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType>& invP): P(P), J(J), invP(invP)
        {
        }

        template<typename ElementType>
        typename ElementTraits<ElementType>::RationalType Algorithms::det(const Matrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(mat.cols == mat.rows, std::length_error);

            typedef typename ElementTraits<ElementType>::RationalType RationalType;

            const unsigned int N = mat.rows;

            Matrix<RationalType> U = mat;
            int perm = 1;

            for (unsigned int i = 0, idxPivot = 0; i < N - 1; ++i, idxPivot += N + 1)
            {
                bool pivotFound = false;
                RationalType pivot = U.elem[idxPivot];
                unsigned int pivotRow = i;
                for (unsigned int j = i, idx = idxPivot; j < N; ++j, idx += N)
                {
                    if (U.elem[idx] == ElementTraits<RationalType>::zero())
                    {
                        continue;
                    }
                    RationalType val = ElementTraits<RationalType>::abs(U.elem[idx]);
                    if (!pivotFound || pivot < val)
                    {
                        pivotFound = true;

                        pivot = val;
                        pivotRow = j;
                    }
                }
                if (!pivotFound)
                {
                    return ElementTraits<RationalType>::zero();
                }
                if (pivotRow != i)
                {
                    unsigned int idxRow1 = i * N + i;
                    unsigned int idxRow2 = pivotRow * N + i;
                    for (unsigned int j = i; j < N; ++j, ++idxRow1, ++idxRow2)
                    {
                        std::swap(U.elem[idxRow1], U.elem[idxRow2]);
                    }
                    perm = -perm;
                }

                for (unsigned int j = i + 1, idxERow = idxPivot + N; j < N; ++j, idxERow += N)
                {
                    RationalType coef = U.elem[idxERow] / U.elem[idxPivot];
                    for (unsigned int k = i, idxP = idxPivot, idxE = idxERow; k < N; ++k, ++idxP, ++idxE)
                    {
                        U.elem[idxE] -= U.elem[idxP] * coef;
                    }
                }
            }

            RationalType prod = ElementTraits<RationalType>::one();
            for (unsigned int i = 0, idxU = 0; i < mat.rows; ++i, idxU += mat.cols + 1)
            {
                prod *= U.elem[idxU];
            }
            return perm * prod;
        }

        template<typename ElementType>
        Matrix<typename ElementTraits<ElementType>::RationalType> Algorithms::invert(const Matrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(mat.cols == mat.rows, std::length_error);

            typedef typename ElementTraits<ElementType>::RationalType RationalType;

            const unsigned int N = mat.rows;
            Matrix<RationalType> A = mat;
            Matrix<RationalType> I = Matrix<RationalType>::identity(N, N);
            for (unsigned int i = 0, idxPivotRow = 0; i < N; ++i, idxPivotRow += N)
            {
                bool pivotFound = false;
                typename ElementTraits<RationalType>::AbsType pivotMax;
                unsigned int pivotRow;
                for (unsigned int j = i, idxA = idxPivotRow + i; j < N; ++j, idxA += N)
                {
                    if (A.elem[idxA] == ElementTraits<RationalType>::zero())
                    {
                        continue;
                    }
                    typename ElementTraits<RationalType>::AbsType val = ElementTraits<RationalType>::abs(A.elem[idxA]);
                    if (!pivotFound || pivotMax < val)
                    {
                        pivotFound = true;
                        pivotMax = val;
                        pivotRow = j;
                    }
                }
                if (!pivotFound)
                {
                    ASSERT_EXCEPTION(!"NOT INVERTIBLE", std::logic_error);
                    break; // NOT INVERTIBLE
                }
                if (pivotRow != i)
                {
                    unsigned int idxP = idxPivotRow;
                    unsigned int idxE = pivotRow * N;
                    for (unsigned int j = 0; j < i; ++j, ++idxP, ++idxE)
                    {
                        std::swap(I.elem[idxP], I.elem[idxE]);
                    }
                    for (unsigned int j = i; j < N; ++j, ++idxP, ++idxE)
                    {
                        std::swap(I.elem[idxP], I.elem[idxE]);
                        std::swap(A.elem[idxP], A.elem[idxE]);
                    }
                }

                for (unsigned int j = i + 1, idxEliminateRow = idxPivotRow + N; j < N; ++j, idxEliminateRow += N)
                {
                    RationalType coef = A.elem[idxEliminateRow + i] / A.elem[idxPivotRow + i];
                    unsigned int idxP = idxPivotRow;
                    unsigned int idxE = idxEliminateRow;
                    for (unsigned int k = 0; k < i; ++k, ++idxP, ++idxE)
                    {
                        I.elem[idxE] -= I.elem[idxP] * coef;
                    }
                    for (unsigned int k = i; k < N; ++k, ++idxP, ++idxE)
                    {
                        A.elem[idxE] -= A.elem[idxP] * coef;
                        I.elem[idxE] -= I.elem[idxP] * coef;
                    }
                }
            }
            for (unsigned int i = N - 1, idxPivotRow = (N - 1) * N; i > 0; --i, idxPivotRow -= N)
            {
                for (unsigned int j = 0, idxP = idxPivotRow; j < N; ++j, ++idxP)
                {
                    I.elem[idxP] /= A.elem[idxPivotRow + i];
                }
                for (unsigned int j = 0, idxEliminateRow = 0; j < i; ++j, idxEliminateRow += N)
                {
                    for (unsigned int k = 0, idxP = idxPivotRow, idxE = idxEliminateRow; k < N; ++k, ++idxP, ++idxE)
                    {
                        I.elem[idxE] -= I.elem[idxP] * A.elem[idxEliminateRow + i];
                    }
                }
            }
            for (unsigned int j = 0, idxP = 0; j < N; ++j, ++idxP)
            {
                I.elem[idxP] /= A.elem[0];
            }

            return I;
        }

#ifdef __unix__

        template<typename ElementType>
        Matrix<ElementType> Algorithms::getAdjoint(const Matrix<ElementType>& mat)
        {
            Matrix<mpz_class> mpz_mat = mat;
            Matrix<mpz_class> mpz_adj = Traits::convertUnsafe<mpq_class, mpz_class>(Algorithms::invert(mpz_mat) * Algorithms::det(mpz_mat));

            Matrix<ElementType> adj(mpz_adj.getRows(), mpz_adj.getCols());
            for (unsigned int i = 0; i < mpz_adj.size(); ++i)
            {
                adj.elem[i] = mpz_get_si(mpz_adj.elem[i].get_mpz_t());
            }

            return adj;
        }

#else

        template<typename ElementType>
        Matrix<ElementType> Algorithms::getAdjoint(const Matrix<ElementType>& mat)
        {
            typedef typename ElementTraits<ElementType>::RationalType RationalType;

            return Traits::convertUnsafe<RationalType, ElementType>(Algorithms::invert(mat) * Algorithms::det(mat));
        }

#endif // __unix__

        template<typename ElementType>
        LU<ElementType> Algorithms::decomposeLU(const Matrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(mat.cols == mat.rows, std::length_error);

            typedef typename ElementTraits<ElementType>::RationalType RationalType;

            const unsigned int N = mat.rows;
            Matrix<RationalType> L = Matrix<RationalType>::identity(N, N);
            Matrix<RationalType> U(N, N);

            for (unsigned int row = 0; row < N; ++row)
            {
                for (unsigned int i = row; i < N; ++i)
                {
                    RationalType sum = ElementTraits<RationalType>::zero();
                    for (unsigned int j = 0; j < row; ++j)
                    {
                        sum += L.elem[row * N + j] * U.elem[j * N + i];
                    }
                    U.elem[row * N + i] = ElementTraits<ElementType>::template asType<RationalType>(mat.elem[row * N + i]) - sum;
                }
                for (unsigned int i = row + 1; i < N; ++i)
                {
                    RationalType sum = ElementTraits<RationalType>::zero();
                    for (unsigned int j = 0; j < row; ++j)
                    {
                        sum += L.elem[i * N + j] * U.elem[j * N + row];
                    }
                    ASSERT_EXCEPTION(U.elem[row * N + row] != ElementTraits<RationalType>::zero(), std::logic_error);
                    L.elem[i * N + row] = (ElementTraits<ElementType>::template asType<RationalType>(mat.elem[i * N + row]) - sum) / U.elem[row * N + row];
                }
            }

            return LU<ElementType>(L, U);
        }

        template<typename ElementType>
        Matrix<typename ElementTraits<ElementType>::RealType> Algorithms::getHouseholderMatrix(const Vector<ElementType>& vct, const unsigned int N)
        {
            typedef typename ElementTraits<ElementType>::RealType RealType;
            typedef typename ElementTraits<RealType>::AbsType AbsRealType;

            const unsigned int M = vct.getLength();

            Vector<RealType> houseVct = vct;
            RealType sgn = ElementTraits<RealType>::sgn(houseVct.elem[0]);
            if (sgn == ElementTraits<RealType>::zero())
            {
                sgn = ElementTraits<RealType>::one();
            }
            houseVct.elem[0] += sgn * PNorm<2>::norm(houseVct);
            RealType normSqr = PNorm<2>::normSqr(houseVct);

            Matrix<RealType> result = Matrix<RealType>::identity(N, N);
            if (ElementTraits<RealType>::abs(normSqr) <= ElementTraits<AbsRealType>::epsilon())
            {
                return result;
            }
            for (unsigned int i = 0, idx = (N - M) * (N + 1); i < M; ++i, idx += N - M)
            {
                for (unsigned int j = 0; j < M; ++j, ++idx)
                {
                    result.elem[idx] -= ElementTraits<int>::template asType<RealType>(2) * houseVct[i] * ElementTraits<RealType>::conj(houseVct[j]) / normSqr;
                }
            }

            return result;
        }

        template<typename ElementType>
        QR<ElementType> Algorithms::decomposeQR(const Matrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(mat.cols == mat.rows, std::length_error);

            typedef typename ElementTraits<ElementType>::RealType RealType;

            const unsigned int N = mat.cols;
            Matrix<RealType> Q = Matrix<RealType>::identity(N, N);
            Matrix<RealType> R = mat;

            for (unsigned int col = 0; col < N - 1; ++col)
            {
                Matrix<RealType> transform = getHouseholderMatrix(Traits::getCol(R, col, col, N), N);

                Q = transform * Q;
                R = transform * R;
            }

            return QR<ElementType>(Q.conjugateTranspose(), R);
        }

        template<typename ElementType>
        SmithNormalForm<ElementType> Algorithms::getSmithNormalForm(const Matrix<ElementType>& mat)
        {
            Matrix<ElementType> v = Matrix<ElementType>::identity(mat.rows, mat.cols);
            Matrix<ElementType> s = mat;
            Matrix<ElementType> u = Matrix<ElementType>::identity(mat.rows, mat.cols);

            for (unsigned int i = 0, idxPivot = 0; i < mat.rows - 1; ++i, idxPivot += mat.cols + 1)
            {
                //
                // STEP 1: PIVOT
                //
                bool pivotFound = false;
                ElementType min = ElementTraits<ElementType>::zero();
                unsigned int pivotRow = i;
                unsigned int pivotCol = i;
                for (unsigned int j = i, idx = idxPivot; j < mat.rows; ++j, idx += i)
                {
                    for (unsigned int k = i; k < mat.cols; ++k, ++idx)
                    {
                        if (s.elem[idx] == ElementTraits<ElementType>::zero())
                        {
                            continue;
                        }
                        ElementType val = ElementTraits<ElementType>::abs(s.elem[idx]);
                        if (!pivotFound || min > val)
                        {
                            pivotFound = true;

                            min = val;
                            pivotRow = j;
                            pivotCol = k;
                        }
                    }
                }
                if (!pivotFound)
                {
                    break;
                }
                if (pivotRow != i)
                {
                    Traits::swapRows(s, pivotRow, i);
                    Traits::swapRows(u, pivotRow, i);
                }
                if (pivotCol != i)
                {
                    Traits::swapCols(s, pivotCol, i);
                    Traits::swapCols(v, pivotCol, i);
                }

                //
                // STEP 2: TRANSFORM COLS / ROWS
                //
                bool flag;
                do
                {
                    flag = false;

                    // ROW TRANSFORM
                    for (unsigned int j = i + 1, idxPCol = idxPivot + mat.cols; j < mat.rows; ++j, idxPCol += mat.cols)
                    {
                        if (!ElementTraits<ElementType>::divisible(s.elem[idxPCol], s.elem[idxPivot]))
                        {
                            ExtendedGCD<ElementType> gcd = ElementTraits<ElementType>::egcd(s.elem[idxPivot], s.elem[idxPCol]);

                            ElementType coef1 = -gcd.b / gcd.gcd;
                            ElementType coef2 = gcd.a / gcd.gcd;

                            unsigned int idxR1 = i * mat.cols;
                            unsigned int idxR2 = j * mat.cols;
                            for (unsigned int k = 0; k < i; ++k, ++idxR1, ++idxR2)
                            {
                                ElementType ui = gcd.cA * u.elem[idxR1] + gcd.cB * u.elem[idxR2];
                                ElementType uj = coef1 * u.elem[idxR1] + coef2 * u.elem[idxR2];
                                u.elem[idxR1] = ui;
                                u.elem[idxR2] = uj;
                            }
                            for (unsigned int k = i; k < mat.cols; ++k, ++idxR1, ++idxR2)
                            {
                                ElementType ui = gcd.cA * u.elem[idxR1] + gcd.cB * u.elem[idxR2];
                                ElementType uj = coef1 * u.elem[idxR1] + coef2 * u.elem[idxR2];
                                u.elem[idxR1] = ui;
                                u.elem[idxR2] = uj;

                                ElementType si = gcd.cA * s.elem[idxR1] + gcd.cB * s.elem[idxR2];
                                ElementType sj = coef1 * s.elem[idxR1] + coef2 * s.elem[idxR2];
                                s.elem[idxR1] = si;
                                s.elem[idxR2] = sj;
                            }

                            flag = true;
                        }
                    }

                    // COL TRANSFORM
                    for (unsigned int j = i + 1, idxPRow = idxPivot + 1; j < mat.cols; ++j, idxPRow += 1)
                    {
                        if (!ElementTraits<ElementType>::divisible(s.elem[idxPRow], s.elem[idxPivot]))
                        {
                            ExtendedGCD<ElementType> gcd = ElementTraits<ElementType>::egcd(s.elem[idxPivot], s.elem[idxPRow]);

                            ElementType coef1 = -gcd.b / gcd.gcd;
                            ElementType coef2 = gcd.a / gcd.gcd;

                            unsigned int idxC1 = i;
                            unsigned int idxC2 = j;
                            for (unsigned int k = 0; k < i; ++k, idxC1 += mat.cols, idxC2 += mat.cols)
                            {
                                ElementType vi = gcd.cA * v.elem[idxC1] + gcd.cB * v.elem[idxC2];
                                ElementType vj = coef1 * v.elem[idxC1] + coef2 * v.elem[idxC2];
                                v.elem[idxC1] = vi;
                                v.elem[idxC2] = vj;
                            }
                            for (unsigned int k = i; k < mat.rows; ++k, idxC1 += mat.cols, idxC2 += mat.cols)
                            {
                                ElementType vi = gcd.cA * v.elem[idxC1] + gcd.cB * v.elem[idxC2];
                                ElementType vj = coef1 * v.elem[idxC1] + coef2 * v.elem[idxC2];
                                v.elem[idxC1] = vi;
                                v.elem[idxC2] = vj;

                                ElementType si = gcd.cA * s.elem[idxC1] + gcd.cB * s.elem[idxC2];
                                ElementType sj = coef1 * s.elem[idxC1] + coef2 * s.elem[idxC2];
                                s.elem[idxC1] = si;
                                s.elem[idxC2] = sj;
                            }

                            flag = true;
                        }
                    }
                }
                while (flag);

                //
                // STEP 3: ELIMINATE COLS / ROWS
                //

                // ROW ELIMINATE
                for (unsigned int k = i + 1, idxCol = idxPivot - i + mat.cols; k < mat.rows; ++k, idxCol += mat.cols)
                {
                    ElementType coef = s.elem[idxCol + i] / s.elem[idxPivot];
                    unsigned int idxPRow = idxPivot - i;
                    unsigned int idxERow = idxCol;
                    for (unsigned int j = 0; j < i; ++j, ++idxPRow, ++idxERow)
                    {
                        u.elem[idxERow] -= u.elem[idxPRow] * coef;
                    }
                    for (unsigned int j = i; j < mat.cols; ++j, ++idxPRow, ++idxERow)
                    {
                        u.elem[idxERow] -= u.elem[idxPRow] * coef;
                        s.elem[idxERow] -= s.elem[idxPRow] * coef;
                    }
                }

                // COL ELIMINATE
                for (unsigned int k = i + 1, idxRow = i + 1; k < mat.rows; ++k, ++idxRow)
                {
                    ElementType coef = s.elem[idxRow + idxPivot - i] / s.elem[idxPivot];
                    unsigned int idxPCol = i;
                    unsigned int idxECol = idxRow;
                    for (unsigned int j = 0; j < i; ++j, idxPCol += mat.cols, idxECol += mat.cols)
                    {
                        v.elem[idxECol] -= v.elem[idxPCol] * coef;
                    }
                    for (unsigned int j = i; j < mat.rows; ++j, idxPCol += mat.cols, idxECol += mat.cols)
                    {
                        v.elem[idxECol] -= v.elem[idxPCol] * coef;
                        s.elem[idxECol] -= s.elem[idxPCol] * coef;
                    }
                }
            }

            //
            // FINAL STEP: CORRECT
            //
            bool corrected;
            do
            {
                corrected = false;

                for (unsigned int i = 0, idxPivot = 0; i < mat.rows - 1; ++i, idxPivot += mat.cols + 1)
                {
                    if (s.elem[idxPivot] == ElementTraits<ElementType>::zero())
                    {
                        break;
                    }

                    if (ElementTraits<ElementType>::divisible(s.elem[idxPivot + mat.cols + 1], s.elem[idxPivot]))
                    {
                        continue;
                    }

                    // ADD COLS (i. := i. + (i + 1).)
                    for (unsigned int j = 0, idxC1 = i, idxC2 = i + 1; j < mat.rows; ++j, idxC1 += mat.cols, idxC2 += mat.cols)
                    {
                        v.elem[idxC1] += v.elem[idxC2];
                    }
                    s.elem[idxPivot + mat.cols] = s.elem[idxPivot + mat.cols + 1];

                    // TRANSFORM
                    bool flag;
                    do
                    {
                        flag = false;

                        // TRANSFORM ROW
                        if (!ElementTraits<ElementType>::divisible(s.elem[idxPivot + mat.cols], s.elem[idxPivot]))
                        {
                            ExtendedGCD<ElementType> gcd = ElementTraits<ElementType>::egcd(s.elem[idxPivot], s.elem[idxPivot + mat.cols]);

                            ElementType coef1 = -gcd.b / gcd.gcd;
                            ElementType coef2 = gcd.a / gcd.gcd;

                            unsigned int idxR1 = i * mat.cols;
                            unsigned int idxR2 = idxR1 + mat.cols;
                            for (unsigned int k = 0; k < mat.cols; ++k, ++idxR1, ++idxR2)
                            {
                                ElementType ui = gcd.cA * u.elem[idxR1] + gcd.cB * u.elem[idxR2];
                                ElementType uj = coef1 * u.elem[idxR1] + coef2 * u.elem[idxR2];
                                u.elem[idxR1] = ui;
                                u.elem[idxR2] = uj;
                            }

                            ElementType si1 = gcd.cA * s.elem[idxPivot] + gcd.cB * s.elem[idxPivot + mat.cols];
                            ElementType sj1 = coef1 * s.elem[idxPivot] + coef2 * s.elem[idxPivot + mat.cols];
                            s.elem[idxPivot] = si1;
                            s.elem[idxPivot + mat.cols] = sj1;

                            ElementType si2 = gcd.cA * s.elem[idxPivot + 1] + gcd.cB * s.elem[idxPivot + mat.cols + 1];
                            ElementType sj2 = coef1 * s.elem[idxPivot + 1] + coef2 * s.elem[idxPivot + mat.cols + 1];
                            s.elem[idxPivot + 1] = si2;
                            s.elem[idxPivot + mat.cols + 1] = sj2;

                            flag = true;
                        }

                        // TRANSFORM COL
                        if (!ElementTraits<ElementType>::divisible(s.elem[idxPivot + 1], s.elem[idxPivot]))
                        {
                            ExtendedGCD<ElementType> gcd = ElementTraits<ElementType>::egcd(s.elem[idxPivot], s.elem[idxPivot + 1]);

                            ElementType coef1 = -gcd.b / gcd.gcd;
                            ElementType coef2 = gcd.a / gcd.gcd;

                            unsigned int idxC1 = i;
                            unsigned int idxC2 = idxC1 + 1;
                            for (unsigned int k = 0; k < mat.rows; ++k, idxC1 += mat.cols, idxC2 += mat.cols)
                            {
                                ElementType vi = gcd.cA * v.elem[idxC1] + gcd.cB * v.elem[idxC2];
                                ElementType vj = coef1 * v.elem[idxC1] + coef2 * v.elem[idxC2];
                                v.elem[idxC1] = vi;
                                v.elem[idxC2] = vj;
                            }

                            ElementType si1 = gcd.cA * s.elem[idxPivot] + gcd.cB * s.elem[idxPivot + 1];
                            ElementType sj1 = coef1 * s.elem[idxPivot] + coef2 * s.elem[idxPivot + 1];
                            s.elem[idxPivot] = si1;
                            s.elem[idxPivot + 1] = sj1;

                            ElementType si2 = gcd.cA * s.elem[idxPivot + mat.cols] + gcd.cB * s.elem[idxPivot + mat.cols + 1];
                            ElementType sj2 = coef1 * s.elem[idxPivot + mat.cols] + coef2 * s.elem[idxPivot + mat.cols + 1];
                            s.elem[idxPivot + mat.cols] = si2;
                            s.elem[idxPivot + mat.cols + 1] = sj2;

                            flag = true;
                        }

                    }
                    while (flag);

                    // ELIMINATE ROW
                    ElementType coefRow = s.elem[idxPivot + mat.cols] / s.elem[idxPivot];
                    unsigned int idxPRow = i * mat.cols;
                    unsigned int idxERow = idxPRow + mat.cols;
                    for (unsigned int j = 0; j < mat.cols; ++j, ++idxPRow, ++idxERow)
                    {
                        u.elem[idxERow] -= u.elem[idxPRow] * coefRow;
                    }
                    s.elem[idxPivot + mat.cols] = 0;
                    s.elem[idxPivot + mat.cols + 1] -= s.elem[idxPivot + 1] * coefRow;

                    // ELIMINATE COL
                    ElementType coefCol = s.elem[idxPivot + 1] / s.elem[idxPivot];
                    unsigned int idxPCol = i;
                    unsigned int idxECol = idxPCol + 1;
                    for (unsigned int j = 0; j < mat.rows; ++j, idxPCol += mat.cols, idxECol += mat.cols)
                    {
                        v.elem[idxECol] -= v.elem[idxPCol] * coefCol;
                    }
                    s.elem[idxPivot + 1] = 0;
                    s.elem[idxPivot + mat.cols + 1] -= s.elem[idxPivot + mat.cols] * coefCol;

                    corrected = true;
                }

            }
            while (corrected);

            for (unsigned int i = 0, idxPivot = 0; i < mat.rows; ++i, idxPivot += mat.cols + 1)
            {
                if (s.elem[idxPivot] < 0)
                {
                    s.elem[idxPivot] = -s.elem[idxPivot];
                    for (unsigned int j = 0, idxU = i * mat.rows; j < mat.rows; ++j, ++idxU)
                    {
                        u.elem[idxU] = -u.elem[idxU];
                    }
                }
            }

            return SmithNormalForm<ElementType>(s, u, v);
        }

        template<typename ElementType>
        HessenbergForm<ElementType> Algorithms::getHessenbergForm(const Matrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(mat.cols == mat.rows, std::length_error);

            typedef typename ElementTraits<ElementType>::RealType RealType;

            const unsigned int N = mat.cols;
            Matrix<RealType> Q = Matrix<RealType>::identity(N, N);
            Matrix<RealType> QT = Matrix<RealType>::identity(N, N);
            Matrix<RealType> A = mat;

            for (unsigned int col = 0; col < N - 1; ++col)
            {
                Matrix<RealType> transform = getHouseholderMatrix(Traits::getCol(A, col, col + 1, N), N);

                Q = transform * Q;
                A = transform * A * transform.conjugateTranspose();
            }

            return HessenbergForm<ElementType>(Q, A, Q.conjugateTranspose());
        }

        template<typename ElementType>
        SchurForm<ElementType> Algorithms::getSchurForm(const Matrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(mat.cols == mat.rows, std::length_error);

            typedef typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType ComplexRealType;
            typedef typename ElementTraits<ComplexRealType>::AbsType AbsComplexRealType;

            HessenbergForm<ElementType> hf = getHessenbergForm(mat);

            Matrix<ComplexRealType> Q = hf.QT;
            Matrix<ComplexRealType> H = hf.A;
            Matrix<ComplexRealType> I = Matrix<ComplexRealType>::identity(H.getRows(), H.getRows());

            bool flag = true;
            do
            {
                unsigned int shiftRow = H.getRows() - 1;
                for (unsigned int round = 0; flag && round < 100; ++round)
                {
                    ////////////////////////////////
                    // Chose o near an eigenvalue //
                    ////////////////////////////////
                    ComplexRealType h11 = H(shiftRow - 1, shiftRow - 1);
                    ComplexRealType h12 = H(shiftRow - 1, shiftRow);
                    ComplexRealType h21 = H(shiftRow, shiftRow - 1);
                    ComplexRealType h22 = H(shiftRow, shiftRow);

                    ComplexRealType sqrt = ElementTraits<ComplexRealType>::sqrt((h11 + h22) * (h11 + h22) - ComplexRealType(4, 0) * (h11 * h22 - h12 * h21));
                    ComplexRealType o = ElementTraits<ComplexRealType>::div((h11 + h22) + sqrt, ComplexRealType(2, 0));
                    ////////////////////////////////

                    Matrix<ComplexRealType> shift = I * o;
                    QR<ComplexRealType> qr = Algorithms::decomposeQR(H - shift);

                    Q = Q * qr.Q;
                    H = qr.R * qr.Q + shift;

                    bool l = true;
                    for (unsigned int j = 0, idxH = H.cols * H.rows - 2; j < H.cols - 1 && l; ++j, idxH -= H.cols + 1)
                    {
                        l = (ElementTraits<ComplexRealType>::abs(H.elem[idxH]) <= ElementTraits<AbsComplexRealType>::epsilon());
                        shiftRow = H.getRows() - (j + 1);
                    }
                    flag = !l;
                }


                for (unsigned int round = 0; flag && round < 10; ++round)
                {
                    // Try a random shift
                    unsigned int r = rand() % (H.rows - 1);
                    Matrix<ComplexRealType> shift = I * H(r, r);
                    QR<ComplexRealType> qr = Algorithms::decomposeQR(H - shift);

                    Q = Q * qr.Q;
                    H = qr.R * qr.Q + shift;

                    bool l = true;
                    for (unsigned int j = 0, idxH = H.cols * H.rows - 2; j < H.cols - 1 && l; ++j, idxH -= H.cols + 1)
                    {
                        l = (ElementTraits<ComplexRealType>::abs(H.elem[idxH]) <= ElementTraits<AbsComplexRealType>::epsilon());
                        shiftRow = H.getRows() - (j + 1);
                    }
                    flag = !l;
                }
            }
            while (flag);

            return SchurForm<ElementType>(Q, H, Q.conjugateTranspose());
        }

        template<typename ElementType>
        unsigned int Algorithms::getRank(const Matrix<ElementType>& mat)
        {
            typedef typename ElementTraits<ElementType>::RationalType RationalType;
            typedef typename ElementTraits<RationalType>::AbsType AbsType;

            const unsigned int N = mat.rows;
            const unsigned int M = mat.cols;

            Matrix<RationalType> A = mat;
            Matrix<RationalType> P = Matrix<RationalType>::identity(N, N);

            unsigned int rank = 0;
            for (unsigned int i = 0, idxPivotRow = 0; i < N; ++i, idxPivotRow += M)
            {
                bool pivotFound = false;
                AbsType pivotMax = ElementTraits<AbsType>::zero();
                unsigned int pivotRow = i;
                unsigned int pivotCol = i;
                for (unsigned int j = i, idxA = idxPivotRow + i; j < N; ++j, idxA += i)
                {
                    for (unsigned int k = i; k < M; ++k, ++idxA)
                    {
                        AbsType val = ElementTraits<RationalType>::abs(A.elem[idxA]);
                        if (val <= ElementTraits<AbsType>::epsilon())
                        {
                            continue;
                        }

                        if (!pivotFound || pivotMax < val)
                        {
                            pivotFound = true;
                            pivotMax = val;
                            pivotRow = j;
                            pivotCol = k;
                        }
                    }
                }
                if (!pivotFound)
                {
                    break;
                }
                if (pivotRow != i)
                {
                    unsigned int idxP = idxPivotRow + i;
                    unsigned int idxE = pivotRow * M + i;
                    for (unsigned int j = i; j < M; ++j, ++idxP, ++idxE)
                    {
                        std::swap(A.elem[idxP], A.elem[idxE]);
                    }
                }
                if (pivotCol != i)
                {
                    unsigned int idxP = i;
                    unsigned int idxE = pivotCol;
                    for (unsigned int j = 0; j < N; ++j, idxP += M, idxE += M)
                    {
                        std::swap(P.elem[idxP], P.elem[idxE]);
                        std::swap(A.elem[idxP], A.elem[idxE]);
                    }
                }

                for (unsigned int j = i + 1, idxEliminateRow = idxPivotRow + M; j < N; ++j, idxEliminateRow += M)
                {
                    RationalType coef = A.elem[idxEliminateRow + i] / A.elem[idxPivotRow + i];
                    unsigned int idxP = idxPivotRow + i;
                    unsigned int idxE = idxEliminateRow + i;
                    for (unsigned int k = i; k < M; ++k, ++idxP, ++idxE)
                    {
                        A.elem[idxE] -= A.elem[idxP] * coef;
                    }
                }

                ++rank;
            }

            return rank;
        }

        template<typename ElementType>
        std::vector<Vector<typename ElementTraits<ElementType>::RationalType>> Algorithms::solveHomogeneous(const Matrix<ElementType>& mat)
        {
            typedef typename ElementTraits<ElementType>::RationalType RationalType;
            typedef typename ElementTraits<RationalType>::AbsType AbsType;

            const unsigned int N = mat.rows;
            const unsigned int M = mat.cols;

            Matrix<RationalType> A = mat;
            Matrix<RationalType> P = Matrix<RationalType>::identity(M, M);

            unsigned int rank = 0;
            for (unsigned int i = 0, idxPivotRow = 0; i < N; ++i, idxPivotRow += M)
            {
                bool pivotFound = false;
                AbsType pivotMax = ElementTraits<AbsType>::zero();
                unsigned int pivotRow = i;
                unsigned int pivotCol = i;
                for (unsigned int j = i, idxA = idxPivotRow + i; j < N; ++j, idxA += i)
                {
                    for (unsigned int k = i; k < M; ++k, ++idxA)
                    {
                        AbsType val = ElementTraits<RationalType>::abs(A.elem[idxA]);
                        if (val <= ElementTraits<AbsType>::epsilon())
                        {
                            continue;
                        }

                        if (!pivotFound || pivotMax < val)
                        {
                            pivotFound = true;
                            pivotMax = val;
                            pivotRow = j;
                            pivotCol = k;
                        }
                    }
                }
                if (!pivotFound)
                {
                    break;
                }
                if (pivotRow != i)
                {
                    unsigned int idxP = idxPivotRow + i;
                    unsigned int idxE = pivotRow * M + i;
                    for (unsigned int j = i; j < M; ++j, ++idxP, ++idxE)
                    {
                        std::swap(A.elem[idxP], A.elem[idxE]);
                    }
                }
                if (pivotCol != i)
                {
                    unsigned int idxP = i;
                    unsigned int idxE = pivotCol;
                    for (unsigned int j = 0; j < N; ++j, idxP += M, idxE += M)
                    {
                        std::swap(P.elem[idxP], P.elem[idxE]);
                        std::swap(A.elem[idxP], A.elem[idxE]);
                    }
                }

                for (unsigned int j = i + 1, idxEliminateRow = idxPivotRow + M; j < N; ++j, idxEliminateRow += M)
                {
                    RationalType coef = A.elem[idxEliminateRow + i] / A.elem[idxPivotRow + i];
                    unsigned int idxP = idxPivotRow + i;
                    unsigned int idxE = idxEliminateRow + i;
                    for (unsigned int k = i; k < M; ++k, ++idxP, ++idxE)
                    {
                        A.elem[idxE] -= A.elem[idxP] * coef;
                    }
                }

                ++rank;
            }

            for (unsigned int i = rank - 1, idxPivotRow = (rank - 1) * M; i > 0; --i, idxPivotRow -= M)
            {
                RationalType coef = A.elem[idxPivotRow + i];
                for (unsigned int j = i, idxA = idxPivotRow + i; j < M; ++j, ++idxA)
                {
                    A.elem[idxA] /= coef;
                }
                for (unsigned int j = 0, idxEliminateRow = 0; j < i; ++j, idxEliminateRow += M)
                {
                    RationalType coef = A.elem[idxEliminateRow + i];
                    for (unsigned int k = i, idxP = idxPivotRow + i, idxE = idxEliminateRow + i; k < M; ++k, ++idxP, ++idxE)
                    {
                        A.elem[idxE] -= A.elem[idxP] * coef;
                    }
                }
            }
            RationalType coef = A.elem[0];
            for (unsigned int idxA = 0; idxA < M; ++idxA)
            {
                A.elem[idxA] /= coef;
            }

            std::vector<Vector<RationalType>> result(N - rank);
            for (unsigned int i = rank; i < N; ++i)
            {
                Vector<RationalType> eigenVector(M, 00);
                for (unsigned int j = 0, idxA = i; j < rank; ++j, idxA += M)
                {
                    eigenVector.elem[j] = -A.elem[idxA];
                }
                for (unsigned int j = rank; j < M; ++j)
                {
                    eigenVector.elem[j] = ElementTraits<RationalType>::zero();
                }
                eigenVector.elem[i] = ElementTraits<RationalType>::one();

                result[i - rank] = P * eigenVector;
            }

            return result;
        }

        template<typename ElementType>
        std::vector<Vector<typename ElementTraits<ElementType>::RationalType>> Algorithms::getEigenVectors(const Matrix<ElementType>& mat, const ElementType& eigenValue)
        {
            return solveHomogeneous(mat -  Matrix<ElementType>::identity(mat.getRows(), mat.getCols()) * eigenValue);
        }

        template<typename ElementType>
        JordanForm<ElementType> Algorithms::getJordanForm(const Matrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(mat.cols == mat.rows, std::length_error);

            typedef typename ElementTraits<ElementType>::RealType RealType;
            typedef typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType ComplexRealType;
            typedef typename ElementTraits<ComplexRealType>::AbsType AbsComplexRealType;

            const unsigned int N = mat.rows;

            SchurForm<ElementType> schurForm = getSchurForm(mat);

            Matrix<ComplexRealType> M = mat;
            Matrix<ComplexRealType> I = Matrix<ComplexRealType>::identity(N, N);

            Matrix<ComplexRealType> P(N, N);
            Matrix<ComplexRealType> J(N, N);

            unsigned int idxJ = 0;
            unsigned int colP = 0;

            std::vector<ComplexRealType> eigenvalues;
            for (unsigned int i = 0; i < N; ++i)
            {
                eigenvalues.push_back(schurForm.U(i, i));
            }
            std::sort(eigenvalues.begin(), eigenvalues.end(), complex_comparator<RealType>);

            ComplexRealType act = eigenvalues[0];
            ComplexRealType actSum = eigenvalues[0];
            unsigned int mul = 1;
            for (unsigned int i = 1; i <= eigenvalues.size(); ++i)
            {
                if (i == eigenvalues.size() || ElementTraits<ComplexRealType>::abs(eigenvalues[i] - act) > 2 * ElementTraits<AbsComplexRealType>::epsilon())
                {
                    act = ElementTraits<ComplexRealType>::div(actSum, ComplexRealType(mul));

                    Matrix<ComplexRealType> base = M - I * act;
                    Matrix<ComplexRealType> pow = I;

                    std::vector<unsigned int> blocks;
                    std::vector<Matrix<ComplexRealType>> gevBases;
                    std::vector<std::vector<Vector<ComplexRealType>>> gevs;

                    do
                    {
                        pow = base * pow;

                        std::vector<Vector<ComplexRealType>> solutions = solveHomogeneous(pow);
                        blocks.push_back(solutions.size());

                        gevBases.push_back(pow);
                        gevs.push_back(solutions);
                    }
                    while (blocks[blocks.size() - 1] < mul);

                    // Calculate block size -> number of blocks
                    for (unsigned int j = blocks.size() - 1; j > 0; --j)
                    {
                        blocks[j] -= blocks[j - 1];
                    }
                    for (unsigned int j = 0; j < blocks.size() - 1; ++j)
                    {
                        blocks[j] -= blocks[j + 1];
                    }

                    // Fill Jordan form with blocks
                    for (unsigned int j = 0; j < blocks.size(); ++j)
                    {
                        for (unsigned int k = 0; k < blocks[j]; ++k)
                        {
                            J.elem[idxJ] = act;
                            ++idxJ;
                            for (unsigned int l = 1; l < j + 1; ++l)
                            {
                                J.elem[idxJ] = ElementTraits<ComplexRealType>::one();
                                idxJ += N;
                                J.elem[idxJ] = act;
                                ++idxJ;
                            }
                            idxJ += N;
                        }
                    }

                    Matrix<ComplexRealType> T(N, mul);
                    unsigned int colT = mul;
                    unsigned int rankT = 0;
                    for (int j = blocks.size() - 1; j >= 0; --j)
                    {
                        for (unsigned int k = 0; k < blocks[j]; ++k)
                        {
                            unsigned int selected;
                            for (unsigned int l = 0; l < gevs[j].size(); ++l)
                            {
                                // NOT ZERO WITH PREVIOUS
                                if (j > 0 && PNorm<00>::norm(gevBases[j - 1] * gevs[j][l]) <= ElementTraits<typename PNorm<00>::template NormType<ComplexRealType>::Type>::epsilon())
                                {
                                    continue;
                                }

                                --colT;
                                for (unsigned int m = 0, idxT = colT; m < N; ++m, idxT += mul)
                                {
                                    T.elem[idxT] = gevs[j][l].elem[m];
                                }
                                if (getRank(T) > rankT)
                                {
                                    selected = l;
                                    break;
                                }
                                else
                                {
                                    ++colT;
                                }
                            }

                            Vector<ComplexRealType> evp = gevs[j][selected];
                            for (int l = 0; l < j; ++l)
                            {
                                --colT;
                                evp = base * evp;
                                for (unsigned int m = 0, idxT = colT; m < N; ++m, idxT += mul)
                                {
                                    T.elem[idxT] = evp.elem[m];
                                }
                            }

                            rankT = getRank(T);
                        }
                    }

                    for (unsigned int j = 0; j < mul; ++j)
                    {
                        for (unsigned int k = 0; k < N; ++k)
                        {
                            P.set(k, colP + j, T(k, j));
                        }
                    }
                    colP += mul;

                    // NEXT
                    if (i != eigenvalues.size())
                    {
                        act = eigenvalues[i];
                        actSum = eigenvalues[i];
                        mul = 1;
                    }
                }
                else
                {
                    ++mul;
                    actSum += eigenvalues[i];
                }
            }

            return JordanForm<ElementType>(invert(P), J, P);
        }

    }
}
