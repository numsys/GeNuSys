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

#include <stdexcept>

#include "linalg_operations.h"
#include "utils.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename ElementType>
        SparseVector<ElementType>::Entry::Entry()
        {
        }

        template<typename ElementType>
        SparseVector<ElementType>::Entry::Entry(unsigned int idx, const ElementType& value): idx(idx), value(value)
        {
        }

        template<typename ElementType>
        SparseVector<ElementType>::SparseVector(): length(0), elem()
        {
        }

        template<typename ElementType>
        SparseVector<ElementType>::SparseVector(unsigned int length, unsigned int size): length(length), elem(size)
        {
        }

        template<typename ElementType>
        SparseVector<ElementType>::SparseVector(unsigned int length): length(length), elem()
        {
        }

        template<typename ElementType>
        SparseVector<ElementType>::SparseVector(const SparseVector<ElementType>& vct): length(vct.length), elem(vct.elem)
        {
        }

        template<typename ElementType>
        template<typename SourceType>
        SparseVector<ElementType>::SparseVector(const SparseVector<SourceType>& vct): length(vct.length), elem(vct.elem.size())
        {
            for (unsigned int i = 0; i < elem.size(); ++i)
            {
                elem[i] = Entry(i, ElementTraits<SourceType>::template asType<ElementType>(vct.elem[i].value));
            }
        }

        template<typename ElementType>
        template<typename SourceType>
        SparseVector<ElementType>::SparseVector(const Vector<SourceType>& vct): length(vct.length)
        {
            unsigned int nnz = 0;
            for (unsigned int i = 0; i < length; ++i)
            {
                if (vct.elem[i] != ElementTraits<SourceType>::zero())
                {
                    ++nnz;
                }
            }

            elem = std::vector<Entry>(nnz);
            unsigned int size = 0;
            for (unsigned int i = 0; i < length; ++i)
            {
                if (vct.elem[i] != ElementTraits<SourceType>::zero())
                {
                    elem[size++] = Entry(i, ElementTraits<SourceType>::template asType<ElementType>(vct.elem[i]));
                }
            }
        }

        template<typename ElementType>
        SparseVector<ElementType>::~SparseVector()
        {
        }

        template<typename ElementType>
        SparseVector<ElementType>& SparseVector<ElementType>::operator =(const SparseVector<ElementType>& vct)
        {
            if (this != &vct)
            {
                length = vct.length;
                elem = vct.elem;
            }

            return *this;
        }

        template<typename ElementType>
        template<typename SourceType>
        SparseVector<ElementType>& SparseVector<ElementType>::operator =(const SparseVector<SourceType>& vct)
        {
            length = vct.length;
            elem = std::vector<Entry>(vct.elem.size());
            for (unsigned int i = 0; i < elem.size(); ++i)
            {
                elem[i] = Entry(i, ElementTraits<SourceType>::template asType<ElementType>(vct.elem[i].value));
            }

            return *this;
        }

        template<typename ElementType>
        template<typename SourceType>
        SparseVector<ElementType>& SparseVector<ElementType>::operator =(const Vector<SourceType>& vct)
        {
            length = vct.length;

            unsigned int nnz = 0;
            for (unsigned int i = 0; i < length; ++i)
            {
                if (vct.elem[i] != ElementTraits<SourceType>::zero())
                {
                    ++nnz;
                }
            }

            elem = std::vector<Entry>(nnz);
            unsigned int size = 0;
            for (unsigned int i = 0; i < length; ++i)
            {
                if (vct.elem[i] != ElementTraits<SourceType>::zero())
                {
                    elem[size++] = Entry(i, ElementTraits<SourceType>::template asType<ElementType>(vct.elem[i]));
                }
            }

            return *this;
        }

        template<typename ElementType>
        unsigned int SparseVector<ElementType>::size() const
        {
            return elem.size();
        }

        template<typename ElementType>
        unsigned int SparseVector<ElementType>::search(unsigned int idx) const
        {
            unsigned int min = 0;
            unsigned int max = size();
            unsigned int mid;
            while (min < max)
            {
                mid = (min + max) / 2;
                if (elem[mid].idx < idx)
                {
                    min = mid + 1;
                }
                else
                {
                    max = mid;
                }
            }
            return min;
        }

        template<typename ElementType>
        void SparseVector<ElementType>::push(unsigned int idx, const ElementType& value)
        {
            elem.push_back(Entry(idx, value));
        }

        template<typename ElementType>
        unsigned int SparseVector<ElementType>::getLength() const
        {
            return length;
        }

        template<typename ElementType>
        void SparseVector<ElementType>::set(unsigned int idx, const ElementType& value)
        {
            ASSERT_EXCEPTION(idx < length, std::out_of_range);

            unsigned int entry_idx = search(idx);
            if (value != ElementTraits<ElementType>::zero())
            {
                if (entry_idx < size() && elem[entry_idx].idx == idx)
                {
                    elem[entry_idx].value = value;
                }
                else
                {
                    elem.insert(elem.begin() + entry_idx, Entry(idx, value));
                }
            }
            else
            {
                if (entry_idx < size() && elem[entry_idx].idx == idx)
                {
                    elem.erase(elem.begin() + entry_idx);
                }
            }
        }

        template<typename ElementType>
        const ElementType& SparseVector<ElementType>::operator [](unsigned int idx) const
        {
            ASSERT_EXCEPTION(idx < length, std::out_of_range);

            unsigned int entry_idx = search(idx);
            if (entry_idx < size() && elem[entry_idx].idx == idx)
            {
                return elem[entry_idx].value;
            }
            else
            {
                return ElementTraits<ElementType>::zero();
            }
        }

        template<typename ElementType>
        SparseVector<ElementType> SparseVector<ElementType>::operator *(const ElementType& value) const
        {
            SparseVector<ElementType> result(length, size());
            Operations::vct_mul(*this, value, result);

            return result;
        }

        template<typename ElementType>
        SparseVector<typename ElementTraits<ElementType>::RationalType> SparseVector<ElementType>::operator /(const ElementType& value) const
        {
            SparseVector<typename ElementTraits<ElementType>::RationalType> result(length, size());
            Operations::vct_div(*this, value, result);

            return result;
        }

        template<typename ElementType>
        ElementType SparseVector<ElementType>::operator *(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_mul(*this, vct);
        }

        template<typename ElementType>
        Vector<ElementType> SparseVector<ElementType>::operator +(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            Vector<ElementType> result(length, 00);
            Operations::vct_add(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        Vector<ElementType> SparseVector<ElementType>::operator -(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            Vector<ElementType> result(length, 00);
            Operations::vct_sub(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        bool SparseVector<ElementType>::operator ==(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_eq(*this, vct);
        }

        template<typename ElementType>
        bool SparseVector<ElementType>::operator !=(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return !Operations::vct_eq(*this, vct);
        }

        template<typename ElementType>
        bool SparseVector<ElementType>::operator <(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_lt(*this, vct);
        }

        template<typename ElementType>
        bool SparseVector<ElementType>::operator <=(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_leq(*this, vct);
        }

        template<typename ElementType>
        bool SparseVector<ElementType>::operator >(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_gt(*this, vct);
        }

        template<typename ElementType>
        bool SparseVector<ElementType>::operator >=(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_geq(*this, vct);
        }

        template<typename ElementType>
        ElementType SparseVector<ElementType>::operator *(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_mul(*this, vct);
        }

        template<typename ElementType>
        SparseVector<ElementType> SparseVector<ElementType>::operator +(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            SparseVector<ElementType> result(length);
            Operations::vct_add(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        SparseVector<ElementType> SparseVector<ElementType>::operator -(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            SparseVector<ElementType> result(length);
            Operations::vct_sub(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        bool SparseVector<ElementType>::operator ==(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_eq(*this, vct);
        }

        template<typename ElementType>
        bool SparseVector<ElementType>::operator !=(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return !Operations::vct_eq(*this, vct);
        }

        template<typename ElementType>
        bool SparseVector<ElementType>::operator <(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_lt(*this, vct);
        }

        template<typename ElementType>
        bool SparseVector<ElementType>::operator <=(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_leq(*this, vct);
        }

        template<typename ElementType>
        bool SparseVector<ElementType>::operator >(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_gt(*this, vct);
        }

        template<typename ElementType>
        bool SparseVector<ElementType>::operator >=(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_geq(*this, vct);
        }

    }
}
