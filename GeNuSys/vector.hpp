#include <stdexcept>

#include "utils.h"
#include "linalg_operations.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename ElementType>
        Vector<ElementType>::Vector(): length(0), elem(length)
        {
        }

        template<typename ElementType>
        Vector<ElementType>::Vector(unsigned int length, int): length(length), elem(length)
        {
        }

        template<typename ElementType>
        Vector<ElementType>::Vector(unsigned int length): length(length), elem(length, ElementTraits<ElementType>::zero())
        {
        }

        template<typename ElementType>
        Vector<ElementType>::Vector(const Vector<ElementType>& vct): length(vct.length), elem(vct.elem)
        {
        }

        template<typename ElementType>
        template<typename SourceType>
        Vector<ElementType>::Vector(const Vector<SourceType>& vct): length(vct.length), elem(vct.length)
        {
            for (unsigned int i = 0; i < length; ++i)
            {
                elem[i] = ElementTraits<SourceType>::template asType<ElementType>(vct.elem[i]);
            }
        }

        template<typename ElementType>
        template<typename SourceType>
        Vector<ElementType>::Vector(const SparseVector<SourceType>& vct): length(vct.length), elem(vct.length, ElementTraits<SourceType>::zero())
        {
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                elem[vct.elem[i].idx] = ElementTraits<SourceType>::template asType<ElementType>(vct.elem[i].value);
            }
        }

        template<typename ElementType>
        Vector<ElementType>::~Vector()
        {
        }

        template<typename ElementType>
        Vector<ElementType>& Vector<ElementType>::operator =(const Vector<ElementType>& vct)
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
        Vector<ElementType>& Vector<ElementType>::operator =(const Vector<SourceType>& vct)
        {
            length = vct.length;
            elem = std::vector<ElementType>(vct.length);
            for (unsigned int i = 0; i < length; ++i)
            {
                elem[i] = ElementTraits<SourceType>::template asType<ElementType>(vct.elem[i]);
            }

            return *this;
        }

        template<typename ElementType>
        template<typename SourceType>
        Vector<ElementType>& Vector<ElementType>::operator =(const SparseVector<SourceType>& vct)
        {
            length = vct.length;
            elem = std::vector<ElementType>(vct.length);
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                elem[vct.elem[i].idx] = ElementTraits<SourceType>::template asType<ElementType>(vct.elem[i].value);
            }

            return *this;
        }

        template<typename ElementType>
        unsigned int Vector<ElementType>::getLength() const
        {
            return length;
        }

        template<typename ElementType>
        const ElementType& Vector<ElementType>::operator [](unsigned int idx) const
        {
            ASSERT_EXCEPTION(idx < length, std::out_of_range);

            return elem[idx];
        }

        template<typename ElementType>
        void Vector<ElementType>::set(unsigned int idx, const ElementType& value)
        {
            ASSERT_EXCEPTION(idx < length, std::out_of_range);

            elem[idx] = value;
        }

        template<typename ElementType>
        Vector<ElementType> Vector<ElementType>::operator *(const ElementType& value) const
        {
            Vector<ElementType> result(length, 00);
            Operations::vct_mul(*this, value, result);

            return result;
        }

        template<typename ElementType>
        Vector<typename ElementTraits<ElementType>::RationalType> Vector<ElementType>::operator /(const ElementType& value) const
        {
            Vector<typename ElementTraits<ElementType>::RationalType> result(length, 00);
            Operations::vct_div(*this, value, result);

            return result;
        }

        template<typename ElementType>
        ElementType Vector<ElementType>::operator *(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_mul(*this, vct);
        }

        template<typename ElementType>
        Vector<ElementType> Vector<ElementType>::operator +(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            Vector<ElementType> result(length, 00);
            Operations::vct_add(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        Vector<ElementType> Vector<ElementType>::operator -(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            Vector<ElementType> result(length, 00);
            Operations::vct_sub(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        bool Vector<ElementType>::operator ==(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_eq(*this, vct);
        }

        template<typename ElementType>
        bool Vector<ElementType>::operator !=(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return !Operations::vct_eq(*this, vct);
        }

        template<typename ElementType>
        bool Vector<ElementType>::operator <(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_lt(*this, vct);
        }

        template<typename ElementType>
        bool Vector<ElementType>::operator <=(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_leq(*this, vct);
        }

        template<typename ElementType>
        bool Vector<ElementType>::operator >(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_gt(*this, vct);
        }

        template<typename ElementType>
        bool Vector<ElementType>::operator >=(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_geq(*this, vct);
        }

        template<typename ElementType>
        ElementType Vector<ElementType>::operator *(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_mul(*this, vct);
        }

        template<typename ElementType>
        Vector<ElementType> Vector<ElementType>::operator +(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            Vector<ElementType> result(length, 00);
            Operations::vct_add(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        Vector<ElementType> Vector<ElementType>::operator -(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            Vector<ElementType> result(length, 00);
            Operations::vct_sub(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        bool Vector<ElementType>::operator ==(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_eq(*this, vct);
        }

        template<typename ElementType>
        bool Vector<ElementType>::operator !=(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return !Operations::vct_eq(*this, vct);
        }

        template<typename ElementType>
        bool Vector<ElementType>::operator <(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_lt(*this, vct);
        }

        template<typename ElementType>
        bool Vector<ElementType>::operator <=(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_leq(*this, vct);
        }

        template<typename ElementType>
        bool Vector<ElementType>::operator >(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_gt(*this, vct);
        }

        template<typename ElementType>
        bool Vector<ElementType>::operator >=(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(length == vct.length, std::length_error);

            return Operations::vct_geq(*this, vct);
        }

    }
}
