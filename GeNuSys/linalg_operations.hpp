#include <stdexcept>
#include <algorithm>
#include <vector>

#include "utils.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename ElementType>
        void Operations::vct_mul(const Vector<ElementType>& vct, const ElementType& value, Vector<ElementType>& result)
        {
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                result.elem[i] = vct.elem[i] * value;
            }
        }

        template<typename ElementType>
        void Operations::vct_mul(Vector<ElementType>& vct, const ElementType& value)
        {
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                vct.elem[i] *= value;
            }
        }

        template<typename ElementType>
        void Operations::vct_mul(const SparseVector<ElementType>& vct, const ElementType& value, SparseVector<ElementType>& result)
        {
            result.elem.resize(vct.size());

            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                result.elem[i] = typename SparseVector<ElementType>::Entry(vct.elem[i].idx, vct.elem[i].value * value);
            }
        }

        template<typename ElementType>
        void Operations::vct_mul(SparseVector<ElementType>& vct, const ElementType& value)
        {
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                vct.elem[i].value *= value;
            }
        }

        template<typename ElementType>
        void Operations::vct_idiv(const Vector<ElementType>& vct, const ElementType& value, Vector<ElementType>& result)
        {
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                result.elem[i] = ElementTraits<ElementType>::idiv(vct.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::vct_idiv(Vector<ElementType>& vct, const ElementType& value)
        {
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                vct.elem[i] = ElementTraits<ElementType>::idiv(vct.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::vct_idiv(const SparseVector<ElementType>& vct, const ElementType& value, SparseVector<ElementType>& result)
        {
            result.elem.resize(vct.size());

            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                result.elem[i] = typename SparseVector<ElementType>::Entry(vct.elem[i].idx, ElementTraits<ElementType>::idiv(vct.elem[i].value, value));
            }
        }

        template<typename ElementType>
        void Operations::vct_idiv(SparseVector<ElementType>& vct, const ElementType& value)
        {
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                vct.elem[i].value = ElementTraits<ElementType>::idiv(vct.elem[i].value, value);
            }
        }

        template<typename ElementType>
        void Operations::vct_mod(const Vector<ElementType>& vct, const ElementType& value, Vector<ElementType>& result)
        {
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                result.elem[i] = ElementTraits<ElementType>::mod(vct.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::vct_mod(Vector<ElementType>& vct, const ElementType& value)
        {
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                vct.elem[i] = ElementTraits<ElementType>::mod(vct.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::vct_mod(const SparseVector<ElementType>& vct, const ElementType& value, SparseVector<ElementType>& result)
        {
            result.elem.resize(vct.size());

            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                result.elem[i] = typename SparseVector<ElementType>::Entry(vct.elem[i].idx, ElementTraits<ElementType>::mod(vct.elem[i].value, value));
            }
        }

        template<typename ElementType>
        void Operations::vct_mod(SparseVector<ElementType>& vct, const ElementType& value)
        {
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                vct.elem[i] = ElementTraits<ElementType>::mod(vct.elem[i].value, value);
            }
        }

        template<typename ElementType>
        void Operations::vct_mods(const Vector<ElementType>& vct, const ElementType& value, Vector<ElementType>& result)
        {
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                result.elem[i] = ElementTraits<ElementType>::mods(vct.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::vct_mods(Vector<ElementType>& vct, const ElementType& value)
        {
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                vct.elem[i] = ElementTraits<ElementType>::mods(vct.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::vct_mods(const SparseVector<ElementType>& vct, const ElementType& value, SparseVector<ElementType>& result)
        {
            result.elem.resize(vct.size());

            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                result.elem[i] = typename SparseVector<ElementType>::Entry(vct.elem[i].idx, ElementTraits<ElementType>::mods(vct.elem[i].value, value));
            }
        }

        template<typename ElementType>
        void Operations::vct_mods(SparseVector<ElementType>& vct, const ElementType& value)
        {
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                vct.elem[i].value = ElementTraits<ElementType>::mods(vct.elem[i].value, value);
            }
        }

        template<typename ElementType>
        void Operations::vct_div(const Vector<ElementType>& vct, const ElementType& value, Vector<typename ElementTraits<ElementType>::RationalType>& result)
        {
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                result.elem[i] = ElementTraits<ElementType>::div(vct.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::vct_div(Vector<ElementType>& vct, const ElementType& value)
        {
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                vct.elem[i] = ElementTraits<ElementType>::div(vct.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::vct_div(const SparseVector<ElementType>& vct, const ElementType& value, SparseVector<typename ElementTraits<ElementType>::RationalType>& result)
        {
            result.elem.resize(vct.size());

            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                result.elem[i] = typename SparseVector<typename ElementTraits<ElementType>::RationalType>::Entry(vct.elem[i].idx, ElementTraits<ElementType>::div(vct.elem[i].value, value));
            }
        }

        template<typename ElementType>
        void Operations::vct_div(SparseVector<ElementType>& vct, const ElementType& value)
        {
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                vct.elem[i].value = ElementTraits<ElementType>::div(vct.elem[i].value, value);
            }
        }

        template<typename ElementType>
        void Operations::vct_add(const Vector<ElementType>& op1, const Vector<ElementType>& op2, Vector<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            for (unsigned int i = 0; i < op1.length; ++i)
            {
                result.elem[i] = op1.elem[i] + op2.elem[i];
            }
        }

        template<typename ElementType>
        void Operations::vct_add(Vector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            for (unsigned int i = 0; i < op1.length; ++i)
            {
                op1.elem[i] += op2.elem[i];
            }
        }

        template<typename ElementType>
        void Operations::vct_add(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2, Vector<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            for (unsigned int i = 0; i < op1.length; ++i)
            {
                result.elem[i] = op1.elem[i];
            }
            for (unsigned int i = 0; i < op2.size(); ++i)
            {
                result.elem[op2.elem[i].idx] += op2.elem[i].value;
            }
        }

        template<typename ElementType>
        void Operations::vct_add(Vector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            for (unsigned int i = 0; i < op2.size(); ++i)
            {
                op1.elem[op2.elem[i].idx] += op2.elem[i].value;
            }
        }

        template<typename ElementType>
        void Operations::vct_add(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2, Vector<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            for (unsigned int i = 0; i < op2.length; ++i)
            {
                result.elem[i] = op2.elem[i];
            }
            for (unsigned int i = 0; i < op1.size(); ++i)
            {
                result.elem[op1.elem[i].idx] += op1.elem[i].value;
            }
        }

        template<typename ElementType>
        void Operations::vct_add(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2, SparseVector<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            result.elem.clear();

            unsigned int i = 0, j = 0;
            while (i < op1.size() && j < op2.size())
            {
                if (op1.elem[i].idx < op2.elem[j].idx)
                {
                    result.push(op1.elem[i].idx, op1.elem[i].value);
                    ++i;
                }
                else if (op1.elem[i].idx > op2.elem[j].idx)
                {
                    result.push(op2.elem[j].idx, op2.elem[j].value);
                    ++j;
                }
                else
                {
                    result.push(op1.elem[i].idx, op1.elem[i].value + op2.elem[j].value);
                    ++i;
                    ++j;
                }
            }
            while (i < op1.size())
            {
                result.push(op1.elem[i].idx, op1.elem[i].value);
                ++i;
            }
            while (j < op2.size())
            {
                result.push(op2.elem[j].idx, op2.elem[j].value);
                ++j;
            }
        }

        template<typename ElementType>
        void Operations::vct_sub(const Vector<ElementType>& op1, const Vector<ElementType>& op2, Vector<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            for (unsigned int i = 0; i < op1.length; ++i)
            {
                result.elem[i] = op1.elem[i] - op2.elem[i];
            }
        }

        template<typename ElementType>
        void Operations::vct_sub(Vector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            for (unsigned int i = 0; i < op1.length; ++i)
            {
                op1.elem[i] -= op2.elem[i];
            }
        }

        template<typename ElementType>
        void Operations::vct_sub(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2, Vector<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            for (unsigned int i = 0; i < op1.length; ++i)
            {
                result.elem[i] = op1.elem[i];
            }
            for (unsigned int i = 0; i < op2.size(); ++i)
            {
                result.elem[op2.elem[i].idx] -= op2.elem[i].value;
            }
        }

        template<typename ElementType>
        void Operations::vct_sub(Vector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            for (unsigned int i = 0; i < op2.size(); ++i)
            {
                op1.elem[op2.elem[i].idx] -= op2.elem[i].value;
            }
        }

        template<typename ElementType>
        void Operations::vct_sub(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2, Vector<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            for (unsigned int i = 0; i < op2.length; ++i)
            {
                result.elem[i] = -op2.elem[i];
            }
            for (unsigned int i = 0; i < op1.size(); ++i)
            {
                result.elem[op1.elem[i].idx] += op1.elem[i].value;
            }
        }

        template<typename ElementType>
        void Operations::vct_sub(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2, SparseVector<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            result.elem.clear();

            unsigned int i = 0, j = 0;
            while (i < op1.size() && j < op2.size())
            {
                if (op1.elem[i].idx < op2.elem[j].idx)
                {
                    result.push(op1.elem[i].idx, op1.elem[i].value);
                    ++i;
                }
                else if (op1.elem[i].idx > op2.elem[j].idx)
                {
                    result.push(op2.elem[j].idx, -op2.elem[j].value);
                    ++j;
                }
                else
                {
                    result.push(op1.elem[i].idx, op1.elem[i].value - op2.elem[j].value);
                    ++i;
                    ++j;
                }
            }
            while (i < op1.size())
            {
                result.push(op1.elem[i].idx, op1.elem[i].value);
                ++i;
            }
            while (j < op2.size())
            {
                result.push(op2.elem[j].idx, -op2.elem[j].value);
                ++j;
            }
        }

        template<typename ElementType>
        ElementType Operations::vct_mul(const Vector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            ElementType result = ElementTraits<ElementType>::zero();
            for (unsigned int i = 0; i < op1.length; ++i)
            {
                result += op1.elem[i] * op2.elem[i];
            }

            return result;
        }

        template<typename ElementType>
        ElementType Operations::vct_mul(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            ElementType result = ElementTraits<ElementType>::zero();
            for (unsigned int i = 0; i < op2.size(); ++i)
            {
                result += op1.elem[op2.elem[i].idx] * op2.elem[i].value;
            }

            return result;
        }

        template<typename ElementType>
        ElementType Operations::vct_mul(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            ElementType result = ElementTraits<ElementType>::zero();
            for (unsigned int i = 0; i < op1.size(); ++i)
            {
                result += op2.elem[op1.elem[i].idx] * op1.elem[i].value;
            }

            return result;
        }

        template<typename ElementType>
        ElementType Operations::vct_mul(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            ElementType result = ElementTraits<ElementType>::zero();
            unsigned int i = 0, j = 0;
            while (i < op1.size() && j < op2.size())
            {
                if (op1.elem[i].idx < op2.elem[j].idx)
                {
                    ++i;
                }
                else if (op1.elem[i].idx > op2.elem[j].idx)
                {
                    ++j;
                }
                else
                {
                    result += op1.elem[i].value * op2.elem[j].value;
                    ++i;
                    ++j;
                }
            }

            return result;
        }

        template<typename ElementType>
        bool Operations::vct_eq(const Vector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            for (unsigned int i = 0; i < op1.length && l; ++i)
            {
                l = (op1.elem[i] == op2.elem[i]);
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_eq(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.length && j < op2.size() && l)
            {
                if (i < op2.elem[j].idx)
                {
                    l = (op1.elem[i] == ElementTraits<ElementType>::zero());
                    ++i;
                }
                else
                {
                    l = (op1.elem[i] == op2.elem[j].value);
                    ++i;
                    ++j;
                }
            }
            while (i < op1.length && l)
            {
                l = (op1.elem[i] == ElementTraits<ElementType>::zero());
                ++i;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_eq(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.size() && j < op2.length && l)
            {
                if (j < op1.elem[i].idx)
                {
                    l = (ElementTraits<ElementType>::zero() == op2.elem[j]);
                    ++j;
                }
                else
                {
                    l = (op1.elem[i].value == op2.elem[j]);
                    ++i;
                    ++j;
                }
            }
            while (j < op2.length && l)
            {
                l = (ElementTraits<ElementType>::zero() == op2.elem[j]);
                ++j;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_eq(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = (op1.size() == op2.size());
            for (unsigned int i = 0; i < op1.size() && l; ++i)
            {
                l = (op1.elem[i].idx == op2.elem[i].idx && op1.elem[i].value == op2.elem[i].value);
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_gt(const Vector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            for (unsigned int i = 0; i < op1.length && l; ++i)
            {
                l = (op1.elem[i] > op2.elem[i]);
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_gt(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.length && j < op2.size() && l)
            {
                if (i < op2.elem[j].idx)
                {
                    l = (op1.elem[i] > ElementTraits<ElementType>::zero());
                    ++i;
                }
                else
                {
                    l = (op1.elem[i] > op2.elem[j].value);
                    ++i;
                    ++j;
                }
            }
            while (i < op1.length && l)
            {
                l = (op1.elem[i] > ElementTraits<ElementType>::zero());
                ++i;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_gt(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.size() && j < op2.length && l)
            {
                if (j < op1.elem[i].idx)
                {
                    l = (ElementTraits<ElementType>::zero() > op2.elem[j]);
                    ++j;
                }
                else
                {
                    l = (op1.elem[i].value > op2.elem[j]);
                    ++i;
                    ++j;
                }
            }
            while (j < op2.length && l)
            {
                l = (ElementTraits<ElementType>::zero() > op2.elem[j]);
                ++j;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_gt(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.size() && j < op2.size() && l)
            {
                if (op1.elem[i].idx < op2.elem[j].idx)
                {
                    l = (op1.elem[i].value > ElementTraits<ElementType>::zero());
                    ++i;
                }
                else if (op1.elem[i].idx > op2.elem[j].idx)
                {
                    l = (ElementTraits<ElementType>::zero() > op2.elem[j].value);
                    ++j;
                }
                else
                {
                    l = (op1.elem[i].value > op2.elem[j].value);
                    ++i;
                    ++j;
                }
            }
            while (i < op1.size() && l)
            {
                l = (op1.elem[i].value > ElementTraits<ElementType>::zero());
                ++i;
            }
            while (j < op2.size() && l)
            {
                l = (ElementTraits<ElementType>::zero() > op2.elem[j].value);
                ++j;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_geq(const Vector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            for (unsigned int i = 0; i < op1.length && l; ++i)
            {
                l = (op1.elem[i] >= op2.elem[i]);
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_geq(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.length && j < op2.size() && l)
            {
                if (i < op2.elem[j].idx)
                {
                    l = (op1.elem[i] >= ElementTraits<ElementType>::zero());
                    ++i;
                }
                else
                {
                    l = (op1.elem[i] >= op2.elem[j].value);
                    ++i;
                    ++j;
                }
            }
            while (i < op1.length && l)
            {
                l = (op1.elem[i] >= ElementTraits<ElementType>::zero());
                ++i;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_geq(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.size() && j < op2.length && l)
            {
                if (j < op1.elem[i].idx)
                {
                    l = (ElementTraits<ElementType>::zero() >= op2.elem[j]);
                    ++j;
                }
                else
                {
                    l = (op1.elem[i].value >= op2.elem[j]);
                    ++i;
                    ++j;
                }
            }
            while (j < op2.length && l)
            {
                l = (ElementTraits<ElementType>::zero() >= op2.elem[j]);
                ++j;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_geq(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.size() && j < op2.size() && l)
            {
                if (op1.elem[i].idx < op2.elem[j].idx)
                {
                    l = (op1.elem[i].value >= ElementTraits<ElementType>::zero());
                    ++i;
                }
                else if (op1.elem[i].idx > op2.elem[j].idx)
                {
                    l = (ElementTraits<ElementType>::zero() >= op2.elem[j].value);
                    ++j;
                }
                else
                {
                    l = (op1.elem[i].value >= op2.elem[j].value);
                    ++i;
                    ++j;
                }
            }
            while (i < op1.size() && l)
            {
                l = (op1.elem[i].value >= ElementTraits<ElementType>::zero());
                ++i;
            }
            while (j < op2.size() && l)
            {
                l = (ElementTraits<ElementType>::zero() >= op2.elem[j].value);
                ++j;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_lt(const Vector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            for (unsigned int i = 0; i < op1.length && l; ++i)
            {
                l = (op1.elem[i] < op2.elem[i]);
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_lt(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.length && j < op2.size() && l)
            {
                if (i < op2.elem[j].idx)
                {
                    l = (op1.elem[i] < ElementTraits<ElementType>::zero());
                    ++i;
                }
                else
                {
                    l = (op1.elem[i] < op2.elem[j].value);
                    ++i;
                    ++j;
                }
            }
            while (i < op1.length && l)
            {
                l = (op1.elem[i] < ElementTraits<ElementType>::zero());
                ++i;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_lt(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.size() && j < op2.length && l)
            {
                if (j < op1.elem[i].idx)
                {
                    l = (ElementTraits<ElementType>::zero() < op2.elem[j]);
                    ++j;
                }
                else
                {
                    l = (op1.elem[i].value < op2.elem[j]);
                    ++i;
                    ++j;
                }
            }
            while (j < op2.length && l)
            {
                l = (ElementTraits<ElementType>::zero() < op2.elem[j]);
                ++j;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_lt(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.size() && j < op2.size() && l)
            {
                if (op1.elem[i].idx < op2.elem[j].idx)
                {
                    l = (op1.elem[i].value < ElementTraits<ElementType>::zero());
                    ++i;
                }
                else if (op1.elem[i].idx > op2.elem[j].idx)
                {
                    l = (ElementTraits<ElementType>::zero() < op2.elem[j].value);
                    ++j;
                }
                else
                {
                    l = (op1.elem[i].value < op2.elem[j].value);
                    ++i;
                    ++j;
                }
            }
            while (i < op1.size() && l)
            {
                l = (op1.elem[i].value < ElementTraits<ElementType>::zero());
                ++i;
            }
            while (j < op2.size() && l)
            {
                l = (ElementTraits<ElementType>::zero() < op2.elem[j].value);
                ++j;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_leq(const Vector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            for (unsigned int i = 0; i < op1.length && l; ++i)
            {
                l = (op1.elem[i] <= op2.elem[i]);
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_leq(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.length && j < op2.size() && l)
            {
                if (i < op2.elem[j].idx)
                {
                    l = (op1.elem[i] <= ElementTraits<ElementType>::zero());
                    ++i;
                }
                else
                {
                    l = (op1.elem[i] <= op2.elem[j].value);
                    ++i;
                    ++j;
                }
            }
            while (i < op1.length && l)
            {
                l = (op1.elem[i] <= ElementTraits<ElementType>::zero());
                ++i;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_leq(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.size() && j < op2.length && l)
            {
                if (j < op1.elem[i].idx)
                {
                    l = (ElementTraits<ElementType>::zero() <= op2.elem[j]);
                    ++j;
                }
                else
                {
                    l = (op1.elem[i].value <= op2.elem[j]);
                    ++i;
                    ++j;
                }
            }
            while (j < op2.length && l)
            {
                l = (ElementTraits<ElementType>::zero() <= op2.elem[j]);
                ++j;
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::vct_leq(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.length == op2.length, std::length_error);

            bool l = true;
            unsigned int i = 0, j = 0;
            while (i < op1.size() && j < op2.size() && l)
            {
                if (op1.elem[i].idx < op2.elem[j].idx)
                {
                    l = (op1.elem[i].value <= ElementTraits<ElementType>::zero());
                    ++i;
                }
                else if (op1.elem[i].idx > op2.elem[j].idx)
                {
                    l = (ElementTraits<ElementType>::zero() <= op2.elem[j].value);
                    ++j;
                }
                else
                {
                    l = (op1.elem[i].value <= op2.elem[j].value);
                    ++i;
                    ++j;
                }
            }
            while (i < op1.size() && l)
            {
                l = (op1.elem[i].value <= ElementTraits<ElementType>::zero());
                ++i;
            }
            while (j < op2.size() && l)
            {
                l = (ElementTraits<ElementType>::zero() <= op2.elem[j].value);
                ++j;
            }

            return l;
        }

        template<typename ElementType>
        void Operations::mat_transpose(const Matrix<ElementType>& mat, Matrix<ElementType>& result)
        {
            for (unsigned int i = 0, idxA = 0; i < mat.rows; ++i)
            {
                for (unsigned int j = 0, idxR = i; j < mat.cols; ++j, ++idxA, idxR += mat.rows)
                {
                    result.elem[idxR] = mat.elem[idxA];
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_transpose(const SparseMatrix<ElementType>& mat, SparseMatrix<ElementType>& result)
        {
            result.elem.resize(mat.size());
            std::fill(result.row_ptr.begin(), result.row_ptr.end(), ElementTraits<ElementType>::zero());

            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                ++result.row_ptr[mat.elem[i].col_idx + 1];
            }
            for (unsigned int i = 1; i <= result.rows; ++i)
            {
                result.row_ptr[i] += result.row_ptr[i - 1];
            }
            std::vector<unsigned int> inserted(result.rows, 0);
            for (unsigned int i = 0; i < mat.rows; ++i)
            {
                for (unsigned int j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; ++j)
                {
                    result.elem[result.row_ptr[mat.elem[j].col_idx] + inserted[mat.elem[j].col_idx]++] =
                        typename SparseMatrix<ElementType>::Entry(i, mat.elem[j].value);
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_conjugate_transpose(const Matrix<ElementType>& mat, Matrix<ElementType>& result)
        {
            for (unsigned int i = 0, idxA = 0; i < mat.rows; ++i)
            {
                for (unsigned int j = 0, idxR = i; j < mat.cols; ++j, ++idxA, idxR += mat.rows)
                {
                    result.elem[idxR] = ElementTraits<ElementType>::conj(mat.elem[idxA]);
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_conjugate_transpose(const SparseMatrix<ElementType>& mat, SparseMatrix<ElementType>& result)
        {
            result.elem.resize(mat.size());
            std::fill(result.row_ptr.begin(), result.row_ptr.end(), ElementTraits<ElementType>::zero());

            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                ++result.row_ptr[mat.elem[i].col_idx + 1];
            }
            for (unsigned int i = 1; i <= result.rows; ++i)
            {
                result.row_ptr[i] += result.row_ptr[i - 1];
            }
            std::vector<unsigned int> inserted(result.rows, 0);
            for (unsigned int i = 0; i < mat.rows; ++i)
            {
                for (unsigned int j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; ++j)
                {
                    result.elem[result.row_ptr[mat.elem[j].col_idx] + inserted[mat.elem[j].col_idx]++] =
                        typename SparseMatrix<ElementType>::Entry(i, ElementTraits<ElementType>::conj(mat.elem[j].value));
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const Matrix<ElementType>& mat, const ElementType& value, Matrix<ElementType>& result)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                result.elem[i] = mat.elem[i] * value;
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(Matrix<ElementType>& mat, const ElementType& value)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                mat.elem[i] *= value;
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<ElementType>& result)
        {
            result.elem.resize(mat.size());

            result.row_ptr = mat.row_ptr;
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                result.elem[i] = typename SparseMatrix<ElementType>::Entry(mat.elem[i].col_idx, mat.elem[i].value * value);
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(SparseMatrix<ElementType>& mat, const ElementType& value)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                mat.elem[i].value *= value;
            }
        }

        template<typename ElementType>
        void Operations::mat_idiv(const Matrix<ElementType>& mat, const ElementType& value, Matrix<ElementType>& result)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                result.elem[i] = ElementTraits<ElementType>::idiv(mat.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::mat_idiv(Matrix<ElementType>& mat, const ElementType& value)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                mat.elem[i] = ElementTraits<ElementType>::idiv(mat.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::mat_idiv(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<ElementType>& result)
        {
            result.elem.resize(mat.size());

            result.row_ptr = mat.row_ptr;
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                result.elem[i] = typename SparseMatrix<ElementType>::Entry(mat.elem[i].col_idx, ElementTraits<ElementType>::idiv(mat.elem[i].value, value));
            }
        }

        template<typename ElementType>
        void Operations::mat_idiv(SparseMatrix<ElementType>& mat, const ElementType& value)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                mat.elem[i].value = ElementTraits<ElementType>::idiv(mat.elem[i].value, value);
            }
        }

        template<typename ElementType>
        void Operations::mat_mod(const Matrix<ElementType>& mat, const ElementType& value, Matrix<ElementType>& result)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                result.elem[i] = ElementTraits<ElementType>::mod(mat.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::mat_mod(Matrix<ElementType>& mat, const ElementType& value)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                mat.elem[i] = ElementTraits<ElementType>::mod(mat.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::mat_mod(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<ElementType>& result)
        {
            result.elem.resize(mat.size());

            result.row_ptr = mat.row_ptr;
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                result.elem[i] = typename SparseMatrix<ElementType>::Entry(mat.elem[i].col_idx, ElementTraits<ElementType>::mod(mat.elem[i].value, value));
            }
        }

        template<typename ElementType>
        void Operations::mat_mod(SparseMatrix<ElementType>& mat, const ElementType& value)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                mat.elem[i].value = ElementTraits<ElementType>::mod(mat.elem[i].value, value);
            }
        }

        template<typename ElementType>
        void Operations::mat_mods(const Matrix<ElementType>& mat, const ElementType& value, Matrix<ElementType>& result)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                result.elem[i] = ElementTraits<ElementType>::mods(mat.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::mat_mods(Matrix<ElementType>& mat, const ElementType& value)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                mat.elem[i] = ElementTraits<ElementType>::mods(mat.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::mat_mods(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<ElementType>& result)
        {
            result.elem.resize(mat.size());

            result.row_ptr = mat.row_ptr;
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                result.elem[i] = typename SparseMatrix<ElementType>::Entry(mat.elem[i].col_idx, ElementTraits<ElementType>::mods(mat.elem[i].value, value));
            }
        }

        template<typename ElementType>
        void Operations::mat_mods(SparseMatrix<ElementType>& mat, const ElementType& value)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                mat.elem[i].value = ElementTraits<ElementType>::mods(mat.elem[i].value, value);
            }
        }

        template<typename ElementType>
        void Operations::mat_div(const Matrix<ElementType>& mat, const ElementType& value, Matrix<typename ElementTraits<ElementType>::RationalType>& result)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                result.elem[i] = ElementTraits<ElementType>::div(mat.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::mat_div(Matrix<ElementType>& mat, const ElementType& value)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                mat.elem[i] = ElementTraits<ElementType>::div(mat.elem[i], value);
            }
        }

        template<typename ElementType>
        void Operations::mat_div(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<typename ElementTraits<ElementType>::RationalType>& result)
        {
            result.elem.resize(mat.size());

            result.row_ptr = mat.row_ptr;
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                result.elem[i] = typename SparseMatrix<ElementType>::Entry(mat.elem[i].col_idx, ElementTraits<ElementType>::div(mat.elem[i].value, value));
            }
        }

        template<typename ElementType>
        void Operations::mat_div(SparseMatrix<ElementType>& mat, const ElementType& value)
        {
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                mat.elem.value[i] = ElementTraits<ElementType>::div(mat.elem[i].value, value);
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const Matrix<ElementType>& mat, const Vector<ElementType>& vct, Vector<ElementType>& result)
        {
            ASSERT_EXCEPTION(mat.cols == vct.length, std::length_error);

            for (unsigned int idxA = 0, idxR = 0; idxA < mat.size(); ++idxR)
            {
                ElementType prod = ElementTraits<ElementType>::zero();
                for (unsigned int idxB = 0; idxB < vct.length; ++idxB, ++idxA)
                {
                    prod += mat.elem[idxA] * vct.elem[idxB];
                }
                result.elem[idxR] = prod;
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const Matrix<ElementType>& mat, const Vector<ElementType>& vct, SparseVector<ElementType>& result)
        {
            ASSERT_EXCEPTION(mat.cols == vct.length, std::length_error);

            result.elem.clear();

            for (unsigned int idxA = 0, idxR = 0; idxA < mat.size(); ++idxR)
            {
                ElementType prod = ElementTraits<ElementType>::zero();
                for (unsigned int idxB = 0; idxB < vct.length; ++idxB, ++idxA)
                {
                    prod += mat.elem[idxA] * vct.elem[idxB];
                }
                result.push(idxR, prod);
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const Matrix<ElementType>& mat, const SparseVector<ElementType>& vct, Vector<ElementType>& result)
        {
            ASSERT_EXCEPTION(mat.cols == vct.length, std::length_error);

            for (unsigned int idxA = 0, idxR = 0; idxA < mat.size(); ++idxR, idxA += mat.cols)
            {
                ElementType prod = ElementTraits<ElementType>::zero();
                for (unsigned int idxB = 0; idxB < vct.size(); ++idxB)
                {
                    prod += mat.elem[idxA + vct.elem[idxB].idx] * vct.elem[idxB].value;
                }
                result.elem[idxR] = prod;
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const Matrix<ElementType>& mat, const SparseVector<ElementType>& vct, SparseVector<ElementType>& result)
        {
            ASSERT_EXCEPTION(mat.cols == vct.length, std::length_error);

            result.elem.clear();

            for (unsigned int idxA = 0, idxR = 0; idxA < mat.size(); ++idxR, idxA += mat.cols)
            {
                ElementType prod = ElementTraits<ElementType>::zero();
                for (unsigned int idxB = 0; idxB < vct.size(); ++idxB)
                {
                    prod += mat.elem[idxA + vct.elem[idxB].idx] * vct.elem[idxB].value;
                }
                result.push(idxR, prod);
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const SparseMatrix<ElementType>& mat, const Vector<ElementType>& vct, Vector<ElementType>& result)
        {
            ASSERT_EXCEPTION(mat.cols == vct.length, std::length_error);

            for (unsigned int row = 0; row < mat.rows; ++row)
            {
                ElementType prod = ElementTraits<ElementType>::zero();
                for (unsigned int idxA = mat.row_ptr[row]; idxA < mat.row_ptr[row + 1]; ++idxA)
                {
                    prod += mat.elem[idxA].value * vct.elem[mat.elem[idxA].col_idx];
                }
                result.elem[row] = prod;
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const SparseMatrix<ElementType>& mat, const Vector<ElementType>& vct, SparseVector<ElementType>& result)
        {
            ASSERT_EXCEPTION(mat.cols == vct.length, std::length_error);

            result.elem.clear();

            for (unsigned int row = 0; row < mat.rows; ++row)
            {
                ElementType prod = ElementTraits<ElementType>::zero();
                for (unsigned int idxA = mat.row_ptr[row]; idxA < mat.row_ptr[row + 1]; ++idxA)
                {
                    prod += mat.elem[idxA].value * vct.elem[mat.elem[idxA].col_idx];
                }
                result.push(row, prod);
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const SparseMatrix<ElementType>& mat, const SparseVector<ElementType>& vct, Vector<ElementType>& result)
        {
            ASSERT_EXCEPTION(mat.cols == vct.length, std::length_error);

            for (unsigned int row = 0; row < mat.rows; ++row)
            {
                ElementType prod = ElementTraits<ElementType>::zero();
                for (unsigned int idxA = mat.row_ptr[row], idxB = 0; idxA < mat.row_ptr[row + 1] && idxB < vct.size();)
                {
                    if (mat.elem[idxA].col_idx < vct.elem[idxB].idx)
                    {
                        ++idxA;
                    }
                    else if (mat.elem[idxA].col_idx > vct.elem[idxB].idx)
                    {
                        ++idxB;
                    }
                    else
                    {
                        prod += mat.elem[idxA].value * vct.elem[idxB].value;
                        ++idxA;
                        ++idxB;
                    }
                }
                result.elem[row] = prod;
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const SparseMatrix<ElementType>& mat, const SparseVector<ElementType>& vct, SparseVector<ElementType>& result)
        {
            ASSERT_EXCEPTION(mat.cols == vct.length, std::length_error);

            result.elem.clear();

            for (unsigned int row = 0; row < mat.rows; ++row)
            {
                ElementType prod = ElementTraits<ElementType>::zero();
                for (unsigned int idxA = mat.row_ptr[row], idxB = 0; idxA < mat.row_ptr[row + 1] && idxB < vct.size();)
                {
                    if (mat.elem[idxA].col_idx < vct.elem[idxB].idx)
                    {
                        ++idxA;
                    }
                    else if (mat.elem[idxA].col_idx > vct.elem[idxB].idx)
                    {
                        ++idxB;
                    }
                    else
                    {
                        prod += mat.elem[idxA].value * vct.elem[idxB].value;
                        ++idxA;
                        ++idxB;
                    }
                }
                result.push(row, prod);
            }
        }

        template<typename ElementType>
        void Operations::mat_add(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            for (unsigned int idxR = 0; idxR < op1.size(); ++idxR)
            {
                result.elem[idxR] = op1.elem[idxR] + op2.elem[idxR];
            }
        }

        template<typename ElementType>
        void Operations::mat_add(Matrix<ElementType>& op1, const Matrix<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            for (unsigned int idxR = 0; idxR < op1.size(); ++idxR)
            {
                op1.elem[idxR] += op2.elem[idxR];
            }
        }

        template<typename ElementType>
        void Operations::mat_add(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, Matrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            for (unsigned int idxR = 0; idxR < op1.size(); ++idxR)
            {
                result.elem[idxR] = op1.elem[idxR];
            }
            for (unsigned int rowB = 0, idxR = 0; rowB < op2.rows; ++rowB, idxR += result.cols)
            {
                for (unsigned int idxB = op2.row_ptr[rowB]; idxB < op2.row_ptr[rowB + 1]; ++idxB)
                {
                    result.elem[idxR + op2.elem[idxB].col_idx] += op2.elem[idxB].value;
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_add(Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            for (unsigned int rowB = 0, idxR = 0; rowB < op2.rows; ++rowB, idxR += op1.cols)
            {
                for (unsigned int idxB = op2.row_ptr[rowB]; idxB < op2.row_ptr[rowB + 1]; ++idxB)
                {
                    op1.elem[idxR + op2.elem[idxB].col_idx] += op2.elem[idxB].value;
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_add(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            for (unsigned int idxR = 0; idxR < op2.size(); ++idxR)
            {
                result.elem[idxR] = op2.elem[idxR];
            }
            for (unsigned int rowA = 0, idxR = 0; rowA < op1.rows; ++rowA, idxR += result.cols)
            {
                for (unsigned int idxA = op1.row_ptr[rowA]; idxA < op1.row_ptr[rowA + 1]; ++idxA)
                {
                    result.elem[idxR + op1.elem[idxA].col_idx] += op1.elem[idxA].value;
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_add(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, SparseMatrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            result.elem.clear();

            for (unsigned int row = 0; row < op1.rows; ++row)
            {
                unsigned int idxA = op1.row_ptr[row];
                unsigned int idxB = op2.row_ptr[row];
                while (idxA < op1.row_ptr[row + 1] && idxB < op2.row_ptr[row + 1])
                {
                    if (op1.elem[idxA].col_idx < op2.elem[idxB].col_idx)
                    {
                        result.push(op1.elem[idxA].col_idx, op1.elem[idxA].value);
                        ++idxA;
                    }
                    else if (op1.elem[idxA].col_idx > op2.elem[idxB].col_idx)
                    {
                        result.push(op2.elem[idxB].col_idx, op2.elem[idxB].value);
                        ++idxB;
                    }
                    else
                    {
                        result.push(op1.elem[idxA].col_idx, op1.elem[idxA].value + op2.elem[idxB].value);
                        ++idxA;
                        ++idxB;
                    }
                }
                while (idxA < op1.row_ptr[row + 1])
                {
                    result.push(op1.elem[idxA].col_idx, op1.elem[idxA].value);
                    ++idxA;
                }
                while (idxB < op2.row_ptr[row + 1])
                {
                    result.push(op2.elem[idxB].col_idx, op2.elem[idxB].value);
                    ++idxB;
                }
                result.row_ptr[row + 1] = result.size();
            }
        }

        template<typename ElementType>
        void Operations::mat_sub(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            for (unsigned int idxR = 0; idxR < op1.size(); ++idxR)
            {
                result.elem[idxR] = op1.elem[idxR] - op2.elem[idxR];
            }
        }

        template<typename ElementType>
        void Operations::mat_sub(Matrix<ElementType>& op1, const Matrix<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            for (unsigned int idxR = 0; idxR < op1.size(); ++idxR)
            {
                op1.elem[idxR] -= op2.elem[idxR];
            }
        }

        template<typename ElementType>
        void Operations::mat_sub(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, Matrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            for (unsigned int idxR = 0; idxR < op1.size(); ++idxR)
            {
                result.elem[idxR] = op1.elem[idxR];
            }
            for (unsigned int rowB = 0, idxR = 0; rowB < op2.rows; ++rowB, idxR += result.cols)
            {
                for (unsigned int idxB = op2.row_ptr[rowB]; idxB < op2.row_ptr[rowB + 1]; ++idxB)
                {
                    result.elem[idxR + op2.elem[idxB].col_idx] -= op2.elem[idxB].value;
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_sub(Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            for (unsigned int rowB = 0, idxR = 0; rowB < op2.rows; ++rowB, idxR += op1.cols)
            {
                for (unsigned int idxB = op2.row_ptr[rowB]; idxB < op2.row_ptr[rowB + 1]; ++idxB)
                {
                    op1.elem[idxR + op2.elem[idxB].col_idx] -= op2.elem[idxB].value;
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_sub(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            for (unsigned int idxR = 0; idxR < op2.size(); ++idxR)
            {
                result.elem[idxR] = -op2.elem[idxR];
            }
            for (unsigned int rowA = 0, idxR = 0; rowA < op1.rows; ++rowA, idxR += result.cols)
            {
                for (unsigned int idxA = op1.row_ptr[rowA]; idxA < op1.row_ptr[rowA + 1]; ++idxA)
                {
                    result.elem[idxR + op1.elem[idxA].col_idx] += op1.elem[idxA].value;
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_sub(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, SparseMatrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            result.elem.clear();

            for (unsigned int row = 0; row < op1.rows; ++row)
            {
                unsigned int idxA = op1.row_ptr[row];
                unsigned int idxB = op2.row_ptr[row];
                while (idxA < op1.row_ptr[row + 1] && idxB < op2.row_ptr[row + 1])
                {
                    if (op1.elem[idxA].col_idx < op2.elem[idxB].col_idx)
                    {
                        result.push(op1.elem[idxA].col_idx, op1.elem[idxA].value);
                        ++idxA;
                    }
                    else if (op1.elem[idxA].col_idx > op2.elem[idxB].col_idx)
                    {
                        result.push(op2.elem[idxB].col_idx, -op2.elem[idxB].value);
                        ++idxB;
                    }
                    else
                    {
                        result.push(op1.elem[idxA].col_idx, op1.elem[idxA].value - op2.elem[idxB].value);
                        ++idxA;
                        ++idxB;
                    }
                }
                while (idxA < op1.row_ptr[row + 1])
                {
                    result.push(op1.elem[idxA].col_idx, op1.elem[idxA].value);
                    ++idxA;
                }
                while (idxB < op2.row_ptr[row + 1])
                {
                    result.push(op2.elem[idxB].col_idx, -op2.elem[idxB].value);
                    ++idxB;
                }
                result.row_ptr[row + 1] = result.size();
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.cols == op2.rows, std::length_error);

            // TODO: Improve indexing
            for (unsigned int rowR = 0, idxR = 0; rowR < op1.rows; ++rowR)
            {
                for (unsigned int colR = 0; colR < op2.cols; ++colR, ++idxR)
                {
                    ElementType prod = ElementTraits<ElementType>::zero();
                    for (unsigned int idx = 0, idxA = rowR * op1.cols, idxB = colR; idx < op1.cols; ++idx, ++idxA, idxB += op2.cols)
                    {
                        prod += op1.elem[idxA] * op2.elem[idxB];
                    }
                    result.elem[idxR] = prod;
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2, SparseMatrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.cols == op2.rows, std::length_error);

            result.elem.clear();

            for (unsigned int rowR = 0; op1.rowR < op1.rows; ++rowR)
            {
                for (unsigned int colR = 0; colR < op2.cols; ++colR)
                {
                    ElementType prod = ElementTraits<ElementType>::zero();
                    for (unsigned int idx = 0, idxA = rowR * op1.cols, idxB = colR; idx < op1.cols; ++idx, ++idxA, idxB += op2.cols)
                    {
                        prod += op1.elem[idxA] * op2.elem[idxB];
                    }
                    result.push(colR, prod);
                }
                result.row_ptr[rowR + 1] = result.size();
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, Matrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.cols == op2.rows, std::length_error);

            result.elem.clear();

            for (unsigned int rowA = 0, idxA = 0, idxR = 0; rowA < op1.rows; ++rowA, idxR += op2.cols)
            {
                for (unsigned int colA = 0; colA < op1.cols; ++colA, ++idxA)
                {
                    for (unsigned int idxB = op2.row_ptr[colA]; idxB < op2.row_ptr[colA + 1]; ++idxB)
                    {
                        result.elem[idxR + op2.elem[idxB].col_idx] += op1.elem[idxA] * op2.elem[idxB].value;
                    }
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, SparseMatrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.cols == op2.rows, std::length_error);

            result.elem.clear();

            std::vector<ElementType> cache(op2.cols);
            for (unsigned int rowA = 0, idxA = 0; rowA < op1.rows; ++rowA)
            {
                std::fill(cache.begin(), cache.end(), ElementTraits<ElementType>::zero());
                for (unsigned int colA = 0; colA < op1.cols; ++colA, ++idxA)
                {
                    for (unsigned int idxB = op2.row_ptr[colA]; idxB < op2.row_ptr[colA + 1]; ++idxB)
                    {
                        cache[op2.elem[idxB].col_idx] += op1.elem[idxA] * op2.elem[idxB].value;
                    }
                }
                for (unsigned int idxC = 0; idxC < op2.cols; ++idxC)
                {
                    result.push(idxC, cache[idxC]);
                }
                result.row_ptr[rowA + 1] = result.size();
            }

        }

        template<typename ElementType>
        void Operations::mat_mul(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.cols == op2.rows, std::length_error);

            std::fill(result.elem.begin(), result.elem.end(), ElementTraits<ElementType>::zero());

            for (unsigned int row = 0, idxR = 0; row < op1.rows; ++row, idxR += op2.cols)
            {
                for (unsigned int idxA = op1.row_ptr[row]; idxA < op1.row_ptr[row + 1]; ++idxA)
                {
                    for (unsigned int idxRR = idxR, idxB = op1.elem[idxA].col_idx * op2.cols; idxRR < idxR + op2.cols; ++idxRR, ++idxB)
                    {
                        result.elem[idxRR] += op1.elem[idxA].value * op2.elem[idxB];
                    }
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2, SparseMatrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.cols == op2.rows, std::length_error);

            result.elem.clear();

            std::vector<ElementType> cache(op2.cols);
            for (unsigned int row = 0; row < op1.rows; ++row)
            {
                std::fill(cache.begin(), cache.end(), ElementTraits<ElementType>::zero());
                for (unsigned int idxA = op1.row_ptr[row]; idxA < op1.row_ptr[row + 1]; ++idxA)
                {
                    unsigned int idxB = op1.elem[idxA].col_idx * op2.cols;
                    for (unsigned int idxC = 0; idxC < op2.cols; ++idxC, ++idxB)
                    {
                        cache[idxC] += op1.elem[idxA].value * op2.elem[idxB];
                    }
                }
                for (unsigned int idxC = 0; idxC < op2.cols; ++idxC)
                {
                    result.push(idxC, cache[idxC]);
                }
                result.row_ptr[row + 1] = result.size();
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, Matrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.cols == op2.rows, std::length_error);

            std::fill(result.elem.begin(), result.elem.end(), ElementTraits<ElementType>::zero());

            for (unsigned int row = 0, idxR = 0; row < op1.rows; ++row, idxR += op2.cols)
            {
                for (unsigned int idxA = op1.row_ptr[row]; idxA < op1.row_ptr[row + 1]; ++idxA)
                {
                    for (unsigned int idxB = op2.row_ptr[op1.elem[idxA].col_idx]; idxB < op2.row_ptr[op1.elem[idxA].col_idx + 1]; ++idxB)
                    {
                        result.elem[idxR + op2.elem[idxB].col_idx] += op1.elem[idxA].value * op2.elem[idxB].value;
                    }
                }
            }
        }

        template<typename ElementType>
        void Operations::mat_mul(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, SparseMatrix<ElementType>& result)
        {
            ASSERT_EXCEPTION(op1.cols == op2.rows, std::length_error);

            std::vector<ElementType> cache(op2.cols);
            for (unsigned int row = 0; row < op1.rows; ++row)
            {
                std::fill(cache.begin(), cache.end(), ElementTraits<ElementType>::zero());
                for (unsigned int idxA = op1.row_ptr[row]; idxA < op1.row_ptr[row + 1]; ++idxA)
                {
                    for (unsigned int idxB = op2.row_ptr[op1.elem[idxA].col_idx]; idxB < op2.row_ptr[op1.elem[idxA].col_idx + 1]; ++idxB)
                    {
                        cache[op2.elem[idxB].col_idx] += op1.elem[idxA].value * op2.elem[idxB].value;
                    }
                }
                for (unsigned int idxC = 0; idxC < op2.cols; ++idxC)
                {
                    result.push(idxC, cache[idxC]);
                }
                result.row_ptr[row + 1] = result.size();
            }
        }

        template<typename ElementType>
        bool Operations::mat_eq(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            bool l = true;
            for (unsigned int i = 0; i < op1.size() && l; ++i)
            {
                l = (op1.elem[i] == op2.elem[i]);
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::mat_eq(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            bool l = true;
            for (unsigned int row = 0, idxA = 0; op1.row < op1.rows && l; ++row)
            {
                unsigned int colA = 0;
                unsigned int idxB = op2.row_ptr[row];
                while (colA < op1.cols && idxB < op2.row_ptr[row + 1] && l)
                {
                    if (colA == op2.elem[idxB].col_idx)
                    {
                        l = (op1.elem[idxA] == op2.elem[idxB].value);
                        ++colA;
                        ++idxA;
                        ++idxB;
                    }
                    else
                    {
                        l = (op1.elem[idxA] == ElementTraits<ElementType>::zero());
                        ++colA;
                        ++idxA;
                    }
                }
                while (colA < op1.cols && l)
                {
                    l = (op1.elem[idxA] == ElementTraits<ElementType>::zero());
                    ++colA;
                    ++idxA;
                }
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::mat_eq(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            bool l = true;
            for (unsigned int row = 0, idxB = 0; row < op1.rows && l; ++row)
            {
                unsigned int idxA = op2.row_ptr[row];
                unsigned int colB = 0;
                while (idxA < op1.row_ptr[row + 1] && colB < op2.cols && l)
                {
                    if (op1.elem[idxA].col_idx == colB)
                    {
                        l = (op1.elem[idxA].value == op2.elem[idxB]);
                        ++idxA;
                        ++colB;
                        ++idxB;
                    }
                    else
                    {
                        l = (ElementTraits<ElementType>::zero() == op2.elem[idxB]);
                        ++colB;
                        ++idxB;
                    }
                }
                while (colB < op2.cols && l)
                {
                    l = (ElementTraits<ElementType>::zero() == op2.elem[idxB]);
                    ++colB;
                    ++idxB;
                }
            }

            return l;
        }

        template<typename ElementType>
        bool Operations::mat_eq(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2)
        {
            ASSERT_EXCEPTION(op1.rows == op2.rows && op1.cols == op2.cols, std::length_error);

            bool l = (op1.row_ptr[0] == op2.row_ptr[0]);
            for (unsigned int row = 0; row < op1.rows && l; ++row)
            {
                l = (op1.row_ptr[row + 1] == op2.row_ptr[row + 1]);
                for (unsigned int idx = op1.row_ptr[row]; idx < op1.row_ptr[row + 1] && l; ++idx)
                {
                    l = ((op1.elem[idx].col_idx == op2.elem[idx].col_idx) && (op1.elem[idx].value == op2.elem[idx].value));
                }
            }

            return l;
        }

    }
}
