#include <vector>

namespace GeNuSys
{
    namespace Algebra
    {

        // http://en.wikipedia.org/wiki/Lehmer%E2%80%93Schur_algorithm

        template<typename ElementType>
        bool lehmerSchur(const Polynomial<ElementType>& poly)
        {
            unsigned int idxPrev = 0, idxAct = 1;
            std::vector<ElementType> ps[2] = { std::vector<ElementType>(poly.getDegree() + 1), std::vector<ElementType>(poly.getDegree() + 1) };
            for (int i = 0; i <= poly.getDegree(); ++i)
            {
                ps[idxAct][i] = poly[i];
            }

            for (unsigned int deg = poly.getDegree(); deg > 0 && ps[idxAct][0] > 0; --deg)
            {
                std::swap(idxPrev, idxAct);

                ElementType lcoef = ps[idxPrev][deg];
                ElementType tcoef = ps[idxPrev][0];
                for (unsigned int k = 0; k < deg; ++k)
                {
                    ps[idxAct][k] = ps[idxPrev][k] * tcoef - ps[idxPrev][deg - k] * lcoef;
                }

                while (deg > 0 && ps[idxAct][deg] == ElementTraits<ElementType>::zero())
                {
                    --deg;
                }
            }

            return ps[idxPrev][0] > 0;
        }

    }
}
