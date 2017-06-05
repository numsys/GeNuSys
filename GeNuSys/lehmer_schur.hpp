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
