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

#ifndef GENUSYS_NUMSYS_BIT_VECTOR_H_
#define GENUSYS_NUMSYS_BIT_VECTOR_H_

#include <vector>

namespace GeNuSys
{
    namespace NumSys
    {

        class BitVector
        {

            private:

                std::vector<unsigned long long> table;

            public:

                BitVector(unsigned long long size): table(size / (sizeof(unsigned long long) * 8) + 1, 0) { }

                bool operator [](unsigned long long idx) const
                {
                    unsigned int d = (int)(idx / (sizeof(unsigned long long) * 8));
                    unsigned int m = (int)(idx % (sizeof(unsigned long long) * 8));

                    return (table[d] & ((unsigned long long) 1 << m)) != 0;
                }

                void set(unsigned long long idx)
                {
                    unsigned int d = (int)(idx / (sizeof(unsigned long long) * 8));
                    unsigned int m = (int)(idx % (sizeof(unsigned long long) * 8));

                    table[d] |= ((unsigned long long) 1 << m);
                }

                void unset(unsigned long long idx)
                {
                    unsigned int d = (int)(idx / (sizeof(unsigned long long) * 8));
                    unsigned int m = (int)(idx % (sizeof(unsigned long long) * 8));

                    table[d] &= ~((unsigned long long) 1 << m);
                }

        };

    }
}

#endif // GENUSYS_NUMSYS_BIT_VECTOR_H_
