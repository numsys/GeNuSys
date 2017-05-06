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
