#ifndef GENUSYS_UTILS_H_
#define GENUSYS_UTILS_H_

#ifndef NDEBUG
#define ASSERT_EXCEPTION(x,y) if(!(x)) throw y{#x}
#else
#define ASSERT_EXCEPTION(x,y)
#endif

#endif // GENUSYS_UTILS_H_
