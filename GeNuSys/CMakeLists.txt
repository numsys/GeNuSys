cmake_minimum_required(VERSION 2.4.0)

project(GeNuSys)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(SRCS
bit_vector.h
digit_set.h
digit_set.hpp
element_traits_complex.hpp
element_traits_double.hpp
element_traits.h
element_traits.hpp
element_traits_int.hpp
element_traits_long_long.hpp
element_traits_mpf.hpp
element_traits_mpq.hpp
element_traits_mpz.hpp
frobenius_norm.h
frobenius_norm.hpp
hash_table.h
hash_table.hpp
lehmer_schur.h
lehmer_schur.hpp
linalg_algorithms.h
linalg_algorithms.hpp
linalg_operations.h
linalg_operations.hpp
linalg_traits.h
linalg_traits.hpp
matrix.h
matrix.hpp
number_system.h
number_system.hpp
numsys_traits.h
numsys_traits.hpp
operator_norm.h
operator_norm.hpp
p_norm.h
p_norm.hpp
polynomial.h
polynomial.hpp
radix_properties.h
radix_properties.hpp
simultan.h
simultan.hpp
smith_hash.h
smith_hash.hpp
sparse_matrix.h
sparse_matrix.hpp
sparse_vector.h
sparse_vector.hpp
test_runner.h
tests.h
test_suite.h
test_utils.h
test_utils.hpp
vector_coder.h
vector_coder.hpp
vector.h
vector.hpp
main.cpp
test_runner.cpp
test_suite.cpp
utils.h
)

add_definitions(-std=c++11 -Wall -Wextra -Woverloaded-virtual -pedantic -rdynamic)

add_executable(GeNuSys ${SRCS})
target_link_libraries(GeNuSys -lgmp)
