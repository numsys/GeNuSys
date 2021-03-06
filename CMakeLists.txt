cmake_minimum_required(VERSION 2.4.0)

project(GeNuSys)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

find_package(GMP)
find_package(GeNuSys REQUIRED)
find_package(Threads) #Not required if threading is disabled by defining the GENUSYS_NO_THREADING macro

include_directories(${GeNuSys_INCLUDE_DIRS})

set(test_SRCS
tests/test_runner.cpp
tests/test_suite.cpp
tests/main.cpp
)

set(example0_SRCS
examples/example0.cpp
)

set(example1_SRCS
examples/example1.cpp
)

set(cyclotomic_SRCS
examples/cyclotomic.cpp
)

add_definitions(-std=c++11 -Wall -Wextra -Woverloaded-virtual -pedantic -rdynamic)

add_executable(GeNuSys-example0 ${example0_SRCS})
add_executable(GeNuSys-example1 ${example1_SRCS})
add_executable(GeNuSys-tests ${test_SRCS})
add_executable(GeNuSys-cyclotomic ${cyclotomic_SRCS})
target_link_libraries(GeNuSys-example0 ${GMP_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(GeNuSys-example1 ${GMP_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(GeNuSys-tests ${GMP_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(GeNuSys-cyclotomic ${GMP_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

install(DIRECTORY ./GeNuSys/ DESTINATION include/GeNuSys FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ./GeNuSys/ DESTINATION include/GeNuSys FILES_MATCHING PATTERN "*.hpp")
