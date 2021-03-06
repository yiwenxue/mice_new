cmake_minimum_required(VERSION 2.8.9)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_FLAGS "-Wall -m64 -O2" )
set(CMAKE_CXX_STANDARD 17)
project(main)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include_directories(src/include)

add_library(mylib STATIC src/utility.cpp src/mathematics.cpp)

add_executable(main src/main.cpp)
add_executable(run src/run.cpp)
add_executable(doc src/doc.cpp)
add_executable(test src/test.cpp)

FIND_PACKAGE(PkgConfig REQUIRED)

FIND_PACKAGE(ROOT REQUIRED)
INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIRS})
LINK_DIRECTORIES(${ROOT_LIBRARIES})

FIND_PACKAGE(GSL REQUIRED)
LINK_DIRECTORIES(${GSL_LIBRARIES})
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})

FIND_PACKAGE(FFTW3 REQUIRED)
include_directories(${FFTW3_INCLUDE_DIR})
LINK_DIRECTORIES(${FFTW3_LIBRARIES})


TARGET_LINK_LIBRARIES(main mylib m ${GSL_LIBRARIES} ${GSL_CLBAS_LIBRARY} ${FFTW3_LIBRARIES} ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(test m ${ROOT_LIBRARIES})
