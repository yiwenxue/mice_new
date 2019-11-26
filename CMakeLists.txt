cmake_minimum_required(VERSION 2.8.9)
set(CMAKE_CXX_COMPILER clang++)
set( CMAKE_CXX_FLAGS "-Wall -m64 -O2" )
set(CXX_STANDARD 11)
project(main)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include_directories(include)

add_library(mylib  src/utility.cpp src/mathematics.cpp)

set(SRC src)

add_executable(main ${SRC}/main.cpp)
add_executable(run ${SRC}/new.cpp)
add_executable(doc ${SRC}/doc.cpp)

FIND_PACKAGE(PkgConfig REQUIRED)

FIND_PACKAGE(GSL REQUIRED)
LINK_DIRECTORIES(${GSL_LIBRARIES})
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})

FIND_PACKAGE(FFTW3 REQUIRED)
include_directories(${FFTW3_INCLUDE_DIR})
LINK_DIRECTORIES(${FFTW3_LIBRARIES})

TARGET_LINK_LIBRARIES(main mylib m ${GSL_LIBRARIES} ${GSL_CLBAS_LIBRARY} ${FFTW3_LIBRARIES} )
