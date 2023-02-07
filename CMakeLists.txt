cmake_minimum_required(VERSION 3.24)
project(carbondate)

set(CMAKE_CXX_STANDARD 11)

add_executable(carbondate
        main.cpp src/WalkerDPMM.cpp src/WalkerDPMM.h src/helpers.cpp src/sort.cpp src/helpers.h src/csv_helpers.cpp src/csv_helpers.h)
target_include_directories(carbondate PUBLIC /usr/local/include)
target_link_libraries(carbondate /Users/sara/Downloads/R-4.2.2/src/nmath/standalone/libRmath.a)
