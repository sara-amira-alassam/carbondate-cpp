cmake_minimum_required(VERSION 3.24)
project(carbondate)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(Rmath)

add_executable(carbondate
        main.cpp src/WalkerDPMM.cpp src/WalkerDPMM.h src/helpers.cpp src/sort.cpp src/helpers.h src/csv_helpers.cpp src/csv_helpers.h src/DensityOutput.cpp src/DensityOutput.h src/read_data.cpp src/read_data.h)
target_include_directories(carbondate PUBLIC Rmath/include)
target_link_libraries(carbondate PUBLIC Rmath)
