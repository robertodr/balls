find_package(Eigen3 3.4 REQUIRED)

add_library(balls SHARED LeopardiPartition.cpp lebedev/Lebedev.cpp lebedev/detail/lebedev_utils.cpp)

target_include_directories(balls PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(balls PUBLIC cxx_std_17)
target_link_libraries(balls PUBLIC Eigen3::Eigen)

add_executable(example example.cpp)

target_link_libraries(example PRIVATE balls)
