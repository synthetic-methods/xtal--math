cmake_minimum_required (VERSION 3.20)
project (test_package)

set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package (xtal-math CONFIG REQUIRED)
find_package (Catch2    CONFIG REQUIRED)
find_package (Eigen3    CONFIG REQUIRED)

cmake_path (GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH xtal-math_SOURCE_DIR)

if (EXISTS ${xtal-math_SOURCE_DIR}/include/FIXME.cc)
	add_executable (test ${xtal-math_SOURCE_DIR}/include/FIXME.cc)
else ()
	add_executable (test ${xtal-math_SOURCE_DIR}/include/xtal/math/all.cc)
endif ()

target_link_libraries (test xtal-math::xtal-math Catch2::Catch2WithMain Eigen3::Eigen)

target_precompile_headers (test INTERFACE xtal-math::xtal-math Catch2::Catch2WithMain Eigen3::Eigen)
